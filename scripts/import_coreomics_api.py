#!/usr/bin/env python3
"""
import_coreomics_api.py — Load Adam Olshen's coreomics submissions directly
from the live REST API into the DE-LIMP PG database.

Supersedes import_coreomics_xlsx.py for routine sync: hitting the API is
real-time, idempotent, doesn't require Adam to email a fresh export, and
captures per-sample condition_name / unique_id directly from
`submission_data.samples`. Use the xlsx importer only for offline one-shot
backfills.

USAGE
-----
    # Token: in env or file at ~/.coreomics_token
    PGPASSWORD=$(cat /Volumes/proteomics-grp/brett/.pgfarm_delimp_token) \\
    COREOMICS_TOKEN=$(cat ~/.coreomics_token) \\
    python3 scripts/import_coreomics_api.py \\
        --host pgfarm.library.ucdavis.edu \\
        --user brettsp \\
        --database uc-davis-genome-center-proteomics-core/delimp \\
        --auto-link-raws

    # Limit run (testing — fetch only N submissions)
    python3 scripts/import_coreomics_api.py --max-submissions 50 --dry-run

ENDPOINTS
---------
    GET <base>/submissions/?lab=PROTEOMICS&page=N&page_size=K   (paginated list)
    GET <base>/submissions/<id>/                                (detail w/ samples)

PG TARGET
---------
Same schema as the xlsx importer:
    coreomics_submissions_cache
    coreomics_samples_cache
    delimp_raw_files_internal (with --auto-link-raws)

Idempotent UPSERTs by (submission_id) and (submission_id, unique_id).
Re-run as often as you like; only changed rows write.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path

import psycopg2
import psycopg2.extras


DEFAULT_BASE_URL = "https://ucdavis.coreomics.com/server/api"
DEFAULT_LAB = "PROTEOMICS"


# ────────────────────────────────────────────────────────────────────────────
# Token + HTTP
# ────────────────────────────────────────────────────────────────────────────

def load_token(token_file: str | None) -> str:
    """Resolve the coreomics token: --token-file > COREOMICS_TOKEN env > ~/.coreomics_token"""
    if token_file:
        return Path(token_file).read_text().strip()
    env = os.environ.get("COREOMICS_TOKEN")
    if env:
        return env.strip()
    default = Path.home() / ".coreomics_token"
    if default.exists():
        return default.read_text().strip()
    sys.exit("FATAL: no coreomics token. Set COREOMICS_TOKEN env, "
             "create ~/.coreomics_token, or pass --token-file.")


def api_get(base_url: str, path: str, token: str, params: dict | None = None,
            timeout: int = 30, retries: int = 3) -> dict:
    """GET a coreomics endpoint with auth. Returns parsed JSON. Retries on transient errors."""
    url = f"{base_url.rstrip('/')}/{path.lstrip('/')}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    last_exc = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Authorization": f"Token {token}",
                "Accept": "application/json",
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code in (401, 403):
                sys.exit(f"FATAL: HTTP {e.code} on {url} — token expired or unauthorized.")
            if e.code == 404:
                return None
            last_exc = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_exc = e
        if attempt < retries - 1:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last_exc}")


def list_all_submissions(base_url: str, token: str, lab: str,
                         page_size: int = 100, max_pages: int = 100,
                         verbose: bool = True) -> list:
    """Walk paginated submissions endpoint. Returns the full list."""
    all_subs = []
    for page in range(1, max_pages + 1):
        d = api_get(base_url, "submissions/", token,
                    params={"lab": lab, "page": page, "page_size": page_size})
        if not d or not d.get("results"):
            break
        all_subs.extend(d["results"])
        if verbose:
            print(f"  page {page}: {len(d['results'])} ({len(all_subs)}/{d.get('count', '?')})")
        if not d.get("next"):
            break
    return all_subs


def fetch_submission_detail(base_url: str, token: str, submission_id: str) -> dict:
    """Detail view contains submission_data.samples (not in the list view)."""
    return api_get(base_url, f"submissions/{submission_id}/", token)


# ────────────────────────────────────────────────────────────────────────────
# Mapping from API → PG rows
# ────────────────────────────────────────────────────────────────────────────

def coerce_ts(v):
    if v is None or v == "":
        return None
    if isinstance(v, datetime):
        return v
    try:
        # API returns ISO-8601 with offset like '2026-05-19T00:06:36.995000-07:00'
        return datetime.fromisoformat(str(v))
    except ValueError:
        return None


def coerce_date(v):
    ts = coerce_ts(v)
    if ts:
        return ts.date()
    if v and isinstance(v, str):
        # Try formats like '2026-05-22'
        for fmt in ("%Y-%m-%d", "%m/%d/%Y"):
            try:
                return datetime.strptime(v.strip(), fmt).date()
            except ValueError:
                pass
    return None


def coerce_int(v):
    if v in (None, ""):
        return None
    try:
        return int(float(str(v).strip()))
    except (ValueError, TypeError):
        return None


# These fields are PII we explicitly drop (not even kept in raw_payload).
DROP_API_FIELDS = {"phone", "pi_phone"}


def submission_to_pg_row(sub: dict, source: str) -> dict:
    """Map an API submission detail dict → coreomics_submissions_cache row."""
    sd = sub.get("submission_data") or {}

    # Pull selected fields out of submission_data into typed columns. Anything
    # we don't pull stays in raw_payload (minus the explicit PII drops).
    record = {
        "submission_id": sub.get("id"),
        "internal_id": sub.get("internal_id"),
        "type": _typename(sub.get("type")),
        "submitted_at": coerce_ts(sub.get("submitted")),
        "status": sub.get("status"),
        "send_date": coerce_date(sd.get("send_date")),

        "submitter_first_name": sub.get("first_name"),
        "submitter_last_name": sub.get("last_name"),
        "submitter_email": sub.get("email"),
        "pi_first_name": sub.get("pi_first_name"),
        "pi_last_name": sub.get("pi_last_name"),
        "pi_email": sub.get("pi_email"),
        "institute": sub.get("institute"),

        "num_samples": coerce_int(len(sd.get("samples") or [])),
        "organism": sd.get("organism"),
        "species": sd.get("species"),
        "prot_or_pep": sd.get("prot_or_pep"),
        "proteomics_type": sd.get("proteomics_type"),
        "mass_spec_wanted": sd.get("mass_spec_wanted"),
        "sample_prep": sd.get("sample_prep"),
        "gradient_length": sd.get("gradient_length"),
        "dia": sd.get("dia"),
        "tmt": sd.get("tmt"),
        "description": sd.get("description"),
        "other_info": sd.get("other_info"),

        "biohazard": sd.get("biohazard"),
        "pathogenic": sd.get("pathogenic"),
        "nih_s10_user": sd.get("nih_s10_user"),
        "is_nih_major_user": sd.get("is_nih_major_user"),
        "transgenic": sd.get("transgenic"),
        "po_account_number": sd.get("po_account_number"),

        "source_export": source,
        "source_export_md5": None,                       # not applicable for API
    }

    # raw_payload = everything left over (less PII)
    keep = {k: v for k, v in sub.items()
            if k not in ("id", "submitted", "status", "first_name", "last_name",
                         "email", "pi_first_name", "pi_last_name", "pi_email",
                         "institute", "submission_data", "internal_id", "type")
               and k not in DROP_API_FIELDS}
    if sd:
        # Keep the nested submission_data too (minus samples — those go to
        # the per-sample table — and minus PII duplicates)
        keep["submission_data"] = {k: v for k, v in sd.items()
                                    if k != "samples" and k not in DROP_API_FIELDS}
    record["raw_payload"] = json.dumps(keep, default=str)
    return record


def _typename(t):
    """`type` in API is a dict like {'id':..., 'name':'Proteomics', ...} — pluck the name."""
    if isinstance(t, dict):
        return t.get("name")
    return t


def samples_for_submission(sub: dict) -> list[dict]:
    """Extract per-sample rows from submission_data.samples."""
    sd = sub.get("submission_data") or {}
    samples = sd.get("samples") or []
    sub_id = sub.get("id")
    out = []
    for s in samples:
        if not isinstance(s, dict):
            continue
        unique_id = s.get("unique_id") or s.get("uniqueId") or s.get("id")
        if not unique_id:
            continue
        out.append({
            "submission_id": sub_id,
            "unique_id": str(unique_id),
            "sample_name": s.get("sample_name") or s.get("name"),
            "condition_name": s.get("condition_name") or s.get("condition"),
            "amt_to_inject": s.get("amt_2_inject") or s.get("amt_to_inject"),
            "internal_id": s.get("internal_id"),
            "internal_notes": s.get("internal_notes") or s.get("notes"),
        })
    return out


# ────────────────────────────────────────────────────────────────────────────
# PG UPSERT (same SQL as the xlsx version — schemas match)
# ────────────────────────────────────────────────────────────────────────────

UPSERT_SUBMISSION_SQL = """
INSERT INTO coreomics_submissions_cache (
    submission_id, internal_id, type, submitted_at, status, send_date,
    submitter_first_name, submitter_last_name, submitter_email,
    pi_first_name, pi_last_name, pi_email,
    institute, num_samples, organism, species, prot_or_pep, proteomics_type,
    mass_spec_wanted, sample_prep, gradient_length, dia, tmt,
    description, other_info,
    biohazard, pathogenic, nih_s10_user, is_nih_major_user, transgenic,
    po_account_number,
    raw_payload, source_export, source_export_md5
) VALUES (
    %(submission_id)s, %(internal_id)s, %(type)s, %(submitted_at)s,
    %(status)s, %(send_date)s,
    %(submitter_first_name)s, %(submitter_last_name)s, %(submitter_email)s,
    %(pi_first_name)s, %(pi_last_name)s, %(pi_email)s,
    %(institute)s, %(num_samples)s, %(organism)s, %(species)s,
    %(prot_or_pep)s, %(proteomics_type)s, %(mass_spec_wanted)s,
    %(sample_prep)s, %(gradient_length)s, %(dia)s, %(tmt)s,
    %(description)s, %(other_info)s,
    %(biohazard)s, %(pathogenic)s, %(nih_s10_user)s, %(is_nih_major_user)s,
    %(transgenic)s, %(po_account_number)s,
    %(raw_payload)s::jsonb, %(source_export)s, %(source_export_md5)s
)
ON CONFLICT (submission_id) DO UPDATE SET
    internal_id = EXCLUDED.internal_id,
    status = EXCLUDED.status,
    send_date = EXCLUDED.send_date,
    num_samples = EXCLUDED.num_samples,
    organism = EXCLUDED.organism,
    species = EXCLUDED.species,
    description = EXCLUDED.description,
    other_info = EXCLUDED.other_info,
    raw_payload = EXCLUDED.raw_payload,
    source_export = EXCLUDED.source_export,
    imported_at = NOW();
"""

UPSERT_SAMPLE_SQL = """
INSERT INTO coreomics_samples_cache (
    submission_id, unique_id, sample_name, condition_name,
    amt_to_inject, internal_id, internal_notes
) VALUES (
    %(submission_id)s, %(unique_id)s, %(sample_name)s, %(condition_name)s,
    %(amt_to_inject)s, %(internal_id)s, %(internal_notes)s
)
ON CONFLICT (submission_id, unique_id) DO UPDATE SET
    sample_name = EXCLUDED.sample_name,
    condition_name = EXCLUDED.condition_name,
    amt_to_inject = EXCLUDED.amt_to_inject,
    internal_id = EXCLUDED.internal_id,
    internal_notes = EXCLUDED.internal_notes,
    imported_at = NOW();
"""

AUTO_LINK_SQL = """
WITH unlinked AS (
    SELECT rf.raw_path, rf.raw_basename
    FROM raw_files rf
    LEFT JOIN delimp_raw_files_internal rfi ON rf.raw_path = rfi.raw_path
    WHERE rfi.coreomics_sample_unique_id IS NULL OR rfi IS NULL
),
candidates AS (
    SELECT u.raw_path, cs.submission_id, cs.unique_id, cs.sample_name,
           cs.condition_name,
           ROW_NUMBER() OVER (PARTITION BY u.raw_path
                              ORDER BY LENGTH(cs.sample_name) DESC) AS rank
    FROM unlinked u
    CROSS JOIN coreomics_samples_cache cs
    WHERE cs.sample_name IS NOT NULL
      AND regexp_replace(u.raw_basename, '[._-]', '', 'g')
          LIKE regexp_replace(cs.sample_name, '[._-]', '', 'g') || '%%'
)
INSERT INTO delimp_raw_files_internal (
    raw_path, coreomics_submission_id, coreomics_sample_unique_id,
    coreomics_sample_name, coreomics_condition_name, coreomics_last_synced_at
)
SELECT raw_path, submission_id, unique_id, sample_name, condition_name, NOW()
FROM candidates
WHERE rank = 1
ON CONFLICT (raw_path) DO UPDATE SET
    coreomics_submission_id = EXCLUDED.coreomics_submission_id,
    coreomics_sample_unique_id = EXCLUDED.coreomics_sample_unique_id,
    coreomics_sample_name = EXCLUDED.coreomics_sample_name,
    coreomics_condition_name = EXCLUDED.coreomics_condition_name,
    coreomics_last_synced_at = NOW()
WHERE delimp_raw_files_internal.coreomics_sample_unique_id IS NULL;
"""


# ────────────────────────────────────────────────────────────────────────────
# Driver
# ────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-url", default=os.environ.get("COREOMICS_BASE_URL", DEFAULT_BASE_URL))
    ap.add_argument("--lab", default=os.environ.get("COREOMICS_LAB", DEFAULT_LAB))
    ap.add_argument("--token-file", help="Override token source (default: COREOMICS_TOKEN env or ~/.coreomics_token)")
    ap.add_argument("--page-size", type=int, default=100)
    ap.add_argument("--max-pages", type=int, default=200)
    ap.add_argument("--max-submissions", type=int, default=None,
                    help="Stop after N submissions (for testing)")
    ap.add_argument("--skip-detail", action="store_true",
                    help="Don't fetch per-submission detail (skips samples!) — list only")

    ap.add_argument("--host", default="pgfarm.library.ucdavis.edu")
    ap.add_argument("--port", type=int, default=5432)
    ap.add_argument("--database", default="uc-davis-genome-center-proteomics-core/delimp")
    ap.add_argument("--user", default=os.environ.get("USER", "brettsp"))
    ap.add_argument("--auto-link-raws", action="store_true",
                    help="After import, run the auto-linkage pass on raw_files")
    ap.add_argument("--dry-run", action="store_true",
                    help="Pull from API + parse but don't write to PG")
    args = ap.parse_args()

    token = load_token(args.token_file)

    # ── Pull list ─────────────────────────────────────────────────────────
    print(f"[+] Listing submissions from {args.base_url} (lab={args.lab})...")
    subs = list_all_submissions(args.base_url, token, args.lab,
                                 page_size=args.page_size,
                                 max_pages=args.max_pages)
    print(f"    {len(subs)} submissions total")
    if args.max_submissions:
        subs = subs[:args.max_submissions]
        print(f"    (capping at --max-submissions {args.max_submissions} = {len(subs)})")

    # ── Pull detail per submission (for samples) ──────────────────────────
    sub_rows = []
    sample_rows = []
    source = f"api:{args.base_url}:{datetime.utcnow().isoformat()}Z"

    print("[+] Fetching detail for each submission (this is where samples live)...")
    for i, sub_summary in enumerate(subs, 1):
        sub_id = sub_summary.get("id")
        if not sub_id:
            continue
        if args.skip_detail:
            detail = sub_summary
        else:
            detail = fetch_submission_detail(args.base_url, token, sub_id) or sub_summary
        sub_rows.append(submission_to_pg_row(detail, source))
        sample_rows.extend(samples_for_submission(detail))
        if i % 50 == 0:
            print(f"  {i}/{len(subs)} detail fetched ({len(sample_rows)} samples so far)")

    print(f"    parsed: {len(sub_rows)} submissions, {len(sample_rows)} samples")

    if args.dry_run:
        print("\n[dry-run] would upsert above counts. Sample submission record:")
        for k, v in list(sub_rows[0].items())[:10]:
            print(f"    {k}: {str(v)[:80]}")
        print("\nSample sample record:")
        if sample_rows:
            for k, v in sample_rows[0].items():
                print(f"    {k}: {v}")
        return

    # ── PG upsert ─────────────────────────────────────────────────────────
    pgpwd = os.environ.get("PGPASSWORD")
    if not pgpwd:
        sys.exit("FATAL: PGPASSWORD env required for PG connection.")

    print(f"[+] Connecting to {args.host}:{args.port}/{args.database} as {args.user}")
    conn = psycopg2.connect(
        host=args.host, port=args.port, database=args.database,
        user=args.user, password=pgpwd, sslmode="require",
    )
    conn.autocommit = False
    try:
        with conn.cursor() as cur:
            print(f"[+] Upserting {len(sub_rows)} submissions...")
            psycopg2.extras.execute_batch(cur, UPSERT_SUBMISSION_SQL, sub_rows, page_size=200)
            print(f"[+] Upserting {len(sample_rows)} samples...")
            psycopg2.extras.execute_batch(cur, UPSERT_SAMPLE_SQL, sample_rows, page_size=500)
            if args.auto_link_raws:
                print("[+] Running auto-link against raw_files...")
                cur.execute(AUTO_LINK_SQL)
                print(f"    linked {cur.rowcount} raw files")
        conn.commit()
        print("[+] Done.")
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
