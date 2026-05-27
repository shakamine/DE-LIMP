#!/usr/bin/env python3
"""
rewrite_transdecoder_headers.py — v1.0

Rewrites TransDecoder peptide FASTA headers into DE-LIMP proteogenomics
format, classifying each ORF by gffcompare class code:

  REF | NOVEL_GENE | NOVEL_ISOFORM | UNPARSED

Output header format:
  >sp|<protein_id>|<gene_symbol>_<TAG> source=<class> ORF_type=<type> \\
    strand=<+|-> len=<aa> coords=<...> parent_gene=<id> transcript=<id>

Requires the proteog_helpers conda env (gffutils, biopython).

USAGE:
  python rewrite_transdecoder_headers.py \\
    --transdecoder TRANSDECODER_PEP \\
    --merged-gtf MERGED_GTF \\
    --gffcompare-tmap GFFCMP_TMAP \\
    --project-tag TAG \\
    --output OUTPUT_FASTA

PROJECT-CONVENTION NOTE (maintainers, please read):

  This script uses a locked CLASS_CODE_MAP. When you encounter a gffcompare
  class code in real data that isn't in the map, the rewriter halts via
  exit code 1 (Rule 4: silent classification failure is banned).

  The right response is NOT to disable the check or to work around the lock.
  The right response is to EXTEND CLASS_CODE_MAP with the new code and a
  one-line comment justifying the bucket choice — then commit. The map
  exists to make extension deliberate, not to prevent it.

  This same pattern applies to the STAR threshold tiers in helpers_rnaseq.R
  and any future locked constants in this pipeline: locked != immutable;
  locked means "expansion requires documented reasoning."

CLAUDE.md Rule 4 alignment: exits non-zero if any UNPARSED entries are
produced. Silent classification failure would be a data-quality bug that
could poison downstream DE-LIMP reports.
"""

import argparse
import re
import sys
from collections import Counter

import gffutils

# -----------------------------------------------------------------------------
# Locked gffcompare class-code → source-class mapping (v1.0).
# See header docstring for the "locked != immutable" project convention.
# -----------------------------------------------------------------------------
CLASS_CODE_MAP = {
    # === Reference-matching ============================================
    "=": "REF",            # exact intron-chain match
    "c": "REF",            # query contained in reference (shorter)

    # === Novel isoforms (alternative splicing of known genes) ==========
    "j": "NOVEL_ISOFORM",  # multi-exon novel splice junction
    "e": "NOVEL_ISOFORM",  # single-exon novel splice junction
    "k": "NOVEL_ISOFORM",  # query CONTAINS reference (UTR extension,
                           #   boundary novel exon — opposite of "c")
    "m": "NOVEL_ISOFORM",  # retained intron, full chain match
    "n": "NOVEL_ISOFORM",  # retained intron, partial chain
    "y": "NOVEL_ISOFORM",  # edge case: contains reference within intron;
                           #   rare; same-strand overlap so treated as
                           #   isoform rather than intergenic

    # === Novel genes (intergenic or intronic) ==========================
    "u": "NOVEL_GENE",     # fully intergenic
    "i": "NOVEL_GENE",     # intronic — uORF/sORF candidates worth surfacing

    # === Low-discovery-value overlaps (reference-related) ==============
    "o": "REF",            # generic exonic overlap, same/other strand
    "x": "REF",            # exonic overlap, opposite strand
    "s": "REF",            # intron match, opposite strand (mapping error usually)
    "p": "REF",            # possible polymerase run-on (no actual overlap)
    "r": "REF",            # repeat-containing

    # === Unmapped (default) ============================================
    # Codes ".", "?", or absent → UNPARSED via dict.get(code, "UNPARSED")
    # Halts via Rule 4. See "PROJECT-CONVENTION NOTE" in the docstring.
}


# TransDecoder .pep header format (validated 2026-05-20):
#   >ENSMUST00000000001.5.p1 GENE.ENSMUST00000000001.5~~ENSMUST00000000001.5.p1  ORF type:complete (+),score=62.30 len:354 ENSMUST00000000001.5:142-1206(+)
TD_HEADER_RE = re.compile(
    r"^>(?P<protein_id>\S+)\s+"
    r"GENE\.(?P<td_gene>\S+?)~~(?P<td_protein>\S+?)\s+"
    r"ORF\s+type:(?P<orf_type>\S+)\s+\((?P<strand>[+\-.])\),"
    r"score=(?P<score>\S+)\s+"
    r"len:(?P<len>\d+)\s+"
    r"(?P<coords>\S+)"
)


def parse_td_header(line):
    """Parse a TransDecoder pep header line. Returns dict or None if unparseable."""
    m = TD_HEADER_RE.match(line)
    if not m:
        return None
    d = m.groupdict()
    # Strip the .pY suffix to recover parent transcript_id.
    d["transcript_id"] = re.sub(r"\.p\d+$", "", d["protein_id"])
    return d


def parse_gffcompare_tmap(tmap_path):
    """
    Parse a gffcompare .tmap file. Standard column order (gffcompare v0.12.6):
      ref_gene_id  ref_id  class_code  qry_gene_id  qry_id  num_exons \\
      FPKM  TPM  cov  len  major_iso_id  ref_match_len

    Returns:
      dict mapping qry_id (StringTie transcript_id) → class_code (single char).
    """
    mapping = {}
    with open(tmap_path) as fh:
        header_line = fh.readline().rstrip("\n")
        if not header_line:
            sys.exit(f"FATAL: empty tmap file: {tmap_path}")
        header = header_line.split("\t")
        try:
            ci_qry_id = header.index("qry_id")
            ci_class_code = header.index("class_code")
        except ValueError:
            sys.exit(
                f"FATAL: tmap header lacks qry_id or class_code column.\n"
                f"  Got columns: {header}"
            )
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(ci_qry_id, ci_class_code):
                continue
            mapping[cols[ci_qry_id]] = cols[ci_class_code]
    return mapping


def preflight_unmapped_codes(class_map):
    """
    Scan the tmap class-code dict for codes not in CLASS_CODE_MAP and emit a
    warning BEFORE classification begins.

    Returns a Counter of unmapped codes (for inclusion in the diagnostic).
    """
    counts = Counter(class_map.values())
    # "." is the documented "unclassified" sentinel — expected, not unknown.
    unmapped = {c: n for c, n in counts.items()
                if c not in CLASS_CODE_MAP and c != "."}
    if unmapped:
        total = sum(unmapped.values())
        sys.stderr.write(
            f"\n⚠ WARNING: Found {total:,} entries with unmapped class codes in tmap.\n"
            f"These will be classified as UNPARSED and will cause exit code 1\n"
            f"if any of them have downstream TransDecoder ORFs.\n\n"
            f"Unmapped codes:\n"
        )
        for code, ct in sorted(unmapped.items(), key=lambda x: -x[1]):
            sys.stderr.write(f"  '{code}': {ct:,} entries\n")
        sys.stderr.write(
            "\nTo resolve: extend CLASS_CODE_MAP in this script with documented\n"
            "reasoning, or confirm these are intentionally excluded.\n"
            "See the PROJECT-CONVENTION NOTE in the script docstring.\n\n"
        )
        sys.stderr.flush()
    return Counter(unmapped)


def build_gtf_db(merged_gtf):
    """Build an in-memory gffutils SQLite DB from the StringTie merged GTF."""
    sys.stderr.write(f"[1/4] Building gffutils DB from {merged_gtf}...\n")
    sys.stderr.flush()
    db = gffutils.create_db(
        merged_gtf,
        dbfn=":memory:",
        force=True,
        keep_order=False,
        merge_strategy="merge",
        sort_attribute_values=False,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
    n = sum(1 for _ in db.features_of_type("transcript"))
    sys.stderr.write(f"  {n:,} transcript records indexed\n")
    sys.stderr.flush()
    return db


def _clean_attr(val):
    """
    Strip stray GTF delimiters from an attribute value.

    gffutils' attribute parser has a known quirk: on GTF lines that don't end
    with trailing whitespace before the newline, the closing `"` and `;` of
    the final attribute get kept in the value. This affects every HAVANA-
    source line in a gencode-derived merged.gtf — about half the entries in
    our validation data. We defensively strip the delimiters so downstream
    consumers never see artifacts regardless of which gffutils version is
    installed.
    """
    if val is None:
        return None
    return val.strip(' "\';\t\r\n') or None


def get_transcript_attrs(db, transcript_id):
    """
    Look up (gene_id, ref_gene_id, gene_name) for a transcript.
    Returns (None, None, None) if not found. All values defensively cleaned
    of stray GTF delimiters (see _clean_attr docstring).
    """
    try:
        t = db[transcript_id]
    except (KeyError, gffutils.exceptions.FeatureNotFoundError):
        return (None, None, None)
    gene_id     = _clean_attr((t.attributes.get("gene_id") or [None])[0])
    ref_gene_id = _clean_attr((t.attributes.get("ref_gene_id") or [None])[0])
    gene_name   = _clean_attr((t.attributes.get("gene_name") or [None])[0])
    return (gene_id, ref_gene_id, gene_name)


def pick_symbol(source, gene_id, ref_gene_id, gene_name):
    """Pick the gene symbol that goes into the sp|ID|SYMBOL_TAG accession."""
    if source in ("REF", "NOVEL_ISOFORM"):
        return gene_name or ref_gene_id or gene_id or "Unknown"
    # NOVEL_GENE — use gene_id (typically MSTRG.X)
    return gene_id or "Unknown"


def rewrite_one(td, class_code, db, project_tag, unknown_codes):
    """Build the new header for one TransDecoder record."""
    tid = td["transcript_id"]
    source = CLASS_CODE_MAP.get(class_code)
    if source is None:
        if class_code:
            unknown_codes[class_code] += 1
        source = "UNPARSED"

    gene_id, ref_gene_id, gene_name = get_transcript_attrs(db, tid)
    sym = pick_symbol(source, gene_id, ref_gene_id, gene_name)
    parent_gene = ref_gene_id or gene_id or "Unknown"

    desc = (
        f"source={source} "
        f"ORF_type={td['orf_type']} "
        f"strand={td['strand']} "
        f"len={td['len']} "
        f"coords={td['coords']} "
        f"parent_gene={parent_gene} "
        f"transcript={tid}"
    )
    return f">sp|{td['protein_id']}|{sym}_{project_tag} {desc}"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--transdecoder", required=True,
                    help="TransDecoder .pep FASTA")
    ap.add_argument("--merged-gtf", required=True,
                    help="StringTie merged.gtf")
    ap.add_argument("--gffcompare-tmap", required=True,
                    help="gffcompare .tmap (e.g. gffcmp.merged.gtf.tmap)")
    ap.add_argument("--project-tag", required=True,
                    help="Project tag glued to symbol with _ (e.g. MM39TEST)")
    ap.add_argument("--output", required=True,
                    help="Output FASTA")
    args = ap.parse_args()

    if not re.match(r"^[A-Za-z0-9_-]+$", args.project_tag):
        sys.exit(
            f"FATAL: project-tag must match [A-Za-z0-9_-]+; got: {args.project_tag}"
        )

    db = build_gtf_db(args.merged_gtf)

    sys.stderr.write(f"[2/4] Loading gffcompare tmap from {args.gffcompare_tmap}...\n")
    sys.stderr.flush()
    class_map = parse_gffcompare_tmap(args.gffcompare_tmap)
    sys.stderr.write(f"  {len(class_map):,} qry_id → class_code mappings\n")
    sys.stderr.flush()

    # Pre-flight: warn about any class codes outside the locked map BEFORE
    # we start classifying. Makes the next "unknown code surfaces in real
    # data" event self-explanatory instead of looking like a mysterious halt.
    preflight_unmapped_codes(class_map)

    sys.stderr.write(
        f"[3/4] Rewriting headers: {args.transdecoder} → {args.output}\n"
    )
    sys.stderr.flush()
    n_in = 0
    n_parsed = 0
    n_unparsed_td = 0
    source_counts = Counter()
    unknown_codes = Counter()
    missing_from_tmap = 0

    with open(args.transdecoder) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                td = parse_td_header(line)
                if td is None:
                    fout.write(
                        line.rstrip("\n")
                        + " source=UNPARSED reason=td_header_parse_failed\n"
                    )
                    n_unparsed_td += 1
                    source_counts["UNPARSED"] += 1
                    continue
                class_code = class_map.get(td["transcript_id"])
                if class_code is None:
                    missing_from_tmap += 1
                    class_code = "."  # → UNPARSED via the lock-in map
                new_header = rewrite_one(td, class_code, db, args.project_tag,
                                         unknown_codes)
                source = CLASS_CODE_MAP.get(class_code, "UNPARSED")
                source_counts[source] += 1
                fout.write(new_header + "\n")
                n_parsed += 1
            else:
                fout.write(line)

    sys.stderr.write("[4/4] Done.\n")
    sys.stderr.write(f"  Headers in:                {n_in:,}\n")
    sys.stderr.write(f"  Parsed cleanly:            {n_parsed:,}\n")
    if n_unparsed_td:
        sys.stderr.write(
            f"  TransDecoder-header-unparseable: {n_unparsed_td:,}\n"
        )
    if missing_from_tmap:
        sys.stderr.write(
            f"  Transcripts missing from tmap (classified UNPARSED): "
            f"{missing_from_tmap:,}\n"
        )
    sys.stderr.write("  Source classification:\n")
    for src in ("REF", "NOVEL_GENE", "NOVEL_ISOFORM",
                "VARIANT", "UNIPROT", "UNPARSED"):
        sys.stderr.write(f"    {src:14s} {source_counts[src]:>10,}\n")
    if unknown_codes:
        sys.stderr.write(
            "  Unknown class codes encountered during rewrite "
            "(escaped pre-flight check):\n"
        )
        for code, ct in sorted(unknown_codes.items(), key=lambda x: -x[1]):
            sys.stderr.write(f"    '{code}': {ct:,}\n")

    # CLAUDE.md Rule 4 — silent classification failure is banned.
    if source_counts["UNPARSED"] > 0:
        sys.stderr.write(
            f"\nFATAL: {source_counts['UNPARSED']:,} UNPARSED entries — "
            "pipeline halted.\n"
            "Likely causes:\n"
            "  1. gffcompare tmap is missing transcripts present in TransDecoder pep\n"
            "  2. TransDecoder header format has drifted (file an issue)\n"
            "  3. gffcompare produced class codes outside the locked map\n"
            "     (see pre-flight WARNING above for the specific codes)\n"
        )
        sys.exit(1)

    sys.stderr.write(f"\nOutput: {args.output}\n")


if __name__ == "__main__":
    main()
