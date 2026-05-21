#!/usr/bin/env python3
"""Rewrite TransDecoder .pep FASTA headers into DE-LIMP proteogenomics format.

Input:
  --gtf   StringTie merged GTF (transcript ↔ gene mapping)
  --pep   TransDecoder peptide FASTA
  --out   output FASTA path
  --project-tag  tag appended to gene symbol (default MM39TEST)

Output header format (UniProt-like, DIA-NN parser-friendly):
  >sp|<protein_id>|<gene_symbol>_<TAG> source=<REF|NOVEL_ISOFORM|NOVEL_GENE> \\
      ORF_type=<...> len=<aa> coords=<...> parent_gene=<id> transcript=<id>

The leading sp|ACC|NAME format lets DE-LIMP's existing gene-symbol parser
recover the symbol; downstream tools that don't recognise the new TAG
still see a parseable accession.
"""

import argparse
import re
import sys
from collections import Counter

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')
TD_HEADER_RE = re.compile(
    r"^>(?P<protein_id>\S+)\s+"
    r"GENE\.(?P<td_gene>\S+?)~~(?P<td_protein>\S+?)\s+"
    r"ORF\s+type:(?P<orf_type>\S+)\s+\((?P<strand>[+\-.])\),"
    r"score=(?P<score>\S+)\s+"
    r"len:(?P<len>\d+)\s+"
    r"(?P<coords>\S+)"
)


def parse_gtf_mapping(gtf_path):
    """Return dict: transcript_id -> {gene_id, ref_gene_id, gene_name}"""
    mapping = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            mapping[tid] = {
                "gene_id": attrs.get("gene_id", ""),
                "ref_gene_id": attrs.get("ref_gene_id", ""),
                "gene_name": attrs.get("gene_name", ""),
            }
    return mapping


def parse_td_header(line):
    m = TD_HEADER_RE.match(line)
    if not m:
        return None
    d = m.groupdict()
    d["transcript_id"] = re.sub(r"\.p\d+$", "", d["protein_id"])
    return d


def classify(td, info):
    """Classify ORF source as REF | NOVEL_ISOFORM | NOVEL_GENE."""
    tid = td["transcript_id"]
    gene_id = info.get("gene_id", "") if info else ""
    ref_gene_id = info.get("ref_gene_id", "") if info else ""

    is_mstrg_transcript = tid.startswith("MSTRG.")
    is_mstrg_gene = gene_id.startswith("MSTRG.") and not ref_gene_id

    if is_mstrg_transcript and is_mstrg_gene:
        return "NOVEL_GENE"
    if is_mstrg_transcript and ref_gene_id:
        return "NOVEL_ISOFORM"
    return "REF"


def pick_symbol(td, info, source):
    info = info or {}
    if source == "REF":
        return info.get("gene_name") or info.get("ref_gene_id") or info.get("gene_id") or "Unknown"
    if source == "NOVEL_ISOFORM":
        return info.get("gene_name") or info.get("ref_gene_id") or info.get("gene_id") or "Unknown"
    # NOVEL_GENE
    return info.get("gene_id") or "Unknown"


def rewrite_header(td, gtf_map, project_tag):
    tid = td["transcript_id"]
    info = gtf_map.get(tid)
    source = classify(td, info)
    sym = pick_symbol(td, info, source)
    parent_gene = (info or {}).get("ref_gene_id") or (info or {}).get("gene_id") or "Unknown"

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


def rewrite_fasta(pep_in, fasta_out, gtf_map, project_tag):
    n_in = 0
    n_parsed = 0
    n_unparsed = 0
    source_counts = Counter()

    with open(pep_in) as fin, open(fasta_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                td = parse_td_header(line)
                if td is None:
                    fout.write(line.rstrip() + f" [unparsed]_{project_tag}\n")
                    n_unparsed += 1
                    source_counts["UNPARSED"] += 1
                    continue
                info = gtf_map.get(td["transcript_id"])
                source = classify(td, info)
                source_counts[source] += 1
                fout.write(rewrite_header(td, gtf_map, project_tag) + "\n")
                n_parsed += 1
            else:
                fout.write(line)
    return n_in, n_parsed, n_unparsed, source_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="StringTie merged GTF")
    ap.add_argument("--pep", required=True, help="TransDecoder peptide FASTA")
    ap.add_argument("--out", required=True, help="Output FASTA path")
    ap.add_argument("--project-tag", default="MM39TEST")
    args = ap.parse_args()

    print(f"[1/3] Parsing GTF: {args.gtf}", file=sys.stderr)
    gtf_map = parse_gtf_mapping(args.gtf)
    print(f"  Transcript→gene mappings: {len(gtf_map):,}", file=sys.stderr)

    print(f"[2/3] Rewriting headers: {args.pep} → {args.out}", file=sys.stderr)
    n_in, n_parsed, n_unparsed, source_counts = rewrite_fasta(
        args.pep, args.out, gtf_map, args.project_tag
    )
    print(f"  Headers in: {n_in:,}", file=sys.stderr)
    print(f"  Parsed: {n_parsed:,}", file=sys.stderr)
    print(f"  Unparsed: {n_unparsed:,}", file=sys.stderr)

    print(f"[3/3] Source classification:", file=sys.stderr)
    for src, cnt in sorted(source_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {src}: {cnt:,}", file=sys.stderr)

    print(f"Done. Output: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
