#!/usr/bin/env python3
"""Per-peptide LCA species attribution from a DIAMOND nr search (staxids col).

Reads NCBI taxonomy from the PRISTINE nodes.dmp.preDmnd (real ranks — the
diamond DB itself was built with ranks neutralized to "no rank", which only
affects diamond's labels, not the tree, so we recover real ranks here).

Usage: lca_attribute.py <blast_results.tsv> <out_prefix>
  blast tsv cols: qseqid sseqid pident len mm gap qs qe ss se evalue bitscore staxids
"""
import sys, collections

NODES = "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/ncbi_nr/nodes.dmp.preDmnd"
NAMES = "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/ncbi_nr/names.dmp"
HITS, OUTP = sys.argv[1], sys.argv[2]

# anchor taxa for kingdom/domain classification
METAZOA, BACTERIA, ARCHAEA, VIRUSES, FUNGI, VIRIDIPLANTAE = 33208, 2, 2157, 10239, 4751, 33090
DIAG_RANKS = {"species", "subspecies", "species group", "species subgroup",
              "genus", "subgenus", "strain", "isolate", "varietas", "forma"}

sys.stderr.write("[lca] loading taxonomy...\n")
parent, rank = {}, {}
for line in open(NODES):
    p = line.split("\t|\t")
    t = int(p[0]); parent[t] = int(p[1]); rank[t] = p[2]
name = {}
for line in open(NAMES):
    p = line.split("\t|\t")
    if len(p) >= 4 and p[3].rstrip("\t|\n") == "scientific name":
        name[int(p[0])] = p[1]
sys.stderr.write(f"[lca] {len(parent)} nodes, {len(name)} names\n")

def lineage(t):
    out, seen = [], set()
    while t and t not in seen:
        seen.add(t); out.append(t)
        pt = parent.get(t)
        if pt is None or pt == t:
            break
        t = pt
    return out  # taxon -> ... -> root

def lca(taxids):
    taxids = [t for t in set(taxids) if t in parent]
    if not taxids:
        return None
    lins = [lineage(t) for t in taxids]
    common = set.intersection(*[set(l) for l in lins])
    for t in lins[0]:          # deepest (lowest) shared node
        if t in common:
            return t
    return 1

# group hits by peptide
hits = collections.defaultdict(list)   # pep -> [(bitscore, [taxids], pident)]
for line in open(HITS):
    f = line.rstrip("\n").split("\t")
    if len(f) < 13:
        continue
    tids = [int(x) for x in f[12].replace(";", " ").split() if x.isdigit() and int(x) > 0]
    hits[f[0]].append((float(f[11]), tids, float(f[2])))

cat = collections.Counter(); host_sp = collections.Counter(); micro = collections.Counter()
n_diag = 0
with open(OUTP + "_peptide_lca.tsv", "w") as o:
    o.write("peptide\tn_hits\ttop_pident\tlca_taxid\tlca_name\tlca_rank\tcategory\tdiagnostic\n")
    for pep, hl in hits.items():
        topbs = max(h[0] for h in hl)
        keep = [h for h in hl if h[0] >= 0.9 * topbs]   # MEGAN-style top-10% window
        l = lca([t for h in keep for t in h[1]])
        if l is None:
            continue
        lin = set(lineage(l))
        c = ("host" if METAZOA in lin else
             "microbiome" if (BACTERIA in lin or ARCHAEA in lin or VIRUSES in lin) else
             "plant/fungal" if (FUNGI in lin or VIRIDIPLANTAE in lin) else "other/conserved")
        rk = rank.get(l, "no rank"); nm = name.get(l, str(l))
        diag = rk in DIAG_RANKS
        n_diag += diag
        o.write(f"{pep}\t{len(hl)}\t{max(h[2] for h in keep):.1f}\t{l}\t{nm}\t{rk}\t{c}\t{int(diag)}\n")
        cat[c] += 1
        if c == "host" and diag:
            host_sp[f"{nm} ({rk})"] += 1
        if c == "microbiome":
            micro[nm] += 1

print(f"=== peptides with LCA: {sum(cat.values())}  (diagnostic species/genus: {n_diag}) ===")
print("--- category ---")
for c, n in cat.most_common():
    print(f"  {n:>6}  {c}")
print("--- top HOST species (diagnostic, species/genus LCA) ---")
for s, n in host_sp.most_common(20):
    print(f"  {n:>6}  {s}")
print("--- top microbiome taxa ---")
for s, n in micro.most_common(12):
    print(f"  {n:>6}  {s}")
