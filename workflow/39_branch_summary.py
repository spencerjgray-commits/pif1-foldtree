#!/usr/bin/env python3
"""
39_branch_summary.py  --  branch-length numbers quoted in Results (review RES-R8 / RES-R15) from an
AA branch-length tree: the two daughter stem branches below the ScPif1/ScRrm3 duplication node, and
mean terminal branches for Saccharomycotina Pif1 / Rrm3 copies vs single- and multi-copy genes
outside the duplication clade. Topology queries use tip sets; branch lengths only (no support labels).

Usage: python workflow/39_branch_summary.py TREE
"""
import csv, sys
from collections import Counter
import dendropy

tree = dendropy.Tree.get(path=sys.argv[1], schema="newick", preserve_underscores=True)
acc = {r["tip_label"]: r["accession"] for r in csv.DictReader(open("data/seqs/tip_map.tsv"), delimiter="\t")}
taxid = {r["tip_label"]: r["taxid"] for r in csv.DictReader(open("data/seqs/tip_map.tsv"), delimiter="\t")}
sub = {r["accession"]: r["subphylum"] for r in csv.DictReader(open("manifest.csv"))}
copies = Counter(taxid.values())
lab = {acc[t]: t for t in acc}
term_len = {l.taxon.label: l.edge.length for l in tree.leaf_node_iter()}   # before rerooting
hum = tree.find_node_with_taxon_label(lab["Q9H611"])
tree.reroot_at_edge(hum.edge, update_bipartitions=True)
dup = tree.mrca(taxon_labels=[lab["P07271"], lab["P38766"]])
kids = dup.child_nodes()
inside = {l.taxon.label for l in dup.leaf_iter()}
for k in kids:
    tips = {l.taxon.label for l in k.leaf_iter()}
    which = "Pif1" if lab["P07271"] in tips else "Rrm3"
    term = [term_len[l.taxon.label] for l in k.leaf_iter()]
    print(f"{which} daughter: {len(tips)} tips, stem branch {k.edge.length:.3f}, mean terminal {sum(term)/len(term):.3f}")
out_single, out_multi = [], []
for l in tree.leaf_node_iter():
    t = l.taxon.label
    if t in inside or t == lab["Q9H611"]:
        continue
    if sub.get(acc[t]) == "Saccharomycotina":
        continue
    (out_single if copies[taxid[t]] == 1 else out_multi).append(term_len[t])
allout = [term_len[t] for t in term_len if t not in inside]
m = lambda v: sum(v) / len(v)
print(f"outside duplication clade, all tips (incl. 8 Sacch + hPIF1): mean terminal {m(allout):.3f} (n={len(allout)})")
print(f"non-Saccharomycotina single-copy: {m(out_single):.3f} (n={len(out_single)}); multi-copy: {m(out_multi):.3f} (n={len(out_multi)})")
