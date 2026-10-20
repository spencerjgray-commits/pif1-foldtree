#!/usr/bin/env python3
"""
43_count_split_in_trees.py  --  in a multi-tree file (e.g. a .ufboot), count how many trees contain a
given tip set as a split, and how often the outgroup-rooted ScPif1/ScRrm3 MRCA is all-Saccharomycotina.

Usage: python workflow/43_count_split_in_trees.py TREES CLADE_TIPS.txt [--root Q9H611]
"""
import argparse, csv, os, sys, tempfile
from collections import Counter
from importlib import import_module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SS = import_module("33_split_supports")
ap = argparse.ArgumentParser()
ap.add_argument("trees"); ap.add_argument("clade"); ap.add_argument("--root", default="Q9H611")
a = ap.parse_args()
clade = frozenset(l.strip() for l in open(a.clade) if l.strip())
sac = {r["accession"] for r in csv.DictReader(open("manifest.csv")) if r["subphylum"] == "Saccharomycotina"}
lines = [l for l in open(a.trees) if l.strip()]
tmp = tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False).name
n_split = n_sac = 0
sizes = Counter()
for l in lines:
    open(tmp, "w").write(l)
    st = SS.SplitTree(tmp, a.root)
    n_split += st.support(clade)[0]
    m = st.mrca(["P07271", "P38766"])
    sizes[len(m)] += 1
    n_sac += all(t.split("_")[0] in sac for t in m)
os.remove(tmp)
print(f"trees: {len(lines)}; clade ({len(clade)} tips) present as a split: {n_split}; "
      f"anchor MRCA all-Saccharomycotina: {n_sac}; MRCA sizes {sizes.most_common(6)}")
