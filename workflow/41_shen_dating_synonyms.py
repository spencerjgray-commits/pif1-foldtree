#!/usr/bin/env python3
"""
41_shen_dating_synonyms.py  --  re-date the duplication clade on Shen et al. 2018 using the NCBI-current
name matching of workflow/37 (results/rerun_2026-10/grafted/shen2018_synonyms.tsv) instead of the
organism-string matching of workflow/17, and report the dated node's composition (which Shen tips it
excludes), crown and stem ages.

Usage: python workflow/41_shen_dating_synonyms.py --gene-tree results/seq_tree/pif1_aa3di_supported.treefile
"""
import argparse, csv, os, sys
from importlib import import_module
import dendropy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SS = import_module("33_split_supports")

ap = argparse.ArgumentParser()
ap.add_argument("--gene-tree", required=True)
ap.add_argument("--syn", default="results/rerun_2026-10/grafted/shen2018_synonyms.tsv")
ap.add_argument("--shen", default="data/species_tree/shen2018_timetree.newick")
a = ap.parse_args()
st = SS.SplitTree(a.gene_tree, "Q9H611")
clade = st.mrca(["P07271", "P38766"])
tax = {r["tip_label"]: r["taxid"] for r in csv.DictReader(open("data/seqs/tip_map.tsv"), delimiter="\t")}
syn = {r["taxid"]: r for r in csv.DictReader(open(a.syn), delimiter="\t")}
taxa = {tax[t] for t in clade}
exact = sorted({syn[t]["shen_tip_or_clade"] for t in taxa if t in syn and syn[t]["placement"] == "exact"})
print(f"duplication clade: {len(clade)} genes, {len(taxa)} taxids; exact Shen matches: {len(exact)} tips")
tt = dendropy.Tree.get(path=a.shen, schema="newick", preserve_underscores=True, rooting="force-rooted")
tt.calc_node_ages(ultrametricity_precision=0.01)
m = tt.mrca(taxon_labels=exact)
under = {l.taxon.label for l in m.leaf_iter()}
allt = {l.taxon.label for l in tt.leaf_node_iter()}
print(f"MRCA subtends {len(under)} of {len(allt)} Shen tips; excluded: {sorted(allt - under)}")
print(f"crown {m.age*100:.1f} Mya; stem {m.parent_node.age*100:.1f} Mya; Shen root {tt.seed_node.age*100:.1f} Mya")
