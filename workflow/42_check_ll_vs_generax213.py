#!/usr/bin/env python3
"""
42_check_ll_vs_generax213.py  --  evaluate workflow/35's reconciliation log-likelihood at the D/L rates
GeneRax 2.1.3 estimated (results/rerun_2026-10/udl_validation/generax213_optrates_reference.tsv) and
compare with the ll GeneRax printed; also name species nodes by NCBI clade.
"""
import csv, math, os, sys
from importlib import import_module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
U = import_module("35_undated_dl")
RC = "results/reconciliation"
SP = {"ncbi": f"{RC}/aa3di/eval/species_trees/inferred_species_tree.newick",
      "graftv1": f"{RC}/aa3di_grafted/run/species_trees/starting_species_tree.newick",
      "graftv2": "results/rerun_2026-10/grafted/grafted_species_v2.binary.nwk"}
GT = {"aa3di50": "results/seq_tree/pif1_aa3di.treefile", "supported": "results/seq_tree/pif1_aa3di_supported.treefile"}
g2s = {l.split()[0]: l.split()[1] for l in open("data/species_tree/gene_species.map") if l.strip()}
out = []
for r in csv.DictReader(open("results/rerun_2026-10/udl_validation/generax213_optrates_reference.tsv"), delimiter="\t"):
    sp = U.Species(SP[r["species_tree"]])
    gt = U.GeneTree(GT[r["gene_tree"]], keep=[g for g in g2s if g2s[g] in sp.leaves])
    ll = U.UndatedDL(sp, gt, g2s, float(r["gx213_D"]), float(r["gx213_L"])).sum_loglik()
    r["udl_ll_at_gx213_rates"] = f"{ll:.3f}"
    r["diff"] = f"{ll - float(r['gx213_ll']):.3f}"
    out.append(r)
    print("\t".join(r.values()), flush=True)
with open("results/rerun_2026-10/udl_validation/ll_check_vs_generax213.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
