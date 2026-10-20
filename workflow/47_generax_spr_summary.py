#!/usr/bin/env python3
"""
47_generax_spr_summary.py  --  summarise GeneRax --strategy SPR runs started from different gene trees.

Per run (directory given with --run LABEL=DIR):
  - the joint / reconciliation / libpll log-likelihood after every search round (from DIR.stdout) and
    the D/L rates GeneRax used;
  - the final gene tree DIR/results/pif1/geneTree.newick: is the 197-gene split present, and is the
    hPIF1-rooted ScPif1/ScRrm3 MRCA all-Saccharomycotina (workflow/43);
  - where the anchor duplication maps (workflow/35, the validated UndatedDL re-implementation, run on the
    final tree with ML rates on the same species tree; written to DIR/recon35/).

Usage (repo root, pif1 env):
  python workflow/47_generax_spr_summary.py --species-tree S.nwk --map M --out OUT.tsv \
      --run deep=results/.../generax/v2fix_aa_deep --run 197=results/.../generax/v2fix_aa_197
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys

CLADE = "results/rerun_2026-10/supports/clade197_supported.txt"


def rounds(stdout):
    out, rates = [], None
    for line in open(stdout):
        m = re.search(r"D=([\d.e+-]+), L=([\d.e+-]+), RecLL= ([-\d.e+]+)", line)
        if m:
            rates = (float(m.group(1)), float(m.group(2)))
        m = re.search(r"JointLL=([-\d.e+]+) RecLL=([-\d.e+]+) LibpllLL=([-\d.e+]+)", line)
        if m:
            out.append((rates, float(m.group(1)), float(m.group(2)), float(m.group(3))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, help="LABEL=DIR")
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--map", default="results/rerun_2026-10/species_tree_aware/species_trees/gene_species_957.map")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rows = []
    for spec in args.run:
        label, d = spec.split("=", 1)
        r = rounds(d + ".stdout")
        tree = os.path.join(d, "results", "pif1", "geneTree.newick")
        row = {"run": label, "dir": d, "rounds": len(r)}
        if r:
            row.update({"start_joint": r[0][1] if len(r) else None})
        # the initial LL is printed by the worker before the first round
        init = None
        log = os.path.join(d, "gene_optimization_0", "per_job_logs", "pif1_out.txt")
        if os.path.exists(log):
            m = re.search(r"joint: ([-\d.e+]+)\s+libpll: ([-\d.e+]+)\s+reconciliation: ([-\d.e+]+)", open(log).read())
            if m:
                init = tuple(float(x) for x in m.groups())
        row.update({"init_joint": init[0] if init else None, "init_libpll": init[1] if init else None,
                    "init_rec": init[2] if init else None})
        if r:
            rates, j, rec, lib = r[-1]
            row.update({"final_joint": j, "final_rec": rec, "final_libpll": lib,
                        "final_D": rates[0] if rates else None, "final_L": rates[1] if rates else None})
        print(f"{label}: init {init}; per round (D, L, joint, rec, libpll):")
        for rr in r:
            print(f"   {rr[0]}  {rr[1]:.2f}  {rr[2]:.2f}  {rr[3]:.2f}")
        if os.path.exists(tree):
            s = subprocess.run([sys.executable, "workflow/43_count_split_in_trees.py", tree, CLADE],
                               capture_output=True, text=True).stdout.strip()
            print("   final tree:", s)
            m = re.search(r"present as a split: (\d+); anchor MRCA all-Saccharomycotina: (\d+); MRCA sizes \[\((\d+)", s)
            row.update({"split197": int(m.group(1)), "anchor_mrca_all_sacch": int(m.group(2)),
                        "anchor_mrca_size": int(m.group(3))})
            od = os.path.join(d, "recon35")
            subprocess.run([sys.executable, "workflow/35_undated_dl.py", "--gene-tree", tree, "--species-tree",
                            args.species_tree, "--map", args.map, "--optimize-rates", "--quiet", "--outdir", od],
                           check=True, capture_output=True)
            j = json.load(open(os.path.join(od, "summary.json")))
            a = j["anchor_duplication"]
            row.update({"recon35_D": j["rates"][0], "recon35_L": j["rates"][1], "recon35_ll_sum": j["sum_loglik"],
                        "anchor_event": a["event"], "anchor_node": a["species_node"],
                        "anchor_node_leaves": a["species_node_leaves"], "anchor_gene_clade": a["gene_clade_size"],
                        "D_total": j["totals"].get("D"), "SL_total": j["totals"].get("SL")})
            print(f"   workflow/35 on final tree: rates {j['rates']}, rec ll (sum over roots) {j['sum_loglik']:.2f}; "
                  f"anchor {a['event']} -> {a['species_node']} ({a['species_node_leaves']} leaves), gene clade {a['gene_clade_size']}")
        rows.append(row)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(args.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
