#!/usr/bin/env python3
"""
50_local_swap_summary.py  --  joint (sequence + reconciliation) comparison of the local-swap trees.

For each tree of a workflow/49 set:
  sequence lnL   from the IQ-TREE -z run (branch lengths and model re-optimised per tree; AU p-value)
  rec lnL        workflow/35 at ML D/L rates, summed over gene-tree roots, on each species tree
                 (results in <dir>/recon/<tree>_<species tree>/summary.json)
  joint          sequence lnL + rec lnL
The pairs that matter are a background tree and its swapped version (same rest-of-tree): the joint
difference within a pair is the species-tree-aware preference for one placement over the other.

Usage (repo root): python workflow/50_local_swap_summary.py --dir results/rerun_2026-10/species_tree_aware/local_swap
       python workflow/50_local_swap_summary.py --dir .../local_swap/generax_r1 --sets aa --species-trees v2fix
"""
import argparse
import csv
import json
import os
import re


def au_table(path, names):
    rows, on = [], False
    for line in open(path):
        if line.startswith("Tree      logL"):
            on = True
            continue
        if on:
            if line.startswith("---"):
                continue
            m = re.match(r"\s*(\d+)\s+(-?[\d.]+)\s+([\d.e-]+)\s+.*?([\d.e-]+)\s*[+-]\s*$", line)
            if not m:
                if rows:
                    break
                continue
            rows.append((int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))))
    return {names[i - 1]: {"seq_lnL": l, "seq_dlnL": d, "p_AU": p} for i, l, d, p in rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="results/rerun_2026-10/species_tree_aware/local_swap")
    ap.add_argument("--species-trees", nargs="+", default=["v2fix", "ncbi_plus"])
    ap.add_argument("--sets", nargs="+", default=["aa", "aa3di"])
    args = ap.parse_args()
    out = []
    for s in args.sets:
        names = [l.strip() for l in open(os.path.join(args.dir, f"{s}_trees.names")) if l.strip()]
        au = au_table(os.path.join(args.dir, f"au_{s}.iqtree"), names)
        for nm in names:
            row = {"set": s, "tree": nm, **au.get(nm, {})}
            for st in args.species_trees:
                f = os.path.join(args.dir, "recon", f"{nm}_{st}", "summary.json")
                if not os.path.exists(f):
                    continue
                j = json.load(open(f))
                a = j["anchor_duplication"]
                row[f"rec_{st}"] = round(j["sum_loglik"], 2)
                row[f"rates_{st}"] = "%.3f/%.3f" % tuple(j["rates"])
                row[f"anchor_{st}"] = f"{a['event']}@{a['species_node']}({a['species_node_leaves']})"
                row[f"joint_{st}"] = round(row["seq_lnL"] + j["sum_loglik"], 2) if "seq_lnL" in row else None
            out.append(row)
    # within-pair differences (swapped minus background)
    by = {r["tree"]: r for r in out}
    for r in out:
        base = re.sub(r"_(to197|todeep)$", "", r["tree"])
        if base != r["tree"] and base in by:
            b = by[base]
            r["pair_dseq"] = round(r["seq_lnL"] - b["seq_lnL"], 3)
            for st in args.species_trees:
                if f"rec_{st}" in r and f"rec_{st}" in b:
                    r[f"pair_drec_{st}"] = round(r[f"rec_{st}"] - b[f"rec_{st}"], 2)
                    r[f"pair_djoint_{st}"] = round(r[f"joint_{st}"] - b[f"joint_{st}"], 2)
    keys = []
    for r in out:
        keys += [k for k in r if k not in keys]
    path = os.path.join(args.dir, "local_swap_summary.tsv")
    with open(path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader(); w.writerows(out)
    for r in out:
        print("\t".join(str(r.get(k, "")) for k in keys))
    print(f"-> {path}")


if __name__ == "__main__":
    main()
