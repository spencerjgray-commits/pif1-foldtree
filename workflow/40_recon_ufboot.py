#!/usr/bin/env python3
"""
40_recon_ufboot.py  --  reconcile a sample of bootstrap (UFBoot) gene trees against species trees
(review §1.5: "a real species-tree-aware test would reconcile the UFBoot sample ... and report how often
the duplication lands on Saccharomycotina").

For every sampled gene tree and species tree, the ML UndatedDL scenario (workflow/35, GeneRax 2.0.4
EVAL semantics, D = L = 0.2) is computed and the species node receiving the ScPif1/ScRrm3 duplication
is recorded, together with whether that node is inside Saccharomycotina (all descendant leaves are
Saccharomycotina taxa). The gene-tree-only readout (size of the outgroup-rooted anchor MRCA and whether
it is all-Saccharomycotina) is recorded as well.

Usage (repo root):
  python workflow/40_recon_ufboot.py --trees results/seq_tree/pif1.ufboot --n 200 \
      --species-trees NCBI.nwk V2.nwk --out results/rerun_2026-10/ufboot_recon/aa_ufboot.tsv
"""
import argparse
import csv
import os
import sys
from collections import Counter
from importlib import import_module

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
U = import_module("35_undated_dl")
SS = import_module("33_split_supports")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True, help="file with one newick tree per line")
    ap.add_argument("--n", type=int, default=200, help="number of evenly spaced trees to use (0 = all)")
    ap.add_argument("--species-trees", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", help="short names for the species trees")
    ap.add_argument("--map", default="data/species_tree/gene_species.map")
    ap.add_argument("--root", default="Q9H611", help="outgroup tip for the gene-tree-only readout")
    ap.add_argument("--anchors", default="P07271_Saccharomyces,P38766_Saccharomyces")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    labels = args.labels or [os.path.basename(s) for s in args.species_trees]
    gene2sp = {}
    for line in open(args.map):
        p = line.split()
        if len(p) >= 2:
            gene2sp[p[0]] = p[1]
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    sps = [U.Species(s) for s in args.species_trees]
    sacch_leaves = {t for sp in sps for t in sp.leaves if 147537 in ncbi.get_lineage(int(t))}
    sac_tips = set()
    for r in csv.DictReader(open("manifest.csv")):
        if r["subphylum"] == "Saccharomycotina":
            sac_tips.add(r["accession"])

    lines = [l.strip() for l in open(args.trees) if l.strip()]
    idx = list(range(len(lines))) if args.n in (0, None) or args.n >= len(lines) else \
        [round(i * (len(lines) - 1) / (args.n - 1)) for i in range(args.n)]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    tmp = args.out + ".tree.tmp"
    a1, a2 = args.anchors.split(",")
    rows = []
    with open(args.out, "w") as fh:
        w = None
        for k, i in enumerate(idx):
            open(tmp, "w").write(lines[i] + "\n")
            st = SS.SplitTree(tmp, args.root)
            m = st.mrca([a1.split("_")[0], a2.split("_")[0]])
            gt_size = len(m) if m else None
            gt_sacch = all(t.split("_")[0] in sac_tips for t in m) if m else None
            row = {"tree_index": i, "gt_anchor_mrca_size": gt_size, "gt_anchor_mrca_all_sacch": gt_sacch}
            for lab, sp in zip(labels, sps):
                gt = U.GeneTree(tmp, keep=[g for g in gene2sp if gene2sp[g] in sp.leaves])
                rec = U.UndatedDL(sp, gt, gene2sp, 0.2, 0.2)
                counts, tot = rec.infer()
                best = None
                for cid, (ev, e) in rec.node_event.items():
                    if ev == "Leaf":
                        continue
                    tips = rec.clade_tips(cid)
                    if a1 in tips and a2 in tips and (best is None or len(tips) < best[0]):
                        best = (len(tips), ev, e)
                node = sp.names[best[2]]
                desc = sp.desc[best[2]]
                row[f"{lab}_anchor_node"] = node
                row[f"{lab}_anchor_node_leaves"] = len(desc)
                row[f"{lab}_anchor_in_sacch"] = desc <= sacch_leaves
                row[f"{lab}_anchor_gene_clade"] = best[0]
                row[f"{lab}_D_total"] = tot["D"]
            if w is None:
                w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
                w.writeheader()
            w.writerow(row); fh.flush()
            rows.append(row)
            if (k + 1) % 20 == 0:
                print(f"[40] {k + 1}/{len(idx)} trees", file=sys.stderr, flush=True)
    os.remove(tmp)
    n = len(rows)
    print(f"trees: {n} from {args.trees}")
    print(f"gene-tree readout (rooted on {args.root}): anchor MRCA all-Saccharomycotina in "
          f"{sum(1 for r in rows if r['gt_anchor_mrca_all_sacch'])}/{n}; MRCA size distribution "
          f"{Counter(r['gt_anchor_mrca_size'] for r in rows).most_common(6)}")
    for lab in labels:
        ins = sum(1 for r in rows if r[f"{lab}_anchor_in_sacch"])
        nodes = Counter(f"{r[f'{lab}_anchor_node']}({r[f'{lab}_anchor_node_leaves']})" for r in rows)
        print(f"{lab}: duplication maps inside Saccharomycotina in {ins}/{n} ({100 * ins / n:.1f}%); "
              f"nodes {nodes.most_common(6)}")


if __name__ == "__main__":
    main()
