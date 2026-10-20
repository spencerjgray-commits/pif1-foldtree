#!/usr/bin/env python3
"""
38_recon_species_trees.py  --  reconcile one or more gene trees against one or more species trees
with workflow/35_undated_dl.py (GeneRax 2.0.4 EVAL / UndatedDL, D = L = 0.2) and tabulate:
  - the species node the ScPif1/ScRrm3 duplication maps to (and its leaf count)
  - D at node_448 (the Saccharomycotina crown in every species tree used here)
  - D at the anchor-mapped node, and its rank among all species nodes
  - total D / SL, and the top-5 nodes by D

Usage (repo root):
  python workflow/38_recon_species_trees.py --gene-trees A.treefile B.treefile \
      --species-trees S1.nwk S2.nwk ... --out results/.../table.tsv [--outdir-prefix DIR]
"""
import argparse
import csv
import os
import sys
from importlib import import_module

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
U = import_module("35_undated_dl")


def run(gene_tree, species_tree, gene2sp, anchors, outdir=None):
    sp = U.Species(species_tree)
    gt = U.GeneTree(gene_tree, keep=[g for g in gene2sp if gene2sp[g] in sp.leaves])
    m = U.UndatedDL(sp, gt, gene2sp, 0.2, 0.2)
    counts, tot = m.infer()
    a1, a2 = anchors
    anchor = None
    for cid, (ev, e) in m.node_event.items():
        if ev == "Leaf":
            continue
        tips = m.clade_tips(cid)
        if a1 in tips and a2 in tips and (anchor is None or len(tips) < anchor[0]):
            anchor = (len(tips), ev, e)
    D = {sp.names[e]: c[1] for e, c in counts.items()}
    ranked = sorted(D.items(), key=lambda kv: -kv[1])
    an = sp.names[anchor[2]]
    rank_of = lambda node: 1 + sum(1 for _, d in ranked if d > D.get(node, 0))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "speciesEventCounts.txt"), "w") as fh:
            for e in range(sp.n):
                c = counts.get(e, [0, 0, 0])
                fh.write(f"{sp.names[e]} {c[0]} {c[1]} {c[2]} 0 \n")
    return {"gene_tree": gene_tree, "species_tree": species_tree, "n_genes": len(gt.leafname),
            "anchor_clade_genes": anchor[0], "anchor_event": anchor[1], "anchor_node": an,
            "anchor_node_leaves": len(sp.desc[anchor[2]]), "D_anchor_node": D.get(an, 0),
            "rank_anchor_node": rank_of(an), "D_node_448": D.get("node_448", 0),
            "rank_node_448": rank_of("node_448"), "D_node_578": D.get("node_578", 0),
            "D_total": tot["D"], "SL_total": tot["SL"],
            "top5": "; ".join(f"{n}({len(sp.desc[sp.names.index(n)])}):{d}" for n, d in ranked[:5])}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-trees", nargs="+", required=True)
    ap.add_argument("--species-trees", nargs="+", required=True)
    ap.add_argument("--map", default="data/species_tree/gene_species.map")
    ap.add_argument("--anchors", default="P07271_Saccharomyces,P38766_Saccharomyces")
    ap.add_argument("--out", required=True)
    ap.add_argument("--outdir-prefix", help="also write per-run speciesEventCounts.txt under this dir")
    args = ap.parse_args()
    gene2sp = {}
    for line in open(args.map):
        p = line.split()
        if len(p) >= 2:
            gene2sp[p[0]] = p[1]
    rows = []
    for g in args.gene_trees:
        for s in args.species_trees:
            od = None
            if args.outdir_prefix:
                od = os.path.join(args.outdir_prefix, os.path.basename(g).split(".")[0] + "__" +
                                  os.path.basename(s).rsplit(".", 1)[0])
            r = run(g, s, gene2sp, args.anchors.split(","), od)
            rows.append(r)
            print(f"{os.path.basename(g)} x {os.path.basename(s)}: anchor {r['anchor_node']}({r['anchor_node_leaves']}) "
                  f"D={r['D_anchor_node']} rank {r['rank_anchor_node']} | node_448 D={r['D_node_448']} | "
                  f"D/SL {r['D_total']}/{r['SL_total']} | {r['top5']}", flush=True)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
