#!/usr/bin/env python3
"""
51_generax_tree_distances.py  --  pairwise unrooted Robinson-Foulds distances between GeneRax gene trees
(e.g. the radius-1 and final trees of the three SPR starts), and the sister group of the Saccharomycotina
Pif1 daughter clade in each tree (rooted on hPIF1; daughter as defined in workflow/49).

Usage (repo root, pif1 env):
  python workflow/51_generax_tree_distances.py --tree NAME=PATH [--tree ...] [--out OUT.txt]
"""
import argparse
import importlib.util
import itertools
import os

import dendropy
from dendropy.calculate import treecompare

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ls49", os.path.join(HERE, "49_local_swap_trees.py"))
ls49 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ls49)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--clade", default="results/rerun_2026-10/supports/clade197_supported.txt")
    ap.add_argument("--out")
    args = ap.parse_args()
    c197 = {l.strip() for l in open(args.clade) if l.strip()}
    ns = dendropy.TaxonNamespace()
    names, trees, lines = [], [], []
    for s in args.tree:
        nm, path = s.split("=", 1)
        t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, taxon_namespace=ns,
                              rooting="force-unrooted")
        t.encode_bipartitions()
        names.append(nm); trees.append(t)
        r = ls49.load(path)
        p = ls49.daughter(r, ls49.PIF1, ls49.RRM3, c197)
        sis = set().union(*(ls49.tips(c) for c in p.parent_node.child_nodes() if c is not p))
        lines.append(f"{nm}: {len(t.leaf_nodes())} tips; Pif1 daughter {len(ls49.tips(p))}, its sister {len(sis)} tips"
                     + (f" ({', '.join(sorted(sis))})" if len(sis) <= 20 else ""))
    n_int = len(trees[0].leaf_nodes()) - 3
    lines.append(f"pairwise unrooted RF (max {2 * n_int}):")
    for (a, ta), (b, tb) in itertools.combinations(zip(names, trees), 2):
        rf = treecompare.symmetric_difference(ta, tb)
        lines.append(f"  {a} vs {b}: RF {rf} (normalised {rf / (2 * n_int):.3f})")
    txt = "\n".join(lines)
    print(txt)
    if args.out:
        open(args.out, "w").write(txt + "\n")


if __name__ == "__main__":
    main()
