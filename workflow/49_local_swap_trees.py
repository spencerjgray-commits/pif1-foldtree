#!/usr/bin/env python3
"""
49_local_swap_trees.py  --  trees that differ from a background tree only in where the Pif1 daughter
clade attaches, so the two duplication placements can be scored with everything else held fixed.

Why: in the deep AA tree the Pif1 and Rrm3 Saccharomycotina daughters are 8 edges apart, beyond
GeneRax's maximum SPR radius (5); a GeneRax search from either basin cannot test the other placement
in a single move.

Definitions (trees rooted on hPIF1, Q9H611_Homo):
  daughter(anchor)  the largest clade containing the anchor (ScPif1 P07271 or ScRrm3 P38766), not the
                    other anchor, whose tips all belong to the 197-gene set
                    (results/rerun_2026-10/supports/clade197_supported.txt)
  to197(T)          prune the Pif1 daughter of T and regraft it as sister to the Rrm3 daughter of T
                    -> Saccharomycotina-specific duplication, rest of T unchanged
  todeep(T, D)      prune the Pif1 daughter of T and regraft it as sister to MRCA_T(S), where S is the
                    tip set of the Pif1 daughter's sister in the deep reference tree D
                    -> D's deep placement, rest of T unchanged
The regrafted clade keeps its own branch lengths; the target edge is split in half. Branch lengths are
re-optimised when the trees are scored (IQ-TREE -z), so they only matter as a starting point.

Usage (repo root, pif1 env):
  python workflow/49_local_swap_trees.py --deep D.nwk --tree name=T.nwk [--tree ...] --out-trees OUT.nwk --out-names OUT.txt
Writes, for every --tree: the tree itself, then its swapped version (to197 for trees whose anchors'
MRCA is larger than the 197 clade, todeep otherwise). Prints the clade sizes and the move distance.
"""
import argparse

import dendropy

ROOT = "Q9H611_Homo"
PIF1, RRM3 = "P07271_Saccharomyces", "P38766_Saccharomyces"


def load(path, ns=None):
    t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, rooting="force-rooted",
                          taxon_namespace=ns)
    t.reroot_at_edge(t.find_node_with_taxon_label(ROOT).edge, update_bipartitions=False)
    for nd in t.preorder_node_iter():
        if not nd.is_leaf():
            nd.label = None
    return t


def tips(n):
    return {l.taxon.label for l in n.leaf_iter()}


def daughter(t, anchor, other, c197):
    n = t.find_node_with_taxon_label(anchor)
    best = n
    while n.parent_node is not None:
        n = n.parent_node
        s = tips(n)
        if s <= c197 and other not in s:
            best = n
        else:
            break
    return best


def mrca(t, labels):
    """MRCA by walking ancestors (dendropy's mrca() uses bipartition encodings that go stale after edits)."""
    nodes = [t.find_node_with_taxon_label(x) for x in labels]
    common = [nodes[0]] + list(nodes[0].ancestor_iter())
    for n in nodes[1:]:
        anc = set([n] + list(n.ancestor_iter()))
        common = [x for x in common if x in anc]
    return common[0]


def path_edges(a, b):
    anc = lambda n: [n] + list(n.ancestor_iter())
    pa, pb = anc(a), anc(b)
    common = next(x for x in pa if x in pb)
    return pa.index(common) + pb.index(common)


def prune(t, node):
    par = node.parent_node
    par.remove_child(node)
    if len(par.child_nodes()) == 1:          # suppress the now-unifurcating parent
        child = par.child_nodes()[0]
        gp = par.parent_node
        if gp is None:
            t.seed_node = child
            child.parent_node = None
        else:
            el = (par.edge.length or 0) + (child.edge.length or 0)
            gp.remove_child(par)
            par.remove_child(child)
            gp.add_child(child)
            child.edge.length = el
    return node


def regraft(target, clade):
    par = target.parent_node
    el = target.edge.length or 0.0
    par.remove_child(target)
    new = par.new_child()
    new.edge.length = el / 2
    new.add_child(target)
    target.edge.length = el / 2
    new.add_child(clade)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deep", required=True, help="deep reference tree D (defines the deep attachment point)")
    ap.add_argument("--tree", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--clade", default="results/rerun_2026-10/supports/clade197_supported.txt")
    ap.add_argument("--out-trees", required=True)
    ap.add_argument("--out-names", required=True)
    args = ap.parse_args()
    c197 = {l.strip() for l in open(args.clade) if l.strip()}

    D = load(args.deep)
    pD = daughter(D, PIF1, RRM3, c197)
    sis = [c for c in pD.parent_node.child_nodes() if c is not pD]
    S = set().union(*(tips(c) for c in sis))
    print(f"deep reference: Pif1 daughter {len(tips(pD))} tips; its sister clade {len(S)} tips")

    out, names = [], []
    for spec in args.tree:
        nm, path = spec.split("=", 1)
        T = load(path)
        mrca_ = mrca(T, [PIF1, RRM3])
        n_mrca = len(tips(mrca_))
        p = daughter(T, PIF1, RRM3, c197)
        r = daughter(T, RRM3, PIF1, c197)
        print(f"{nm}: anchor MRCA {n_mrca} tips; Pif1 daughter {len(tips(p))}; Rrm3 daughter {len(tips(r))}; "
              f"edges between their stems {path_edges(p, r)}")
        out.append(T.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True).strip())
        names.append(nm)
        T2 = load(path)
        p2 = daughter(T2, PIF1, RRM3, c197)
        if n_mrca > len(tips(p)) + len(tips(r)):          # deep background -> make the 197 configuration
            r2 = daughter(T2, RRM3, PIF1, c197)
            prune(T2, p2)
            regraft(r2, p2)
            tag = "to197"
        else:                                              # 197 background -> put Pif1 where D has it
            prune(T2, p2)
            tgt = mrca(T2, sorted(S & tips(T2.seed_node)))
            if tgt is T2.seed_node:                        # S = everything but the root tip: use its edge
                tgt = next(c for c in T2.seed_node.child_nodes() if ROOT not in tips(c))
            regraft(tgt, p2)
            tag = "todeep"
        m2 = mrca(T2, [PIF1, RRM3])
        print(f"   {nm}_{tag}: anchor MRCA now {len(tips(m2))} tips")
        out.append(T2.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True).strip())
        names.append(f"{nm}_{tag}")
    open(args.out_trees, "w").write("\n".join(out) + "\n")
    open(args.out_names, "w").write("\n".join(names) + "\n")


if __name__ == "__main__":
    main()
