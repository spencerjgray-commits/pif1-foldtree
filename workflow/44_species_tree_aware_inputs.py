#!/usr/bin/env python3
"""
44_species_tree_aware_inputs.py  --  species trees and a complete gene->species map for the
species-tree-aware tests (GeneRax SPR, ALE) of where the Pif1/Rrm3 duplication sits.

Fixes two input problems found in the 2026-10 re-run:

  1. Map. data/species_tree/gene_species.map covers 941 of the 957 genes. The 16 missing genes carry
     a taxid that is not a leaf of the species trees:
       - 14 carry a species-level taxid (e.g. Rhizoctonia solani 456999) while the leaf is a strain
         of the same species (R. solani AG-1 IB 1108050). They are mapped to the leaf of the same
         NCBI species; where a species has more than one leaf (R. solani: 1108050, 1423351) the
         lowest taxid is used.
       - 2 (Wickerhamomyces ciferrii K0KM86, K0KZ61, both inside the 197-gene clade) carry the merged
         taxid 1206466 -> 1041607, a species with no leaf at all. A leaf 1041607 is added.
  2. Grafted tree. workflow/37 (with its genus-fallback fix) is re-run on a v1 backbone that carries
     the added leaf, so W. ciferrii is placed on the Shen tree like every other budding yeast.

An added leaf goes in as sister to the MRCA of the existing leaves that share its lowest NCBI
ancestor (genus, else family, ...). The new parent node is labelled add_<taxid>.

Outputs (--outdir):
  gene_species_957.map               gene  species-leaf  (all 957 genes)
  map_changes.tsv                    the 16 remapped genes and how
  ncbi_plus.nwk                      NCBI binary species tree (GeneRax inferred_species_tree) + added leaves
  v1_plus.nwk                        v1 grafted binary tree + added leaves (input to workflow/37)
  v2fix/grafted_species_v2.binary.nwk  rebuilt graft (fixed fallback, W. ciferrii placed via Shen)
  species_trees_report.txt

Run (repo root, pif1 env):
  python workflow/44_species_tree_aware_inputs.py
"""
import argparse
import csv
import os
import subprocess
import sys
import warnings
from collections import defaultdict

import dendropy
from ete3 import NCBITaxa

warnings.filterwarnings("ignore", message="taxid .* was translated")


def lineage(ncbi, taxid):
    return ncbi.get_lineage(int(taxid))          # resolves merged taxids


def species_of(ncbi, taxid):
    lin = lineage(ncbi, taxid)
    r = ncbi.get_rank(lin)
    sp = [x for x in lin if r[x] == "species"]
    return sp[-1] if sp else lin[-1]


def add_leaf(tree, new_label, ncbi, report):
    """Attach new_label as sister to the MRCA of the leaves sharing its lowest NCBI ancestor."""
    leaves = [l.taxon.label for l in tree.leaf_node_iter()]
    lin_new = lineage(ncbi, new_label)
    lins = {l: set(lineage(ncbi, l)) for l in leaves}
    for anc in reversed(lin_new[:-1]):
        hits = [l for l in leaves if anc in lins[l]]
        if hits:
            break
    node = tree.mrca(taxon_labels=hits) if len(hits) > 1 else tree.find_node_with_taxon_label(hits[0])
    par = node.parent_node
    par.remove_child(node)
    holder = par.new_child(label=f"add_{new_label}")
    holder.edge.length = 1.0
    holder.add_child(node)
    leaf = holder.new_child(taxon=tree.taxon_namespace.require_taxon(label=new_label))
    leaf.edge.length = 1.0
    name = ncbi.get_taxid_translator([anc]).get(anc, anc)
    rank = ncbi.get_rank([anc]).get(anc, "")
    report.append(f"   added leaf {new_label} as sister to the MRCA of {len(hits)} leaves of {rank} {name} ({anc})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--old-map", default="data/species_tree/gene_species.map")
    ap.add_argument("--ncbi", default="results/reconciliation/aa3di/eval/species_trees/inferred_species_tree.newick")
    ap.add_argument("--v1", default="results/reconciliation/aa3di_grafted/run/species_trees/starting_species_tree.newick")
    ap.add_argument("--outdir", default="results/rerun_2026-10/species_tree_aware/species_trees")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ncbi = NCBITaxa()
    rep = []
    log = lambda s: (rep.append(s), print(s))

    trees = {}
    for key, path in (("ncbi", args.ncbi), ("v1", args.v1)):
        trees[key] = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, rooting="force-rooted")
    leaves = sorted(l.taxon.label for l in trees["ncbi"].leaf_node_iter())
    assert leaves == sorted(l.taxon.label for l in trees["v1"].leaf_node_iter()), "NCBI and v1 leaf sets differ"
    sp2leaves = defaultdict(list)
    for l in leaves:
        sp2leaves[species_of(ncbi, l)].append(l)

    old = {}
    for line in open(args.old_map):
        p = line.split()
        if len(p) >= 2:
            old[p[0]] = p[1]
    tips = list(csv.DictReader(open(args.tip_map), delimiter="\t"))
    gmap, changes, to_add = {}, [], set()
    for r in tips:
        g, tx = r["tip_label"], r["taxid"]
        if g in old:
            gmap[g] = old[g]
            continue
        cur = str(lineage(ncbi, tx)[-1])
        sp = species_of(ncbi, tx)
        if cur in leaves:
            leaf, how = cur, "merged taxid is a leaf"
        elif sp2leaves.get(sp):
            cands = sorted(sp2leaves[sp], key=int)
            leaf = cands[0]
            how = "same NCBI species" + (f" (1 of {len(cands)} leaves: {','.join(cands)})" if len(cands) > 1 else "")
        else:
            leaf, how = str(sp), "no leaf for this species: leaf added"
            to_add.add(leaf)
        gmap[g] = leaf
        changes.append({"gene": g, "taxid": tx, "current_taxid": cur, "organism": r["organism"][:60],
                        "species": sp, "species_name": ncbi.get_taxid_translator([sp]).get(sp, ""),
                        "leaf": leaf, "leaf_name": ncbi.get_taxid_translator([int(leaf)]).get(int(leaf), ""), "how": how})
    log(f"genes: {len(tips)}; in old map: {sum(1 for r in tips if r['tip_label'] in old)}; remapped: {len(changes)}; leaves to add: {sorted(to_add)}")
    for c in changes:
        log(f"   {c['gene']:28s} {c['taxid']:>8s} -> {c['leaf']:>8s} {c['leaf_name'][:45]:45s} [{c['how']}]")
    with open(os.path.join(args.outdir, "map_changes.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(changes[0]), delimiter="\t")
        w.writeheader(); w.writerows(changes)
    with open(os.path.join(args.outdir, "gene_species_957.map"), "w") as fh:
        for r in tips:
            fh.write(f"{r['tip_label']} {gmap[r['tip_label']]}\n")

    for key in ("ncbi", "v1"):
        log(f"{key}:")
        for leaf in sorted(to_add):
            add_leaf(trees[key], leaf, ncbi, rep)
        t = trees[key]
        assert all(len(n.child_nodes()) in (0, 2) for n in t.preorder_node_iter()), "not binary"
        out = os.path.join(args.outdir, f"{key}_plus.nwk")
        t.write(path=out, schema="newick", suppress_rooting=True, unquoted_underscores=True)
        log(f"   wrote {out}: {len(list(t.leaf_node_iter()))} leaves")

    # rebuilt graft: workflow/37 (fixed fallback) on the v1 backbone with the added leaf
    v2dir = os.path.join(args.outdir, "v2fix")
    cmd = [sys.executable, "workflow/37_graft_species_tree_v2.py", "--v1-binary", os.path.join(args.outdir, "v1_plus.nwk"),
           "--outdir", v2dir]
    log("running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    v2 = dendropy.Tree.get(path=os.path.join(v2dir, "grafted_species_v2.binary.nwk"), schema="newick",
                           preserve_underscores=True, rooting="force-rooted")
    v2_leaves = {l.taxon.label for l in v2.leaf_node_iter()}
    missing = sorted({s for s in gmap.values() if s not in v2_leaves})
    log(f"v2fix: {len(v2_leaves)} leaves; genes whose species leaf is missing: {missing}")
    open(os.path.join(args.outdir, "species_trees_report.txt"), "w").write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
