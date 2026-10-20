#!/usr/bin/env python3
"""
37_graft_species_tree_v2.py  --  rebuild the grafted species tree (2026-10-15 review §1.5, MET-M9, RES-R3).

The v1 graft (workflow/16) matched Shen et al. 2018 tips by the organism string only, so six species
were missed through synonymy (Candida auris = Candidozyma auris, ...), and every unmatched budding
yeast without a *sampled-and-matched* congener was attached at the base of the grafted clade: a 27-way
polytomy that was then arbitrarily binarised for GeneRax.

v2:
  1. Every Shen tip and every sampled Saccharomycotina taxid is mapped to its CURRENT NCBI species,
     genus and family (local ete3 NCBI taxonomy, which resolves synonyms such as Candida_auris ->
     Candidozyma auris). The resulting name table is written out (shen2018_synonyms.tsv).
  2. Each sampled taxon is placed on the FULL 332-tip Shen tree:
       exact   same NCBI species as a Shen tip            -> takes that tip's place
       genus   genus present in Shen                       -> sister to the Shen genus MRCA
       family  genus absent, family present in Shen        -> sister to the Shen family MRCA
       root    neither                                     -> child of the Shen root
     then the Shen tree is pruned to the sampled taxa.
  3. Polytomies that remain (several taxa attached at one point, or strains of one species) are
     resolved (a) deterministically (caterpillar in taxid order) for the main tree and (b) at random
     for --n-random replicate trees, so the sensitivity of D to the resolution can be measured.
  4. The new Saccharomycotina subtree replaces the Saccharomycotina node of the v1 binary species tree
     (GeneRax's starting_species_tree.newick), so every non-Saccharomycotina node and label is unchanged.
     The new subtree's root keeps the label node_448.

Outputs (--outdir):
  shen2018_synonyms.tsv            per sampled taxid: names, Shen tip, placement level
  grafted_species_v2.binary.nwk    main binary species tree
  grafted_species_v2.rand<i>.nwk   random polytomy resolutions
  graft_v2_report.txt

Run (repo root): python workflow/37_graft_species_tree_v2.py
"""
import argparse
import csv
import os
import random
import sys
from collections import defaultdict

import dendropy
from ete3 import NCBITaxa

SACCHAROMYCOTINA = 147537
# Shen labels that the local NCBI names table does not resolve, mapped by hand to the current name
MANUAL_SHEN_NAMES = {"Candida_parapsilosis": "Lodderomyces parapsilosis"}


def ranks_of(ncbi, taxid):
    lin = ncbi.get_lineage(taxid)
    r = ncbi.get_rank(lin)
    pick = lambda rank: next((x for x in reversed(lin) if r[x] == rank), None)
    return pick("species"), pick("genus"), pick("family"), lin


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v1-binary", default="results/reconciliation/aa3di_grafted/run/species_trees/starting_species_tree.newick")
    ap.add_argument("--shen", default="data/species_tree/shen2018_timetree.newick")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--outdir", default="results/rerun_2026-10/grafted")
    ap.add_argument("--n-random", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ncbi = NCBITaxa()
    rep = []
    log = lambda s: (rep.append(s), print(s))

    org = {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        org.setdefault(r["taxid"], r["organism"])

    v1 = dendropy.Tree.get(path=args.v1_binary, schema="newick", preserve_underscores=True, rooting="force-rooted")
    v1_leaves = [l.taxon.label for l in v1.leaf_node_iter()]
    sac = [t for t in v1_leaves if SACCHAROMYCOTINA in ncbi.get_lineage(int(t))]
    sac_node = v1.mrca(taxon_labels=sac)
    assert len(list(sac_node.leaf_iter())) == len(sac), "Saccharomycotina not monophyletic in v1 tree"
    log(f"v1 binary tree: {len(v1_leaves)} leaves; Saccharomycotina node {sac_node.label} with {len(sac)} leaves")

    # ---- Shen tips -> NCBI
    shen = dendropy.Tree.get(path=args.shen, schema="newick", preserve_underscores=True, rooting="force-rooted")
    shen_tips = [l.taxon.label for l in shen.leaf_node_iter()]
    qname = lambda s: MANUAL_SHEN_NAMES.get(s, s.replace("_", " "))
    tr = ncbi.get_name_translator([qname(s) for s in shen_tips])
    shen_info = {}
    for s in shen_tips:
        t = tr.get(qname(s))
        shen_info[s] = ranks_of(ncbi, t[0])[:3] if t else (None, None, None)
    # Untranslated tip: fall back to the genus name in the label, but only for a genus that no translated
    # tip represents. Form genera such as Candida are polyphyletic in Shen, so two untranslated labels
    # (Candida_corydali, Candida_carpophila) used to stretch the 'Candida' target to an 88-tip MRCA.
    # (Fix of 2026-10-16; the outputs in results/rerun_2026-10/grafted/ predate it.)
    translated_genera = {v[1] for v in shen_info.values() if v[0]}
    for s in shen_tips:
        if shen_info[s][0]:
            continue
        g = ncbi.get_name_translator([s.split("_")[0]]).get(s.split("_")[0])
        if g:
            _, gg, ff, _ = ranks_of(ncbi, g[0])
            if gg not in translated_genera:
                shen_info[s] = (None, gg, ff)
    by_sp, by_gen, by_fam = defaultdict(list), defaultdict(list), defaultdict(list)
    for s, (spc, gen, fam) in shen_info.items():
        if spc: by_sp[spc].append(s)
        if gen: by_gen[gen].append(s)
        if fam: by_fam[fam].append(s)
    log(f"Shen tips: {len(shen_tips)}; translated to an NCBI species: {sum(1 for v in shen_info.values() if v[0])}")

    # ---- place each sampled taxon
    name = lambda x: ncbi.get_taxid_translator([x]).get(x, "") if x else ""
    placements = []      # (taxid, level, target node, shen_label)
    rows = []
    for t in sac:
        spc, gen, fam, _ = ranks_of(ncbi, int(t))
        if spc in by_sp:
            lvl, labs = "exact", by_sp[spc][:1]
        elif gen in by_gen:
            lvl, labs = "genus", by_gen[gen]
        elif fam in by_fam:
            lvl, labs = "family", by_fam[fam]
        else:
            lvl, labs = "root", []
        if labs:
            node = shen.mrca(taxon_labels=labs) if len(labs) > 1 else shen.find_node_with_taxon_label(labs[0])
        else:
            node = shen.seed_node
        placements.append((t, lvl, node))
        our_norm = "_".join(org.get(t, "").split()[:2])
        rows.append({"taxid": t, "organism_in_dataset": org.get(t, ""), "ncbi_current_species": name(spc),
                     "ncbi_genus": name(gen), "ncbi_family": name(fam), "placement": lvl,
                     "shen_tip_or_clade": labs[0] if lvl == "exact" else (f"MRCA of {len(labs)} {name(gen) if lvl=='genus' else name(fam)} tips" if labs else "Shen root"),
                     "synonym_hit": "yes" if lvl == "exact" and labs[0].lower() != our_norm.lower() else ""})
    with open(os.path.join(args.outdir, "shen2018_synonyms.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    cnt = defaultdict(int)
    for _, lvl, _ in placements: cnt[lvl] += 1
    log(f"placements: {dict(cnt)}; exact matches through a synonym/current name: "
        f"{sum(1 for r in rows if r['synonym_hit'])}")
    for r in rows:
        if r["synonym_hit"]:
            log(f"   synonym: {r['organism_in_dataset'][:45]} -> {r['shen_tip_or_clade']}")

    # ---- attach to the Shen tree (all targets resolved on the original tree first)
    tax_ns = shen.taxon_namespace
    groups = defaultdict(list)
    for t, lvl, node in placements:
        groups[id(node)].append((t, lvl, node))
    for key, items in groups.items():
        node = items[0][2]
        new = [dendropy.Node(taxon=tax_ns.require_taxon(label=t)) for t, _, _ in items]
        if node is shen.seed_node and items[0][1] == "root":
            for n in new: node.add_child(n)
            continue
        par = node.parent_node
        if par is None:
            for n in new: node.add_child(n)
            continue
        par.remove_child(node)
        holder = par.new_child()
        holder.add_child(node)
        for n in new: holder.add_child(n)
    sub = shen.extract_tree_with_taxa_labels(labels=sac)
    sub.suppress_unifurcations()
    polys = [nd for nd in sub.preorder_node_iter() if len(nd.child_nodes()) > 2]
    log(f"pruned Saccharomycotina subtree: {len(list(sub.leaf_node_iter()))} leaves; "
        f"{len(polys)} polytomies to resolve (sizes {sorted(len(p.child_nodes()) for p in polys)})")
    for p in polys:
        kids = [c.taxon.label if c.is_leaf() else f"[{len(list(c.leaf_iter()))}-leaf clade]" for c in p.child_nodes()]
        log(f"   polytomy: {kids}")

    def resolve(tree, rng):
        for nd in list(tree.postorder_node_iter()):
            kids = nd.child_nodes()
            if len(kids) <= 2: continue
            if rng is None:
                kids = sorted(kids, key=lambda c: (not c.is_leaf(), c.taxon.label if c.is_leaf() else ""))
            else:
                kids = kids[:]; rng.shuffle(kids)
            for c in kids: nd.remove_child(c)
            cur = kids[0]
            for c in kids[1:-1]:
                m = dendropy.Node(); m.add_child(cur); m.add_child(c); cur = m
            nd.add_child(cur); nd.add_child(kids[-1])
        return tree

    def splice(subtree):
        v = dendropy.Tree.get(path=args.v1_binary, schema="newick", preserve_underscores=True, rooting="force-rooted")
        leaves = [l.taxon.label for l in v.leaf_node_iter()]
        sn = v.mrca(taxon_labels=sac)
        s = dendropy.Tree.get(data=subtree.as_string(schema="newick"), schema="newick",
                              preserve_underscores=True, rooting="force-rooted", taxon_namespace=v.taxon_namespace)
        i = 0
        for nd in s.preorder_node_iter():
            if not nd.is_leaf():
                nd.label = "node_448" if nd is s.seed_node else f"sac_{i}"
                i += 1
            nd.edge.length = 1.0
        par = sn.parent_node
        par.remove_child(sn)
        par.add_child(s.seed_node)
        out = sorted(l.taxon.label for l in v.leaf_node_iter())
        assert out == sorted(leaves), "leaf set changed"
        assert all(len(n.child_nodes()) in (0, 2) for n in v.preorder_node_iter()), "not binary"
        return v

    main_sub = resolve(dendropy.Tree(sub), None)
    t_main = splice(main_sub)
    t_main.write(path=os.path.join(args.outdir, "grafted_species_v2.binary.nwk"), schema="newick",
                 suppress_rooting=True, unquoted_underscores=True)
    rng = random.Random(args.seed)
    for i in range(args.n_random):
        ti = splice(resolve(dendropy.Tree(sub), rng))
        ti.write(path=os.path.join(args.outdir, f"grafted_species_v2.rand{i:02d}.nwk"), schema="newick",
                 suppress_rooting=True, unquoted_underscores=True)
    # Saccharomycotina subtree structure: what sits at the base
    kids = main_sub.seed_node.child_nodes()
    desc = lambda n: [org.get(l.taxon.label, l.taxon.label)[:30] for l in n.leaf_iter()]
    log("v2 Saccharomycotina root split: " + " | ".join(
        f"{len(list(k.leaf_iter()))} leaves ({', '.join(desc(k)[:3])}{'...' if len(desc(k)) > 3 else ''})" for k in kids))
    open(os.path.join(args.outdir, "graft_v2_report.txt"), "w").write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
