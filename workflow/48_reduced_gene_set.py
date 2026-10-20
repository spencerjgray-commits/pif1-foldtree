#!/usr/bin/env python3
"""
48_reduced_gene_set.py  --  a reduced but representative gene set for GeneRax SPR searches.

A GeneRax 2.0.4 SPR search on all 957 genes / 720 species applies one move per round at ~15 min a
round on this Mac (2026-10-16), i.e. days per run. The reduced set:
  - keeps every Saccharomycotina gene (205 genes, 105 species: the whole 197-gene clade and the 8
    Saccharomycotina genes outside it);
  - keeps, for every other NCBI order (class/phylum where a species has no order), up to --per-order
    species drawn at random (--seed), with ALL genes of each kept species, so within-species paralogy
    is not broken up; hPIF1 (Homo sapiens) is the only Primates species and is always kept.
The species trees are pruned to the kept species (internal labels such as node_0, node_448 and sac_1
survive pruning because those nodes keep two children); the starting gene trees are pruned to the
kept genes; the alignment is subset by row.

Outputs (--outdir): genes.txt, composition.tsv, aln.fasta, gene_species.map, <name>.species.nwk for
every --species-tree NAME=PATH, start_<name>.nwk for every --start NAME=PATH, report.txt

Usage (repo root, pif1 env):
  python workflow/48_reduced_gene_set.py --outdir OUT --species-tree v2fix=S1 --species-tree ncbi=S2 \
     --start deep=results/seq_tree/pif1.treefile --start 197=...
"""
import argparse
import csv
import os
import random
from collections import Counter, defaultdict

import dendropy
from ete3 import NCBITaxa

SACCHAROMYCOTINA = 147537


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", default="results/rerun_2026-10/species_tree_aware/species_trees/gene_species_957.map")
    ap.add_argument("--aln", default="results/seq_tree/aln.trim.fasta")
    ap.add_argument("--species-tree", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--start", action="append", default=[], help="NAME=PATH")
    ap.add_argument("--per-order", type=int, default=3)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ncbi = NCBITaxa()
    rep = []
    log = lambda s: (rep.append(s), print(s))

    g2s = {}
    for line in open(args.map):
        p = line.split()
        if len(p) >= 2:
            g2s[p[0]] = p[1]
    sp_genes = defaultdict(list)
    for g, s in g2s.items():
        sp_genes[s].append(g)

    def group(s):
        lin = ncbi.get_lineage(int(s))
        ranks = ncbi.get_rank(lin)
        for r in ("order", "class", "phylum"):
            hit = [x for x in lin if ranks[x] == r]
            if hit:
                return hit[-1], r
        return lin[-1], "none"

    rng = random.Random(args.seed)
    keep_sp, comp = set(), []
    groups = defaultdict(list)
    for s in sorted(sp_genes, key=int):
        if SACCHAROMYCOTINA in ncbi.get_lineage(int(s)):
            keep_sp.add(s)
        else:
            groups[group(s)].append(s)
    for (gid, rank), sps in sorted(groups.items(), key=lambda kv: kv[0][0]):
        pick = sorted(rng.sample(sps, min(args.per_order, len(sps))), key=int)
        keep_sp.update(pick)
        name = ncbi.get_taxid_translator([gid]).get(gid, str(gid))
        comp.append({"group": name, "rank": rank, "species_total": len(sps),
                     "genes_total": sum(len(sp_genes[x]) for x in sps), "species_kept": len(pick),
                     "genes_kept": sum(len(sp_genes[x]) for x in pick), "kept_taxids": ",".join(pick)})
    genes = sorted(g for g, s in g2s.items() if s in keep_sp)
    n_sac = sum(1 for g in genes if SACCHAROMYCOTINA in ncbi.get_lineage(int(g2s[g])))
    log(f"kept {len(genes)} of {len(g2s)} genes ({n_sac} Saccharomycotina) from {len(keep_sp)} of {len(sp_genes)} species; "
        f"{len(groups)} outgroup orders/classes, up to {args.per_order} species each (seed {args.seed})")
    with open(os.path.join(args.outdir, "composition.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comp[0]), delimiter="\t")
        w.writeheader(); w.writerows(comp)
    open(os.path.join(args.outdir, "genes.txt"), "w").write("\n".join(genes) + "\n")
    with open(os.path.join(args.outdir, "gene_species.map"), "w") as fh:
        for g in genes:
            fh.write(f"{g} {g2s[g]}\n")

    seqs, name = {}, None
    for line in open(args.aln):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    with open(os.path.join(args.outdir, "aln.fasta"), "w") as fh:
        for g in genes:
            fh.write(f">{g}\n{''.join(seqs[g])}\n")

    for spec in args.species_tree:
        nm, path = spec.split("=", 1)
        t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, rooting="force-rooted")
        t.retain_taxa_with_labels(keep_sp)
        t.suppress_unifurcations()
        labels = {nd.label for nd in t.preorder_node_iter() if not nd.is_leaf()}
        out = os.path.join(args.outdir, f"{nm}.species.nwk")
        t.write(path=out, schema="newick", suppress_rooting=True, unquoted_underscores=True)
        log(f"species tree {nm}: {len(list(t.leaf_node_iter()))} leaves; node_0 {'node_0' in labels}, "
            f"node_448 {'node_448' in labels}, sac_1 {'sac_1' in labels}")
    for spec in args.start:
        nm, path = spec.split("=", 1)
        t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, rooting="force-unrooted")
        t.retain_taxa_with_labels(set(genes))
        t.suppress_unifurcations()
        t.collapse_basal_bifurcation()
        for nd in t.preorder_node_iter():
            if not nd.is_leaf():
                nd.label = None
        out = os.path.join(args.outdir, f"start_{nm}.nwk")
        t.write(path=out, schema="newick", suppress_rooting=True, unquoted_underscores=True,
                suppress_internal_node_labels=True)
        log(f"start tree {nm}: {len(list(t.leaf_node_iter()))} leaves -> {out}")
    open(os.path.join(args.outdir, "report.txt"), "w").write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
