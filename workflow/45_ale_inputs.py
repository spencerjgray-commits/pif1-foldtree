#!/usr/bin/env python3
"""
45_ale_inputs.py  --  write ALE inputs: gene-tree samples relabelled to ALE's naming and a species tree
without internal labels.

ALE reads the species of a gene from its name: with separators=_ ALEml_undated takes the first
'_'-separated field (species_field = 0). Gene labels ACC_Genus are therefore rewritten as
<species leaf>_ACC using the 957-gene map from workflow/44. ALE parses internal species-tree labels as
bootstrap values, so they are dropped, and it needs branch lengths, so every branch gets length 1; workflow/46 maps ALE's branch numbers back to node names by
leaf set.

  --trees      one or more newick tree files (UFBoot samples); concatenated in the order given
  --every      keep every k-th tree of each file (default 1 = all)
  --prune      prune the species tree to the species present in the sample (for reduced gene sets)
  --genes      keep only these genes (one per line, e.g. workflow/48 genes.txt): every tree is pruned

Usage (repo root, pif1 env):
  python workflow/45_ale_inputs.py --trees results/seq_tree/pif1.ufboot --out DIR/aa_LG.treelist \
      --species-tree S.nwk --species-out DIR/S.ale.nwk
"""
import argparse
import os
import re

import dendropy

LEAF = re.compile(r"(?<=[(,])([^(),:;]+)(?=[:),;])")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", nargs="+", required=True)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--map", default="results/rerun_2026-10/species_tree_aware/species_trees/gene_species_957.map")
    ap.add_argument("--out", required=True)
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--species-out", required=True)
    ap.add_argument("--prune", action="store_true")
    ap.add_argument("--genes", help="file with the genes to keep; trees are pruned to them")
    args = ap.parse_args()

    g2s = {}
    for line in open(args.map):
        p = line.split()
        if len(p) >= 2:
            g2s[p[0]] = p[1]
    rename = {g: f"{s}_{g.split('_')[0]}" for g, s in g2s.items()}
    assert len(set(rename.values())) == len(rename), "accessions not unique"

    keep = {l.strip() for l in open(args.genes) if l.strip()} if args.genes else None
    present, n_out = set(), 0
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        for path in args.trees:
            for i, line in enumerate(open(path)):
                line = line.strip()
                if not line or i % args.every:
                    continue
                if keep is not None:
                    t = dendropy.Tree.get(data=line, schema="newick", preserve_underscores=True, rooting="force-unrooted")
                    t.retain_taxa_with_labels(keep)
                    t.suppress_unifurcations()
                    line = t.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True,
                                       suppress_internal_node_labels=True).strip()
                names = LEAF.findall(line)
                missing = [n for n in names if n not in rename]
                if missing:
                    raise SystemExit(f"{path} tree {i}: unmapped leaves {missing[:5]}")
                present.update(g2s[n] for n in names)
                fh.write(LEAF.sub(lambda m: rename[m.group(1)], line) + "\n")
                n_out += 1
    print(f"wrote {n_out} trees to {args.out}; species with genes: {len(present)}")

    st = dendropy.Tree.get(path=args.species_tree, schema="newick", preserve_underscores=True, rooting="force-rooted")
    if args.prune:
        st.retain_taxa_with_labels(present)
        st.suppress_unifurcations()
    leaves = {l.taxon.label for l in st.leaf_node_iter()}
    absent = leaves - present
    for nd in st.preorder_node_iter():
        if not nd.is_leaf():
            nd.label = None
        nd.edge.length = 1.0
    st.write(path=args.species_out, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True,
             unquoted_underscores=True)
    print(f"wrote {args.species_out}: {len(leaves)} leaves ({len(absent)} without genes in this sample)")


if __name__ == "__main__":
    main()
