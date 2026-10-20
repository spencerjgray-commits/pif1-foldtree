#!/bin/bash
# workflow/35 (ML D/L rates) on every local-swap tree, v2fix and NCBI(+W. ciferrii) backbones; 3 at a time
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel 2>/dev/null || (cd "$HERE/../../../.." && pwd))
cd "$REPO"
PY=${PY:-python}   # python from the pif1 conda env (ete3 3.1.3, dendropy 5.0.8)
L=results/rerun_2026-10/species_tree_aware/local_swap; ST=results/rerun_2026-10/species_tree_aware/species_trees
n=0
for t in $L/trees/*.nwk; do b=$(basename $t .nwk)
  for st in v2fix ncbi_plus; do S=$ST/$st.nwk; [ $st = v2fix ] && S=$ST/v2fix/grafted_species_v2.binary.nwk
    $PY workflow/35_undated_dl.py --gene-tree $t --species-tree $S --map $ST/gene_species_957.map --optimize-rates --quiet --outdir $L/recon/${b}_$st > /dev/null 2>&1 &
    n=$((n+1)); [ $((n % 3)) -eq 0 ] && wait
  done
done
wait; echo done
