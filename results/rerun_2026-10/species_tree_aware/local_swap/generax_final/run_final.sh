#!/bin/bash
# Local-swap test on the final (radius-5) GeneRax trees, as in ../generax_r1:
#   IQ-TREE -z (LG+I+G4, reduced 408-gene alignment) + workflow/35 (ML D/L) on reduced v2fix.
# Two deep references for the todeep move (workflow/49 --deep):
#   .        results/seq_tree/pif1.treefile        (published deep AA tree; Pif1 sister to 855 tips -> node_0)
#   r1ref/   ../generax_r1/gr_r1_aa_deep.nwk       (GeneRax radius-1 deep tree, the reference used in generax_r1)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel 2>/dev/null || (cd "$HERE/../../../../.." && pwd))
cd "$REPO"
PY=${PY:-python}   # python from the pif1 conda env (ete3 3.1.3, dendropy 5.0.8)
IQ=tools/iqtree-3.1.2-macOS-arm/bin/iqtree3
B=results/rerun_2026-10/species_tree_aware; F=$B/local_swap/generax_final
for d in $F $F/r1ref; do
  $IQ -s $B/reduced/aln.fasta -z $d/aa_trees.nwk -m LG+I+G4 -n 0 -zb 10000 -au -zw -T 1 -seed 42 -pre $d/au_aa > $d/au_aa.stdout 2>&1 &
done
n=0
for d in $F $F/r1ref; do
  for t in $d/trees/*.nwk; do b=$(basename $t .nwk)
    [ -f $d/recon/${b}_v2fix/summary.json ] && continue
    $PY workflow/35_undated_dl.py --gene-tree $t --species-tree $B/reduced/v2fix.species.nwk --map $B/reduced/gene_species.map \
        --optimize-rates --quiet --outdir $d/recon/${b}_v2fix > /dev/null 2>&1 &
    n=$((n+1)); [ $((n % 3)) -eq 0 ] && wait
  done
done
wait
for d in $F $F/r1ref; do $PY workflow/50_local_swap_summary.py --dir $d --sets aa --species-trees v2fix; done
echo done
