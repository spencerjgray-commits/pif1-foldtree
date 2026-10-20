#!/bin/bash
# fixed-gene-tree baseline: workflow/35 (validated UndatedDL re-implementation), ML D/L rates, 957-gene map
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel 2>/dev/null || (cd "$HERE/../../../.." && pwd))
cd "$REPO"
PY=${PY:-python}   # python from the pif1 conda env (ete3 3.1.3, dendropy 5.0.8)
run() { $PY workflow/35_undated_dl.py --gene-tree $1 --species-tree $2 --map results/rerun_2026-10/species_tree_aware/species_trees/gene_species_957.map --optimize-rates --quiet --outdir results/rerun_2026-10/species_tree_aware/fixed_tree_baseline/$3 > results/rerun_2026-10/species_tree_aware/fixed_tree_baseline/$3.log 2>&1; }
for st in v2fix ncbi_plus; do
  S=results/rerun_2026-10/species_tree_aware/species_trees/$st.nwk; [ $st = v2fix ] && S=results/rerun_2026-10/species_tree_aware/species_trees/v2fix/grafted_species_v2.binary.nwk
  run results/seq_tree/pif1.treefile $S aa_LG_deep_$st &
  run results/rerun_2026-10/au/aa_constr197_s12345.treefile $S aa_LG_197_$st &
  run results/rerun_2026-10/ufboot_recon/aa3di_full_rep_mac.treefile $S aa3di_rep_deep_$st &
  run results/seq_tree/pif1_aa3di_supported.treefile $S aa3di_supported_197_$st &
  wait
done
