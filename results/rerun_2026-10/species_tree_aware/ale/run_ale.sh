#!/bin/bash
# usage: run_ale.sh <dir> <ale species tree, relative to ale/> <named species tree for workflow/46> [delta lambda]
#   <dir>/<dir>.treelist must exist (workflow/45). With delta and lambda given, the rates are fixed (one
#   likelihood evaluation); otherwise ALEml_undated optimises them (Nelder-Mead).
# ALE = ssolo/ALE @ b59dfaf3 (2025-05-14) built x86_64 (run under Rosetta: 80-bit long double is needed for
# 957 genes), with one patched line in src/undated.cpp (sample_undated):
#     if (r*uq_sum<uq_resum or no_T)   ->   if (r*uq_sum<uq_resum)
# Unpatched, a DL-only run (tau=0) samples a duplication whenever the speciation draw fails, so the sampled
# reconciliations are biased toward duplications. The likelihood code is untouched (logL identical).
# tau=0: DL only. (A small positive tau is not equivalent: logL rises with tau, -1432.53 at 0, -1428.22 at 1e-9.)
set -e
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel 2>/dev/null || (cd "$HERE/../../../.." && pwd))
D=$REPO/results/rerun_2026-10/species_tree_aware/ale/$1
cd $D
AB=${ALE_BIN:-$REPO/tools/ale-b59dfaf-x86_64-patched}   # patched x86_64 build; see ../software/BUILD.md
PY=${PY:-python}   # python from the pif1 conda env (ete3 3.1.3, dendropy 5.0.8)
[ -s $1.treelist.ale ] || $AB/ALEobserve $1.treelist burnin=0 > /dev/null 2>&1
cp ../$2 ./S.nwk
RATES=""; [ -n "$4" ] && RATES="delta=$4 lambda=$5"
/usr/bin/time -p $AB/ALEml_undated S.nwk $1.treelist.ale sample=1000 separators=_ tau=0 $RATES seed=42 > ALEml_undated.raw 2>&1
tr '\r' '\n' < ALEml_undated.raw | grep -a -v -E '^\s*$|nodes loaded|^[-/\\|0-9 *%]+$' > ALEml_undated.log
mv S.nwk_$1.treelist.ale.uml_rec $1.uml_rec; mv S.nwk_$1.treelist.ale.uTs $1.uTs 2>/dev/null || true
cd "$REPO"
$PY workflow/46_ale_anchor_posterior.py --uml-rec $D/$1.uml_rec --species-tree $3 --label $1 --out $D/anchor_posterior.tsv > $D/anchor_posterior.txt
gzip -f $D/$1.uml_rec
