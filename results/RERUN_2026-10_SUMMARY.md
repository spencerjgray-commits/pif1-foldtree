# Re-run of the analyses flagged by the 2026-10-15 manuscript review

*Produced 2026-10-15/16; last updated 2026-10-16 05:30 (PMSF status updated 2026-10-20). Section codes (§1.1, MET-M1, RES-R3, …) refer to an internal manuscript review (2026-10-15) that is not part of this repository.
New outputs are all under `results/rerun_2026-10/`, and new scripts are `workflow/33_*.py`–`workflow/43_*.py`. Existing outputs were not modified.*

**Software.**
- IQ-TREE 3.1.2 macOS ARM: `tools/iqtree-3.1.2-macOS-arm/bin/iqtree3`, written `$IQ` below.
- Python 3.11 in the `pif1` conda env: ete3 3.1.3 (with the local NCBI taxonomy DB of 2026-06-16) and dendropy 5.0.8.
- All commands run from the repo root.

**Status legend.** ✅ = final. ⏳ = still running when this file was last updated (see §9).

---

## Bottom line for the revision

1. **Topology tests.** No test rejects either placement.
   - AA data: p-AU ≥ 0.30 for every tree under both LG+I+G4 and the BIC-selected Q.yeast+R10.
   - AA+3Di data: the 197-clade tree has p-AU 0.39 against the best deep tree.
   - The draft's p = 0.033 is gone.
2. **The 197-gene clade is not a stable feature of the AA+3Di ML tree.**
   - A replicate full AA+3Di search finds a tree 32 lnL better than `pif1_aa3di_supported` that places the duplication deep. None of its 1,000 UFBoot trees contains the clade.
   - The clade's own support in `pif1_aa3di_supported` is 0/13.
   - Conversely, the best AA-only trees found (under both models) **do** contain the clade.
   - The likelihood surface at this node is flat and multi-peaked under every dataset. Search outcome, not data type, decides the topology.
3. **All three structure-robustness variants (R1–R3) place the duplication deep** when run as full UFBoot searches (MRCA 763/828, 712/750, 876/957). The draft's 50-iteration tree (`pif1_aa3di`), which all figures, reconciliations and dates used, is rejected by its own AA+3Di data (p-AU 0.002).
4. **Reconciliation when the clade is present:**
   - The duplication maps to the Saccharomycotina crown on NCBI.
   - On the rebuilt, synonym-corrected graft it maps to **crown Saccharomycotina after *Lipomyces* diverged** (D = 7, not 43). That agrees with the 326–383 Mya Shen node, which also excludes *Lipomyces*.
   - The D = 43 "largest in the tree" was a species-tree artifact.
5. **Bootstrap reconciliation:** the duplication maps inside Saccharomycotina in 1.5 % (AA), 34.5 % (LBA control) and 0 % (AA+3Di replicate) of UFBoot trees.

---

## 0. Where the reconciliation numbers come from (read this first)

**No number in this file comes from GeneRax.** Every reconciliation number (D counts, the node the anchor duplication maps to, D/SL totals, likelihoods) comes from `workflow/35_undated_dl.py`, a Python re-implementation of GeneRax 2.0.4 `--strategy EVAL -r UndatedDL --unrooted-gene-tree`.
- *Why:* the only GeneRax on this Mac is 2.1.3, an Intel build under Rosetta. It segfaults in `Routines::inferReconciliation` before writing any event counts, and the archived runs used 2.0.4.
- *Fixed rates (D = L = 0.2, the draft's setting), validated exactly.* The script was run on the archived inputs of all seven GeneRax 2.0.4 runs in the repo:
  - aa3di, aa3di_supp and lgig on the NCBI tree;
  - aa3di, R1, R2 and R3 on the v1 grafted tree.
  
  For every run, all 1,437 species nodes have identical (S+SL+leaf, D, SL) counts, and the totals are identical. Outputs are in `results/rerun_2026-10/udl_validation/<run>/`; the comparison was made with `workflow/36_compare_event_counts.py`. For example:

  | archived run | totals D / SL (2.0.4 → re-impl.) | D(node_448) | rows differing |
  |---|---|---|---|
  | aa3di / NCBI | 464/3226 → 464/3226 | 9 → 9 | 0 of 1437 |
  | aa3di_supp / NCBI | 477/3348 → 477/3348 | 8 → 8 | 0 |
  | lgig / NCBI | 472/3236 → 472/3236 | 8 → 8 | 0 |
  | aa3di / grafted v1 | 476/3618 → 476/3618 | 43 → 43 | 0 |
  | R1 / grafted v1 | 402/3245 → 402/3245 | 37 → 37 | 0 |
  | R2 / grafted v1 | 386/3352 → 386/3352 | 39 → 39 | 0 |
  | R3 / grafted v1 | 489/3920 → 489/3920 | 43 → 43 | 0 |

- *Optimized rates, a separate check.* The 2.0.4 validation covers only fixed rates on a fixed gene tree. For optimized rates I compared against GeneRax 2.1.3, which prints its ML rates and log-likelihood before it crashes (`results/rerun_2026-10/udl_validation/generax213_optrates_reference.tsv`).
  - **Likelihood function** (`workflow/42_check_ll_vs_generax213.py` → `udl_validation/ll_check_vs_generax213.tsv`). Evaluated at GeneRax 2.1.3's own rates, my log-likelihood equals 2.1.3's printed value to ≤ 0.005 log units for the three `pif1_aa3di` configurations. For the three `pif1_aa3di_supported` configurations it differs by a constant +0.43–0.44. I did not find the cause: the input topologies are identical, and the fixed-rate scenarios match 2.0.4 exactly for both trees.
  - **Optimizer:** mine is Nelder–Mead on log rates. It finds rates within 1–4 % of 2.1.3's and log-likelihoods 0.0–0.9 units higher on every configuration (flat surface; table in §4.1).
- **Implication for the paper:** either describe and deposit `35_undated_dl.py` (with this validation as a supplementary note), or rerun the final configuration in GeneRax 2.0.4 on the Linux host before submission. The fixed-rate validation means 2.0.4 will reproduce the fixed-rate numbers here.

---

## 1. Old vs new for every number the review flagged

| # | Quantity (review §) | Draft / old value | New value | Status |
|---|---|---|---|---|
| 1a | AU test, AA tree vs converged **205-gene** Saccharomycotina constraint (§1.1, MET-M1) | ΔlnL 142.78, p-AU 0.0328 (unconverged tree) | **ΔlnL 2.58, p-AU 0.476** (the review got 0.481) | ✅ |
| 1b | AU test, AA tree vs **197-gene AA+3Di duplication clade** constraint (CRIT-K1) | not run | The best constrained tree scores **higher** than the AA ML tree: the AA tree has ΔlnL = 2.50, p-AU = 0.482; the constrained tree has p-AU = 1 | ✅ |
| 1c | ML score under the 197 constraint, 3 seeds | — | −95,216.05 (seed 12345), −95,226.17 (42), −95,234.75 (2026); the unconstrained AA search scored −95,217.42. Raw scores are only approximately comparable (constrained runs keep the 40 identical sequences); use the AU-evaluated lnL | ✅ |
| 2a | Duplication-node (197-clade) support, `pif1_aa3di_supported` (§1.2, MET-M6) | "SH-aLRT 98 / UFBoot 29" | **0 / 13** | ✅ |
| 2b | Pif1 daughter (97 tips) | — | **98.2 / 29** (this is where the draft's value came from) | ✅ |
| 2c | Rrm3 daughter (100 tips) | — | **88 / 16** | ✅ |
| 2d | LBA-control 198-tip clade (RES-R7) | 0 / 39 (earlier internal note, not included) | **71.7 / 54**. The 0/39 belongs to its 95-tip Pif1 side; the Rrm3 side is 90.1/82 | ✅ |
| 3a | Main AA+3Di tree (§1.2) | `pif1_aa3di` (50-iteration seeded search) | `pif1_aa3di_supported` (full search) | ✅ |
| 3b | Anchor duplication mapping, NCBI and v1 graft | node_448 | **node_448 (unchanged)** on both | ✅ |
| 3c | D(node_448), NCBI backbone (§1.5) | 9 (tied 3rd) | **8, rank 9**. Top: node_578 29, node_183 22, node_582 12, node_63 11; inside Saccharomycotina node_515 has 9 | ✅ |
| 3d | D(node_448), v1 grafted backbone | 43 (top) | **47 (top)**; node_578 29, node_183 22 | ✅ |
| 3e | Family composition (RES-R11, RES-R14) | Pif1 98 genes / 9 families, Rrm3 99 / 8, "Dipodascaceae Pif1-only" | **Pif1 97 / 8, Rrm3 100 / 8.** Dipodascaceae appear only in Pif1 (4 genes) and Trichomonascaceae only in Rrm3 (5 genes) | ✅ |
| 3f | Dated node on the Shen 2018 chronogram (§1.4) | crown 326.0, stem 383.5 Mya | **326.0 / 383.5 Mya (unchanged)**. Same 197 genes; the node subtends 317/332 tips and excludes 9 *Lipomyces* + 6 Trigonopsidales. With the synonym-aware matching, 59 instead of 54 species match, giving the same node and ages | ✅ |
| 3g | Duplication-node-to-tip AA path length, Pif1 vs Rrm3 (RES-R8) | 0.877 vs 1.073 (ratio 1.22) | **0.896 vs 1.084 (ratio 1.21)**; permutation p still reported as 0.0000 (anti-conservative) | ✅ |
| 3h | Daughter stem branches | Rrm3 0.148, Pif1 0.011 | **Rrm3 0.082, Pif1 0.035** | ✅ |
| 3i | Mean terminal branches (RES-R15) | Pif1 0.150, Rrm3 0.191, outside 0.082 | **Pif1 0.142, Rrm3 0.198, outside 0.086**; non-Sacch single-copy 0.048, multi-copy 0.167 (old 0.046 / 0.163) | ✅ |
| 3j | Tajima relative-rate test (MET-M7) | χ² 0.35 (not reported) | χ² 0.35, n.s. (alignment-based, so it cannot change) | ✅ |
| 3k | **Reproducibility of the main AA+3Di tree** (new check) | `pif1_aa3di_supported` (Linux, full search) lnL −185,669.93, contains the 197 clade | A replicate full search here (same alignment, partition model, `-B 1000 -bnni`, seed 42; ARM build) finds a **better tree, lnL −185,637.52 (+32.4), that does not contain the 197 clade**: the anchor MRCA is 880 tips (79.5/56); Pif1 daughter 97 tips (99.2/86); the Rrm3 copies nest among other fungi. Anchor dup → node_0 on NCBI, v1 and v2. AU under the partition model: supported tree ΔlnL 34.1, **p-AU 0.392 (not rejected)**; the draft's `pif1_aa3di` ΔlnL 576.9, **p-AU 0.0022 (rejected)** | ✅ |
| 4a | R2 tip count (§2) | 748 | **750**. Two *Candidozyma* core-only models were 72.81 → **88.73**. A0A9W9M2N2 = 86.998, which rounds to "87.0", explaining the 749-vs-748 discrepancy | ✅ |
| 4b | R1 (AFDB-only) human-rooted MRCA | 765 tips, SH-aLRT 79.7 (no UFBoot); anchor → node_0 | **763 of 828 tips (deep), 0 / 21**; anchor → **node_0** on NCBI, v1 and v2 | ✅ |
| 4c | R2 (pLDDT ≥ 87) MRCA | 157 tips of 748, SH-aLRT 48.9; anchor → node_448 | **712 of 750 tips (deep), SH-aLRT/UFBoot 92.3 / 76**; Pif1 daughter (78 tips, all Sacch) 100/92; anchor → **node_0** on NCBI, v1 and v2. The change comes from the search, not the 2 added tips (see §5) | ✅ |
| 4d | R3 (GTR20 3Di) MRCA | 876 tips, SH-aLRT 91.8; anchor → node_0 | **876 tips (deep), 84.2 / 32**; anchor → **node_0** on NCBI, v1 and v2 | ✅ |
| 5a | AA model (MET-M4) | "LG+I+G4 selected by ModelFinder" (never run) | **Q.yeast+R10** by BIC (199,722.5; LG+I+G4 200,962.1, ΔBIC 1,240). R10 is the `-cmax` ceiling. The unconstrained Q.yeast+R10 search from a parsimony start gives a deep tree (lnL −94,554.21; anchor MRCA 956/957). But a 197-constrained Q.yeast search beats it by 21.8 lnL (AU-evaluated), and an unconstrained search started from that tree **keeps the 197 clade** (lnL −94,534.78). Under the best-fit model, **the best AA tree found contains the 197 clade**; no tree in the AU set is rejected (all p-AU ≥ 0.30) | ✅ |
| 5b | LG+C20 PMSF final tree (MET-M5) | unfinished at iteration 300; intermediate MRCA 916/917 | Resumed from the checkpoint. The best tree (lnL −93,386.09, found at iteration 236) did not improve through iteration 330, and the run finished on 10-16 at 16:52 (final lnL −93,385.89). **Anchor MRCA = all tips but hPIF1 (956 of 957 with identical sequences kept), deep; that is the trivial split, so it has no support value.** The 197 split is absent from the ML tree and present in 2 of 1,000 UFBoot trees. The Saccharomycotina Pif1 clade (94 tips) has SH-aLRT/UFBoot 97.1/92 and the Rrm3 clade (102 tips) 97.6/85, but they are not sisters | ✅ |
| 6a | Grafted species tree (§1.5, MET-M9, RES-R3) | 57 matched, 18 next to congeners, 29 at the subphylum root (27-way polytomy, arbitrary binarization) | **62 exact (6 via synonyms), 36 by genus, 6 by family, 0 at the root**; 9 small polytomies (sizes 3–5) | ✅ |
| 6b | D(node_448) on the rebuilt grafted tree | 43, "the largest in the tree" | **0.** The duplication now maps to **sac_1 = Saccharomycotina minus *Lipomyces*** (103 leaves): **D = 7, rank 7** (`pif1_aa3di`) / **rank 9** (`_supported`). Identical across 20 random polytomy resolutions | ✅ |
| 7 | UFBoot-sample reconciliation (§1.5) | not run | Duplication maps inside Saccharomycotina in **AA 3/200 (1.5 %)**, **LBA control 69/200 (34.5 %)** and **AA+3Di (replicate run) 0/200 (0 %)**. The 197 split is present in 0 of the replicate's 1,000 UFBoot trees | ✅ |

---

## 2. Task 1 – AU tests (§1.1, CRIT-K1)

**Constraint files.**
- 205-gene: `results/seq_tree/constraint_sacch.nwk` (unchanged; converged tree `results/seq_tree/aa_constrained.treefile`).
- 197-gene: `results/rerun_2026-10/au/constraint_clade197.nwk` = ((197 tips of the hPIF1-rooted ScPif1/ScRrm3 MRCA in `pif1_aa3di_supported`), all 760 other tips). The 197-tip list is `results/rerun_2026-10/supports/clade197_supported.txt`; it is the identical set in `pif1_aa3di`.

**Constrained searches** (to IQ-TREE's default convergence, 100 unsuccessful iterations):
```
$IQ -s results/seq_tree/aln.trim.fasta -g results/rerun_2026-10/au/constraint_clade197.nwk -m LG+I+G4 -T 1 -seed S \
    -pre results/rerun_2026-10/au/aa_constr197_sS          # S = 42, 12345, 2026
```
| seed | best lnL | wall time | 197-clade present |
|---|---|---|---|
| 12345 | **−95,216.054** | 3 h 05 m | yes |
| 42 | −95,226.167 | 6 h 07 m (machine slept) | yes |
| 2026 | −95,234.749 | 3 h 31 m | yes |

For reference, the unconstrained AA search (`pif1.log`) scored −95,217.417, and the 205-constraint search (`aa_constrained.log`) −95,220.378.
- *Caveat on raw scores:* constrained (`-g`) runs keep all 957 sequences, while unconstrained runs drop 40 identical ones, so raw search scores are only approximately comparable. The AU tables below evaluate every tree identically and are the authoritative comparison.
- **There, the best 197-constrained tree has a higher likelihood than the published unconstrained AA ML tree (by 2.50).** The unconstrained search had not found the ML tree: the 197 clade is not just compatible with the AA data, it is in the best LG+I+G4 tree found.

**AU tests.** All use `-m LG+I+G4 -n 0 -zb 10000 -au -zw -T 2 -seed 42`, the same flags as the draft and the review.
```
$IQ -s results/seq_tree/aln.trim.fasta -z results/rerun_2026-10/au/au_trees_205conv.nwk ... -pre results/rerun_2026-10/au/au_205conv
$IQ -s results/seq_tree/aln.trim.fasta -z results/rerun_2026-10/au/au_trees_197best.nwk ... -pre results/rerun_2026-10/au/au_197best
$IQ -s results/seq_tree/aln.trim.fasta -z results/rerun_2026-10/au/au_trees_5set.nwk   ... -pre results/rerun_2026-10/au/au_5set
```
| test | tree | lnL (re-optimised) | ΔlnL | p-AU |
|---|---|---|---|---|
| `au_205conv` | `pif1.treefile` (AA ML) | −95,220.662 | 0 | 1 |
| | `aa_constrained.treefile` (205 constraint, converged) | −95,223.240 | **2.578** | **0.476** |
| `au_197best` | `pif1.treefile` | −95,220.662 | **2.501** | **0.482** |
| | `aa_constr197_s12345.treefile` | −95,218.161 | 0 | 1 |
| `au_5set` | `pif1.treefile` | −95,220.662 | 2.501 | 0.499 |
| | `aa_constrained` (205) | −95,223.240 | 5.079 | 0.472 |
| | 197-constrained seed 12345 | −95,218.161 | 0 | 1 |
| | 197-constrained seed 42 | −95,228.333 | 10.172 | 0.435 |
| | 197-constrained seed 2026 | −95,238.245 | 20.084 | 0.392 |

Notes:
- The review reports p-AU 0.481 for the 205 test. I get 0.476 with the same seed; RELL resampling differs slightly between thread counts.
- The draft's p = 0.033 / ΔlnL 143 is superseded either way. **The AA alignment does not reject either constraint.** For the hypothesis the paper actually argues (the 197-gene clade), the AA data slightly *favour* it.

**Sensitivity to the selected model (Q.yeast+R10).** Same five trees, `-m Q.YEAST+R10` (`au_5set_QyeastR10`):

| tree | lnL | ΔlnL | p-AU |
|---|---|---|---|
| `pif1.treefile` | −94,585.452 | 0 | 1 |
| 205-constrained | −94,591.977 | 6.52 | 0.490 |
| 197-constrained s12345 | −94,593.172 | 7.72 | 0.510 |
| 197-constrained s42 | −94,607.590 | 22.14 | 0.377 |
| 197-constrained s2026 | −94,610.448 | 25.00 | 0.369 |

These trees were all optimised under LG+I+G4. The set below adds the unconstrained Q.yeast+R10 ML tree (`modelfinder/pif1_mfp7.treefile`) → `au_6set_QyeastR10`:

| tree | lnL | ΔlnL | p-AU |
|---|---|---|---|
| Q.yeast+R10 ML tree (`pif1_mfp7`) | −94,574.739 | 0 | 1 |
| `pif1.treefile` (LG+I+G4 ML) | −94,585.394 | 10.66 | 0.512 |
| 205-constrained | −94,591.932 | 17.19 | **0.462** |
| 197-constrained s12345 | −94,593.127 | 18.39 | **0.445** |
| 197-constrained s42 | −94,607.541 | 32.80 | 0.355 |
| 197-constrained s2026 | −94,610.406 | 35.67 | 0.338 |

Neither constraint is rejected under the selected model either.

**Constrained searches under Q.yeast+R10**, started from the LG-constrained trees:
```
$IQ -s results/seq_tree/aln.trim.fasta -g results/rerun_2026-10/au/constraint_clade197.nwk -t results/rerun_2026-10/au/aa_constr197_s12345.treefile -m Q.YEAST+R10 -T 2 -seed 12345 -pre results/rerun_2026-10/au/aa_constr197_QyeastR10   # search lnL −94,536.28
$IQ -s results/seq_tree/aln.trim.fasta -g results/seq_tree/constraint_sacch.nwk -t results/seq_tree/aa_constrained.treefile -m Q.YEAST+R10 -T 2 -seed 42 -pre results/rerun_2026-10/au/aa_constr205_QyeastR10                 # search lnL −94,555.61
```
Final AU under Q.yeast+R10 (`au_Qfinal`, same flags as above):

| tree | lnL | ΔlnL | p-AU |
|---|---|---|---|
| **197-constrained, Q.yeast+R10 search** | **−94,552.950** | **0** | 1 |
| 205-constrained, Q.yeast+R10 search | −94,573.926 | 20.98 | 0.453 |
| unconstrained Q.yeast+R10 ML tree (`pif1_mfp7`) | −94,574.784 | 21.83 | 0.451 |
| `pif1.treefile` (LG+I+G4 ML) | −94,585.444 | 32.50 | 0.345 |
| 197-constrained, LG search (s12345) | −94,593.170 | 40.22 | 0.304 |

**Under the best-fit model, the highest-likelihood AA tree found contains the 197-gene clade.** It beats the unconstrained Q.yeast+R10 search by 21.8 lnL under identical evaluation, but nothing is rejected.
- The AA likelihood surface around this node is flat and multi-peaked. Unconstrained searches from random or parsimony starts do not find the clade; searches that start in it find higher scores.
- An **unconstrained** Q.yeast+R10 search started from the 197-constrained tree **stays in the clade** and improves to lnL −94,534.78 (`-keep-ident`, so all 957 tips). Its anchor MRCA is exactly the 197 genes, and it reconciles to node_448 (NCBI, D = 8) / sac_1 (v2, D = 6). This is the best AA tree found under the best-fit model.
```
$IQ -s results/seq_tree/aln.trim.fasta -t results/rerun_2026-10/au/aa_constr197_QyeastR10.treefile -m Q.YEAST+R10 -keep-ident -T 3 -seed 7 -pre results/rerun_2026-10/au/aa_free_from197_QyeastR10
```

- The old `au_test.iqtree` result still stands: the AA data strongly reject the *whole* `pif1_aa3di` topology (ΔlnL 1,891). They do not reject the duplication clade by itself, so the conflict lies elsewhere in the tree.

---

## 3. Task 2 – supports by bipartition (§1.2, RES-R7)

`workflow/33_split_supports.py` turns every IQ-TREE label into an unrooted split at parse time. A "clade" under outgroup rooting is then looked up as the side of a split not containing the outgroup tip, so no label can move with a reroot. Outputs: `results/rerun_2026-10/supports/*.json`.
```
python workflow/33_split_supports.py --tree results/seq_tree/pif1_aa3di_supported.treefile --write-clade results/rerun_2026-10/supports/clade197_supported.txt
python workflow/33_split_supports.py --tree results/lba/lba_noEDF.treefile --root A0A074S981 --json results/rerun_2026-10/supports/lba_noEDF.json
```
| tree | node | tips | SH-aLRT / UFBoot |
|---|---|---|---|
| `pif1_aa3di_supported` (rooted on hPIF1) | duplication (197 clade) | 197 | **0 / 13** |
| | Pif1 daughter | 97 | **98.2 / 29** |
| | Rrm3 daughter | 100 | **88 / 16** |
| `lba_noEDF` (rooted on *Rhizoctonia* A0A074S981) | 198-tip clade (197 Sacch + *Austropuccinia*) | 198 | **71.7 / 54** |
| | Pif1 side (incl. *Austropuccinia*) | 95 | 0 / 39 |
| | Rrm3 side | 103 | 90.1 / 82 |
| `pif1_aa3di` (50-iter; labels inherited from the AA tree, **not valid supports**) | 197 clade | 197 | (98.4/75: stale label) |
| `pif1.treefile` (AA) | anchor MRCA | 950 | 28.3 / 19 |

---

## 4. Task 3 – main AA+3Di tree = `pif1_aa3di_supported` (§1.2, CRIT-K5)

**Branch lengths.** The review asked for `-keep-ident`, so identical sequences stay where the input topology puts them.
```
$IQ -s results/seq_tree/aln.trim.fasta -te results/seq_tree/pif1_aa3di_supported.treefile -m LG+I+G4 -keep-ident -T 1 -seed 42 \
    -pre results/rerun_2026-10/main_supported/aa_bl_supported            # lnL −97,704.43
python workflow/15_rate_check.py --tree results/rerun_2026-10/main_supported/aa_bl_supported.treefile   # -> rate_check_aa_bl_supported.txt
python workflow/39_branch_summary.py results/rerun_2026-10/main_supported/aa_bl_supported.treefile        # -> branch_summary_old_vs_supported.txt
```
The numbers are in table rows 3g–3j. The old tree's values reproduce the draft and review exactly (`rate_check_old_aa_bl.txt`, `branch_summary_old_vs_supported.txt`). On the supported topology, the Rrm3-stem share of the path-length difference falls from about 70 % to about 25 % (stems 0.082 vs 0.035). The mean path-length ratio is unchanged (1.21).

**Reconciliation**, fixed D = L = 0.2:
```
python workflow/35_undated_dl.py --gene-tree results/seq_tree/pif1_aa3di_supported.treefile \
    --species-tree results/reconciliation/aa3di/eval/species_trees/inferred_species_tree.newick --outdir results/rerun_2026-10/main_supported/recon_ncbi
python workflow/35_undated_dl.py --gene-tree results/seq_tree/pif1_aa3di_supported.treefile \
    --species-tree results/reconciliation/aa3di_grafted/run/species_trees/starting_species_tree.newick --outdir results/rerun_2026-10/main_supported/recon_grafted_old
```
| backbone | anchor dup → | D(node_448) | node_578 | node_183 | totals D / SL |
|---|---|---|---|---|---|
| NCBI, old tree | node_448 (104) | 9 | 23 | 16 | 464 / 3226 |
| NCBI, **supported** | node_448 (104) | **8** (rank 9; node_515 = 9) | 29 | 22 | 477 / 3348 |
| v1 graft, old tree | node_448 | 43 | 23 | 16 | 476 / 3618 |
| v1 graft, **supported** | node_448 | **47** | 29 | 22 | 492 / 3803 |

(`results/rerun_2026-10/grafted/recon_main_trees.tsv` has all six tree × backbone combinations, including v2.)

**Shen chronogram placement:**
```
python workflow/17_place_on_timetree.py --timetree data/species_tree/shen2018_timetree.newick --gene-tree results/seq_tree/pif1_aa3di_supported.treefile
python workflow/41_shen_dating_synonyms.py --gene-tree results/seq_tree/pif1_aa3di_supported.treefile
```
Output is byte-identical to the old tree (`main_supported/shen_placement_*.txt`): 326.0 / 383.5 Mya, and the Shen root is 404.2 Mya. The synonym-aware version (59 exact matches) gives the same node.

### 4.1 Optimized D/L rates (review MET-M10 asks for sensitivity to the fixed rates)
```
python workflow/35_undated_dl.py --gene-tree G --species-tree S --optimize-rates --outdir results/rerun_2026-10/optrates/<tree>_<backbone>
```
| gene tree / backbone | ML D, L (re-impl.) | ll | GeneRax 2.1.3 D, L, ll | anchor dup → | D there | totals D / SL |
|---|---|---|---|---|---|---|
| aa3di50 / NCBI | 0.410, 1.084 | −6073.43 | 0.404, 1.081, −6073.48 | node_448 (104) | 9 | 464 / 3226 |
| supported / NCBI | 0.425, 1.110 | −6227.21 | 0.408, 1.085, −6227.99 | node_448 (104) | 8 | 477 / 3348 |
| aa3di50 / graft v1 | 0.416, 1.116 | −6619.98 | 0.407, 1.084, −6620.82 | **node_3 = Saccharomyceta (549)** | node_448 still has D = 42 | 476 / 3619 |
| supported / graft v1 | 0.433, 1.146 | −6866.39 | 0.416, 1.116, −6867.31 | **node_3 (549)** | node_448 D = 46 | 492 / 3804 |
| aa3di50 / graft v2 | 0.405, 1.069 | −5872.39 | 0.407, 1.073, −5872.39 | sac_1 (103) | 7 | 459 / 3084 |
| supported / graft v2 | 0.423, 1.099 | −6057.96 | 0.408, 1.085, −6058.58 | sac_1 (103) | 7 | 474 / 3223 |

- The ML scenario at the optimized rates is the same as at 0.2/0.2 on NCBI and on v2.
- On the v1 graft, the higher loss rate moves the anchor duplication up to Saccharomyceta. That is another symptom of v1's artificial basal resolution.
- The reconciliation likelihood ranks the backbones **v2 ≫ NCBI ≫ v1**: −6058 vs −6227 vs −6866 for the supported tree, the same ordering in both implementations.

---

## 5. Task 4 – R2 slicing fix and full R1–R3 searches (§1.3, §2, CRIT-K10)

**Core pLDDT fixed** (`workflow/34_core_plddt_fixed.py` → `results/rerun_2026-10/robustness/core_plddt_fixed.tsv`, `aln_plddt87_fixed.fasta`):
- AFDB values are recomputed from the local `.cif` CA B-factors over core_from..core_to.
- ColabFold full-length models are sliced the same way. The 3 core-only models (`data/structures/af3/_cores_note.txt`) are averaged over the whole model.
- Only 2 of 957 values change: A0ABX8I2X5 and A0A2V1AIT0, 72.81 → 88.73. The R2 set is 750 tips.
- Ten ColabFold models are 1–9 residues shorter than the UniProt length. They are sliced by residue number exactly as before, and their values are unchanged.

**Searches** (same partition files as the original runs; `-B 1000 -bnni` added as the review asked):
```
$IQ -s results/reviewer/aln_afdb.fasta                       -p results/seq_tree/aa3di_partition.nex -B 1000 -bnni -alrt 1000 -T 1 -seed 42 -pre results/rerun_2026-10/robustness/R1_afdb_B1000
$IQ -s results/rerun_2026-10/robustness/aln_plddt87_fixed.fasta -p results/seq_tree/aa3di_partition.nex -B 1000 -bnni -alrt 1000 -T 1 -seed 42 -pre results/rerun_2026-10/robustness/R2_plddt87fixed_B1000
$IQ -s results/seq_tree/alnAA_3di.fasta                      -p results/reviewer/part_gtr20.nex      -B 1000 -bnni -alrt 1000 -T 1 -seed 42 -pre results/rerun_2026-10/robustness/R3_gtr20_B1000
```
Results (human-rooted ScPif1/ScRrm3 MRCA by bipartition, `workflow/33`; reconciliation `workflow/38` → `robustness/recon_R*_B1000.tsv`, `recon_R2fixed.tsv`):

| run | tips | final lnL | anchor MRCA (tips, SH-aLRT/UFBoot) | anchor dup → (NCBI / v1 / v2) | archived run (alrt only) |
|---|---|---|---|---|---|
| R1 AFDB-only | 828 | −163,924.06 | **763 (0 / 21)** | node_0 / node_0 / node_0 | 765 tips (79.7), node_0; lnL −163,945.2 |
| R2 core pLDDT ≥ 87 (fixed) | **750** | −142,993.90 | **712 (92.3 / 76)**; Pif1 daughter 78 tips, all Sacch (100/92); Rrm3 side 634 tips | node_0 / node_0 / node_0 | 157 of 748 tips (48.9), node_448; lnL −142,729.8 |
| R3 GTR20 3Di | 957 | −192,430.43 | **876 (84.2 / 32)** | node_0 / node_0 / node_0 | 876 (91.8), node_0; lnL −192,393.2 |

**Why R2 changed.** I pruned the two added *Candidozyma* tips from the new R2 tree and compared it with the archived R2 tree on the archived 748-tip alignment (`robustness/au_R2_748_old_vs_newpruned.iqtree`; same partition model, `-n 0 -zb 10000 -au`).
- The new, deep topology scores **15.9 lnL higher** than the archived tree, whose 157-tip clade gets p-AU 0.439.
- So the archived R2 result came from a worse local optimum of an alrt-only search, not from excluding the two models. Neither topology is rejected.

**Tally.** With full UFBoot searches, **R1, R2 and R3 all place the duplication deep**: the anchor duplication maps to node_0 (all Fungi) on every backbone. Of the full-search AA+3Di runs, only `pif1_aa3di_supported` recovers the 197-gene clade, at 0/13 support. A replicate full search of the same data finds a better-scoring tree without it (§8). The new R3 search ended 37 lnL below the archived R3 tree (−192,430 vs −192,393); both are deep. D(node_448) on the v1 graft stays high in every run (R2 37, R3 42) even when the anchor maps to node_0, so it does not diagnose where the duplication maps.

---

## 6. Task 5 – ModelFinder and PMSF (MET-M4, MET-M5)

**ModelFinder.**
- IQ-TREE 3's default protein set is 1,232 models; at ~4 min per model on this machine that is about 80 h. That run was stopped after 9 models: `modelfinder/aborted_default_set/`. Its partial BIC ranking already had LG+R6 (200,353.4) ahead of LG+I+G4 (200,962.1).
- Run instead (182 models; 7 matrices that include the yeast- and Pfam-trained Q matrices):
```
$IQ -s results/seq_tree/aln.trim.fasta -m MFP -mset LG,WAG,JTT,Q.pfam,Q.yeast,Q.insect,VT -mfreq FU,F -mrate E,I,G,I+G,R -cmax 10 \
    -alrt 1000 -T 2 -seed 42 -pre results/rerun_2026-10/modelfinder/pif1_mfp7
```
Result (`modelfinder/pif1_mfp7.iqtree`, `.log`): **Q.yeast+R10 is selected by BIC (199,722.5)**. Runners-up: Q.yeast+R9 199,737.8, Q.insect+R10 199,738.6. LG+R10 is 200,105.4 and **LG+I+G4 200,962.1 (ΔBIC 1,240)**.
- R10 is the upper bound allowed by `-cmax 10`, and R9 is only 15 BIC units worse, so further categories would add little.
- The MFP tree search under Q.yeast+R10 (with `-alrt 1000`, no UFBoot) gives lnL −94,554.21. The hPIF1-rooted anchor MRCA is **956 of 957 tips (everything but hPIF1)**, so the AA placement stays deep under the better-fitting model. The 197 clade is not a split in this tree.
- AU tests under Q.yeast+R10 are in §2.

**PMSF**, resumed from the archived checkpoint in a mirror directory so the archived files are untouched. The command line is identical to the original minus `-redo`, with `-T 2`. The run restored the model parameters, the 1,000 UFBoot trees and the candidate set, and resumed at iteration 302 (`pmsf/work/results/seq_tree/`, run from `results/rerun_2026-10/pmsf/work/`):
```
iqtree3 -s results/seq_tree/aln.trim.fasta -m LG+C20+F+G -ft results/seq_tree/pif1.treefile -B 1000 -alrt 1000 -bnni -T 2 -seed 12345 -pre results/seq_tree/pif1_pmsf
```
Status at 05:25 on 10-16:
- Iterations 302→330 ran without improving on the best tree from iteration 236 (lnL −93,386.090). IQ-TREE's default stopping rule is 100 unsuccessful iterations; it is now in a UFBoot-refinement pass, about 46 % done after ~10 h, mostly asleep or contended.
- The best tree (`results/rerun_2026-10/pmsf/work/results/seq_tree/pif1_pmsf.treefile`), read by bipartition: **anchor MRCA 916 of 917 tips (all but hPIF1)**. The 197 clade is absent.
- **Finished 10-16 at 16:52** (final lnL −93,385.89; `pif1_pmsf.iqtree`/`.contree` in the same folder). Read by
  bipartition on the final tree: the anchor MRCA is all tips but hPIF1 (deep; the trivial split, so no support value).
  The 197 split is absent from the ML tree and present in 2 of the 1,000 UFBoot trees. The Saccharomycotina Pif1
  clade (94 tips) is 97.1/92 and the Rrm3 clade (102 tips) 97.6/85 (SH-aLRT/UFBoot), but they are not sisters.

---

## 7. Task 6 – rebuilt grafted species tree (§1.5, MET-M9, RES-R3)
```
python workflow/37_graft_species_tree_v2.py          # -> results/rerun_2026-10/grafted/  (NOTE: grafted/ was made by the pre-fix
                                                      # version; workflow/37 now has the genus-fallback fix of species_tree_aware/SUMMARY.md §1.1,
                                                      # so pass --outdir elsewhere to avoid overwriting grafted/)
python workflow/38_recon_species_trees.py --gene-trees results/seq_tree/pif1_aa3di.treefile results/seq_tree/pif1_aa3di_supported.treefile \
   --species-trees <NCBI> <v1> results/rerun_2026-10/grafted/grafted_species_v2.binary.nwk --out results/rerun_2026-10/grafted/recon_main_trees.tsv
python workflow/38_recon_species_trees.py --gene-trees <same two> --species-trees results/rerun_2026-10/grafted/grafted_species_v2.rand*.nwk \
   --out results/rerun_2026-10/grafted/recon_v2_random_resolutions.tsv
```
**Name matching.** Every Shen tip and every sampled taxid is mapped to its current NCBI species/genus/family; NCBI resolves the synonyms. One Shen label NCBI doesn't resolve is mapped by hand (`Candida_parapsilosis` → *Lodderomyces parapsilosis*). The resulting table is `grafted/shen2018_synonyms.tsv`. The six synonym hits are the six the review named: *Candidozyma auris*, *Sungouiella intermedia*, *Dekkera bruxellensis*, *Huiozyma naganishii*, *Henningerozyma blattae*, *Nakaseomyces bracarensis*.

**Placement** on the full 332-tip Shen tree (exact tip / sister to the genus MRCA / sister to the family MRCA), then pruned to the 104 sampled taxa:
- 62 exact, 36 genus, 6 family (*Diutina rugosa*, *Insectozyma verbasci*, *Australozyma saopauloensis*, *Waltiozyma* ×2, *Saccharomycodes ludwigii*), 0 at the root.
- 9 residual polytomies (sizes 3,3,3,3,3,3,4,4,5), all within a genus or family.
- The new subtree replaces node_448 in the v1 binary tree, so every non-Saccharomycotina node and label is unchanged. The new subtree's root keeps the name node_448, and its 103-leaf child is sac_1.
- *Lipomyces orientalis* is now sister to all other sampled budding yeasts, as in Shen et al.

**Result.**
- **The anchor duplication maps to sac_1 (Saccharomycotina minus *Lipomyces*), D = 7, and D(node_448) = 0.** This holds for both AA+3Di trees, on the deterministic resolution and on all 20 random resolutions of the 9 polytomies. Under those, total D varies 455–462 / 470–477, and the rank of sac_1 is 7 (`pif1_aa3di`) or 9–10 (`_supported`).
- The top nodes are now node_578 (Agaricomycetes, 23/29) and node_183 (16/22).
- The rebuilt graft also needs far fewer losses (SL 3084 vs 3618 for `pif1_aa3di`) and has a much higher reconciliation likelihood (§4.1).
- **The v1 D = 43 "largest in the tree" was a species-tree artifact.** The mapping to sac_1 also resolves the §1.4 conflict: the reconciliation and the dated Shen node now agree on "crown Saccharomycotina after *Lipomyces* diverged".

---

## 8. Task 7 – reconciling bootstrap tree samples (§1.5)
```
python workflow/40_recon_ufboot.py --trees results/seq_tree/pif1.ufboot --n 200 --species-trees <NCBI> <v2> --labels ncbi graftv2 \
   --out results/rerun_2026-10/ufboot_recon/aa_LGIG4_ufboot200.tsv
python workflow/40_recon_ufboot.py --trees results/lba/lba_noEDF.ufboot --n 200 --root A0A074S981 --species-trees <NCBI> <v2> --labels ncbi graftv2 \
   --out results/rerun_2026-10/ufboot_recon/lba_noEDF_ufboot200.tsv
```
(200 evenly spaced trees of the 1,000; fixed D = L = 0.2.)

| UFBoot sample | dup maps inside Saccharomycotina (NCBI / v2) | other mappings |
|---|---|---|
| AA, LG+I+G4 (`pif1.ufboot`) | 3/200 (1.5 %) / 3/200 | node_0 all Fungi 196, node_1 Dikarya 1 |
| LBA control (`lba_noEDF.ufboot`) | 69/200 (34.5 %) / 69/200 | node_1 Dikarya 131 (the clade picks up *Austropuccinia*) |
| AA+3Di, replicate full search (`aa3di_full_rep_mac.ufboot`) | **0/200 (0 %) / 0/200** | node_0 198, node_1 1, node_2 1 |

```
python workflow/40_recon_ufboot.py --trees results/rerun_2026-10/ufboot_recon/aa3di_full_rep_mac.ufboot --n 200 --species-trees <NCBI> <v2> --labels ncbi graftv2 \
   --out results/rerun_2026-10/ufboot_recon/aa3di_rep_ufboot200.tsv
python workflow/43_count_split_in_trees.py results/rerun_2026-10/ufboot_recon/aa3di_full_rep_mac.ufboot results/rerun_2026-10/supports/clade197_supported.txt
```
Read directly from all 1,000 replicate UFBoot trees, the 197-gene split is present in **0**, and the anchor MRCA is all-Saccharomycotina in 4 (`ufboot_recon/aa3di_rep_ufboot1000_clade197_count.txt`).

The `pif1_aa3di_supported` run wrote its `.ufboot` on the Linux host only, so it isn't in this repo. A replicate full AA+3Di search with `-wbt` was therefore run here:
```
$IQ -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition.nex -B 1000 -bnni -alrt 1000 -wbt -T 1 -seed 42 \
    -pre results/rerun_2026-10/ufboot_recon/aa3di_full_rep_mac
```
**Replicate result.** Best lnL −185,637.52 after 14 h 43 m, which is **32.4 better than the Linux run** (−185,669.93), and it **does not recover the 197 clade**. AU test on the AA+3Di data:
```
$IQ -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition.nex -z results/rerun_2026-10/ufboot_recon/au_aa3di_trees.nwk -n 0 -zb 10000 -au -zw -T 2 -seed 42 \
    -pre results/rerun_2026-10/ufboot_recon/au_aa3di_supported_vs_rep
```
| tree | lnL | ΔlnL | p-AU |
|---|---|---|---|
| replicate (Mac) ML tree, deep placement | −185,649.805 | 0 | 1 |
| `pif1_aa3di_supported` (Linux), 197 clade | −185,683.903 | 34.10 | 0.392 |
| `pif1_aa3di` (draft's tree: 50-iter, AA-seeded) | −186,226.698 | 576.89 | **0.0022** |

If `results/seq_tree/pif1_aa3di_supported.ufboot` can be copied from the Linux host, run `40_recon_ufboot.py` on it as well. It will give the UFBoot picture around the 197-clade optimum (the split itself has UFBoot 13).

---

## 9. Running jobs / caveats / follow-ups

- **Machine sleep.** From about 14:15 on 10-15 to about 00:00 on 10-16 the Mac was lid-closed on battery, in maintenance sleep (about 8 s awake every 2–6 min). Jobs barely advanced. IQ-TREE wall times include this sleep and are not CPU benchmarks.
- PMSF finished on 10-16 at 16:52; its result is in row 5b. Nothing is still running.
- Not in scope, still open from the review:
  - the 16 genes (957 → 941) dropped from every reconciliation (MET-M8);
  - the 5 taxa with all-gap 3Di rows (MET-M15/RES-R20);
  - the FoldTree alntmscore/lddt trees (MET-M2).
- ModelFinder selected Q.yeast+R10, and the AU tests and constrained searches were repeated under it (§2).
- Suggested next steps:
  - Several more independent full AA+3Di searches (different seeds and starting trees, including one started from `pif1_aa3di_supported`), to find the actual AA+3Di ML tree.
  - A GeneRax 2.0.4 run of the final configuration on the Linux host, so the paper can cite GeneRax directly (§0).
