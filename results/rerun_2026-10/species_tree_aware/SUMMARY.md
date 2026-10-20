# Species-tree-aware test of the Pif1/Rrm3 duplication (ALE, GeneRax SPR, local swaps)

*2026-10-16. Question: once the species tree is taken into account, do the data favour the 197-gene
Saccharomycotina clade (anchor duplication → sac_1 = crown Saccharomycotina after Lipomyces) or the deep placement
(→ node_0, all Fungi)? All outputs are under `results/rerun_2026-10/species_tree_aware/`. New scripts are
`workflow/44_*.py`–`51_*.py`, and `workflow/37` got a one-block bug fix. Logs and raw outputs record the original
machine's paths, including temporary build directories; the binaries used are now in the git-ignored `tools/` (see `software/BUILD.md`).
Updated 2026-10-20 with the finished GeneRax SPR runs (§4, Bottom line item 4). Nothing is still running (§7).*

## Bottom line

**With the species tree included, every analysis that can see both topologies favours a Saccharomycotina-specific
duplication over a fungus-wide one. This holds for the AA data and for the AA+3Di gene-tree samples. The preference
comes from a small, consistent reconciliation term, about 4–6 log-likelihood units for the duplication placement
on the rebuilt graft. It decides the question because the sequences are close to indifferent at that node. Where
on the stem the duplication sits (crown Saccharomycotina after *Lipomyces*, or earlier with losses in
Pezizomycotina) depends on the species tree.**

1. **The reconciliation term, isolated.** Take any background tree and move only the Pif1 clade between the two
   positions (§5). Every time, the reconciliation favours the Saccharomycotina placement:
   - by **5.2–5.7 lnL** on the v2fix graft, **3.6–3.7** on NCBI;
   - in AA and AA+3Di trees, in deep and 197 backgrounds alike;
   - a deep duplication needs only a handful of extra losses (the ML loss rate is high, L ≈ 1.1), not hundreds, so
     the penalty is modest.
2. **AA data, joint likelihood: Saccharomycotina.**
   - The best LG+I+G4 tree found is a new search started inside the 197 clade (§3.1). It has lnL −95,198.8,
     18.6 above the published deep AA tree, and contains the exact 197-gene split.
   - In that tree, moving Pif1 to the deep position costs 4.8 (sequence) + 5.4 (reconciliation) = **10.2 joint lnL**.
   - All 197-clade AA trees beat the deep AA tree by **25–66 joint lnL** on v2fix (14–55 on NCBI).
   - The sequence data alone still do not reject the deep tree (p-AU 0.52).
3. **ALE posteriors (§3).**
   - Samples that contain both topologies:
     - the new from-197 UFBoot samples (AA and AA+3Di);
     - pooled deep + from-197 samples, where two-thirds (AA) to three-quarters (AA+3Di) of the input trees are deep;
     - the LBA control.
   - In those (408 and 957 genes), the duplicated clade is
     Saccharomycotina-only in **95–100 %** of sampled reconciliations.
   - On v2fix the duplication maps to **sac_1 in 56–77 %**, node_448 in 16–27 %, the Saccharomyceta or Ascomycota
     stem in 2–13 %, and anything deeper in ≤ 5 %.
   - On the NCBI backbone (pooled AA), the mapping is node_448 50 %, Saccharomyceta 34 %, Ascomycota 12 %,
     deeper 3 %.
   - Samples that contain only deep trees stay deep: 0–0.4 % Saccharomycotina.
   - In the pooled samples, prior odds of about 1:2 (AA) and 1:3 (AA+3Di) in favour of deep trees become
     ≥ 300:1 (AA) and ≈ 25–125:1 (AA+3Di) in favour of a Saccharomycotina-only clade.
   - Both results hold for 408 genes with ML rates and for 957 genes with fixed rates.
4. **GeneRax SPR (§4): all three starts, the deep one included, end with a Saccharomycotina duplication**
   (408-gene subset; the 957-gene search was intractable).
   - **The search from the deep tree crosses on its own.** After radius 1 it was still deep. In the final
     (radius-5) tree it has the exact 197-gene split, and the anchor duplication maps to **sac_1**, as in the two
     197 starts. In the radius-1 tree the two clades were 8 edges apart, more than the maximum SPR radius of 5,
     so the crossing took several moves.
   - **Final joint lnL:** deep start −62,144.0, 197 −62,144.4, 197-Q.yeast −62,141.0. The spread is 3.4; after
     radius 1 it was 75. The final trees are closer to one another (RF 42–60 of 810) than to their own radius-1
     trees.
   - **Local swap on the final trees:** moving Pif1 back to the fungus-wide position costs **15.3–16.6 joint
     lnL** (sequence 8.4–9.5, reconciliation 6.8–7.1) in all three. On sequence alone the swapped trees are still
     not rejected in the 197 backgrounds (p-AU 0.14).
5. **AA+3Di: the duplication node favours Saccharomycotina; the best whole tree found is still deep.**
   - The ALE results for AA+3Di:
     - from-197 sample: 100 % Saccharomycotina-only clade (sac_1 77 %);
     - pooled sample (only 26 % of input trees Saccharomycotina): 99.2 % (408 genes) and 96.2 % (957 genes);
     - deep-only sample: 0 %.
   - On the `_supported` background, the duplication-local swap favours the Saccharomycotina placement by 4.9 joint
     lnL (3.1 on NCBI).
   - The best AA+3Di tree found (the deep replicate, re-run §8) still beats the best 197-clade AA+3Di tree by 34 in
     sequence lnL and 58 in reconciliation lnL. That gap comes from the rest of the tree: the deep replicate splits
     the Saccharomycotina Rrm3 genes, so no single local move connects the two.
   - GeneRax cannot fit the AA+3Di partition model.
6. **How strong.**
   - For the AA data the support is **consistent but moderate**: a likelihood ratio of e^3–e^10 (≈ 20–25,000) for
     the duplication placement given the rest of the tree, and an ALE posterior near 1 once the gene-tree sample
     contains both topologies.
   - None of this is a significance test against sequence-level uncertainty. The deep AA tree is still not rejected
     by the alignment (p-AU 0.19–0.52), and the AA+3Di searches keep finding deep trees with better
     sequence likelihood.
   - **Recommended wording:** "species-tree-aware reconciliation favours a duplication early in crown
     Saccharomycotina (after Lipomycetales); sequence data alone cannot distinguish this from a fungus-wide
     duplication". This is a qualified version of the Saccharomycotina claim, not the original "structure resolves
     it" framing.

---

## 1. Inputs (fixed first)

### 1.1 Grafted species tree: genus-fallback bug fixed (`workflow/37`, run through `workflow/44`)
- **Bug.** Shen labels that NCBI does not translate fell back to the genus in the label. `Candida_corydali` and
  `Candida_carpophila` thereby joined the 4 translated *Candida* tips, so the "Candida" target became the MRCA of
  6 labels, an 88-tip clade. *C. maltosa* (1245528) and *C. viswanathii* (5486) were attached as its sister.
- **Fix** (in place in `workflow/37`, lines ~90–105): an untranslated label falls back to its genus only if no
  translated Shen tip represents that genus.
  - *C. maltosa* and *C. viswanathii* now sit sister to the 4-tip *C. albicans/dubliniensis/tropicalis/sojae* clade.
  - The fix also stops untranslated labels from widening four other targets (Debaryomycetaceae 56→52 tips,
    *Ogataea* 19→17, *Blastobotrys* 9→7, Wickerhamomycetaceae 24→23). The polytomies are otherwise the same.
  - The older outputs in `results/rerun_2026-10/grafted/` predate the fix and were left untouched.
- **Rebuilt tree:** `species_trees/v2fix/grafted_species_v2.binary.nwk`, 720 leaves.
  - Produced by `workflow/37 --v1-binary species_trees/v1_plus.nwk --outdir species_trees/v2fix`, called from
    `workflow/44`.
  - sac_1 = Saccharomycotina minus *Lipomyces*, now 104 leaves (with *W. ciferrii*); node_448 = 105 leaves.
- **Secondary tree:** NCBI, i.e. the GeneRax `inferred_species_tree` + *W. ciferrii*, as `species_trees/ncbi_plus.nwk`.

### 1.2 Gene → species map: all 957 genes (`workflow/44` → `species_trees/gene_species_957.map`, `map_changes.tsv`)
941 genes keep their mapping. The 16 that were dropped:
- **14 genes** carry a species-level taxid whose leaf is a strain of the same species. They are mapped to that leaf:
  - *Clonostachys rosea* 29856→1349383; *Akanthomyces lecanii* 2714763→1081108;
  - *Rhodotorula toruloides* 5286→1130832;
  - *Puccinia graminis* f. sp. *tritici* 56615→418459 (3 genes); *P. triticina* 208348→630390;
    *P. striiformis* 27350→1165861;
  - *Rhizopus delemar* 936053→246409;
  - *Rhizoctonia solani* 456999→1108050 (5 genes). *R. solani* has two leaves; the lower taxid was used.
- **2 genes:** *Wickerhamomyces ciferrii* K0KM86 and K0KZ61, both in the 197 clade.
  - They carry the merged taxid 1206466 → 1041607, a species with no leaf, so leaf 1041607 was added.
  - On NCBI it is sister to the 2 sampled Wickerhamomycetaceae.
  - On the graft it is placed through its own Shen tip (`Wickerhamomyces_ciferrii`, exact match).

### 1.3 Reduced gene set for GeneRax (`workflow/48` → `reduced/`)
- **Why reduce.** At 957 genes and 720 species, a GeneRax 2.0.4 SPR round evaluates every move and applies the single
  best one. That took ~15 min per round here, i.e. days per search (aborted runs: `generax/full957_aborted/`).
- **What is kept** (408 genes, 264 species):
  - **all 205 Saccharomycotina genes** (105 species, the whole 197 clade);
  - for each of the 74 other NCBI orders (class where a species has no order), up to **3 species drawn at random
    (seed 42), with all their genes**. That is 159 of 615 outgroup species and 203 of 752 outgroup genes.
    Per-order counts are in `reduced/composition.tsv`.
  - hPIF1 (*Homo sapiens*) is kept.
- **Trees.** Species trees are pruned to the kept species; node_0, node_448 and sac_1 survive. Starting gene trees
  are pruned to the kept genes and keep their placement:
  - deep start: anchor MRCA 406 of 408;
  - 197 starts: exactly the 197 clade.
- **ALE was run on both sets:** 408 genes with ML rates, and 957 genes at fixed rates.

---

## 2. Software (patches and the build recipe: `software/BUILD.md`)

| tool | version / build | validation and notes |
|---|---|---|
| GeneRax | **2.0.4** (e4fab40, 2021-08-03), built from source for **arm64**, upstream build type (asserts on) | Reproduces the archived 2.0.4 `--strategy EVAL` event counts **byte for byte**. On the new 957-gene/v2fix input it equals `workflow/35` on all 1,439 species nodes (`validation/compare_generax204_vs_35_aa_deep_v2fix.txt`). |
| ALE | ssolo/ALE master **b59dfaf** (2025-05-14), **x86_64** build run under Rosetta, **one-line sampler patch** | arm64 builds fail at 957 genes: on arm64 macOS `long double` = `double`, so ALEobserve overflows (tgamma) and likelihoods underflow. |
| IQ-TREE | 3.1.2 macOS arm (`tools/iqtree-3.1.2-macOS-arm`) | UFBoot searches; tree scoring with `-z -n 0 -zb 10000 -au` |
| Python | 3.11, `pif1` env: ete3 3.1.3, dendropy 5.0.8 | `workflow/35` (the validated UndatedDL re-implementation) for every fixed-tree reconciliation |

Software problems found and handled:
- **GeneRax Release builds cannot do SPR.** `src/core/search/Moves.cpp:103` and `Rollbacks.cpp:10` call
  `pllmod_utree_spr` / `pllmod_tree_rollback` *inside* `assert()`. With `-DNDEBUG` the moves are never applied and the
  SPR worker segfaults at the first move. This was reproduced:
  - with a Release 2.0.4 on arm64 and on x86_64;
  - with the bioconda 2.1.3 arm64 build (bus error at its first SPR move; in EVAL it also crashes in
    `Scenario::savePerSpeciesEventsCounts`, as the old Intel 2.1.3 did).

  Building without a build type, as GeneRax's own `install.sh` does, fixes it. This likely explains the Mac 2.1.3
  crashes. The Linux host's 2.0.4 builds worked.
- **GeneRax cannot fit AA+3Di.** One substitution model per family: `PLLTreeInfo::buildModel` asserts a single
  submodel, and only `partitions[0]` is used. The LG+G (AA) + 3Di-matrix+G (3Di) partition model cannot be run.
  AA+3Di is therefore covered by ALE and by the local-swap test, which uses IQ-TREE's partition model.
- **ALE's DL-only sampler is broken (τ = 0).** In `sample_undated`, `if (r*uq_sum<uq_resum or no_T)` returns a
  duplication whenever the first speciation draw fails. Unpatched, on the reduced AA sample with D=0.36, L=0.90:
  - all 1,000 samples were one identical scenario with 406 D and 5,395 L each;
  - patched: 120 D and 489 L per sample on average, at **identical logL** (−1432.53).

  See `validation/ale_sampler_patch_check.txt`. The patch removes `or no_T`; the likelihood code is untouched.
- **τ must be exactly 0 for a DL model.** ALE's logL rises with τ: −1432.53 at 0 (and 1e-12), −1428.22 at 1e-9,
  −1413.16 at 1e-6. Even tiny transfer rates absorb discordance.
- **Environment.** Two new conda envs hold the toolchains: `strecon` (arm64) and `strecon_x86` (osx-64). No existing
  env or global config was changed.
  - The PMSF job (PID 45461) finished on its own at 16:52 and was never touched.
  - The Mac slept lid-closed on battery from about 15:00 to 16:40, which stalled all jobs for that time.

---

## 3. ALE: posterior of the anchor-duplication mapping (`ale/`, `workflow/45`, `workflow/46`)

**Command per sample** (driver `ale/run_ale.sh`; outputs in `ale/<sample>/`: `anchor_posterior.txt/.tsv`,
`ALEml_undated.log`, `<sample>.uml_rec.gz`):
```
python workflow/45_ale_inputs.py --trees <ufboot files> [--genes reduced/genes.txt] --out ale/<s>/<s>.treelist \
    --species-tree <named species tree> --species-out ale/<ale species tree>
ALEobserve <s>.treelist burnin=0
ALEml_undated S.nwk <s>.treelist.ale sample=1000 separators=_ tau=0 [delta=D lambda=L] seed=42
python workflow/46_ale_anchor_posterior.py --uml-rec <s>.uml_rec --species-tree <named species tree> --out ...
```
- Tree lists and CCP files are stored gzipped (`*.treelist.gz`, `*.treelist.ale.gz`); gunzip them before
  re-running `run_ale.sh`.
- Gene labels are rewritten as `<species leaf>_<accession>`.
- ALE's species tree is the named tree with internal labels dropped and unit branch lengths. `workflow/46` maps ALE
  branch ids back to node names by leaf set.
- The anchor event is the MRCA of ScPif1 (P07271) and ScRrm3 (P38766) in each sampled reconciled gene tree.
- Node meanings:

  | node | contents (v2fix) |
  |---|---|
  | sac_1 | Saccharomycotina minus *Lipomyces* |
  | node_448 | crown Saccharomycotina |
  | node_3 | Saccharomyceta: Pezizomycotina + Saccharomycotina |
  | node_2 | Ascomycota minus Taphrinomycotina |
  | node_1 | Dikarya |
  | node_0 | Fungi |
  | root (node_717) | Fungi + human |

### 3.1 Tree samples
| sample | source | trees | raw trees with all-Sacch anchor MRCA (197 split) |
|---|---|---|---|
| AA deep | `results/seq_tree/pif1.ufboot` (LG+I+G4, published AA search) | 1000 | 7 (0) |
| AA from-197 | **new**: unconstrained LG+I+G4 UFBoot search started from `aa_constr197_s12345` (`ufboot197/aa_LG_from197.*`) | 1000 | 647 (546) |
| AA pooled | AA deep + AA from-197 | 2000 | 654 (546) |
| AA+3Di deep | `ufboot_recon/aa3di_full_rep_mac.ufboot` (replicate full search) | 1000 | 4 (0) |
| AA+3Di from-197 | **new**: partition-model UFBoot search started from `pif1_aa3di_supported` (`ufboot197/aa3di_from197.*`) | 1000 | 522 (435) |
| AA+3Di pooled | AA+3Di deep + AA+3Di from-197 | 2000 | 526 (435) |
| LBA control | `results/lba/lba_noEDF.ufboot` (no early-diverging fungi; species tree pruned, root = Dikarya) | 1000 | 305 (273) |

The new searches:
```
$IQ -s results/seq_tree/aln.trim.fasta -m LG+I+G4 -t results/rerun_2026-10/au/aa_constr197_s12345.treefile \
    -B 1000 -bnni -wbt -T 2 -seed 12345 -pre results/rerun_2026-10/species_tree_aware/ufboot197/aa_LG_from197
$IQ -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition.nex -t results/seq_tree/pif1_aa3di_supported.treefile \
    -B 1000 -bnni -wbt -T 2 -seed 42 -pre results/rerun_2026-10/species_tree_aware/ufboot197/aa3di_from197
```
- **AA from-197:** best lnL **−95,198.81**, 18.6 above the published AA ML search (−95,217.42). Its tree contains the
  exact 197-gene split, with UFBoot 55. **This is now the best AA LG+I+G4 tree found.**
- **AA+3Di from-197:** best lnL −185,667.95, no improvement on its start in 100 iterations. That is 30.4 below the
  deep replicate (−185,637.52).

### 3.2 Results
**Reduced gene set (408 genes, 264 species), D and L optimised by ALE, v2fix:**

| sample | ALE D / L | logL | Sacch-only anchor gene clade | sac_1 | node_448 | node_3 | node_2 | deep (node_1 / node_0 / root) |
|---|---|---|---|---|---|---|---|---|
| AA deep | 0.174 / 0.501 | −1379.78 | 0.4 % | 0.2 % | 0.1 % | 0.1 % | – | **99.6 %** (– / 65.7 / 33.9) |
| AA from-197 | 0.173 / 0.485 | −1356.32 | 99.5 % | **66.5 %** | 23.4 % | 7.1 % | 2.2 % | 0.8 % |
| **AA pooled** | 0.161 / 0.465 | −1332.33 | **100 %** | **64.8 %** | 24.5 % | 8.8 % | 1.7 % | 0.2 % (0.2 / 0 / 0) |
| AA pooled, **NCBI tree** | 0.193 / 0.567 | −1515.27 | 99.0 % | (no sac_1 on NCBI) | **50.4 %** | **34.2 %** | 12.1 % | 3.3 % (2.2 / 0.6 / 0.5) |
| AA+3Di deep | 0.163 / 0.482 | −1380.00 | 0 % | 0 | 0 | 0 | 0 | **100 %** (– / 81.2 / 18.8) |
| AA+3Di from-197 | 0.172 / 0.508 | −1410.38 | **100 %** | **76.9 %** | 19.3 % | 3.1 % | 0.5 % | 0.2 % (0.2 / 0 / 0) |
| **AA+3Di pooled** | 0.152 / 0.450 | −1360.10 | **99.2 %** | **62.9 %** | 24.8 % | 8.7 % | 2.1 % | 1.5 % (0.7 / 0.7 / 0.1) |
| LBA control (244 species) | 0.168 / 0.488 | −1281.85 | 95.0 % | **76.4 %** | 16.4 % | 2.1 % | 0.2 % | 4.9 % (Dikarya = root) |

**Full gene set (957 genes, 720 species), D and L fixed at the single-tree ML of each data type** (`workflow/35` on
the ML tree: AA 0.418 / 1.079, AA+3Di 0.423 / 1.095). One ALE likelihood evaluation takes ~20 min at this size, so
the rates were not optimised.

| sample | logL | Sacch-only anchor gene clade | sac_1 | node_448 | node_3 | node_2 | deep |
|---|---|---|---|---|---|---|---|
| AA deep | −4691.67 | 0.4 % | 0.4 % | 0 | 0 | 0 | **99.6 %** (node_0 65.1, root 34.5) |
| AA from-197 | −4627.52 | 99.7 % | **62.1 %** | 26.4 % | 7.9 % | 2.7 % | 0.9 % |
| **AA pooled** | −4566.67 | **100 %** | **64.3 %** | 25.1 % | 7.2 % | 2.9 % | 0.5 % |
| AA+3Di deep | −4794.25 | 0 % | 0 | 0 | 0 | 0 | **100 %** (node_0 81.6, root 18.4) |
| **AA+3Di pooled** | −4728.54 | **96.2 %** | **55.7 %** | 27.3 % | 9.8 % | 2.8 % | 4.4 % (node_0 2.8, root 1.1, node_1 0.5) |

**Reading:**
- **Single-optimum samples reproduce their optimum.** The deep samples contain the 197 split in 0 of 1000 trees,
  and ALE can only amalgamate clades it has seen.
- **Mixed samples decide in favour of Saccharomycotina.** Where the sample contains both topologies, the
  reconciliation moves the posterior strongly toward the Saccharomycotina placement:
  - **AA from-197:** 647 → 995 of 1000;
  - **AA pooled:** 654 of 2000 (≈ 1:2) → 1000 of 1000;
  - **AA+3Di from-197:** 522 → 1000 of 1000;
  - **AA+3Di pooled:** 526 of 2000 (≈ 1:3) → 992 (408 genes) / 962 (957 genes) of 1000;
  - **LBA control:** 305 → 950 of 1000.
- **Size of the effect.** The implied reconciliation factor is ≳ e^5 to e^7, consistent with the ~5.5-lnL
  duplication-local reconciliation term in §5.
- **Where the duplication maps.** On v2fix, about 90 % of the mass is on sac_1 or node_448: early in crown
  Saccharomycotina, mostly *after* Lipomyces. The rest (7–11 %) puts the same Saccharomycotina-only duplication on
  the Saccharomyceta/Ascomycota stem, with losses in Pezizomycotina.
- **On the NCBI backbone the mass moves up the stem:** node_448 50 %, Saccharomyceta 34 %, Ascomycota (node_2)
  12 %, deeper 3 %. NCBI's arbitrary resolution of the budding-yeast backbone makes a crown placement costlier.
  The gene content of the duplicated clade is still Saccharomycotina-only in 99 % of samples, and NCBI fits far
  worse (logL −1515 vs −1332). The fungus-wide placement stays ≤ 3 % on both backbones; *where on the
  Saccharomycotina stem* the duplication sits depends on the species tree.
- **The ALE logL difference between the from-197 and deep samples** measures the whole tree distribution,
  including the rest of the tree, not just the duplication node: AA +23.5 (408 genes) / +64.2 (957); AA+3Di −30.4
  (408).

---

## 4. GeneRax SPR from the deep and the 197 starts (`generax/`)

**Command** (408-gene set, v2fix species tree; `generax/inputs/families_red_<start>.txt` with `subst_model = LG+I+G4`):
```
generax -f generax/inputs/families_red_<start>.txt -s reduced/v2fix.species.nwk -r UndatedDL --unrooted-gene-tree \
        --strategy SPR --max-spr-radius 5 --seed 42 -p generax/red_v2fix_<start>
```
- D and L are optimised before every radius phase (GeneRax default).
- Starts, pruned to 408 genes:
  - **aa_deep:** `results/seq_tree/pif1.treefile`;
  - **aa_197:** `results/rerun_2026-10/au/aa_constr197_s12345.treefile` (best LG 197-constrained tree);
  - **aa_197Q:** `results/rerun_2026-10/au/aa_free_from197_QyeastR10.treefile` (best Q.yeast+R10 AA tree; it
    contains the 197 clade).

All three runs finished (radius 5; "End of GeneRax execution"): aa_197Q on 2026-10-17 after 23.6 h, aa_197 on
2026-10-17 after 25.4 h, aa_deep on 2026-10-19 after 77.0 h. Summary: `generax/spr_summary.tsv` (`workflow/47`,
command in §7); final trees `generax/red_v2fix_<start>/results/pif1/geneTree.newick`.

| start | initial joint (seq + rec) | after radius 1: joint / rec / seq | 197 split / anchor MRCA after radius 1 | **final (radius 5): joint / rec / seq** | final D, L (GeneRax) | **final 197 split / anchor MRCA** | final anchor duplication (`workflow/35`, ML D/L, reduced v2fix) |
|---|---|---|---|---|---|---|---|
| aa_deep | −62,742.7 (−60,947.5 − 1,795.2) | −62,360.1 / −1,332.4 / −61,027.7 | absent / 407 of 408 (deep) | **−62,144.0** / −1,043.9 / −61,100.1 | 0.134, 0.456 | **present / 197** | **D → sac_1** (104 leaves); rec −1,037.5 at 0.127 / 0.426; 83 D, 354 SL |
| aa_197 | −62,695.0 (−60,940.2 − 1,754.9) | −62,307.0 / −1,318.6 / −60,988.4 | present / 197 | **−62,144.4** / −1,096.9 / −61,047.5 | 0.141, 0.468 | present / 197 | **D → sac_1** (104); rec −1,091.0 at 0.139 / 0.462; 89 D, 384 SL |
| aa_197Q | −62,696.9 (−60,926.3 − 1,770.6) | −62,285.4 / −1,286.8 / −60,998.6 | present / 197 | **−62,141.0** / −1,098.6 / −61,042.4 | 0.142, 0.457 | present / 197 | **D → sac_1** (104); rec −1,092.7 at 0.144 / 0.463; 92 D, 381 SL |

Joint lnL after each radius (1 → 5): aa_deep −62,360.1 / −62,223.0 / −62,179.2 / −62,168.9 / −62,144.0;
aa_197 −62,307.0 / −62,188.2 / −62,150.0 / −62,145.6 / −62,144.4; aa_197Q −62,285.4 / −62,175.0 / −62,145.7 /
−62,142.0 / −62,141.0.

**Convergence.**
- **The deep start crosses to the Saccharomycotina placement.** After radius 1 its tree was still deep (anchor MRCA
  407 of 408). In the final tree the Saccharomycotina Pif1 daughter (93 genes) is sister to the Rrm3 daughter
  (104 genes): the exact 197-gene split, with the anchor duplication on sac_1. In the radius-1 tree the two
  daughters were **8 edges apart** (`workflow/49` prints this), more than the maximum SPR radius of 5, so the
  crossing took several moves. GeneRax overwrites `results/pif1/geneTree.newick` after each radius, so the radius
  at which the tree crossed is not recorded.
- **All three starts end in the same basin, at nearly the same joint lnL.** The spread is **3.4 joint lnL**
  (aa_197Q −62,141.0, aa_deep −62,144.0, aa_197 −62,144.4). After radius 1 it was 74.7.
  - The deep start ends with a better reconciliation (−1,043.9 vs −1,097) and a worse sequence likelihood (by
    53–58) than the 197 starts. The joint likelihood is about the same.
  - The final trees are closer to one another (unrooted RF 42–60 of 810) than each is to its own radius-1 tree
    (RF 92–146) (`workflow/51` → `generax/final_tree_distances.txt`).
- **Rates.** As the gene trees become more concordant, the re-optimised rates fall:
  - D≈0.36, L≈0.90 on the 408-gene starting trees;
  - D≈0.18, L≈0.58 after radius 1;
  - D≈0.13–0.14, L≈0.46 at the end.

  On the 957-gene starting trees (`generax/full957_aborted/`), GeneRax 2.0.4's gradient search stalled at
  L = 1.0. `workflow/35` and GeneRax 2.1.3 give L≈1.07 there.
- **Swap on the radius-1 trees.** The local-swap test (§5) applied to these GeneRax-optimised trees
  (`local_swap/generax_r1/local_swap_summary.tsv`; reduced alignment, `workflow/35` on reduced v2fix):

  | background (after radius 1) | move | Δ seq | Δ rec | **Δ joint** |
  |---|---|---|---|---|
  | GeneRax deep tree | Pif1 → sister of Rrm3 (197 placement) | +2.68 | +0.88 | **+3.55** |
  | GeneRax 197 tree | Pif1 → deep position | −11.85 | −4.76 | **−16.61** |
  | GeneRax 197Q tree | Pif1 → deep position (MRCA 219, partial) | −16.95 | −8.98 | **−25.93** |

  In the GeneRax-optimised deep tree, the Saccharomycotina Pif1 clade has become sister to a 14-gene Agaricomycotina
  group: Agaricomycetes + Tremellomycetes, the second Agaricomycotina paralog lineage of review §1.6. Even there,
  the 197 placement is better on **both** sequence and reconciliation likelihood. At radius 1 the deep start
  stayed deep only because of the search radius; by radius 5 it had crossed (above).
- **Joint ranking of the six radius-1 trees:**
  - gr_197Q −62,184.8;
  - gr_197Q_todeep −62,210.8;
  - gr_197 −62,213.2;
  - gr_197_todeep −62,229.8;
  - gr_deep_to197 −62,264.1;
  - gr_deep −62,267.6.

  Sequence lnL alone (408-gene alignment): gr_deep is 39.8 below gr_197, p-AU 0.195.
- **Swap on the final (radius-5) trees** (`local_swap/generax_final/`, driver `run_final.sh`; same scoring as
  `generax_r1`). All three final trees have the 197 placement, so each is paired with a *todeep* copy. Two deep
  references were used:
  - **`pif1.treefile`** (main): Pif1 goes back to its position in the published deep AA tree, and the anchor
    MRCA becomes 407 of 408 genes (node_0);
  - **the radius-1 GeneRax deep tree** (`r1ref/`, the reference used in `generax_r1`): Pif1 becomes sister to
    the 14-gene Agaricomycotina group, giving anchor MRCA 219 and a duplication on node_1 (Dikarya).

  | background (final tree) | seq lnL, IQ-TREE re-optimised (p-AU) | rec lnL | joint | → node_0 (`pif1.treefile` ref): Δ seq / Δ rec / **Δ joint** | → node_1 (r1 ref): Δ seq / Δ rec / **Δ joint** |
  |---|---|---|---|---|---|
  | from deep start (gr5_deep) | −61,052.59 (0.043) | −1,037.47 | **−62,090.06** | −9.48 / −7.09 / **−16.57** | −10.66 / −5.75 / **−16.41** |
  | from 197 start (gr5_197) | −61,005.71 (0.552) | −1,091.04 | −62,096.75 | −8.44 / −6.83 / **−15.27** | −8.88 / −5.53 / **−14.41** |
  | from 197Q start (gr5_197Q) | **−60,999.10** (1.000) | −1,092.66 | −62,091.76 | −9.18 / −6.84 / **−16.02** | −10.27 / −5.54 / **−15.81** |

  The AU column is from the `pif1.treefile`-reference set. The `r1ref/` set gives the same lnL values for the
  three final trees.
  - **Every final tree favours the Saccharomycotina placement on both terms.** Moving Pif1 to the fungus-wide
    position costs 8.4–9.5 in sequence and 6.8–7.1 in reconciliation, **15.3–16.6 joint lnL**. The Dikarya
    position costs 14.4–16.4.
  - **On sequence alone, the todeep copies are not rejected in the 197 backgrounds:** p-AU 0.14 (gr5_197_todeep
    and gr5_197Q_todeep, `pif1.treefile` reference) and 0.048–0.050 (r1 reference). gr5_deep_todeep has
    p-AU 0.009 / 0.002.
  - **The three final trees are 6.7 joint lnL apart once branch lengths are re-optimised.** The deep start's
    final tree has the best joint lnL (−62,090.1) and the best reconciliation, but it is 53.5 below gr5_197Q in
    sequence lnL (p-AU 0.04). It trades sequence fit for concordance elsewhere in the tree.
  - IQ-TREE's re-optimised sequence lnL is 42–48 higher than GeneRax's libpll value for the same trees. This is
    probably because GeneRax does not fully re-optimise branch lengths and model parameters during SPR. The
    GeneRax joint values in the table above are therefore not comparable with these.

**AA+3Di in GeneRax: not possible** (one model per family; §2).

---

## 5. Local swap: the two placements on the same background tree (`local_swap/`, `workflow/49`, `workflow/50`)
- **Construction.** `workflow/49` takes a tree rooted on hPIF1 and makes a copy that differs **only** in where the
  Saccharomycotina Pif1 daughter attaches:
  - **to197:** sister to the Rrm3 daughter;
  - **todeep:** where the deep reference tree has it.
- **Scoring.**
  - Sequence lnL: IQ-TREE `-z ... -n 0 -zb 10000 -au`, with branch lengths and model re-optimised per tree.
    AA uses LG+I+G4; AA+3Di uses the partition model.
  - Reconciliation lnL: `workflow/35` at ML D/L rates, summed over gene-tree roots, all 957 genes.
  - Joint = sequence + reconciliation.

```
python workflow/49_local_swap_trees.py --deep results/seq_tree/pif1.treefile --tree aa_deep=... --tree aa_197LG=... \
    --tree aa_197Q=... --tree aa_from197LG=... --out-trees local_swap/aa_trees.nwk --out-names local_swap/aa_trees.names
$IQ -s results/seq_tree/aln.trim.fasta -z local_swap/aa_trees.nwk -m LG+I+G4 -n 0 -zb 10000 -au -zw -T 1 -seed 42 -pre local_swap/au_aa
$IQ -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition.nex -z local_swap/aa3di_trees.nwk -n 0 -zb 10000 -au -zw -T 1 -seed 42 -pre local_swap/au_aa3di
local_swap/run_recon.sh      # workflow/35 --optimize-rates on every tree x {v2fix, ncbi_plus}
python workflow/50_local_swap_summary.py
```

| background | anchor MRCA → after swap | seq lnL (p-AU) | Δ seq | Δ rec v2fix / NCBI | **Δ joint v2fix / NCBI** |
|---|---|---|---|---|---|
| AA deep (`pif1.treefile`) | 950 → 196 (to197) | −95,220.66 (0.517) | −11.46 | +5.20 / +3.62 | **−6.25 / −7.83** |
| AA 197, LG (`aa_constr197_s12345`) | 197 → 956 (todeep) | −95,218.16 (0.456) | −0.28 | −5.45 / −3.70 | **−5.74 / −3.98** |
| AA 197, Q.yeast tree (`aa_free_from197_QyeastR10`) | 197 → 956 | −95,231.48 (0.383) | −0.40 | −5.45 / −3.71 | **−5.85 / −4.11** |
| **AA from-197, best LG tree** | 197 → 956 | **−95,201.72** (1.000) | −4.75 | −5.40 / −3.65 | **−10.16 / −8.41** |
| AA+3Di deep (replicate) | 880 → 147 (to197; Rrm3 is split there) | −185,649.80 (1.000) | −57.80 | +5.70 / +3.72 | **−52.11 / −54.07** |
| AA+3Di 197 (`pif1_aa3di_supported`) | 197 → 946 | −185,683.90 (0.428) | +0.56 | −5.44 / −3.63 | **−4.87 / −3.07** |

A negative Δ for a todeep move favours the 197 placement; a negative Δ for a to197 move favours deep.
- **The reconciliation term** always favours the Saccharomycotina placement: 5.2–5.7 on v2fix, 3.6–3.7 on NCBI.
- **Joint, 197 backgrounds:** the 197 placement wins (3.1–10.2).
- **Joint, deep backgrounds:** deep wins, because those trees are arranged around the deep placement. The AA+3Di
  deep tree splits the Saccharomycotina Rrm3 genes, so its swap is only a partial 147-gene clade.
- **Across backgrounds:** every AA tree with the 197 clade has a higher joint lnL than the deep AA tree:
  - from-197 −101,184.1 / −101,378.1;
  - 197 LG −101,173.2 / −101,387.2;
  - 197 Q.yeast −101,143.1 / −101,346.1;
  - deep −101,209.6 / −101,401.6 (v2fix / NCBI).

  For AA+3Di it is the reverse: deep −191,687.2 vs 197 −191,779.8 on v2fix. That gap comes from the rest of the
  tree (sequence −34.1, reconciliation −58.5), not from the duplication node.

---

## 6. Fixed-tree baseline (`fixed_tree_baseline/`, `workflow/35`, ML rates, all 957 genes)
| gene tree | species tree | D / L | rec lnL (sum over roots) | anchor dup → |
|---|---|---|---|---|
| AA deep (`pif1.treefile`) | v2fix | 0.418 / 1.079 | −5988.93 | node_0 |
| AA 197 (`aa_constr197_s12345`) | v2fix | 0.414 / 1.073 | **−5954.99** | sac_1 |
| AA+3Di deep (replicate) | v2fix | 0.423 / 1.095 | **−6037.41** | node_0 |
| AA+3Di 197 (`_supported`) | v2fix | 0.432 / 1.104 | −6095.93 | sac_1 |
| AA deep | NCBI+ | 0.425 / 1.096 | −6180.88 | node_0 |
| AA 197 | NCBI+ | 0.424 / 1.095 | **−6169.05** | node_448 |
| AA+3Di deep | NCBI+ | 0.425 / 1.108 | **−6230.39** | node_0 |
| AA+3Di 197 | NCBI+ | 0.433 / 1.115 | −6284.39 | node_448 |

- With all 957 genes mapped, the duplication maps to **sac_1 (v2fix) / node_448 (NCBI)** whenever the 197 clade is
  present, and to node_0 otherwise.
- v2fix fits better than NCBI by 188–214 lnL for every tree.

---

## 7. Running jobs, caveats, follow-ups
- **Nothing is still running.** The three reduced GeneRax runs finished radius 5 on 2026-10-17 (aa_197Q, aa_197)
  and 2026-10-19 (aa_deep); §4 has the final numbers. They were summarised (2026-10-20) with:
  ```
  python workflow/47_generax_spr_summary.py --species-tree results/rerun_2026-10/species_tree_aware/reduced/v2fix.species.nwk \
    --map results/rerun_2026-10/species_tree_aware/reduced/gene_species.map --out results/rerun_2026-10/species_tree_aware/generax/spr_summary.tsv \
    --run deep=results/rerun_2026-10/species_tree_aware/generax/red_v2fix_aa_deep --run 197=results/rerun_2026-10/species_tree_aware/generax/red_v2fix_aa_197 \
    --run 197Q=results/rerun_2026-10/species_tree_aware/generax/red_v2fix_aa_197Q
  python workflow/51_generax_tree_distances.py --tree r1_deep=... --tree final_deep=... --out generax/final_tree_distances.txt
  local_swap/generax_final/run_final.sh      # workflow/49 trees made first; see §4
  ```
  Three iCloud conflict copies with a " 2" suffix in `red_v2fix_aa_deep/` (from an aborted 3-rank start, dated
  2026-10-16 14:54–14:55, before this run started) were not part of the result and were left out of this repository.
- **Caveats.**
  - The reconciliation term is model-dependent: undated DL, no transfers, uniform origination. With high ML loss
    rates (single-tree L ≈ 1.1; CCP-based L ≈ 0.5), extra losses are cheap. That is why the deep placement costs
    only ~5 lnL.
  - The joint "sequence + reconciliation" sums in §5 are not a calibrated test. They are likelihoods under a joint
    model, as in GeneRax.
  - ALE posteriors depend on which topologies the UFBoot sample contains. That is why the from-197 and pooled
    samples were built.
  - The deep AA and AA+3Di trees are still not rejected by their alignments.
  - The ALE runs use the patched sampler (§2). Tables from an unpatched ALE with τ=0 would be wrong.
- **For the paper.**
  - Cite GeneRax 2.0.4 and ALE with the build notes, or rerun §3–§5 on the Linux host (GeneRax 2.0.4 as built
    there; ALE built x86_64 with the one-line patch).
  - The new best AA tree (from-197, with the 197 clade) should replace `pif1.treefile` in any AA-vs-AA+3Di
    comparison.
  - The PMSF run finished at 16:52; its supports are in `results/rerun_2026-10/pmsf/work/results/seq_tree/pif1_pmsf.iqtree`
    and have not been read here.
