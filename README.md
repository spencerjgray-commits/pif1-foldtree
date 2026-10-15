# PIF1/RRM3 duplication in budding yeasts — code and data

Code, data and outputs for a manuscript in preparation (Gray, Kumcu & Bochman) that dates the gene duplication
producing the budding-yeast Pif1 and Rrm3 helicases.

**Status:** pre-submission. Analyses are being revised, and the tree, reconciliation and dating outputs here will
be replaced by the revised runs before a tagged, Zenodo-archived release.

## What the pipeline does
1. Collects the cellular PIF1-family helicases (InterPro **IPR048293**, not the Helitron-contaminated Pfam PF05970)
   from UniProt reference proteomes across Fungi, plus human PIF1. That gives 957 proteins from 727 fungal species and human.
2. Gathers structures: 828 AlphaFold DB (v6) models, plus 129 ColabFold (AlphaFold2) models for proteins without one.
3. Trims sequences and structures to the PF05970 helicase core.
4. Builds gene trees: amino-acid ML (IQ-TREE), FoldTree structural distance trees (Foldseek → FastME), and a
   partitioned amino-acid + 3Di ML tree.
5. Reconciles the gene trees with species trees (GeneRax 2.0.4, UndatedDL), using an NCBI-taxonomy backbone and a
   backbone grafted onto the Shen et al. (2018) budding-yeast chronogram.
6. Places the duplication on the Shen et al. (2018) time tree.

## Layout
```
manifest.csv                 one row per protein: accession, taxonomy, structure source, pLDDT, structure path
environment.yml, env.lock.txt   conda environment (GeneRax 2.0.4 was run from a separate env; see below)
workflow/                    numbered pipeline scripts (01–19) and _aps_* analysis/plotting helpers
data/seqs/                   IPR048293 tables per phylum, selected.tsv/.faa (957), cores.faa, tip_map.tsv;
                             _superseded_PF05970/ keeps the Pfam-based candidate sets behind Fig. S1
data/structures/af3/         the 129 ColabFold models (AF2, despite the folder name) + inputs and how-to
data/species_tree/           NCBI and grafted species trees, gene→species map, Shen et al. 2018 time tree
data/3di/                    3Di substitution matrices
data/expansion238/           238-structure expansion check (inputs)
results/seq_tree/            alignments, AA / AA+3Di / PMSF trees, topology tests, branch-length re-estimation
results/struct_tree/         FoldTree distance trees (fident, alntmscore, lddt)
results/reconciliation/      GeneRax runs (AA, AA+3Di on both backbones, robustness variants)
results/reviewer/            robustness re-runs (R1 AFDB-only, R2 core pLDDT ≥ 87, R3 GTR20 3Di matrix)
results/lba/                 outgroup-removal long-branch-attraction control
results/expansion/           expansion-check outputs
docs/figures/                figure files (vector SVG + font-outlined PDF) and their generating scripts
```

Not included because they are large and regenerable:
- the 828 AlphaFold DB models (`workflow/03_fetch_structures.py` re-downloads them);
- core-trimmed structures (`workflow/07_structcorecut.py`);
- whole proteomes.

## Figures
Generating scripts are in `docs/figures/scripts/` and run from the repository root. File names predate the current
manuscript numbering:

| Manuscript | File |
|---|---|
| Fig. 1 | `fig1_pipeline` (Graphviz `.dot`) |
| Fig. 2 | `fig2_copynumber` |
| Fig. 3 | `fig3_combined` (panels also as `fig3_aa_vs_aa3di_tree`, `fig3b_branchlength`) |
| Fig. 4 | `fig4_reconciliation` |
| Fig. 5 | `fig5_timetree` |
| Fig. S1 | `figS1_family_filter` |
| Fig. S2 | `figS2_plddt` |
| Fig. S3 | `figS3_reconciliation_thirdkind` |
| Fig. S4 | `figS5_robustness.png` (raster; a scripted vector version is still to be written) |
| Fig. S5 | `figS4_duplication_counts` |
| — | `Q1_saccharomycotina_paralog_retention`, `Q2_pif1_vs_rrm3_identity` (`workflow/_aps_figures.py`) |

## Reproducing
`README_REPRODUCTION_PACKAGE.md` gives the stage-by-stage commands, inputs and outputs. Software:
- MAFFT 7.526, trimAl 1.5.1, HMMER 3.4, IQ-TREE 3.1.2, Foldseek 10, FastME 2.1.6.3, GeneRax 2.0.4
  (2.1.3, which is pinned in `env.lock.txt`, crashes during reconciliation), ete3 3.1.3;
- R with ape, ggtree, ggplot2, svglite and deeptime for the figures.

The 129 structures were predicted with ColabFold on an external GPU machine. Those models are included because
their exact software versions were not recorded.

## License
MIT (see `LICENSE`).
