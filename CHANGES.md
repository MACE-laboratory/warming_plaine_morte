# Changes that can alter numbers or figures (work in progress)

Only changes that affect results are listed here. Scripts 6_b and 7–11 are not yet cleaned up.

| # | Change | Where | Effect |
|---|---|---|---|
| 1 | **CLR now centred per sample across features** (textbook CLR). Before, `clr()` received a features x samples table and centred each *feature* across samples, i.e. plain log relative abundance shifted per feature. Verified on a toy table. | `0_` `clr_with_min_replacement()`; used for MAG, KO, PFAM, CAZy and 16S tables | Scripts 7 and 6_b results change (GAM fits of MAG/KO CLR vs fluxes; Spearman metaG–metaP correlations). **Re-check the MAG IDs hard-coded in script 7** (`sub_assembly_1_bin_11719`, `sub_assembly_3_bin_752717`): they were picked with the old transform. |
| 2 | **Each sample pair counted once** in pairwise dissimilarity tables (was A-B and B-A). | `4_` `pairwise_dissimilarity()` | `n_pairs` halves; regression / correlation p-values and the Wilcoxon / Kruskal tests of SFig 9 change. Pairs still share samples, so p-values remain descriptive. |
| 3 | 16S CLR computed from relative abundances with the same transform (was counts, wrong axis). | `0_` | 16S CLR table is not used by any numbered script. |
| 4 | 16S rarefaction seeded (`GLOBAL_SEED = 23`, RNG state restored afterwards). | `0_` `load_16sr_data()` | Rarefied-table results (SFig 6) change slightly vs unseeded runs; now reproducible. |
| 5 | PCA axis labels computed from the data (were typed; time and aspect labels were identical). | `1_` | Labels may differ. |
| 6 | PCAtest in the aspect section tests `pca_mat_aspect` (was `pca_mat_time`); PCAtest, spec.gen, RAC bootstrap seeded. | `1_`, `3_` | Re-check the hard-coded `sig_vars_*` indices against `stats/01_pcatest_*.txt`. |
| 7 | Duplicate rank-abundance recomputation removed (identical result). | `3_` | none |

Left unchanged on purpose: the specialist analysis in `3_` is restricted to 2024 samples (see the NOTE in the script).
