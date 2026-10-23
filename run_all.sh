#!/usr/bin/env bash
# Runs the numbered scripts in order, switching conda env for the DEP2 step.
# Usage (from the repo root):  bash run_all.sh
# Requires: conda envs `warming_plaine_morte` and `warming_plaine_morte_dep2`
# (see README) and ./data from Zenodo.
set -euo pipefail

MAIN=warming_plaine_morte
DEP2=warming_plaine_morte_dep2

run_main() { echo "=== [$MAIN] $1"; conda run --no-capture-output -n "$MAIN" Rscript "$1"; }
run_dep2() { echo "=== [$DEP2] $1"; conda run --no-capture-output -n "$DEP2" Rscript "$1"; }

run_main 1_PCA_temperature.R
run_main 2_Abundance_Barplots_EVO.R
run_main 3_Microbiome_alpha_diversity.R
run_main 4_Microbiome_beta_diversity.R
run_main 5_Differential_abundance_MAGs_EVO.R   # slow (ANCOM-BC2)
run_main 6_a1_extract_mag_genes.R
run_dep2 6_a2_impute_normalise_proteins.R      # <- the only script in the DEP2 env
run_main 6_b_functional_differential_abundances.R
run_main 7_modelling_gas_fluxes.R
run_main 8_pmoas_trees.R
run_main 9_Proteomics_EVO.R
run_main 10_eukaryotes.R
run_main 11_metaG_metaP.R
echo "Done. Figures in figures/, statistics in stats/."
