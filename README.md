# Warming Plaine Morte — reproduction code

This repository contains R scripts and supporting configuration to reproduce the study analyses.

## 1) Get the code

```bash
git clone https://github.com/MACE-laboratory/warming_plaine_morte.git
cd warming_plaine_morte
```

## 2) Add the data from Zenodo

The repository expects a `data/` directory that is **not versioned in Git**.

1. Download the Zenodo archive for this study: **<ZENODO DOI/LINK HERE>**
2. Unzip it
3. Copy/merge the provided `data/` folder into the repository root so you have:

```text
warming_plaine_morte/
  data/
  (R scripts ...)
  (conda env files ...)
```

If your Zenodo download contains a different folder name, rename it to `data`.

## 3) Install conda environments

This project uses conda environments for reproducible execution.

### Prerequisites
- Conda (Miniconda or Anaconda), or Mambaforge recommended.
- (Optional but faster) `mamba` installed.

### Create the envs
Environment YAML files (`*.yml` / `*.yaml`) are included in the repository. Create the environments from those files (replace the filenames below with the ones present in this repo):

```bash
# Example (replace with the actual env file paths in this repo)
conda env create -f envs/warming_plaine_morte.yml
conda env create -f envs/warming_plaine_morte_dep2.yml
```

If you use mamba:

```bash
mamba env create -f envs/warming_plaine_morte.yml
mamba env create -f envs/warming_plaine_morte_dep2.yml
```

### Activate an env

For most scripts use:
```bash
conda activate warming_plaine_morte
```

For the 6_a_imput_normalise_proteins.R script use:
```bash
conda activate warming_plaine_morte_dep2
```

## 4) Important: `dep2` must be installed with `devtools` (proteins normalization/imputation) after conda installation

The proteins normalization/imputation step requires the R package **`dep2`**, installed from source (GitHub) using `devtools`.

After activating the conda environment you can run the `6_a_imput_normalise_proteins.R` script that installs **`dep2`** with devtools before loading it to process the proteomics data using dep2.

Notes:
- On some systems, `devtools` requires a compiler toolchain (Rtools on Windows, Xcode Command Line Tools on macOS, `build-essential` on Linux).

## 5) Run the analysis (R scripts)

All R scripts should be run in order, from `1_...` through `7_...`.

The R script `0_functions_and_packages.R` is loaded by others, it serves as a library loader and includes all custom functions.

Activate the `warmin_plaine_morte` conda environment, then run:

```bash
# Example: run all the following numbered scripts in that order
Rscript 1_PCA_temperature.R
Rscript 2_Abundance_Barplots_EVO.R
Rscript 3_Microbiome_alpha_diversity.R
Rscript 4_Microbiome_beta_diversity.R
Rscript 5_Differential_abundance_MAGs_EVO.R
Rscript 6_b_functional_differential_abundances.R
Rscript 7_modelling_gas_fluxes.R
```

## 6) Outputs

By default, scripts may write outputs to locations such as:
- `stats/` will contain tables and models' summaries in txt or tsv/csv formats.
- `figures/` will contain named figures in pdf format.

## Troubleshooting

- **`data/` missing**: ensure the Zenodo download has been placed at `./data`.
- **`dep2` install fails**: install system build tools and retry; then verify `library(dep2)` works.
