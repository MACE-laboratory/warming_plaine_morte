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
conda env create -f <ENV_FILE_1>.yml
conda env create -f <ENV_FILE_2>.yml
```

If you use mamba:

```bash
mamba env create -f <ENV_FILE_1>.yml
mamba env create -f <ENV_FILE_2>.yml
```

### Activate an env

```bash
conda activate <ENV_NAME>
```

## 4) Important: `dep2` must be installed with `devtools` (proteins normalization/imputation)

The proteins normalization/imputation step requires the R package **`dep2`**, installed from source (GitHub) using `devtools`.

After activating the conda environment you use for the proteomics step, open R and run:

```r
install.packages("devtools")
devtools::install_github("vitek-lab/dep2")  # update if the script points to a different dep2 repo
```

Notes:
- On some systems, `devtools` requires a compiler toolchain (Rtools on Windows, Xcode Command Line Tools on macOS, `build-essential` on Linux).

## 5) Run the analysis (R scripts)

All R scripts should be run in order, from `0_...` through `7_...`.

### Option A — run from the command line

Activate the appropriate conda environment, then run:

```bash
# Example: run all numbered scripts in order (adjust if your naming differs)
Rscript 0_*.R
Rscript 1_*.R
Rscript 2_*.R
Rscript 3_*.R
Rscript 4_*.R
Rscript 5_*.R
Rscript 6_*.R
Rscript 7_*.R
```

If there are multiple scripts per number, run them in the intended order (e.g., by filename).

### Option B — run interactively

Open R / RStudio (inside the active conda env) and source the scripts in order from `0_...` to `7_...`.

## 6) Outputs

By default, scripts may write outputs to locations such as:
- `results/`
- `figures/`
- `tables/`

If these folders do not exist, create them:

```bash
mkdir -p results figures tables
```

## Troubleshooting

- **`data/` missing**: ensure the Zenodo download has been placed at `./data`.
- **`dep2` install fails**: install system build tools and retry; then verify `library(dep2)` works.
