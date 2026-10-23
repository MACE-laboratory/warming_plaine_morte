# Warming Plaine Morte — reproduction code

R scripts to reproduce the figures and statistics of the *warming Plaine Morte* study.

## 1) Get the code

```bash
git clone https://github.com/MACE-laboratory/warming_plaine_morte.git
cd warming_plaine_morte
```

## 2) Add the data (Zenodo)

Download the archive **<ZENODO DOI/LINK HERE>** and unpack it so that `data/` sits in the repository root (it is not versioned in Git):

```text
warming_plaine_morte/
  data/
  envs/
  figures/   (created by the scripts)
  stats/     (created by the scripts)
  *.R
```

## 3) Software environments

Two conda environments built from conda-forge + bioconda (do **not** use the `r`/`defaults` channels). If your conda adds `defaults` from its own config, run once `conda config --set channel_priority strict`. Use `mamba`/`micromamba` if you have it; the commands are identical. Both env files were solved for **osx-arm64** (Apple Silicon) and linux-64.

| env | used by | contents |
|---|---|---|
| `warming_plaine_morte` | every script **except** 6_a2 | R 4.4, tidyverse, vegan, mgcv, phyloseq, ANCOMBC, … |
| `warming_plaine_morte_dep2` | `6_a2_impute_normalise_proteins.R` only | R 4.4 + DEP2 and its Bioconductor dependencies |

```bash
# main environment
conda env create -f envs/warming_plaine_morte.yml
conda activate warming_plaine_morte
Rscript envs/install_github_packages.R main     # ANCOMBC, PCAtest, EcolUtils (GitHub) + gratia, khroma, spaa, metacoder, taxonomizr (CRAN)
conda deactivate

# DEP2 environment
conda env create -f envs/warming_plaine_morte_dep2.yml
conda activate warming_plaine_morte_dep2
Rscript envs/install_github_packages.R dep2     # DEP2 (GitHub) + fdrtool (CRAN)
conda deactivate
```

Notes
* GitHub and CRAN-only packages (`gratia`, `khroma`, `spaa`, `metacoder`, `taxonomizr`; `fdrtool` for dep2 — several have no conda build on Apple Silicon) are installed by `install_github_packages.R` at pinned commits / a dated CRAN snapshot. This needs internet and a compiler toolchain (the env provides one).
* The ANCOMBC pin is the last `@bugfix` commit that supports R 4.4. **Check it matches the version used for the paper** (`packageDescription("ANCOMBC")`), and edit the pin if not.
* Bioconductor 3.20 is used on purpose: ggtree ≥ 4 requires ggplot2 ≥ 4, whereas the scripts were written for ggplot2 3.5.

### RStudio
RStudio is not distributed through conda. Install it from posit.co, then make it use the conda R:

* macOS/Linux: `conda activate warming_plaine_morte && export RSTUDIO_WHICH_R=$(which R) && open -a RStudio warming_plaine_morte.Rproj` (Linux: `rstudio warming_plaine_morte.Rproj`)
* Windows: Tools ▸ Global Options ▸ General ▸ R version ▸ choose `<conda>/envs/warming_plaine_morte/bin/R`

Open `warming_plaine_morte.Rproj` so the working directory is the repository root. For `6_a2_…` do the same with the `warming_plaine_morte_dep2` env.

## 4) Run the analysis

All scripts must be run **from the repository root**. `0_functions_and_packages.R` is sourced by the others (packages, loaders, output helpers); do not run it alone.

```bash
bash run_all.sh        # runs everything in order and switches env for the DEP2 step
```

Order and environment:

| script | env |
|---|---|
| `1_PCA_temperature.R` … `5_Differential_abundance_MAGs_EVO.R` | main |
| `6_a1_extract_mag_genes.R` (MAG genes + annotation for the proteomics) | main |
| `6_a2_impute_normalise_proteins.R` (DEP2 filtering / VSN / KNN) | **dep2** |
| `6_b_…`, `7_…`, `8_…`, `9_…`, `10_…`, `11_…` | main |

## 5) Outputs

* `figures/` — `Fig_N.pdf/.png` (main), `SFig_N.pdf/.png` (supplementary). Names are defined in one place: `FIGURE_REGISTRY` in `0_functions_and_packages.R`.
* `stats/` — `NN_<analysis>_<detail>.csv|txt`, `NN` = script number; `NN_sessionInfo.txt` records package versions.

## Troubleshooting
* **`data/` missing** — unpack the Zenodo archive into `./data`.
* **"no package called GenomeInfoDbData"** — the bioconda post-link download did not run; `Rscript envs/install_github_packages.R <main|dep2>` retries it (needs internet).
* **Rplots.pdf appears** — a plot was drawn to the default device; delete it, it is git-ignored.
