#!/usr/bin/env Rscript
################################################################################
# Installs the packages that are NOT available from conda-forge / bioconda,
# pinned to exact versions so every install is identical.
#
#   conda activate warming_plaine_morte        && Rscript envs/install_github_packages.R main
#   conda activate warming_plaine_morte_dep2   && Rscript envs/install_github_packages.R dep2
#
# main : ANCOMBC, PCAtest, EcolUtils (GitHub) + gratia, khroma, spaa, metacoder,
#        taxonomizr (CRAN snapshot)
# dep2 : fdrtool (CRAN snapshot) + DEP2 (GitHub)
################################################################################
target <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(target) || !target %in% c("main", "dep2")) {
  stop("Usage: Rscript envs/install_github_packages.R <main|dep2>", call. = FALSE)
}

# Date-pinned CRAN snapshot (Posit Package Manager) for CRAN-only packages
CRAN_SNAPSHOT <- "https://packagemanager.posit.co/cran/2026-10-01"

# ---- GitHub pins (commit SHAs) ------------------------------------------------
# ANCOMBC: the project installed "FrederickHuangLin/ANCOMBC@bugfix". That branch
# moves; its HEAD (2.15.x) needs R >= 4.5. This is the LAST @bugfix commit that
# supports R 4.4 (v2.9.2, 2025-03-16).
# >>> VERIFY against the version used for the paper:
# >>>   packageDescription("ANCOMBC")[c("Version", "RemoteSha")]
GITHUB_PINS <- list(
  main = c("FrederickHuangLin/ANCOMBC@6d46f622ebc8a785dcb89ec0f6b18b464376532c",
           "arleyc/PCAtest@7e10abe7c40db30241d2e5eb4de1c6c73e7a361d",
           "GuillemSalazar/EcolUtils@2b57a8a8b76eba0852c56ce30f93d92c30320161"),
  dep2 = c("mildpiggy/DEP2@fbc71b5dc499fd6114fd8705a2f73e4b429e60f5")
)

# bioconda ships GenomeInfoDbData as a placeholder; a conda post-link script
# (installBiocDataPackage.sh) downloads the real data. If that step was skipped
# (offline install, --no-... flags), every Bioconductor package importing
# GenomeInfoDb fails to load with "no package called GenomeInfoDbData".
if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
  message("== Re-running the GenomeInfoDbData post-link installer")
  system2("bash", c("-c", shQuote("installBiocDataPackage.sh genomeinfodbdata-1.2.13")))
  if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
    stop("GenomeInfoDbData still missing. Try: conda install -c bioconda ",
         "--force-reinstall bioconductor-genomeinfodbdata  (needs internet)", call. = FALSE)
  }
}

# ---- CRAN-only packages (installed first; compiled from source if needed) -------
CRAN_PKGS <- list(main = c("gratia", "khroma", "spaa", "metacoder", "taxonomizr"),
                  dep2 = c("fdrtool"))
for (pkg in CRAN_PKGS[[target]]) {
  if (requireNamespace(pkg, quietly = TRUE)) next
  message("== Installing ", pkg, " from ", CRAN_SNAPSHOT)
  install.packages(pkg, repos = c(CRAN = CRAN_SNAPSHOT))
}

# Download each pinned commit as a plain tarball (codeload) instead of using
# remotes::install_github(): no GitHub API call, so a stale token in your git
# credential store ("HTTP error 401 Bad credentials") cannot break the install.
# conda already provides all dependencies -> never pull anything else.
github_tarball_url <- function(ref) {
  repo <- sub("@.*$", "", ref); sha <- sub("^.*@", "", ref)
  paste0("https://github.com/", repo, "/archive/", sha, ".tar.gz")
}
for (ref in GITHUB_PINS[[target]]) {
  pkg <- sub("@.*$", "", basename(ref))
  message("== Installing ", ref)
  remotes::install_url(github_tarball_url(ref), dependencies = FALSE, upgrade = "never")
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " failed to install", call. = FALSE)
}

pk <- c(sub("@.*$", "", basename(GITHUB_PINS[[target]])), CRAN_PKGS[[target]])
ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = pk, installed = ok,
                 version = vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, "")))
if (!all(ok)) quit(status = 1)
