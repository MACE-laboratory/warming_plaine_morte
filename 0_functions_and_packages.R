################################################################################
# 0_functions_and_packages.R
#
# Shared setup for the "Warming Plaine Morte" analysis. This file is sourced by
# the numbered scripts (1_ ... 11_); it is not meant to be run on its own.
#
#   1. Package check + loading   (packages come from the conda env, see README)
#   2. Output helpers            (save_fig / save_stat_* + the FIGURE REGISTRY)
#   3. Data loaders              (MAGs, metagenome KO/PFAM/CAZy, 16S, metadata)
#   4. Small utilities           (normality-driven transform + scaling,
#                                 taxonomy aggregation helpers)
#
# Run every script from the repository root (open warming_plaine_morte.Rproj in
# RStudio, or `cd` into the repo before calling Rscript).
################################################################################


# ==============================================================================
# 1. PACKAGES
# ==============================================================================
# Nothing is installed from inside the scripts any more. Everything comes from
#   conda env create -f envs/warming_plaine_morte.yml
#   Rscript envs/install_github_packages.R        (ANCOMBC, PCAtest, EcolUtils, DEP2)
# If a package is missing we stop early with an explicit message instead of
# failing half-way through an analysis.

required_pkgs <- c(
  "tidyverse", "vegan", "ANCOMBC", "reshape2", "mgcv", "ggpubr", "ggrepel",
  "KEGGREST", "performance", "gratia", "ggstance", "cowplot", "PCAtest",
  "EcolUtils", "paletteer", "scales", "phyloseq", "compositions", "microbiome",
  "withr"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
       "  1) conda env create -f envs/warming_plaine_morte.yml\n",
       "  2) conda activate warming_plaine_morte\n",
       "  3) Rscript envs/install_github_packages.R", call. = FALSE)
}

suppressPackageStartupMessages({
  library(tidyverse)     # dplyr, tidyr, readr, purrr, stringr, ggplot2, ...
  library(vegan)
  library(ANCOMBC)
  library(reshape2)
  library(mgcv)
  library(ggpubr)
  library(ggrepel)
  library(KEGGREST)
  library(performance)
  library(gratia)
  library(ggstance)
  library(cowplot)
  library(PCAtest)
  library(EcolUtils)
  library(paletteer)
  library(scales)
  library(phyloseq)
  library(compositions)
  library(microbiome)
})

# Several loaded packages (MASS, plyr-style packages, stats) mask dplyr verbs.
# Pin the dplyr versions so the order of library() calls never matters.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# One global seed for everything stochastic (NMDS, permutations, bootstraps...).
# Individual scripts re-seed where a block must be reproducible on its own.
GLOBAL_SEED <- 23
set.seed(GLOBAL_SEED)


# ==============================================================================
# 2. OUTPUT HELPERS
# ==============================================================================
# Convention
#   figures/<Fig_N | SFig_N>.pdf + .png      (names come from FIGURE_REGISTRY)
#   stats/<NN>_<analysis>_<detail>.csv|txt   (NN = number of the script)

FIG_DIR   <- "figures"
STATS_DIR <- "stats"
for (d in c(FIG_DIR, STATS_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists("data")) {
  stop("`data/` not found. Run from the repository root and make sure the ",
       "Zenodo archive has been unpacked into ./data (see README).", call. = FALSE)
}

# ---- FIGURE REGISTRY ---------------------------------------------------------
# Single place that maps a figure's internal key to its manuscript name.
# To renumber figures for the paper, edit ONLY the right-hand side here.
#   Main figures          -> Fig_N
#   Supplementary figures -> SFig_N   (panels of one SFig: SFig_10a, SFig_10b...)
# Numbers marked (?) are provisional and need to be confirmed against the manuscript.
FIGURE_REGISTRY <- c(
  # ---- main --------------------------------------------------------------
  fig_pca_temp_fluxes      = "Fig_1",   # 1_  PCA + soil temperature + CH4/CO2 fluxes
  fig_alpha_mags           = "Fig_2",   # 3_  MAG taxonomic trees + alpha diversity vs soil temp
  fig_beta_mags            = "Fig_3",   # 4_  NMDS / PERMANOVA / dissimilarity vs temp (MAGs)
  fig_diff_abundance_mags  = "Fig_4",   # 5_  ANCOM-BC2 MAGs + class ridgelines
  fig_functional_ordin_da  = "Fig_5",   # 6_b NMDS + volcano plots (metaG vs metaP)
  fig_flux_gam_partial     = "Fig_6",   # 7_  abiotic GAM partial effects (CH4 / CO2)
  fig_flux_biotic          = "Fig_7",   # 7_  CH4 ~ air temp + water + MAG; CO2 ~ soil temp + water
  fig_proteomics_network   = "Fig_8",   # 9_  (?) was "Figure_5_Proteomics_EVO" - collides with Fig_5
  fig_ko_category_panels   = "Fig_9",   # 11_ (?) KO category panels metaG vs metaP

  # ---- supplementary -----------------------------------------------------
  sfig_pca_time            = "SFig_1",  # 1_  (was SuppFigure1)
  sfig_pca_aspect          = "SFig_2",  # 1_  (was SuppFigure2)
  sfig_yearly_temp         = "SFig_3",  # 1_  (was SuppFigure3_yearly_temp)
  sfig_rac_gam             = "SFig_4",  # 3_  (was SuppFigure_RAC)
  sfig_euk_phylum          = "SFig_5",  # 10_ (kept: was already SFig_5)
  sfig_alpha_16s           = "SFig_6",  # 3_  (was Figure_S??_alpha_16s)
  sfig_alpha_seasonal      = "SFig_7",  # 3_  (previously not saved)
  sfig_rac_bootstrap_rho   = "SFig_8",  # 3_  (was SuppFigure3_RAC - name clashed with yearly temp)
  sfig_beta_condition      = "SFig_9",  # 4_  (was Figure_S??_bray_jaccard_condition)
  sfig_beta_16s            = "SFig_10", # 4_  (was Figure_S??_16s_betadiv)
  sfig_barplot_class_16s   = "SFig_11a",# 2_  (was Figure2_Class_16S.png)
  sfig_barplot_genus_16s   = "SFig_11b",# 2_  (was Figure2_Genus_16S.png)
  sfig_barplot_class_mag   = "SFig_11c",# 2_  (was Figure2_Class_MAG.png)
  sfig_barplot_genus_mag   = "SFig_11d",# 2_  (was Figure2_Genus_MAG.png)
  sfig_da_forest_mags_top  = "SFig_12", # 5_  (was forest_plot_MAG_FDR_soiltemp_EVO.png)
  sfig_da_forest_class_top = "SFig_13", # 5_  (was forest_plot_CLASS_FDR_soiltemp_EVO.png)
  sfig_ko_prot_correlation = "SFig_14", # 6_b (previously not saved)
  sfig_prot_qc             = "SFig_15", # 6_a (previously not saved) DEP2 QC panels
  sfig_pmoa_threshold      = "SFig_16", # 8_  (was SI_Figure_pmoa_threshold_phydist)
  sfig_pmoa_tree_liebner   = "SFig_17", # 8_  (previously not saved)
  sfig_pmoa_tree_metag     = "SFig_18", # 8_  (previously not saved)
  sfig_pmoa_dist_boxplot   = "SFig_19", # 8_  (previously not saved)
  sfig_pmoa_tree_final     = "SFig_20", # 8_  (previously not saved)
  sfig_pmoa_scatter_liebner= "SFig_21", # 8_  (previously not saved)
  sfig_pcatest_condition   = "SFig_22", # 1_  PCAtest diagnostics (previously only drawn on screen)
  sfig_pcatest_time        = "SFig_23", # 1_  PCAtest diagnostics (previously only drawn on screen)
  sfig_pcatest_aspect      = "SFig_24"  # 1_  PCAtest diagnostics (previously only drawn on screen)
)

# Save a ggplot (or ggarrange/patchwork object) as <FIG_DIR>/<name>.pdf and .png.
save_fig <- function(plot, key, width = 8, height = 6, dpi = 300, ...) {
  if (!key %in% names(FIGURE_REGISTRY)) {
    stop("Unknown figure key '", key, "'. Add it to FIGURE_REGISTRY in ",
         "0_functions_and_packages.R", call. = FALSE)
  }
  stem <- file.path(FIG_DIR, FIGURE_REGISTRY[[key]])
  # cairo_pdf embeds fonts and handles the degree/superscript glyphs used in labels
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height,
                  device = pdf_device, limitsize = FALSE, ...)
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white", limitsize = FALSE, ...)
  message("Saved figure: ", stem, ".pdf / .png")
  invisible(stem)
}

# Save BASE-R graphics (e.g. PCAtest, vegan plots). The expression is evaluated
# ONCE (inside the pdf device), the plot is recorded and replayed into the png.
save_base_fig <- function(expr, key, width = 8, height = 6, dpi = 300) {
  if (!key %in% names(FIGURE_REGISTRY)) {
    stop("Unknown figure key '", key, "'. Add it to FIGURE_REGISTRY.", call. = FALSE)
  }
  stem <- file.path(FIG_DIR, FIGURE_REGISTRY[[key]])
  expr <- substitute(expr); env <- parent.frame()
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  pdf_device(paste0(stem, ".pdf"), width = width, height = height)
  grDevices::dev.control("enable")
  rec <- NULL
  tryCatch({ eval(expr, env); rec <- grDevices::recordPlot() },
           finally = grDevices::dev.off())
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi)
  tryCatch(grDevices::replayPlot(rec), finally = grDevices::dev.off())
  message("Saved figure: ", stem, ".pdf / .png")
  invisible(stem)
}

# Save SEVERAL plots (list of ggplots) as one multi-page PDF under the key's name
# (used for QC panels that do not need a PNG each).
save_fig_multipage <- function(plots, key, width = 8, height = 6) {
  if (!key %in% names(FIGURE_REGISTRY)) {
    stop("Unknown figure key '", key, "'. Add it to FIGURE_REGISTRY.", call. = FALSE)
  }
  path <- file.path(FIG_DIR, paste0(FIGURE_REGISTRY[[key]], ".pdf"))
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  pdf_device(path, width = width, height = height, onefile = TRUE)
  tryCatch(for (p in plots) print(p), finally = grDevices::dev.off())
  message("Saved figure: ", path)
  invisible(path)
}

# Write a data.frame/tibble to stats/<name>.csv   (name e.g. "03_alpha_mags_gam_shannon")
save_stat_table <- function(x, name, row.names = FALSE, ...) {
  path <- file.path(STATS_DIR, paste0(name, ".csv"))
  utils::write.csv(x, path, row.names = row.names, ...)
  message("Saved stats: ", path)
  invisible(path)
}

# Write the printed output of one or more R objects (e.g. summary(model), an
# adonis2 table, an envfit object) to stats/<name>.txt. Replaces sink(), which
# silently swallows all output if the script errors before the closing sink().
save_stat_text <- function(name, ...) {
  path <- file.path(STATS_DIR, paste0(name, ".txt"))
  out <- unlist(lapply(list(...), function(o) c(utils::capture.output(print(o)), "")))
  writeLines(out, path)
  message("Saved stats: ", path)
  invisible(path)
}

# Same, for hand-written lines of text.
save_stat_lines <- function(name, lines) {
  path <- file.path(STATS_DIR, paste0(name, ".txt"))
  writeLines(lines, path)
  message("Saved stats: ", path)
  invisible(path)
}

# Record the R / package versions used by a script (stats/<NN>_sessionInfo.txt).
save_session_info <- function(script_id) {
  path <- file.path(STATS_DIR, paste0(script_id, "_sessionInfo.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), path)
  invisible(path)
}


# ==============================================================================
# 3. LABELS / CONSTANTS
# ==============================================================================

# Plot-label (plotmath) lookup for the scaled environmental variables.
PCA_label_conversion <- c(
  "Temp_soil_scaled" = "Temp[soil]",
  "Temp_air_scaled" = "Temp[air]",
  "Luminosity_scaled" = "Luminosity",
  "CH4_Flux_scaled" = "CH[4]~Flux",
  "CO2_Flux_scaled" = "CO[2]~Flux",
  "TC_scaled" = "TC","TOC_scaled" = "TOC","pH_scaled" = "pH",
  "Water_Content_Lab_scaled" = "Water~Content[Lab]",
  "Water_Content_Field_scaled" = "Water~Content[Field]","F._scaled" = "F^'-'",
  "NO3._scaled" = "NO[3]^'-'","Cl._scaled" = "Cl^'-'",
  "NO2._scaled" = "NO[2]^'-'","PO4._scaled" = "PO[4]^'3-'",
  "SO4._scaled" = "SO[4]^'2-'","Br._scaled" = "Br^'-'",
  "Li._scaled" = "Li^'+'","NH4._scaled" = "NH[4]^'+'",
  "Ca2._scaled" = "Ca^'2+'","Sr2._scaled" = "Sr^'2+'",
  "Na._scaled" = "Na^'+'", "Mg2._scaled" = "Mg^'2+'",
  "K._scaled" = "K^'+'", "Formate_scaled" = "Formate",
  "Malate_scaled" = "Malate", "Propionate_scaled" = "Propionate",
  "Lactate_scaled" = "Lactate", "Butyrate_scaled" = "Butyrate",
  "Oxalate_scaled" = "Oxalate", "Acetate_scaled" = "Acetate",
  "Citrate_scaled" = 'Citrate')   # (key was "Citrate_scale" - typo fixed)

# Environmental variables that go through check_normality_change() in every script
ENV_VARS_TO_TRANSFORM <- c(
  "Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
  "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",
  "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",
  "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
  "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
  "Other_Cations","Grouped_OC")

# Condition colours used across figures
COND_COLOURS <- c(control = '#012169', warming = '#C8102E')


# ==============================================================================
# 4. DATA LOADERS
# ==============================================================================

# Centred log-ratio (CLR) transform of a features x samples relative-abundance
# table. Zeros are replaced by the smallest non-zero relative abundance of the
# table, then each SAMPLE is centred across its features: log(x) - mean(log(x)).
# compositions::clr() treats ROWS as observations, so the table is transposed in
# and out. (Until this version the table was passed untransposed, which centred
# each FEATURE across samples = plain log(relative abundance) shifted per feature.)
# Returns a features x samples data.frame with the original dimnames.
clr_with_min_replacement <- function(relab_tab) {
  relab_tab[relab_tab == 0] <- min(relab_tab[relab_tab > 0])
  out <- t(unclass(clr(t(relab_tab))))
  as.data.frame(out, row.names = rownames(relab_tab))
}

# MAG abundance (CoverM trimmed mean), taxonomy (GTDB-Tk) and quality (CheckM2).
# Sample names are rewritten to 'L_<subsite>_<layer>_<dd>_<mm>_<yy>'.
load_mags_data <- function(){
  ab_tab = read.csv('data/mags/mags_trimmed_mean_abundances.tsv', sep='\t')
  tx_tab = read.csv('data/mags/gtdbtk.bac120.summary.tsv', sep='\t')
  ql_tab = read.csv('data/mags/quality_report.tsv', sep='\t')
  ql_tab = ql_tab %>% filter(Name %in% ab_tab$Genome)

  rownames(ab_tab) = ab_tab$Genome
  ab_tab$Genome = NULL

  # raw file names -> harmonised sample IDs
  colnames(ab_tab) = map_chr(colnames(ab_tab), function(x) strsplit(x, split = '_S')[[1]][1])
  new_colnames = c()
  for (i in 1:length(colnames(ab_tab))){
    colname = colnames(ab_tab)[i]
    subsite = substr(colname, 10, 10)
    layer = substr(colname, 12, 12)
    year = substr(colname, 2, 3)
    month = substr(colname, 4, 5)
    day = substr(colname, 6, 7)
    new_colname = paste(c('L', subsite, layer, day, month, year), collapse = '_')
    new_colnames = c(new_colnames, new_colname)
  }
  colnames(ab_tab) = new_colnames

  relab_tab = sweep(ab_tab, 2, colSums(ab_tab), FUN = "/")
  clr_tab = clr_with_min_replacement(relab_tab)

  # MAG info: quality + taxonomy split into ranks
  mag_info = data.frame(MAG=ql_tab$Name, Completeness=ql_tab$Completeness,
                        Contamination=ql_tab$Contamination, Length=ql_tab$Genome_Size)
  mag_info$Taxonomy = map_chr(mag_info$MAG, function(x) tx_tab$classification[tx_tab$user_genome == x])
  mag_info = mag_info %>% extract(Taxonomy,
                                  into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                  regex = "([^;]*);?([^;]*);?([^;]*);([^;]*);?([^;]*);?([^;]*);?([^;]*)",
                                  remove = FALSE)

  mag_info <- mag_info[match(rownames(ab_tab), mag_info$MAG), ]

  return(list(raw_table=ab_tab, relab_table=relab_tab, clr_table=clr_tab, mag_info=mag_info))
}

# Build raw / relative-abundance / CLR versions of a feature table.
.expand_feature_table <- function(feature_tab) {
  relab <- feature_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
  list(raw = feature_tab, relab = relab, clr = clr_with_min_replacement(relab))
}

.pack_functional_tables <- function(ko_tab, pfam_tab, cazy_tab) {
  ko   <- .expand_feature_table(ko_tab)
  pfam <- .expand_feature_table(pfam_tab)
  cazy <- .expand_feature_table(cazy_tab)
  list(ko_raw_table = ko$raw,   pfam_raw_table = pfam$raw,   cazy_raw_table = cazy$raw,
       ko_relab_table = ko$relab, pfam_relab_table = pfam$relab, cazy_relab_table = cazy$relab,
       ko_clr_table = ko$clr,   pfam_clr_table = pfam$clr,   cazy_clr_table = cazy$clr)
}

# Metagenome functional tables (KO / PFAM / CAZy x sample).
# First call builds them from the contig coverage + eggNOG annotations and caches
# them in data/{ko,pfam,cazy}_table.tsv; later calls just read the cache.
load_meta_data <- function() {
  if (file.exists("data/ko_table.tsv")) {
    message("Loading existing KO abundance table from data/ko_table.tsv")
    ko_tab   <- read_tsv("data/ko_table.tsv",   show_col_types = FALSE) %>% column_to_rownames("KEGG_ko")
    pfam_tab <- read_tsv("data/pfam_table.tsv", show_col_types = FALSE) %>% column_to_rownames("PFAMs")
    cazy_tab <- read_tsv("data/cazy_table.tsv", show_col_types = FALSE) %>% column_to_rownames("CAZy")
    return(.pack_functional_tables(ko_tab, pfam_tab, cazy_tab))
  }

  message("Creating KO abundance table from input files...")

  # ---- contig abundance (CoverM trimmed mean) --------------------------------
  ab_tab <- read_tsv("data/functional/all_contigs_trimmed_mean.tsv", show_col_types = FALSE)
  ab_tab[is.na(ab_tab)] <- 0
  colnames(ab_tab)[1] <- "contig"
  colnames(ab_tab)[-1] <- map_chr(colnames(ab_tab)[-1], ~strsplit(.x, "\\.\\.\\.")[[1]][1])
  colnames(ab_tab)[-1] <- gsub("^final\\.contigs\\.renamed\\.filtered\\.fa\\.", "", colnames(ab_tab)[-1])

  # raw file names -> harmonised sample IDs 'L_<subsite>_<layer>_<dd>_<mm>_<yy>'
  new_colnames <- c("contig")
  for (i in 2:ncol(ab_tab)) {
    colname <- colnames(ab_tab)[i]
    subsite <- substr(colname, 9, 9)
    layer <- substr(colname, 11, 11)
    year <- substr(colname, 1, 2)
    month <- substr(colname, 3, 4)
    day <- substr(colname, 5, 6)
    new_colnames <- c(new_colnames, paste(c('L', subsite, layer, day, month, year), collapse = '_'))
  }
  colnames(ab_tab) <- new_colnames

  # sanity checks on sequencing depth per sub-assembly (printed for inspection)
  print('Samples sums in the abundance file:')
  print(quantile(colSums(ab_tab[,colnames(ab_tab) != 'contig'])))
  for (i in 0:3) {
    print(paste0('Samples sums for assembly ', i, ':'))
    print(quantile(colSums(ab_tab[startsWith(ab_tab$contig, paste0('sub_assembly_', i)),
                                  colnames(ab_tab) != 'contig'])))
  }
  print('Range of contigs abundances (please check carefully there are no zeros):')
  print(quantile(rowSums(ab_tab[,colnames(ab_tab) != 'contig'])))

  # ---- eggNOG annotations of the 4 sub-assemblies ------------------------------
  ann_files <- list.files("data/functional", pattern = "assembly_\\d\\.emapper.annotations$", full.names = TRUE)
  ann_list <- lapply(ann_files, function(file) {
    prefix <- str_extract(basename(file), "assembly_\\d+")
    prefix <- paste0('sub_', prefix)
    lines <- read_lines(file)
    lines <- lines[!startsWith(lines, "##")]

    tmp <- tempfile()
    write_lines(lines, tmp)

    df <- read_tsv(tmp, show_col_types = FALSE)
    colnames(df)[1] = 'query'

    df %>%
      mutate(query = paste0(prefix, "_", query)) %>%
      mutate(
        KEGG_ko = gsub("ko:", "", KEGG_ko),
        contig = str_remove(query, "_\\d+$")
      )
  })
  annotations <- bind_rows(ann_list)
  # (The MAG-only annotation subset used for the metaproteomics inference is now
  #  written by 6_a_impute_normalise_proteins.R, where the MAG gene list is built.)

  # ---- merge coverage of identical sample columns (sub-assemblies) -------------
  ab_tab_filtered = ab_tab[ab_tab$contig %in% annotations$contig,]
  rowsums_before = rowSums(ab_tab_filtered[,colnames(ab_tab_filtered) != 'contig'])
  ab_tab_filtered = ab_tab_filtered %>% column_to_rownames('contig')
  ab_tab_filtered = sapply(unique(colnames(ab_tab_filtered)), function(x) rowSums(ab_tab_filtered[,grepl(x, colnames(ab_tab_filtered))]))
  ab_tab_filtered = as.data.frame(ab_tab_filtered) %>% rownames_to_column(var='contig')
  rowsums_after = rowSums(ab_tab_filtered[,colnames(ab_tab_filtered) != 'contig'])
  print(mean(rowsums_after == rowsums_before)) # check that no counts were increased/reduced (should be 1)

  # ---- one row per (contig, function) then sum over contigs --------------------
  ko_annot   = annotations %>% dplyr::select(KEGG_ko, contig) %>% separate_longer_delim(KEGG_ko, ',') %>% filter(KEGG_ko != '-')
  pfam_annot = annotations %>% dplyr::select(PFAMs, contig)   %>% separate_longer_delim(PFAMs, ',')   %>% filter(PFAMs != '-')
  cazy_annot = annotations %>% dplyr::select(CAZy, contig)    %>% separate_longer_delim(CAZy, ',')    %>% filter(CAZy != '-')

  ko_merged   = full_join(ab_tab_filtered, ko_annot,   by='contig') %>% filter(!(is.na(KEGG_ko)))
  pfam_merged = full_join(ab_tab_filtered, pfam_annot, by='contig') %>% filter(!(is.na(PFAMs)))
  cazy_merged = full_join(ab_tab_filtered, cazy_annot, by='contig') %>% filter(!(is.na(CAZy)))

  ko_tab = ko_merged %>% dplyr::select(-contig) %>% group_by(KEGG_ko) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))
  pfam_tab = pfam_merged %>% dplyr::select(-contig) %>% group_by(PFAMs) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))
  cazy_tab = cazy_merged %>% dplyr::select(-contig) %>% group_by(CAZy) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))

  # cache
  write_tsv(ko_tab, "data/ko_table.tsv")
  write_tsv(pfam_tab, "data/pfam_table.tsv")
  write_tsv(cazy_tab, "data/cazy_table.tsv")

  .pack_functional_tables(column_to_rownames(ko_tab, "KEGG_ko"),
                          column_to_rownames(pfam_tab, "PFAMs"),
                          column_to_rownames(cazy_tab, "CAZy"))
}

# Near-full-length 16S OTUs (2024 samples only: names ending in '_24').
# Filters: bacteria only, no chloroplast/mitochondria, no OTU in the negative
# controls, >1 read overall. Also returns a rarefied table (8000 reads).
load_16sr_data <- function(){
  tax_table = read.table("data/16S/otu_taxonomy_centroid.tsv", sep = ',', header = T)
  tax_table = tax_table %>% extract(Taxonomy,
                                  into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                  regex = "([^;]*);?([^;]*);?([^;]*);([^;]*);?([^;]*);?([^;]*);?([^;]*)",
                                  remove = FALSE)
  colnames(tax_table)[1] = 'OTU'

  tax_table = tax_table %>% filter(Order!="o__Chloroplast", Family!="f__Mitochondria", Domain=="d__Bacteria")

  otu_table = read.table("data/16S/otu_table.tsv", sep = '\t', header = T)
  rownames(otu_table) = otu_table$OTU.ID
  otu_table$OTU.ID = NULL
  otu_table = otu_table[rownames(otu_table) %in% tax_table$OTU,]

  # drop OTUs seen in negative controls, then the controls themselves
  otu_table = otu_table[rowSums(otu_table[,c('negative_control', 'negative_control_24')]) == 0,]
  otu_table = otu_table %>% dplyr::select(-negative_control, -negative_control_24)
  otu_table = otu_table[rowSums(otu_table) > 1,]
  tax_table = tax_table[tax_table$OTU %in% rownames(otu_table),]

  relab_tab = sweep(otu_table, 2, colSums(otu_table), FUN = "/")

  # CLR on the relative abundances, same transform as the MAG / KO tables
  # (previously computed from the count table and with the wrong orientation;
  # the 16S CLR table is not used by any of the numbered scripts).
  clr_tab = clr_with_min_replacement(relab_tab)

  # rarefaction by multinomial resampling to 8000 reads (samples below depth -> NA)
  # The seed makes the rarefied table reproducible; with_seed() restores the
  # global RNG state afterwards so downstream random numbers are unaffected.
  rarefy_depth <- 8000

  rar_table <- withr::with_seed(GLOBAL_SEED, apply(otu_table, 2, function(x){
    if(sum(x) < rarefy_depth){
      return(rep(NA, length(x)))
    } else {
      probs <- x / sum(x)
      as.vector(rmultinom(1, size = rarefy_depth, prob = probs))
    }
  }))

  rownames(rar_table) <- rownames(otu_table)
  colnames(rar_table) <- colnames(otu_table)

  # keep the 2024 samples only
  to_keep = colnames(rar_table)[endsWith(colnames(rar_table), suffix = '_24')]
  otu_table = otu_table[,to_keep]
  relab_tab = relab_tab[,to_keep]
  clr_tab = clr_tab[,to_keep]
  rar_table = rar_table[,to_keep]

  return(list(raw_table=otu_table, relab_table=relab_tab, clr_table=clr_tab, rar_table=rar_table, tax_table=tax_table))
}

# Environmental metadata: sample IDs harmonised to 'L_<...>' with '_' separators.
load_metadata <- function(){
  metadata = read.csv('data/SiteL23_24_env_variables_final.csv')
  metadata$Sample = paste0('L_', metadata$Sample)
  metadata$Sample = gsub('\\.','_',metadata$Sample)
  return(metadata)
}

# load_metadata() + normality-driven transform/scale of all environmental
# variables (adds the '<var>_scaled' columns) + design factors:
#   Aspect: replicates C, D, E = South, A, B = North
#   OTC:    replicates A, B, C = open-top chamber (warming) 'Yes', otherwise 'No'
prepare_metadata <- function(){
  metadata <- load_metadata()
  metadata <- check_normality_change(metadata, ENV_VARS_TO_TRANSFORM)
  metadata$Aspect <- ifelse(metadata$Replicate %in% c('C','D','E'), 'South', 'North')
  metadata$OTC    <- ifelse(metadata$Replicate %in% c('A','B','C'), 'Yes', 'No')
  metadata
}

# Samples >= 90 days after the start (warming vs control comparison) and the
# control-only seasonal time series. Both come back with Aspect/OTC/Replicate/
# Condition as factors, as needed by the GAMs.
subset_metadata <- function(metadata){
  as_factors <- function(df) {
    df$Aspect <- as.factor(df$Aspect)
    df$OTC <- as.factor(df$OTC)
    df$Replicate <- as.factor(df$Replicate)
    df$Condition <- as.factor(df$Condition)
    df
  }
  list(
    condition  = as_factors(metadata %>% filter(Days_Since_Start >= 90)),
    timeseries = as_factors(metadata %>% filter(Condition == 'control'))
  )
}


# ==============================================================================
# 5. TAXONOMY AGGREGATION HELPERS
# ==============================================================================

# Sum a features x samples table over a taxonomic string built from `rank_cols`
# of `info` (matched on `id_col`). Strings ending in ';' / '; ' (unresolved
# lower ranks) are lumped into 'Others'.
.aggregate_by_taxonomy <- function(tab, info, id_col, rank_cols, sep, group_name) {
  info[[group_name]] = apply(info[, rank_cols], 1, paste, collapse = sep)
  info[[group_name]][endsWith(info[[group_name]], '; ')] = 'Others'   # same test as the original code
  tab[[group_name]] = map_chr(rownames(tab), function(x) info[[group_name]][info[[id_col]] == x])
  tab_sum <- tab %>%
    group_by(.data[[group_name]]) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  tab_sum %>% column_to_rownames(var = group_name)
}

# 16S relative abundance -> family level (Phylum|Class|Order|Family)
a16s_to_family <- function(loaded_16sr){
  tab = as.data.frame(loaded_16sr$relab_table)
  family_info = loaded_16sr$tax_table
  family_info$family_full = apply(family_info[, c('Phylum','Class','Order','Family')], 1, paste, collapse = "|")
  tab$family_full = map_chr(rownames(tab), function(x) family_info$family_full[family_info$OTU == x])
  tab_sum <- tab %>%
    group_by(family_full) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  tab_sum = tab_sum %>% column_to_rownames(var = 'family_full')
  return(tab_sum)
}

# Counts / relative abundances aggregated to genus, order, class
a16s_to_genus <- function(loaded_16sr){
  .aggregate_by_taxonomy(as.data.frame(loaded_16sr$raw_table), loaded_16sr$tax_table, 'OTU',
                         c('Phylum','Class','Order','Family','Genus'), ";", 'genus_full')
}
a16s_to_class <- function(loaded_16sr){
  .aggregate_by_taxonomy(as.data.frame(loaded_16sr$raw_table), loaded_16sr$tax_table, 'OTU',
                         c('Phylum','Class'), "; ", 'class_full')
}
mags_to_genus <- function(loaded_mags){
  .aggregate_by_taxonomy(as.data.frame(loaded_mags$relab_table), loaded_mags$mag_info, 'MAG',
                         c('Phylum','Class','Order','Family','Genus'), ";", 'genus_full')
}
mags_to_order <- function(loaded_mags){
  .aggregate_by_taxonomy(as.data.frame(loaded_mags$relab_table), loaded_mags$mag_info, 'MAG',
                         c('Phylum','Class','Order'), "; ", 'order_full')
}
mags_to_class <- function(loaded_mags){
  .aggregate_by_taxonomy(as.data.frame(loaded_mags$relab_table), loaded_mags$mag_info, 'MAG',
                         c('Phylum','Class'), "; ", 'class_full')
}


# ==============================================================================
# 6. NUMERIC UTILITIES
# ==============================================================================

# log(x + half of the smallest non-zero value)
log_half_min <- function(x){
  half_min = min(x[x > 0], na.rm = T) / 2
  return(log(x + half_min))
}

# shift so the minimum is 0 (makes log/sqrt transforms valid)
align_to_zero <- function(x){
  x = x - min(x, na.rm=T)
  return(x)
}

test_shapiro <- function(x){
  t = shapiro.test(x)
  return(t$p.value)
}

center_scale <- function(x){
  x = scale(x)
}

# For each variable: shift to >= 0, then choose among {none, log, sqrt} the
# transform giving the highest Shapiro-Wilk p-value (only if the raw variable
# fails normality at p < 0.05), centre/scale it and store it as '<var>_scaled'.
# The decision made for each variable is printed.
check_normality_change <- function(metadata, variables){
  new_col = c()
  for (var in variables){
    pval = test_shapiro(metadata[,var])

    if (pval < 0.05){
      new_col = align_to_zero(metadata[,var])
      log_col = log_half_min(new_col)
      pvlog = test_shapiro(log_col)
      sqrt_col = sqrt(new_col)
      pvsqrt = test_shapiro(sqrt_col)

      if ((pval > pvlog) & (pval > pvsqrt)){
        print(paste(var, ' was not transformed'))
        new_col = align_to_zero(metadata[,var])
      }

      if ((pvlog > pval) & (pvlog > pvsqrt)){
        print(paste(var, ' was log transformed'))
        new_col = log_col
      }

      if ((pvsqrt > pval) & (pvsqrt > pvlog)){
        print(paste(var, ' was sqrt transformed'))
        new_col = sqrt_col
      }

    }

    else {
      print(paste(var, ' was not transformed'))
      new_col = align_to_zero(metadata[,var])
    }
    centered_scaled = center_scale(new_col)
    metadata[,paste0(var, '_scaled', collapse = '')] = centered_scaled
  }
  return(metadata)
}
