################################################################################
# 4_Microbiome_beta_diversity.R
#
# Beta diversity of the soil microbiome vs soil temperature, for the MAGs
# (all samples) and the 16S amplicons (2024 samples):
#   - NMDS (Bray-Curtis) with envfit vectors of the environmental variables
#   - PERMANOVA (adonis2) of Bray-Curtis / Jaccard dissimilarity ~ soil temperature,
#     permutations restricted within season x year x replicate
#   - pairwise dissimilarity vs difference in soil temperature (same year & season)
#   - (MAGs only) dissimilarity within / across warming and control plots
#
# Inputs : data/mags/*, data/16S/*, data/SiteL23_24_env_variables_final.csv
# Outputs: figures/Fig_3, SFig_9 (MAGs), SFig_10 (16S)   (.pdf + .png)
#          stats/04_beta_<mags|16s>_*.csv|txt
#
# Methods notes
#  * vegdist(method = "jaccard") on RELATIVE abundances is the quantitative
#    (Ruzicka) form, not presence/absence Jaccard.
#  * Each sample pair is counted once in the pairwise tables. (Earlier versions
#    counted every pair twice, A-B and B-A, which inflated n and the p-values.)
#    The pairs still share samples, so the regression / correlation p-values on
#    the plots are descriptive and not corrected for that non-independence.
################################################################################
source('0_functions_and_packages.R')   # run from the repository root

loaded_mags = load_mags_data()
loaded_16s = load_16sr_data()
metadata = prepare_metadata()

# Human-readable labels for the envfit vectors shown on the NMDS
ENVFIT_LABELS <- c(TIN_scaled = 'TIN', TOC_scaled = 'TOC', TC_scaled = 'TC',
                   Water_Content_Lab_scaled = 'Soil water content',
                   pH_scaled = 'pH', Temp_soil_scaled = 'Soil temperature')

# Pairwise table (long format) of a dissimilarity matrix with the metadata of the
# two samples; keeps pairs from the same year AND season. Each unordered pair
# appears ONCE (the original kept both A-B and B-A, i.e. every pair twice).
pairwise_dissimilarity <- function(dmat, metadata, value_name) {
  pairs <- expand_grid(Sample1 = metadata$Sample, Sample2 = metadata$Sample) %>%
    filter(match(Sample1, metadata$Sample) < match(Sample2, metadata$Sample))
  i1 <- match(pairs$Sample1, metadata$Sample)
  i2 <- match(pairs$Sample2, metadata$Sample)
  pairs[[value_name]] <- dmat[cbind(pairs$Sample1, pairs$Sample2)]
  pairs$TempDiff   <- abs(metadata$Temp_soil_scaled[i1] - metadata$Temp_soil_scaled[i2])
  pairs$Replicate1 <- metadata$Replicate[i1]
  pairs$Replicate2 <- metadata$Replicate[i2]
  pairs$Year1      <- metadata$Year[i1]
  pairs$Year2      <- metadata$Year[i2]
  pairs$Season1    <- metadata$Season[i1]
  pairs$Season2    <- metadata$Season[i2]
  pairs %>% filter(Year1 == Year2, Season1 == Season2)
}

# Scatter of dissimilarity vs temperature difference, with lm line and the
# regression equation / Pearson p-value; returns the plot.
plot_dissim_vs_tempdiff <- function(pairs, value_name, ylab, ylim = NULL) {
  p <- ggplot(pairs, aes(x = TempDiff, y = .data[[value_name]])) +
    geom_point(alpha = 0.2, colour='slategray4') +
    geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
    stat_regline_equation(
      aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
      label.x.npc = "left",
      label.y.npc = 0.95) +
    stat_cor(
      aes(label = after_stat(p.label)),
      label.x.npc = "left",
      label.y.npc = 0.85) +
    theme_classic() +
    labs(x = "Difference in soil temperature (scaled)", y = ylab)
  if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])
  p
}

# Linear-regression / correlation summary of dissimilarity ~ temperature difference
dissim_lm_summary <- function(pairs, value_name, metric) {
  fit <- lm(reformulate("TempDiff", value_name), data = pairs)
  ct  <- cor.test(pairs$TempDiff, pairs[[value_name]])
  tibble(metric = metric, n_pairs = nrow(pairs),
         intercept = coef(fit)[1], slope = coef(fit)[2],
         r_squared = summary(fit)$r.squared,
         pearson_r = unname(ct$estimate), pearson_p = ct$p.value)
}

# Full beta-diversity analysis of one dataset.
#   relab_table : taxa x samples relative abundances
#   metadata    : samples to analyse (must contain 'Sample')
#   dataset     : short label used in file names ('mags' / '16s')
# Returns the assembled figure (A = NMDS, B/C = dissimilarity vs temperature
# difference) and the pairwise tables.
run_beta_diversity <- function(relab_table, metadata, dataset, ylim_bc = NULL, ylim_jc = NULL) {

  set.seed(GLOBAL_SEED)   # NMDS starts and PERMANOVA permutations
  relab_table <- relab_table[, colnames(relab_table) %in% metadata$Sample]
  relab_table <- relab_table[, match(metadata$Sample, colnames(relab_table))]
  stopifnot(all(metadata$Sample == colnames(relab_table)))

  # ---- NMDS + envfit -----------------------------------------------------------
  nmds <- metaMDS(t(relab_table), distance = "bray", k = 2, trymax = 1000, try = 500)

  scores_df <- as.data.frame(scores(nmds, display = "sites"))
  scores_df$Temp_soil <- metadata$Temp_soil
  scores_df$Replicate <- metadata$Replicate
  scores_df$Aspect    <- metadata$Aspect
  scores_df$OTC       <- metadata$OTC

  env_vars <- metadata[, c("Temp_soil_scaled", "TC_scaled", "TOC_scaled", "TIN_scaled",
                           "pH_scaled", "Water_Content_Lab_scaled")]

  # 9999 permutations, restricted within replicate
  ef <- envfit(nmds, env_vars, permutations = 9999, strata = metadata$Replicate)
  print(ef)
  save_stat_text(paste0("04_beta_", dataset, "_envfit"), ef)

  ef_scores <- scores(ef, display = "vectors")
  ef_df <- data.frame(var = rownames(ef_scores), NMDS1 = ef_scores[,1], NMDS2 = ef_scores[,2],
                      r2 = ef$vectors$r, p = ef$vectors$pvals)
  save_stat_table(ef_df, paste0("04_beta_", dataset, "_envfit_vectors"))

  # vectors with p < 0.1 are drawn (dashed if 0.05 <= p < 0.1)
  ef_sig <- ef_df[ef_df$p < 0.1, ]
  ef_sig$Significance <- ifelse(ef_sig$p < 0.05, '< 0.05', '< 0.1')
  ef_sig$var_renamed  <- ifelse(ef_sig$var %in% names(ENVFIT_LABELS),
                                ENVFIT_LABELS[ef_sig$var], ef_sig$var)

  nmds_env <- ggplot(scores_df, aes(NMDS1, NMDS2, color = Temp_soil)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_gradient2(low = '#012169', high = '#C8102E', name = 'Soil temperature [°C]') +
    geom_segment(
      data = ef_sig,
      aes(x = 0, y = 0, linetype = Significance, xend = NMDS1, yend = NMDS2),
      arrow = arrow(length = unit(0.25, "cm")),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "black"
    ) +
    geom_text(
      data = ef_sig,
      aes(x = NMDS1*1.2, y = NMDS2*1.1, label = var_renamed),
      inherit.aes = FALSE,
      size = 4,
      vjust = -0.6
    ) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste0("Stress = ", round(nmds$stress, 3)),
             hjust = 1.1, vjust = 1.1, size = 4) +
    theme(legend.position = 'top')

  # ---- PERMANOVA (Bray-Curtis and Jaccard) ---------------------------------------
  bray    <- vegdist(t(relab_table), method = "bray")
  jaccard <- vegdist(t(relab_table), method = "jaccard")
  # permutations are restricted within Season x Year x Replicate blocks
  metadata$Season_year <- paste0(metadata$Season, metadata$Year, metadata$Replicate)

  perm_bray <- adonis2(bray ~ Temp_soil_scaled, data = metadata, permutations = 999,
                       strata = metadata$Season_year, by = 'margin')
  perm_jacc <- adonis2(jaccard ~ Temp_soil_scaled, data = metadata, permutations = 999,
                       strata = metadata$Season_year, by = 'margin')
  print(perm_bray); print(perm_jacc)
  save_stat_text(paste0("04_beta_", dataset, "_permanova"),
                 "PERMANOVA - Bray-Curtis ~ Temp_soil_scaled", perm_bray,
                 "PERMANOVA - Jaccard (quantitative) ~ Temp_soil_scaled", perm_jacc,
                 paste0("NMDS stress: ", nmds$stress))
  save_stat_table(rbind(cbind(metric = "Bray-Curtis", as.data.frame(perm_bray) %>% rownames_to_column("term")),
                        cbind(metric = "Jaccard",     as.data.frame(perm_jacc) %>% rownames_to_column("term"))),
                  paste0("04_beta_", dataset, "_permanova_table"))

  # ---- pairwise dissimilarity vs temperature difference --------------------------
  pairs_bc <- pairwise_dissimilarity(as.matrix(bray),    metadata, "Bray")
  pairs_jc <- pairwise_dissimilarity(as.matrix(jaccard), metadata, "Jaccard")

  save_stat_table(bind_rows(dissim_lm_summary(pairs_bc, "Bray", "Bray-Curtis"),
                            dissim_lm_summary(pairs_jc, "Jaccard", "Jaccard")),
                  paste0("04_beta_", dataset, "_dissimilarity_vs_tempdiff_lm"))

  bc_temp <- plot_dissim_vs_tempdiff(pairs_bc, "Bray", "Bray\u2013Curtis dissimilarity", ylim_bc)
  jc_temp <- plot_dissim_vs_tempdiff(pairs_jc, "Jaccard", "Jaccard dissimilarity", ylim_jc)

  splot <- ggarrange(bc_temp, jc_temp, labels = c('B', 'C'), nrow = 2, ncol = 1)
  fig <- ggarrange(nmds_env, splot, ncol = 2, labels = c('A', ''), nrow = 1, widths = c(0.65, 0.35))

  list(figure = fig, pairs_bc = pairs_bc, pairs_jc = pairs_jc, nmds = nmds)
}


################################################################################
# MAGs (all samples)  -> Figure 3
################################################################################
beta_mags <- run_beta_diversity(loaded_mags$relab_table, metadata, "mags",
                                ylim_bc = c(0.13, 0.8), ylim_jc = c(0.1, 0.95))
beta_mags$figure
save_fig(beta_mags$figure, 'fig_beta_mags', width = 9, height = 7)

# ---- Dissimilarity within warming / within control / across (SFig 9) -------------
# Warming plots = replicates A, B, C.
add_condition_group <- function(pairs) {
  pairs$Condition1 <- ifelse(pairs$Replicate1 %in% c('A', 'B', 'C'), 'Warming', 'Control')
  pairs$Condition2 <- ifelse(pairs$Replicate2 %in% c('A', 'B', 'C'), 'Warming', 'Control')
  pairs$ConditionGroup <- 'Across'
  pairs$ConditionGroup[(pairs$Condition1 == pairs$Condition2) & (pairs$Condition1 == 'Control')] <- 'Within Control'
  pairs$ConditionGroup[(pairs$Condition1 == pairs$Condition2) & (pairs$Condition1 == 'Warming')] <- 'Within Warming'
  pairs
}
pairs_jc <- add_condition_group(beta_mags$pairs_jc)
pairs_bc <- add_condition_group(beta_mags$pairs_bc)

comparisons = list(c('Within Warming', 'Within Control'),
                   c('Across', 'Within Warming'),
                   c('Across', 'Within Control'))

plot_condition_group <- function(pairs, value_name, ylab, kruskal_y) {
  ggplot(pairs, aes(x = ConditionGroup, fill = ConditionGroup, y = .data[[value_name]])) +
    geom_boxplot() +
    stat_compare_means(label.y = kruskal_y) +                    # Kruskal-Wallis across the 3 groups
    stat_compare_means(method = "wilcox.test", comparisons = comparisons,
                       p.adjust.method = "BH") +                 # pairwise Wilcoxon
    theme_classic() +
    labs(x = "", y = ylab) +
    scale_fill_manual(values = c('violet', '#012169', '#C8102E'))
}
sp1 <- plot_condition_group(pairs_jc, "Jaccard", "Jaccard dissimilarity", 1.1)
sp2 <- plot_condition_group(pairs_bc, "Bray",    "Bray-Curtis dissimilarity", 1)

# The same tests as numbers (Kruskal-Wallis + pairwise Wilcoxon, BH-adjusted)
for (metric in c("Jaccard", "Bray")) {
  d <- if (metric == "Jaccard") pairs_jc else pairs_bc
  save_stat_table(compare_means(reformulate("ConditionGroup", metric), data = d, method = "kruskal.test"),
                  paste0("04_beta_mags_condition_groups_kruskal_", tolower(metric)))
  save_stat_table(compare_means(reformulate("ConditionGroup", metric), data = d, method = "wilcox.test",
                                p.adjust.method = "BH"),
                  paste0("04_beta_mags_condition_groups_wilcoxon_", tolower(metric)))
}

sfig_condition <- ggarrange(sp1 + theme(legend.position = 'none'), sp2 + theme(legend.position = 'none'),
                            ncol = 1, nrow = 2, labels = c('A', 'B'))
sfig_condition
save_fig(sfig_condition, 'sfig_beta_condition', width = 5, height = 7)


################################################################################
# 16S amplicons (2024 samples only)  -> SFig 10
################################################################################
metadata_2024 <- metadata %>% filter(Year == 2024)
beta_16s <- run_beta_diversity(loaded_16s$relab_table, metadata_2024, "16s")
beta_16s$figure
save_fig(beta_16s$figure, 'sfig_beta_16s', width = 9, height = 7)

save_session_info("04")
