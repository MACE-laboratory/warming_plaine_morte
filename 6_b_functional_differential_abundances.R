setwd('~/Documents/MACE/EVO_study/warming_plaine_morte')

source('0_functions_and_packages.R')

loaded_mags = load_mags_data()
loaded_meta = load_meta_data()

metadata = load_metadata()
metadata = check_normality_change(metadata, c("Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
                                              "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",              
                                              "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",          
                                              "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
                                              "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
                                              "Other_Cations","Grouped_OC"))

# NMDSs and plotting
set.seed(23)
relab_kos = loaded_meta$ko_relab_table[,colnames(loaded_mags$relab_table) %in% metadata$Sample]
relab_kos <- relab_kos[, match(metadata$Sample, colnames(relab_kos))]
all(metadata$Sample == colnames(relab_kos))

protein_kos = read.table('data/proteomics/metaP_norm_imp.tsv', header = T)
protein_kos = protein_kos %>% column_to_rownames('name')
rownames(protein_kos) = map_chr(rownames(protein_kos), function(x) gsub('ko:','',x))
colnames(protein_kos) = map_chr(colnames(protein_kos), function(x) gsub('Control_','',x))
colnames(protein_kos) = map_chr(colnames(protein_kos), function(x) gsub('Warming_','',x))

metadata_prot = metadata %>% filter(Sample %in% colnames(protein_kos))
protein_kos <- protein_kos[, match(metadata_prot$Sample, colnames(protein_kos))]
protein_kos <- protein_kos[, match(metadata_prot$Sample, colnames(protein_kos))]
all(metadata_prot$Sample == colnames(protein_kos))

# ============================================================================
create_nmds_with_envfit <- function(abundance_table, metadata, title_name, top_n = 10) {
  
  # NMDS ordination
  nmds <- metaMDS(
    t(abundance_table),
    distance = "bray",
    k = 2, 
    trymax = 500, 
    try = 200
  )
  
  # Extract NMDS scores for sites
  scores_df <- as.data.frame(scores(nmds, display = "sites"))
  scores_df$Sample <- rownames(scores_df)
  scores_df$Temperature <- metadata$Temp_soil
  
  # PERMANOVA
  perm <- adonis2(
    t(abundance_table) ~ Temp_soil_scaled,
    data = metadata,
    method = "bray",
    permutations = 999
  )
  
  R2  <- perm$R2[1]
  p   <- perm$`Pr(>F)`[1]
  Fv  <- perm$F[1]
  
  # envfit: fit each feature (row) as environmental variable
  envfit_result <- envfit(
    nmds,
    env = t(abundance_table),
    perm = 999,
    na.rm = TRUE
  )
  
  # Extract envfit vectors and p-values
  envfit_vectors <- as.data.frame(scores(envfit_result, display = "vectors"))
  envfit_vectors$pval <- envfit_result$vectors$pvals
  envfit_vectors$r_squared <- envfit_result$vectors$r
  envfit_vectors$feature <- rownames(envfit_vectors)
  
  # Filter for significant features first, then sort by effect size
  sig_envfit <- envfit_vectors[envfit_vectors$pval < 0.05, ]
  sig_envfit <- sig_envfit[order(sig_envfit$r_squared, decreasing = TRUE), ]
  top_envfit <- sig_envfit[1:min(top_n, nrow(sig_envfit)), ]
  
  # Annotation text with PERMANOVA stats
  annot_text <- paste0(
    "Stress = ", round(nmds$stress, 3), "\n",
    "PERMANOVA:\n",
    "F = ", round(Fv, 2), 
    ", R² = ", round(R2, 3), 
    ", p = ", signif(p, 3), "\n",
    "Top ", nrow(top_envfit), " features shown"
  )
  
  return(list(
    nmds = nmds,
    scores = scores_df,
    envfit = envfit_result,
    envfit_vectors = envfit_vectors,
    top_envfit = top_envfit,
    stats = list(stress = nmds$stress, R2 = R2, p = p, Fv = Fv),
    annot_text = annot_text
  ))
}

plot_nmds_with_envfit <- function(nmds_result, title, arrow_scaling) {
  
  scores_df <- nmds_result$scores
  top_envfit <- nmds_result$top_envfit
  annot_text <- nmds_result$annot_text
  
  # Base NMDS plot with site points
  p <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = Temperature), size = 3, alpha = 0.7) +
    scale_color_viridis_c(name = "Temperature (°C)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right"
    ) +
    labs(title = title) +
    coord_fixed(ratio = 1)
  
  # Add envfit arrows and labels for top features
  if (nrow(top_envfit) > 0) {
    # Scale arrows for visibility (multiply by factor for readability)
    arrow_scale <- 0.8
    top_envfit$NMDS1_scaled <- top_envfit$NMDS1 * arrow_scale
    top_envfit$NMDS2_scaled <- top_envfit$NMDS2 * arrow_scale
    
    # Add arrows
    p <- p +
      geom_segment(
        data = top_envfit,
        aes(x = 0, y = 0, xend = NMDS1_scaled * arrow_scaling, yend = NMDS2_scaled * arrow_scaling),
        arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
        color = "darkred",
        alpha = 0.6,
        size = 0.8
      )
    
    # Add labels with repulsion to avoid overlap
    p <- p +
      geom_text_repel(
        data = top_envfit,
        aes(x = NMDS1_scaled * arrow_scaling, y = NMDS2_scaled * arrow_scaling, label = feature),
        size = 3,
        color = "darkred",
        box.padding = 0.3,
        point.padding = 0.2,
        max.overlaps = Inf
      )
  }
  
  # Add annotation box with stats
  #p <- p +
  #  annotate(
  #    "text",
  #    x = Inf, y = Inf,
  #    label = annot_text,
  #    hjust = 1.05, vjust = 1.05,
  #    size = 3,
  #    family = "mono",
  #    bbox = list(boxcolors = "white", fill = "white", alpha = 0.8)
  #  )
  
  return(p)
}

# ============================================================================
# NMDS
# ============================================================================

# KO analysis
ko_result <- create_nmds_with_envfit(
  log(relab_kos + 1), 
  metadata, 
  "KOs"
)

# Protein analysis
protein_result <- create_nmds_with_envfit(
  protein_kos, 
  metadata_prot, 
  "Proteins"
)

# Create plots
p_ko <- plot_nmds_with_envfit(ko_result, "Metagenomes", 0.15)
p_protein <- plot_nmds_with_envfit(protein_result, "Metaproteomes", 0.05)
nmds_plots = ggarrange(p_ko, p_protein, ncol = 2, labels = c('A', 'B'), common.legend = T, align = 'h')

#comparison_data = melt(as.matrix(log(relab_kos+1)))
#colnames(comparison_data) = c('KO', 'Sample', 'metaG_abundance')
#comparison_data = comparison_data %>% filter(KO %in% rownames(protein_kos)) %>% filter(Sample %in% colnames(protein_kos))
#comparison_data$metaP_normvalue = map_dbl(1:nrow(comparison_data), function(i) protein_kos %>% 
#                                            filter(row.names(protein_kos) == comparison_data$KO[i]) %>% 
#                                            select(as.character(comparison_data$Sample[i])) %>% pull())
#ggplot(comparison_data, aes(x=log(metaG_abundance), y=metaP_normvalue)) + geom_point() + geom_smooth(method='lm')

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

sink('stats/6_b_NMDS_envfits_functional.txt')
cat("\n=== KO ANALYSIS ===\n")
cat("NMDS Stress:", ko_result$stats$stress, "\n")
cat("PERMANOVA - F:", ko_result$stats$Fv, "R²:", ko_result$stats$R2, "p-value:", ko_result$stats$p, "\n")
cat("Significant features (p < 0.05):", sum(ko_result$envfit_vectors$pval < 0.05), "\n")
print(ko_result$top_envfit)

cat("\n=== PROTEIN ANALYSIS ===\n")
cat("NMDS Stress:", protein_result$stats$stress, "\n")
cat("PERMANOVA - F:", protein_result$stats$Fv, "R²:", protein_result$stats$R2, "p-value:", protein_result$stats$p, "\n")
cat("Significant features (p < 0.05):", sum(protein_result$envfit_vectors$pval < 0.05), "\n")
print(protein_result$top_envfit)
sink()



################################################################################
# differential abundance
rownames(metadata) = metadata$Sample
metadata$year_season = map_chr(1:nrow(metadata), function(i) paste0(metadata$Year[i],
                                                                    metadata$Season[i], 
                                                                    collapse = '_'))

# DA MAGs
min_val = min(loaded_mags$raw_table[loaded_mags$raw_table > 0])
mags_table = as.data.frame(lapply(loaded_mags$raw_table , function(x) round(x / min_val)))
mags_table[is.na(mags_table)] = 0
mags_table = mags_table[,colnames(mags_table) %in% metadata$Sample]
mags_table <- mags_table[, match(metadata$Sample, colnames(mags_table))]
rownames(mags_table) = loaded_mags$mag_info$MAG
all(metadata$Sample == colnames(mags_table))

out_mags = ancombc2(data = mags_table, 
                    meta_data = metadata,
                    fix_formula = "Temp_soil_scaled", rand_formula = "(1|Replicate) + (1 | Replicate/year_season)",
                    p_adj_method = 'fdr', struc_zero = FALSE, 
                    prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                    alpha = 0.05, n_cl = 5, verbose = TRUE)
table(out_mags$res$diff_Temp_soil_scaled)
table(out_mags$res$passed_ss_Temp_soil_scaled)
table(out_mags$res$passed_ss_Temp_soil_scaled, out_mags$res$diff_Temp_soil_scaled)
write.table(out_mags$res, file = 'stats/mag_diff_abundance_bysoiltemp_fdr.csv', 
            quote = F, row.names = T, sep = ',')

# DA KOs
min_val = min(loaded_meta$ko_raw_table[loaded_meta$ko_raw_table > 0])
meta_table = as.data.frame(lapply(loaded_meta$ko_raw_table , function(x) round(x / min_val)))
meta_table[is.na(meta_table)] = 0
meta_table = meta_table[,colnames(meta_table) %in% metadata$Sample]
meta_table <- meta_table[, match(metadata$Sample, colnames(meta_table))]
all(metadata$Sample == colnames(meta_table))

meta_table = meta_table[,colnames(meta_table) %in% metadata$Sample]
meta_table = meta_table[, match(metadata$Sample, colnames(meta_table))]
rownames(meta_table) = rownames(loaded_meta$ko_raw_table)
all(colnames(meta_table) == metadata$Sample)

out_kos = ancombc2(data = meta_table, 
                   meta_data = metadata,
                   fix_formula = "Temp_soil_scaled", rand_formula = "(1|Replicate) + (1 | Replicate/year_season)",
                   p_adj_method = 'fdr', struc_zero = FALSE, 
                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                   alpha = 0.05, n_cl = 5, verbose = TRUE)
table(out_kos$res$diff_Temp_soil_scaled)
table(out_kos$res$passed_ss_Temp_soil_scaled)
table(out_kos$res$diff_Temp_soil_scaled, out_kos$res$passed_ss_Temp_soil_scaled)
write.table(out_kos$res, file = 'stats/ko_diff_abundance_bysoiltemp_fdr.csv', 
            quote = F, row.names = T, sep = ',')


# Differential abundance proteins
out_data = read.table('data/proteomics/metaP_norm_imp.tsv', header = T)
colnames(out_data) = map_chr(colnames(out_data), function(x) gsub('Control_','',x))
colnames(out_data) = map_chr(colnames(out_data), function(x) gsub('Warming_','',x))

metadata_prot = metadata[metadata$Sample %in% colnames(out_data),]
prot_data = out_data %>% column_to_rownames('name')
prot_data <- prot_data[, match(metadata_prot$Sample, colnames(prot_data))]
all(metadata_prot$Sample == colnames(prot_data))

fit_one_gam_bf_perm <- function(
    feature_values,
    feature_name,
    metadata,
    n_perm = 200
) {
  set.seed(23)
  
  df <- data.frame(
    Sample = names(feature_values),
    value = as.numeric(feature_values),
    stringsAsFactors = FALSE
  ) %>%
    left_join(metadata, by = "Sample") %>%
    drop_na(value, Temp_soil_scaled, Replicate, Days_Since_Start)
  
  df$Replicate <- factor(df$Replicate)
  
  # Full model (with temperature)
  m_full <- gam(
    value ~ s(Days_Since_Start, bs = "ts", k = 5) +
      s(Replicate, bs = "re") +
      Temp_soil_scaled,
    data = df,
    method = "REML"
  )
  
  # Null model (without temperature)
  m_null <- gam(
    value ~ s(Days_Since_Start, bs = "ts", k = 5) +
      s(Replicate, bs = "re"),
    data = df,
    method = "REML"
  )
  
  # Bayes factor for observed data
  bf_real <- performance::test_bf(m_full, m_null)$BF[2]
  
  # Extract estimate and SE for Temp_soil_scaled
  coef_table <- summary(m_full)$p.table
  estimate <- coef_table["Temp_soil_scaled", "Estimate"]
  se <- coef_table["Temp_soil_scaled", "Std. Error"]
  
  # Permutations
  bf_perm <- replicate(n_perm, {
    df$Temp_perm <- sample(df$Temp_soil_scaled)
    
    m_perm <- gam(
      value ~ s(Days_Since_Start, bs = "ts", k = 5) +
        s(Replicate, bs = "re") +
        Temp_perm,
      data = df,
      method = "REML"
    )
    
    performance::test_bf(m_perm, m_null)$BF[2]
  })
  
  # Empirical p-value: proportion of permuted BF >= observed BF
  perm_p <- mean(bf_perm >= bf_real, na.rm = TRUE)
  
  # Return results including estimate and SE
  data.frame(
    feature = feature_name,
    BF10 = bf_real,
    perm_p = perm_p,
    estimate = estimate,
    se = se
  )
}

fit_one_gam_pvalue <- function(
    feature_values,
    feature_name,
    metadata,
    n_perm = 200
) {
  set.seed(23)
  
  df <- data.frame(
    Sample = names(feature_values),
    value = as.numeric(feature_values),
    stringsAsFactors = FALSE
  ) %>%
    left_join(metadata, by = "Sample") %>%
    drop_na(value, Temp_soil_scaled, Replicate, Days_Since_Start)
  
  df$Replicate <- factor(df$Replicate)
  
  # Full model (with temperature)
  m_full <- gam(
    value ~ s(Days_Since_Start, bs = "ts", k = 5) +
      s(Replicate, bs = "re") +
      Temp_soil_scaled,
    data = df,
    method = "REML"
  )
  
  # Extract p-value and other stats for Temp_soil_scaled from full model
  coef_table <- summary(m_full)$p.table
  pval_real <- coef_table["Temp_soil_scaled", "Pr(>|t|)"]
  estimate <- coef_table["Temp_soil_scaled", "Estimate"]
  se <- coef_table["Temp_soil_scaled", "Std. Error"]
  t_stat <- coef_table["Temp_soil_scaled", "t value"]

  # Permutations: get p-value for temperature effect in each permuted dataset
  pval_perm <- replicate(n_perm, {
    df$Temp_perm <- sample(df$Temp_soil_scaled)
    
    m_perm <- gam(
      value ~ s(Days_Since_Start, bs = "ts", k = 5) +
        s(Replicate, bs = "re") +
        Temp_perm,
      data = df,
      method = "REML"
    )
    
    coef_table_perm <- summary(m_perm)$p.table
    coef_table_perm["Temp_perm", "Pr(>|t|)"]
  })
  
  # Empirical p-value: proportion of permuted p-values <= observed p-value
  perm_p <- mean(pval_perm <= pval_real, na.rm = TRUE)
  
  # Return results
  data.frame(
    feature = feature_name,
    pval = pval_real,
    perm_p = perm_p,
    estimate = estimate,
    se = se,
    t_stat = t_stat
  )
}

# Apply to all rows
results = data.frame()
for (i in 1:nrow(prot_data)){
  res_i = fit_one_gam_pvalue(prot_data[i,],rownames(prot_data)[i],metadata_prot)
  results = rbind(results, res_i)
}

results %>% arrange(-pval)
results$padj <- p.adjust(results$pval, method = "fdr")

results %>% filter(padj < 0.05)
write.table(results, file = 'stats/proteins_diff_abundance_bysoiltemp_fdr.csv', 
            quote = F, row.names = T, sep = ',')


# Plotting diff ab
results_kos = out_kos$res
results_prots = results

results_kos = read.table('stats/ko_diff_abundance_bysoiltemp_fdr.csv', sep = ',')
results_prots = read.table('stats/proteins_diff_abundance_bysoiltemp_fdr.csv', sep = ',')


library(ggplot2)
library(ggpubr)

# Create volcano plot
volcano_plot_kos <- results_kos %>%
  mutate(
    neg_log_pval = -log10(q_Temp_soil_scaled),
    significant = (passed_ss_Temp_soil_scaled  == TRUE) & (diff_Temp_soil_scaled == TRUE) & (abs(lfc_Temp_soil_scaled) > 0.5),
    label = taxon
  ) %>%
  ggplot(aes(x = lfc_Temp_soil_scaled, y = neg_log_pval, color = significant, label = label)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("FALSE" = "darkgrey", "TRUE" = "darkred"),
    guide = "none"
  ) +
  geom_text_repel(
    data = . %>% filter(significant),
    color = "darkred",
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Metagenomes",
    x = "Log Fold Change (LFC)",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


volcano_plot_prots <- results_prots %>%
  mutate(
    neg_log_pval = -log10(pval),
    significant = (padj < 0.05),
    label = feature
  ) %>%
  ggplot(aes(x = estimate, y = neg_log_pval, color = significant, label = label)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("FALSE" = "darkgrey", "TRUE" = "darkred"),
    guide = "none"
  ) +
  geom_text_repel(
    data = . %>% filter(significant),
    color = "darkred",
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Metaproteomes",
    x = "Estimate effect of Soil Temperature",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

volcanos = ggarrange(volcano_plot_kos, volcano_plot_prots, labels = c('C', 'D'), common.legend = T, align = 'h')

full_p = ggarrange(nmds_plots, volcanos, nrow = 2, ncol = 1)
ggsave('figures/Figure5.pdf', full_p, width = 10, height = 8)

