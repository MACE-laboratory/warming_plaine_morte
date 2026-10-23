################################################################################
# 5_Differential_abundance_MAGs_EVO.R
#
# Differential abundance of MAGs, classes and genera with ANCOM-BC2, by
#   (a) soil temperature (continuous, scaled)   and   (b) Condition (warming vs control)
# Random effects: (1 | Replicate) + (1 | Replicate/year_season)
# Then: forest plots + class-level ridgeline plots of the log-fold changes.
#
# Inputs : data/mags/*  (+ data/mags/mags_class.csv written by 2_Abundance_Barplots_EVO.R)
# Outputs: figures/Fig_4, SFig_12, SFig_13   (.pdf + .png)
#          stats/05_da_<mags|classes|genera>_<soiltemp|condition>_ancombc2.csv
#          stats/05_ridgeline_*.csv, stats/05_forest_*.csv
#
# The ANCOM-BC2 part is slow (hours). The plotting part ("Plots") only needs the
# stats/05_da_*_soiltemp_ancombc2.csv files, so it can be re-run on its own.
################################################################################
source('0_functions_and_packages.R')   # run from the repository root
library(ggridges)

N_CORES <- 5   # parallel workers for ancombc2()

loaded_mags = load_mags_data()
metadata = prepare_metadata()
rownames(metadata) = metadata$Sample
# nesting factor for the random effects
metadata$year_season = paste0(metadata$Year, metadata$Season)

################################################################################
# ANCOM-BC2
################################################################################

# MAG pseudo-count table: CoverM trimmed means / smallest non-zero value, rounded
# to integers, columns ordered as in the metadata.
min_val = min(loaded_mags$raw_table[loaded_mags$raw_table > 0])
mags_table = as.data.frame(lapply(loaded_mags$raw_table , function(x) round(x / min_val)))
mags_table[is.na(mags_table)] = 0
mags_table = mags_table[,colnames(mags_table) %in% metadata$Sample]
mags_table <- mags_table[, match(metadata$Sample, colnames(mags_table))]
rownames(mags_table) = loaded_mags$mag_info$MAG
stopifnot(all(metadata$Sample == colnames(mags_table)))

# phyloseq object (counts + taxonomy + metadata) for the class / genus-level runs
mags_tax_mat = loaded_mags$mag_info[,c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
rownames(mags_tax_mat) = loaded_mags$mag_info$MAG
physeq = phyloseq(otu_table(mags_table, taxa_are_rows = TRUE), tax_table(as.matrix(mags_tax_mat)))
physeq = merge_phyloseq(physeq, sample_data(metadata))

# Run ANCOM-BC2 for one fixed effect, print the diff / sensitivity-score tables,
# and save the full result table. `level`: NULL = MAGs, or 'Class' / 'Genus'.
run_ancombc2 <- function(variable, level_name, level = NULL, iter_control = NULL) {
  args <- list(fix_formula = variable,
               rand_formula = "(1|Replicate) + (1 | Replicate/year_season)",
               p_adj_method = 'fdr', struc_zero = FALSE,
               prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
               alpha = 0.05, n_cl = N_CORES)
  if (is.null(level)) {
    args <- c(list(data = mags_table, meta_data = metadata), args)
  } else {
    args <- c(list(data = physeq, tax_level = level), args)
  }
  if (is.null(iter_control)) args$verbose <- TRUE else args$iter_control <- iter_control

  out <- do.call(ancombc2, args)

  diff_col <- paste0("diff_", variable); ss_col <- paste0("passed_ss_", variable)
  print(table(out$res[[diff_col]]))
  print(table(out$res[[ss_col]]))
  print(table(out$res[[ss_col]], out$res[[diff_col]]))
  save_stat_table(out$res, paste0("05_da_", level_name, "_",
                                  ifelse(variable == "Temp_soil_scaled", "soiltemp", "condition"),
                                  "_ancombc2"))
  out
}

# (a) soil temperature
out_mags_temp    <- run_ancombc2("Temp_soil_scaled", "mags")
out_class_temp   <- run_ancombc2("Temp_soil_scaled", "classes", level = "Class")
out_genera_temp  <- run_ancombc2("Temp_soil_scaled", "genera",  level = "Genus")

# (b) Condition (warming vs control); looser convergence tolerance for the class /
#     genus models, as in the original analysis
out_mags_cond    <- run_ancombc2("Condition", "mags",
                                 iter_control = list(tol = 1e-5, max_iter = 50, verbose = TRUE))
out_class_cond   <- run_ancombc2("Condition", "classes", level = "Class",
                                 iter_control = list(tol = 0.005, max_iter = 50, verbose = TRUE))
out_genera_cond  <- run_ancombc2("Condition", "genera",  level = "Genus",
                                 iter_control = list(tol = 0.005, max_iter = 50, verbose = TRUE))


################################################################################
# Plots
################################################################################
# ANCOM-BC2 results for soil temperature (re-read from stats/ so this part can
# be run without repeating the models above)
diff_ab_MAG_FDR_temp = read.csv("stats/05_da_mags_soiltemp_ancombc2.csv", header = TRUE)
diff_ab_CLASS_FDR_temp = read.csv("stats/05_da_classes_soiltemp_ancombc2.csv", header = TRUE)

diff_ab_CLASS_FDR_temp_named <- diff_ab_CLASS_FDR_temp
diff_ab_CLASS_FDR_temp_named$taxon <- gsub("c__", "", diff_ab_CLASS_FDR_temp$taxon)

# Taxonomic information
mag_taxonomy = read_tsv("data/mags/gtdbtk.bac120.summary.tsv")

# Class abundances from MAG data
class_ab = read.csv("data/mags/mags_class.csv", row.names=1, header=TRUE)


# Custom colours (one per top class + "Others")
custom_palette <- c("#0072B2", "#E69F00", "#009E73", "#CC3333", "#4EAEE8",
                    "#9A7DCC", "#D55E00", "#009CAD", "#FF7F6E", "#7E9A3C", "#999999")

#########################
### General Use Items ###
#########################

# Clean taxonomy data at class level
taxonomy_prepared <- mag_taxonomy %>%
  mutate(
    class = str_extract(classification, "c__[^;]+") %>%
      str_replace("c__", "") %>%
      str_trim() %>%
      replace_na("unknown")
  ) %>%
  select(user_genome, class)

# Top10 classes by mean relative abundance
top_taxa <- class_ab %>%
  as.data.frame() %>%
  rownames_to_column("class_names") %>%
  pivot_longer(cols = -class_names, names_to = "Sample_ID", values_to = "Abundance") %>%
  group_by(Sample_ID) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  group_by(class_names) %>%
  summarise(mean_relative_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(
    class_clean = str_extract(class_names, "c__[^;]+") %>%
      str_replace("c__", "") %>%
      str_trim() %>%
      replace_na("unknown")) %>%
  slice_max(order_by = mean_relative_abundance, n = 10) %>%
  pull(class_clean)

# Associate colour to top taxa
taxa_order <- c(top_taxa, "Others")
colour_mapping <- setNames(custom_palette[1:length(taxa_order)], taxa_order)

########################################################################
### Add taxonomy data to MAG IDs from differential abundance results ###
########################################################################
add_clean_taxonomy <- function(diff_ab_data, taxonomy_data) {
  taxonomy_info <- taxonomy_data %>%
    select(user_genome, classification) %>%
    mutate(
      # Extract and format each taxonomic level
      domain = str_extract(classification, "d__[^;]+") %>% str_replace("d__", "d ") %>% str_replace_na("d unknown"),
      phylum = str_extract(classification, "p__[^;]+") %>% str_replace("p__", "p ") %>% str_replace_na("p unknown"),
      class = str_extract(classification, "c__[^;]+") %>% str_replace("c__", "c ") %>% str_replace_na("c unknown"),
      order = str_extract(classification, "o__[^;]+") %>% str_replace("o__", "o ") %>% str_replace_na("o unknown"),
      family = str_extract(classification, "f__[^;]+") %>% str_replace("f__", "f ") %>% str_replace_na("f unknown"),
      genus = str_extract(classification, "g__[^;]+") %>% str_replace("g__", "g ") %>% str_replace_na("g unknown"),
      
      # Build the taxonomy string based on the last known taxonomic level
      MAG_ID_taxonomy = case_when(
        genus != "g unknown" ~ str_c(class, "; ", order, "; ", family, "; ", genus, "; ", user_genome),
        family != "f unknown" ~ str_c(class, "; ", order, "; ", family, "; g unknown; ", user_genome),
        order != "o unknown" ~ str_c(class, "; ",  order,  "; f unknown; g unknown; ", user_genome),
        class != "c unknown" ~ str_c(class, "; o unknown; f unknown; g unknown; ", user_genome),
        phylum != "p unknown" ~ str_c(phylum, "; c unknown; o unknown; f unknown; g unknown; ", user_genome),
        domain != "d unknown" ~ str_c(domain, "; p unknown; c unknown; o unknown; f unknown; g unknown; ", user_genome)
      )
    ) %>%
    select(user_genome, MAG_ID_taxonomy)
  
  diff_ab_data <- diff_ab_data %>%
    left_join(taxonomy_info, by = c("taxon" = "user_genome")) %>%
    relocate(MAG_ID_taxonomy, .before = taxon)
}

diff_ab_MAG_FDR_temp_named <- add_clean_taxonomy(diff_ab_MAG_FDR_temp, mag_taxonomy)

###############################################
### Ridgeline plot at class level for MAGs ####
###############################################
plot_ridgeline_abundance <- function(diff_ab_data, taxonomy_data, class_ab, n_cutoff = 10) {
# Determine column name for log fold change  
  if ("lfc_Conditionwarming" %in% names(diff_ab_data)) {
    lfc_col <- "lfc_Conditionwarming"
    x_label <- "Log Fold Change - Warming"
  } else if ("lfc_Temp_soil_scaled" %in% names(diff_ab_data)) {
    lfc_col <- "lfc_Temp_soil_scaled"
    x_label <- "Log Fold Change - Soil Temperature"
  } else {
    stop("Relevant LFC column name not found")
  }
  # Merge taxonomy with differential abundance data
  diff_ab_data <- diff_ab_data %>%
    left_join(taxonomy_prepared, by = c("taxon" = "user_genome")) 
  
  # Create display categories based on top abundance taxa
  diff_ab_data <- diff_ab_data %>%
    mutate(display_taxon = if_else(class %in% top_taxa, class, "Others"))
  
  # Order taxa by abundance (most abundant at top) with Others at bottom
  taxa_order <- c(top_taxa, "Others")
  diff_ab_data <- diff_ab_data %>%
    mutate(display_taxon = factor(display_taxon, levels = rev(taxa_order)))
  
  # Calculate median LFC and perform Wilcoxon test for each class
  stats_summary <- diff_ab_data %>%
    group_by(display_taxon) %>%
    summarise(
      median_lfc = median(!!sym(lfc_col), na.rm = TRUE),
      mean_lfc = mean(!!sym(lfc_col), na.rm = TRUE),
      sd_lfc = sd(!!sym(lfc_col), na.rm = TRUE),
      count = n(),
      p_value = wilcox.test(!!sym(lfc_col), mu = 0)$p.value,
      .groups = "drop"
    ) %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ ""))
  
  # Display and create summary statistics
  print(stats_summary)
  # per class: median / mean / sd of the MAG-level LFCs and a Wilcoxon test against 0
  save_stat_table(stats_summary, paste0("05_ridgeline_", lfc_col, "_statistics"))
  
  # Asterisk positioning
  x_range <- range(diff_ab_data[[lfc_col]], na.rm = TRUE)
  asterisk_x <- max(x_range) * 1.2
  
  # Create the ridgeline plot
  p <- ggplot() +
    # Vertical line at x = 0
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    # Ridgelines 
    geom_density_ridges(data = diff_ab_data, aes(x = !!sym(lfc_col), y = display_taxon, fill = display_taxon),
                        alpha = 0.9, scale = 1.5, rel_min_height = 0.01) +
    # Add median lines
    geom_segment(data = stats_summary, aes(x = median_lfc, xend = median_lfc, 
                                           y = display_taxon, 
                                           yend = as.numeric(display_taxon) + 0.6), 
                 linetype = "longdash", color = "black", linewidth = 0.6, alpha = 0.7,
                 inherit.aes = FALSE) +
    # Add significance asterisks
    geom_text(data = stats_summary %>% filter(significance != ""),
              aes(x = asterisk_x, y = display_taxon, label = significance),
              inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 6) +
    scale_fill_manual(values = colour_mapping) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 10), 
  expand = expansion(mult = c(0, 0.1))) +
    #    scale_x_continuous(breaks = seq(-3, 3, by = 0.5)) + #Old scaling
        scale_y_discrete(expand = expansion(add = c(0.1, 0.8))) +
    labs(x = x_label, y = "Class") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_text(size = 8))
  p
}

# Generate the ridgeline plot (panel B of Figure 4)
ridgeline_class_FDR_temp <- plot_ridgeline_abundance(diff_ab_MAG_FDR_temp_named, mag_taxonomy, class_ab, n_cutoff = 10)
print(ridgeline_class_FDR_temp)

##############################################
### MAG Differential Abundance Forest Plot ###
##############################################
# n_cutoff highest and n_cutoff lowest LFCs; `label` ('mags' / 'classes') names the
# saved table (the original script overwrote one file for both levels)
plot_lfc_forest <- function(diff_ab_data, colour_mapping, n_cutoff = 10, label = "mags") {
    # Determine column names for log fold change and q-value
    if ("lfc_Conditionwarming" %in% names(diff_ab_data)) {
      lfc_col <- "lfc_Conditionwarming"
      q_col <- "q_Conditionwarming"
      se_col <- "se_Conditionwarming"
      x_label <- "Log Fold Change - Warming"
    } else if ("lfc_Temp_soil_scaled" %in% names(diff_ab_data)) {
      lfc_col <- "lfc_Temp_soil_scaled"
      q_col <- "q_Temp_soil_scaled"
      se_col <- "se_Temp_soil_scaled"
      x_label <- "Log Fold Change - Soil Temperature"
    } else {
      stop("Relevant LFC column name not found")
    }
  plot_data <- diff_ab_data %>%
    # Select top 10 highest LFC (already in descending order)
    slice_max(order_by = !!sym(lfc_col), n = n_cutoff) %>%
    bind_rows(
      # Select bottom 10 lowest LFC (already in ascending order)
      diff_ab_data %>%
        slice_min(order_by = !!sym(lfc_col), n = n_cutoff)
    ) %>%
    # Add significance and direction
    mutate(
      significant = ifelse(!!sym(q_col) < 0.05, TRUE, FALSE),
      direction = ifelse(!!sym(lfc_col) > 0, "Positive", "Negative")
    )
# Check if the data is MAG-level or class-level
    if ("MAG_ID_taxonomy" %in% names(plot_data)) {
      # For MAG-level data, join with taxonomy_prepared to get class information
      plot_data <- plot_data %>%
        left_join(taxonomy_prepared, by = c("taxon" = "user_genome")) %>%
        mutate(
          display_taxon = if_else(class %in% top_taxa, class, "Others"),
          display_taxon = factor(display_taxon, levels = c(top_taxa, "Others")))
    } else {
      # For class-level data, use the taxon column directly
      plot_data <- plot_data %>%
        mutate(
          display_taxon = if_else(taxon %in% top_taxa, taxon, "Others"),
          display_taxon = factor(display_taxon, levels = c(top_taxa, "Others")))
    }
    # Order (descending) and set taxa factor levels for correct y-axis order
  plot_data <- plot_data %>%  
  arrange(!!sym(lfc_col)) %>%
    mutate(
      new_ID = if ("MAG_ID_taxonomy" %in% names(.)) MAG_ID_taxonomy else taxon,
      new_ID = factor(new_ID, levels = new_ID)
    )
  
  # Display important lfc statistics for top/bottom 10 taxa
  print(plot_data %>%
          select(new_ID, !!sym(lfc_col), !!sym(se_col), !!sym(q_col), significant) )
  save_stat_table(plot_data %>%
                    select(new_ID, !!sym(lfc_col), !!sym(se_col), !!sym(q_col), significant),
                  paste0("05_forest_", label, "_", lfc_col, "_values"))
  
  # Create the plot
  lfc <- ggplot(plot_data, aes(x = !!sym(lfc_col), y = new_ID, fill = display_taxon)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +
    geom_errorbar(aes(xmin = !!sym(lfc_col) - !!sym(se_col), xmax = !!sym(lfc_col) + !!sym(se_col)),
                  width = 0.3, linewidth = 0.6) +
    geom_text(aes(x = max(!!sym(lfc_col)) * 1.38, label = ifelse(significant, "*", "")),
              hjust = 0, vjust = 0.5, size = 6) +
    scale_fill_manual(values = colour_mapping) +
    labs(x = x_label, y = "MAG Taxa", fill = "Class") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_text(size = 8))
  lfc
}

#############################
### Generate Forest Plots ###
#############################
forest_plot_MAG_FDR_temp <- plot_lfc_forest(diff_ab_MAG_FDR_temp_named, colour_mapping, n_cutoff = 10, label = "mags")
print(forest_plot_MAG_FDR_temp)
save_fig(forest_plot_MAG_FDR_temp, 'sfig_da_forest_mags_top', width = 8, height = 10)

forest_plot_CLASS_FDR_temp <- plot_lfc_forest(diff_ab_CLASS_FDR_temp_named, colour_mapping, n_cutoff = 10, label = "classes")
print(forest_plot_CLASS_FDR_temp)
save_fig(forest_plot_CLASS_FDR_temp, 'sfig_da_forest_class_top', width = 8, height = 10)

############################################
### Forest plot of only significant MAGs ###
############################################
diff_ab_MAG_FDR_temp_named_sig <- diff_ab_MAG_FDR_temp_named %>%
  filter(diff_Temp_soil_scaled == TRUE, passed_ss_Temp_soil_scaled == TRUE)
diff_ab_MAG_FDR_temp_named_sig <- arrange(diff_ab_MAG_FDR_temp_named_sig, desc(lfc_Temp_soil_scaled)) %>%
  mutate(MAG_ID_taxonomy = factor(MAG_ID_taxonomy, levels = rev(MAG_ID_taxonomy)))
diff_ab_MAG_FDR_temp_named_sig <- diff_ab_MAG_FDR_temp_named_sig %>%
  left_join(taxonomy_prepared, by = c("taxon" = "user_genome")) %>%
  mutate(display_taxon = if_else(class %in% top_taxa, class, "Others"),
    display_taxon = factor(display_taxon, levels = c(top_taxa, "Others")))

# (panel A of Figure 4)
forest_plot_MAG_temp_significant <- ggplot(diff_ab_MAG_FDR_temp_named_sig,
                                           aes(x = lfc_Temp_soil_scaled,
                                               y = MAG_ID_taxonomy, 
                                                   fill = display_taxon)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5, alpha = 0.25) +
  #  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +
  geom_errorbar(aes(xmin = lfc_Temp_soil_scaled - se_Temp_soil_scaled,
                    xmax = lfc_Temp_soil_scaled + se_Temp_soil_scaled),
                width = 0.3, linewidth = 0.6) +
  geom_point(size = 5, shape = 23) +
  scale_fill_manual(values = colour_mapping) +
  labs(x = "Log Fold Change - Soil Temperature", y = "MAG Taxa", fill = "Class") +
  theme_minimal() +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.142),
        axis.text.y = element_text(size = 8))
forest_plot_MAG_temp_significant

######################################
### Combine plots for final figure ###
######################################

Figure4_diff_ab_CLASS_MAGs <- ggarrange(forest_plot_MAG_temp_significant, ridgeline_class_FDR_temp,
                                        labels = c("A", "B"), font.label = "bold",
                                        nrow = 2, widths = c(1, 1))
save_fig(Figure4_diff_ab_CLASS_MAGs, 'fig_diff_abundance_mags', width = 10, height = 12)

save_session_info("05")
