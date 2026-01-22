library(ANCOMBC)
library(tidyverse)
library(paletteer)
library(ggridges)
library(ggpubr)
library(ggthemes)

#####################
### Load the data ###
#####################

# Differential abundance
diff_ab_MAG_FDR_cond = read.csv("mag_diff_abundance_bycondition_fdr.csv", header = TRUE)
diff_ab_MAG_FDR_temp = read.csv("mag_diff_abundance_bysoiltemp_fdr.csv", header = TRUE)
# diff_ab files below aren't necessary if we only want to plot a ridgeline of the top10 classes
diff_ab_CLASS_FDR_cond = read.csv("classes_diff_abundance_bycondition_fdr.csv", header = TRUE)
diff_ab_CLASS_FDR_temp = read.csv("classes_diff_abundance_bysoiltemp_fdr.csv", header = TRUE)

diff_ab_CLASS_FDR_cond_named <- diff_ab_CLASS_FDR_cond
diff_ab_CLASS_FDR_cond_named$taxon <- gsub("c__", "", diff_ab_CLASS_FDR_cond_named$taxon)
diff_ab_CLASS_FDR_temp_named <- diff_ab_CLASS_FDR_temp
diff_ab_CLASS_FDR_temp_named$taxon <- gsub("c__", "", diff_ab_CLASS_FDR_temp$taxon)

# Taxonomic information
mag_taxonomy = read_tsv("gtdbtk.bac120.summary_MAG_ID.tsv")

# Class abundances from MAG data
class_ab = read.csv("mags_class_abundance.csv", row.names = 1, header = TRUE)

# Custom colours
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
    select(user_genome, MAG_ID, classification) %>%
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
        genus != "g unknown" ~ str_c(class, "; ", order, "; ", family, "; ", genus, "; ", MAG_ID),
        family != "f unknown" ~ str_c(class, "; ", order, "; ", family, "; g unknown; ", MAG_ID),
        order != "o unknown" ~ str_c(class, "; ",  order,  "; f unknown; g unknown; ", MAG_ID),
        class != "c unknown" ~ str_c(class, "; o unknown; f unknown; g unknown; ", MAG_ID),
        phylum != "p unknown" ~ str_c(phylum, "; c unknown; o unknown; f unknown; g unknown; ", MAG_ID),
        domain != "d unknown" ~ str_c(domain, "; p unknown; c unknown; o unknown; f unknown; g unknown; ", MAG_ID)
      )
    ) %>%
    select(user_genome, MAG_ID_taxonomy)
  
  diff_ab_data <- diff_ab_data %>%
    left_join(taxonomy_info, by = c("taxon" = "user_genome")) %>%
    relocate(MAG_ID_taxonomy, .before = taxon)
}

diff_ab_MAG_FDR_cond_named <- add_clean_taxonomy(diff_ab_MAG_FDR_cond, mag_taxonomy)
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
  write_csv(stats_summary, paste0("ridgeline_", lfc_col, "_statistics_EVO.csv"))
  
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

# Generate the ridgeline plot
ridgeline_class_FDR_temp <- plot_ridgeline_abundance(diff_ab_MAG_FDR_temp_named, mag_taxonomy, class_ab, n_cutoff = 10)
print(ridgeline_class_FDR_temp)
ggsave(ridgeline_class_FDR_temp, filename = "ridgeline_class_FDR_soiltemp_EVO.png", width = 6, height = 8, dpi = 300)

ridgeline_class_FDR_cond <- plot_ridgeline_abundance(diff_ab_MAG_FDR_cond_named, mag_taxonomy, class_ab, n_cutoff = 10)
print(ridgeline_class_FDR_cond)
ggsave(ridgeline_class_FDR_cond, filename = "ridgeline_class_FDR_condition_EVO.png", width = 6, height = 8, dpi = 300)

##############################################
### MAG Differential Abundance Forest Plot ###
##############################################
plot_lfc_forest <- function(diff_ab_data, colour_mapping, n_cutoff = 10) {
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
  write_csv(plot_data %>%
              select(new_ID, !!sym(lfc_col), !!sym(se_col), !!sym(q_col), significant),
            paste0("forest_plot_", lfc_col, "_values_EVO.csv"))
  
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
lfc_forest_MAG_FDR_condition <- plot_lfc_forest(diff_ab_MAG_FDR_cond_named, colour_mapping)
print(lfc_forest_MAG_FDR_condition)
ggsave(lfc_forest_MAG_FDR_condition, filename = "forest_plot_MAG_FDR_condition_EVO.png", width = 8, height = 6, dpi = 300)

forest_plot_MAG_FDR_temp <- plot_lfc_forest(diff_ab_MAG_FDR_temp_named, colour_mapping, n_cutoff = 10)
print(forest_plot_MAG_FDR_temp)
ggsave(forest_plot_MAG_FDR_temp, filename = "forest_plot_MAG_FDR_soiltemp_EVO.png", width = 8, height = 10, dpi = 300)

forest_plot_CLASS_FDR_condition <- plot_lfc_forest(diff_ab_CLASS_FDR_cond_named, colour_mapping, n_cutoff = 10)
print(forest_plot_CLASS_FDR_condition)
ggsave(forest_plot_CLASS_FDR_condition, filename = "forest_plot_CLASS_FDR_condition_EVO.png", width = 8, height = 6, dpi = 300)

forest_plot_CLASS_FDR_temp <- plot_lfc_forest(diff_ab_CLASS_FDR_temp_named, colour_mapping, n_cutoff = 10)
print(forest_plot_CLASS_FDR_temp)
ggsave(forest_plot_CLASS_FDR_temp, filename = "forest_plot_CLASS_FDR_soiltemp_EVO.png", width = 8, height = 10, dpi = 300)

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
  theme(legend.position = c(0.85,
                            #0.745,
    0.142), axis.text.y = element_text(size = 8))
forest_plot_MAG_temp_significant
ggsave(forest_plot_MAG_temp_significant, 
       filename = "forest_plot_MAG_FDR_soiltemp_significant_EVO.png", 
       width = 8, height = 6, dpi = 300)

######################################
### Combine plots for final figure ###
######################################

Figure3_diff_ab_CLASS_MAGs_EVO <- ggarrange(ridgeline_class_FDR_temp, forest_plot_MAG_temp_significant, 
                                            labels = c("A", "B"), font.label = "bold",
                                            nrow = 2, widths = c(1, 1))
print(Figure3_diff_ab_CLASS_MAGs_EVO)
ggsave("Figure3_diff_ab_EVO.png", Figure3_diff_ab_CLASS_MAGs_EVO, width = 10, height = 12)

Supp_Fig_diff_ab_CLASS_MAGs_condition_EVO <- ggarrange(ridgeline_class_FDR_cond, lfc_forest_MAG_FDR_condition, 
                                                    labels = c("A", "B"), font.label = "bold",
                                                    nrow = 2, widths = c(1, 1))
print(Supp_Fig_diff_ab_CLASS_MAGs_condition_EVO)
ggsave("Supplementary_Figure7_diff_ab_CLASS_MAGs_condition_EVO.png", Supp_Fig_diff_ab_CLASS_MAGs_condition_EVO, width = 10, height = 12)
