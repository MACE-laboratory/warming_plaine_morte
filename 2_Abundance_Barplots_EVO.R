################################################################################
# 2_Abundance_Barplots_EVO.R
#
# Stacked barplots of the 10 most abundant classes / genera (rest = "Others"),
# by Condition and Year, for the 16S amplicons and the MAGs.
#
# Inputs : data/16S/*, data/mags/*, data/SiteL23_24_env_variables_final.csv
# Outputs: data/mags/mags_{genus,class}.csv, data/16S/16s_{genus,class}.csv
#            (aggregated tables, also read by 5_Differential_abundance_MAGs_EVO.R)
#          figures/SFig_11a-d  (.pdf + .png)
#          stats/02_barplot_mean_relab_*.csv  (mean relative abundance per group)
################################################################################
source('0_functions_and_packages.R')   # run from the repository root

#####################
### Load the data ###
#####################
loaded_mags = load_mags_data()
loaded_meta = load_meta_data()
loaded_16sr = load_16sr_data()

mags_genus = mags_to_genus(loaded_mags)
a16s_genus = a16s_to_genus(loaded_16sr)
mags_class = mags_to_class(loaded_mags)
a16s_class = a16s_to_class(loaded_16sr)
mags_order = mags_to_order(loaded_mags)

# Aggregate to class / genus and cache the tables (re-read below, and by script 5)
write.table(mags_genus, file = 'data/mags/mags_genus.csv', quote = F, row.names = T, sep = ',')
write.table(a16s_genus, file = 'data/16S/16s_genus.csv', quote = F, row.names = T, sep = ',')
write.table(mags_class, file = 'data/mags/mags_class.csv', quote = F, row.names = T, sep = ',')
write.table(a16s_class, file = 'data/16S/16s_class.csv', quote = F, row.names = T, sep = ',')

Class_16S = read.csv("data/16S/16s_class.csv", row.names=1, header=TRUE)
Genus_16S = read.csv("data/16S/16s_genus.csv", row.names=1, header=TRUE)
Class_MAG = read.csv("data/mags/mags_class.csv", row.names=1, header=TRUE)
Genus_MAG = read.csv("data/mags/mags_genus.csv", row.names=1, header=TRUE)

Metadata = read.csv("data/SiteL23_24_env_variables_final.csv", header=TRUE)
Metadata$Sample = gsub('\\.', '_', Metadata$Sample)
Metadata$Sample = map_chr(Metadata$Sample, function(x) paste0('L_', x))

custom_palette <- c(  "#0072B2", "#E69F00", "#009E73", "#CC3333", "#4EAEE8",
  "#9A7DCC", "#D55E00", "#009CAD", "#FF7F6E", "#7E9A3C", "#999999")
inverted_custom_palette = rev(custom_palette)


########################
### Stacked barplots ###
########################

# data    : taxa x samples table; its OBJECT NAME must start with the taxonomic level
#           (e.g. Class_16S -> 'Class'): the first 5 characters name the taxon column
# stat_name: if given, the mean relative abundance per group is saved to stats/
create_stacked_barplot_simple = function(data, metadata, group_by = "Condition",
                                          facet_by = "Year", top_n = 10, stat_name = NULL) {
  # Extract first 5 characters from dataframe name
  data_name = deparse(substitute(data))
  tax_level = substr(data_name, 1, 5)
  
  plot_data = data %>%
    as.data.frame() %>%
    rownames_to_column(tax_level) %>%
    pivot_longer(-all_of(tax_level), names_to = "Sample", values_to = "Abundance") %>%
    group_by(Sample) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    left_join(metadata, by = "Sample") %>%
    # Clean taxonomy labels
    mutate(!!sym(tax_level) := str_replace_all(!!sym(tax_level), "([pcofg])__", " \\1 ")) %>%
    mutate(!!sym(tax_level) := if_else(
      dense_rank(-ave(Abundance, !!sym(tax_level), FUN = sum)) <= top_n,
      !!sym(tax_level),
      "Others"
    )) %>%
    group_by(Sample, .data[[group_by]], .data[[facet_by]], .data[[tax_level]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Order taxa and put "Others" last
  taxa_order = plot_data %>%
    count(.data[[tax_level]], wt = Abundance, sort = TRUE) %>%
    pull(.data[[tax_level]])
  taxa_order = c(taxa_order[taxa_order != "Others"], "Others")
  # Convert to factor with "Others" last
  plot_data[[tax_level]] = factor(plot_data[[tax_level]], levels = rev(taxa_order))
  
  # Print relative abundance for each group in each condition
  print_abundance = plot_data %>%
    group_by(.data[[group_by]], .data[[facet_by]], .data[[tax_level]]) %>%
    summarise(Abundance = mean(Abundance), .groups = "drop") %>%
    pivot_wider(names_from = c(group_by, facet_by), values_from = Abundance)
  print(print_abundance)
  if (!is.null(stat_name)) save_stat_table(print_abundance, stat_name)

  ggplot(plot_data, aes(x = .data[[group_by]], y = Abundance, fill = .data[[tax_level]])) +
    geom_col(position = "fill") +
    facet_wrap(~ .data[[facet_by]]) +
    scale_y_continuous(breaks = breaks_pretty(n = 10), 
                       expand = expansion(mult = c(0, 0))) +
    labs(x = group_by, y = "Relative Abundance", fill = tax_level) +
    scale_fill_manual(values = inverted_custom_palette) +
    theme_minimal()
}

# One call per dataset; each is saved as a supplementary figure panel (SFig_11a-d)
p_class_16s = create_stacked_barplot_simple(Class_16S, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10,
                                            stat_name = "02_barplot_mean_relab_class_16s")
save_fig(p_class_16s, 'sfig_barplot_class_16s', width = 8, height = 5)

p_genus_16s = create_stacked_barplot_simple(Genus_16S, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10,
                                            stat_name = "02_barplot_mean_relab_genus_16s")
save_fig(p_genus_16s, 'sfig_barplot_genus_16s', width = 8, height = 5)

p_class_mag = create_stacked_barplot_simple(Class_MAG, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10,
                                            stat_name = "02_barplot_mean_relab_class_mag")
save_fig(p_class_mag, 'sfig_barplot_class_mag', width = 8, height = 5)

p_genus_mag = create_stacked_barplot_simple(Genus_MAG, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10,
                                            stat_name = "02_barplot_mean_relab_genus_mag")
save_fig(p_genus_mag, 'sfig_barplot_genus_mag', width = 8, height = 5)

save_session_info("02")
