library(tidyverse)
library(paletteer)
library(scales)
library(phyloseq)
library(compositions)
library(microbiome)
library(ANCOMBC)

#####################
### Load the data ###
#####################

Class_16S = read.csv("16s_class.csv", row.names=1, header=TRUE)
Genus_16S = read.csv("16s_genus.csv", row.names=1, header=TRUE)
Class_MAG = read.csv("mags_class.csv", row.names=1, header=TRUE)
Genus_MAG = read.csv("mags_genus.csv", row.names=1, header=TRUE)

Metadata = read.csv("sampleID_metadata.csv", header=TRUE)

custom_palette <- c(  "#0072B2", "#E69F00", "#009E73", "#CC3333", "#4EAEE8",
  "#9A7DCC", "#D55E00", "#009CAD", "#FF7F6E", "#7E9A3C", "#999999")
inverted_custom_palette = rev(custom_palette)


########################
### Stacked barplots ###
########################

create_stacked_barplot_simple = function(data, metadata, group_by = "Condition",
                                          facet_by = "Year", top_n = 10) {
  # Extract first 5 characters from dataframe name
  data_name = deparse(substitute(data))
  tax_level = substr(data_name, 1, 5)
  
  plot_data = data %>%
    as.data.frame() %>%
    rownames_to_column(tax_level) %>%
    pivot_longer(-all_of(tax_level), names_to = "Sample_ID", values_to = "Abundance") %>%
    group_by(Sample_ID) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    left_join(metadata, by = "Sample_ID") %>%
    # Clean taxonomy labels
    mutate(!!sym(tax_level) := str_replace_all(!!sym(tax_level), "([pcofg])__", " \\1 ")) %>%
    mutate(!!sym(tax_level) := if_else(
      dense_rank(-ave(Abundance, !!sym(tax_level), FUN = sum)) <= top_n,
      !!sym(tax_level),
      "Others"
    )) %>%
    group_by(Sample_ID, .data[[group_by]], .data[[facet_by]], .data[[tax_level]]) %>%
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
  
  ggplot(plot_data, aes(x = .data[[group_by]], y = Abundance, fill = .data[[tax_level]])) +
    geom_col(position = "fill") +
    facet_wrap(~ .data[[facet_by]]) +
    scale_y_continuous(breaks = breaks_pretty(n = 10), 
                       expand = expansion(mult = c(0, 0))) +
    labs(x = group_by, y = "Relative Abundance", fill = tax_level) +
    scale_fill_manual(values = inverted_custom_palette) +
    theme_minimal()
}

create_stacked_barplot_simple(Class_16S, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10)
ggsave("Figure2_Class_16S.png", width = 8, height = 5, bg = "white")

create_stacked_barplot_simple(Genus_16S, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10)
ggsave("Figure2_Genus_16S.png", width = 8, height = 5, bg = "white")

create_stacked_barplot_simple(Class_MAG, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10)
ggsave("Figure2_Class_MAG.png", width = 8, height = 5, bg = "white")

create_stacked_barplot_simple(Genus_MAG, Metadata, group_by = "Condition", facet_by = "Year", top_n = 10)
ggsave("Figure2_Genus_MAG.png", width = 8, height = 5, bg = "white")


###########################################################
### Statistical tests for changes in relative abundance ###
###########################################################

# Massimo Bourquin Script (slightly edited to fit purpose)
load_mags_data <- function(){
  ab_tab = read.csv('mags_abundances.tsv', sep='\t')
  tx_tab = read.csv('gtdbtk.bac120.summary.tsv', sep='\t')
  ql_tab = read.csv('quality_report.tsv', sep='\t')
  ql_tab = ql_tab %>% filter(Name %in% ab_tab$Genome)
  
  rownames(ab_tab) = ab_tab$Genome
  ab_tab$Genome = NULL
  
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
  
  relab_tab_with_scalar = relab_tab
  relab_tab_with_scalar[relab_tab_with_scalar == 0] <- min(relab_tab_with_scalar[relab_tab_with_scalar > 0])
  clr_tab = clr(relab_tab_with_scalar)
  
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

# Load metadata function
load_metadata <- function(){
  metadata_full = read.csv('SiteL23_24_env_variables_final.csv')
  metadata_full$Sample = paste0('L_', metadata_full$Sample)
  metadata_full$Sample = gsub('\\.','_',metadata_full$Sample)
  return(metadata_full)
}
# Load MAG and metadata
loaded_mags = load_mags_data()
metadata_full = load_metadata()

# Prepare MAG data
mags_table = loaded_mags$raw_table
minval = min(mags_table[mags_table > 0])
mags_table = mags_table / minval
mags_table[] <- lapply(mags_table, as.integer)
mags_table <- mags_table[, metadata_full$Sample]

# Prepare Class data
mags_tax_mat = loaded_mags$mag_info[,c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
rownames(mags_tax_mat) = loaded_mags$mag_info$MAG

rownames(metadata_full) = metadata_full$Sample
metadata_full$year_season = paste0(metadata_full$Year, metadata_full$Season)
metadata_full$condition_year = paste0(metadata_full$Condition, '_', metadata_full$Year)
META = sample_data(metadata_full)

mags_table <- mags_table[ , metadata_full$Sample]

OTU = otu_table(mags_table, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(mags_tax_mat))
physeq = phyloseq(OTU, TAX)
physeq = merge_phyloseq(physeq, META)

print('Differential abundance of Class...')
out_class = ancombc2(data = physeq, tax_level = "Class",
                      group = 'condition_year',
                      fix_formula = "condition_year", rand_formula = "(1 | Replicate) + (1 | Replicate:year_season)",
                      p_adj_method = "BH",
                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                      struc_zero = FALSE, neg_lb = TRUE,
                      alpha = 0.05, n_cl = 2, verbose = TRUE)

res_class = out_class$res
write.csv(res_class, file = "ANCOMBC2_Class_results_EVO.csv")


# Repeat test using different condition_year groups as reference.

# Ensure that condition_year is a factor
metadata_full$condition_year <- as.factor(metadata_full$condition_year)

# Check current factor levels
levels(metadata_full$condition_year)

# Reorder factor levels
warming23_ref_level <- c("warming_2023", "control_2023", "warming_2024", "control_2024")
control24_ref_level <- c("control_2024", "control_2023", "warming_2024", "warming_2023")

############################################
# Set new reference level to "warming_2023"#
############################################

metadata_full$condition_year <- factor(metadata_full$condition_year, levels = warming23_ref_level)

# Update the phyloseq object and rerun ANCOMBC2
META = sample_data(metadata_full)
physeq = phyloseq(OTU, TAX, META)

# Check the updated factor levels within the phyloseq object
levels(sample_data(physeq)$condition_year)

print('Differential abundance of Class with warming_2023 as reference...')
out_class_warming23_ref = ancombc2(data = physeq, tax_level = "Class",
                      group = 'condition_year',
                      fix_formula = "condition_year", rand_formula = "(1 | Replicate) + (1 | Replicate:year_season)",
                      p_adj_method = "BH",
                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                      struc_zero = FALSE, neg_lb = TRUE,
                      alpha = 0.05, n_cl = 2, verbose = TRUE)

res_class_warming23 = out_class_warming23_ref$res
write.csv(res_class_warming23, file = "ANCOMBC2_Class_results_warming23_ref_EVO.csv")

############################################
# Set new reference level to "control_2024"#
############################################

metadata_full$condition_year <- factor(metadata_full$condition_year, levels = control24_ref_level)

# Update the phyloseq object and rerun ANCOMBC2
META = sample_data(metadata_full)
physeq = phyloseq(OTU, TAX, META)

# Check the updated factor levels within the phyloseq object
levels(sample_data(physeq)$condition_year)

print('Differential abundance of Class with control_2024 as reference...')
out_class_control24_ref = ancombc2(data = physeq, tax_level = "Class",
                                   group = 'condition_year',
                                   fix_formula = "condition_year", rand_formula = "(1 | Replicate) + (1 | Replicate:year_season)",
                                   p_adj_method = "BH",
                                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                   struc_zero = FALSE, neg_lb = TRUE,
                                   alpha = 0.05, n_cl = 2, verbose = TRUE)

res_class_control24 = out_class_control24_ref$res
write.csv(res_class_control24, file = "ANCOMBC2_Class_results_control24_ref_EVO.csv")
