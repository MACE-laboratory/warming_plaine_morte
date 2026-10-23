setwd("/PATH/TO/FILES")
library(tidyverse)
library(KEGGREST)
library(vegan)
library(igraph)
library(ggraph)
library(paletteer)
library(ggpubr)

# Load list of protein-derived KOs
KO_proteomics <- read.delim("metaP_norm_imp.tsv", sep = "\t", header = TRUE)
KO_proteomics <- KO_proteomics %>% relocate(name)
colnames(KO_proteomics)


# Load list of all KOs found in metagenome
KO_complete_list <- read.table("annotations_mags_EDU_2025.emapper.annotations_edited",
                               sep = "\t",
                               header = TRUE,
                               quote = "",           # Don't interpret quotes
                               comment.char = "",    # Don't interpret # as comments
                               fill = TRUE,          # Fill missing values
                               check.names = FALSE)  # Keep original column names

## NA look up
# Check for any rows with NA values
print(paste("Total number of NAs:", sum(is.na(KO_complete_list))))
# See which columns have NAs
colSums(is.na(KO_complete_list))
# Find rows with any NAs
rows_with_na <- which(rowSums(is.na(KO_complete_list)) > 0)
length(rows_with_na)

# "Explode" complete list to have one row for each KO:
KO_complete_list_expanded <- KO_complete_list %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = trimws(KEGG_ko))

# Check how many times each KO appears in KO_complete_list
KO_counts <- KO_complete_list_expanded %>% 
  filter(KEGG_ko %in% KO_proteomics$name) %>% 
  count(KEGG_ko, name = "num_of_times_found") %>% 
  complete(KEGG_ko = KO_proteomics$name, fill = list(num_of_times_found = 0)) %>%
  arrange(desc(num_of_times_found))

p_KO_count_histo <- ggplot(KO_counts, aes(x = reorder(KEGG_ko, -num_of_times_found), y = num_of_times_found)) + 
  geom_col() + theme(axis.text.x = element_text(angle = -90))
p_KO_count_histo

# Subset KO_complete_list with KOs from KO_proteomics
KO_list <- KO_complete_list_expanded %>% filter(KEGG_ko %in% KO_proteomics$name)
nrow(KO_list)

# Collapse list by unique KOs
KO_list_collapsed <- KO_list %>%
  group_by(KEGG_ko) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = " | ")),
            n_occurrences = as.numeric(n()),
    .groups = "drop") %>% relocate(n_occurrences, .after = KEGG_ko)

##########################################################################

# KOs to pathways
kegg_pathway_links <- keggLink("pathway", "ko")

# Create pathway dataset
pathway_mapping <- data.frame(KO = names(kegg_pathway_links),
  pathway = sub("path:", "", kegg_pathway_links))

# Subset pathways with KOs from KO_proteomics
proteomics_pathways <- pathway_mapping %>%
  filter(KO %in% KO_proteomics$name) %>%
  filter(!grepl("map", pathway))

message("Found ", length(unique(proteomics_pathways$KO)), " of ", 
        nrow(KO_proteomics), " KOs in KEGG pathways. ",
        nrow(KO_proteomics) - length(unique(proteomics_pathways$KO)),
        " KOs without assigned pathway")

# Get pathway names
pathway_names <- data.frame(pathway = names(keggList("pathway")), 
                            pathway_description = keggList("pathway")) %>%
  mutate(pathway = sub("map","ko", pathway))

proteomics_pathways <- proteomics_pathways %>%
  left_join(pathway_names, by = "pathway")

# Add KOs without assigned pathway
proteomics_pathways <- proteomics_pathways %>%
  bind_rows(KO_proteomics %>%
              filter(!name %in% proteomics_pathways$KO) %>%
              mutate(pathway = NA, pathway_description = NA) %>%
              select(KO = name, pathway, pathway_description))

# Sanity check
length(unique(proteomics_pathways$KO)) == nrow(KO_proteomics)

###########################################################################

# Function to extract BRITE level 1 and 2 categories
get_brite_categories <- function(pathway_id) {
  
  pathway_info <- tryCatch(
    keggGet(pathway_id)[[1]],
    error = function(e) NULL
  )
  
  if (is.null(pathway_info) || is.null(pathway_info$CLASS)) {
    return(data.frame(
      pathway = pathway_id,
      brite_lvl1 = NA,
      brite_lvl2 = NA
    ))
  }
  
  class_vec <- pathway_info$CLASS
  
  # Split on ";"
  split_class <- str_split(class_vec, ";\\s*")[[1]]
  
  data.frame(
    pathway = pathway_id,
    brite_lvl1 = split_class[1],
    brite_lvl2 = ifelse(length(split_class) >= 2,
                        split_class[2],
                        NA)
  )
}

# Get unique pathways
unique_pathways <- unique(proteomics_pathways$pathway)

# Retrieve BRITE categories
brite_mapping <- map_dfr(unique_pathways, get_brite_categories)

# Add to main table
proteomics_pathways <- proteomics_pathways %>%
  left_join(brite_mapping, by = "pathway")

# Add number of times each KO is found in metagenome
proteomics_pathways <- proteomics_pathways %>%
  left_join(KO_counts, by = c("KO" = "KEGG_ko")) %>% 
  relocate(num_of_times_found, .after = KO)

# New category for plotting
proteomics_pathways <- proteomics_pathways %>%
  mutate(category = case_when(!is.na(brite_lvl2) ~ brite_lvl2,
                              !is.na(pathway_description) ~ pathway_description,
                              TRUE ~ "KOs without assigned pathways"))

# Filter out human-related categories
human_categories <- c(
  "Endocrine system", 
  "Aging", 
  "Circulatory system", 
  "Immune system", 
  "Endocrine and metabolic disease", 
  "Neurodegenerative disease", 
  "Infectious disease: bacterial", 
  "Cancer: overview", 
  "Cancer: specific types", 
  "Cardiovascular disease")

proteomics_pathways_filtered <- proteomics_pathways %>%
  filter(!category %in% human_categories)

# Plot: 
plot_data <- proteomics_pathways_filtered %>%
  group_by(category) %>%
  summarise(total_count = sum(num_of_times_found),
            num_of_prot_KOs = n_distinct(KO), .groups = "drop") %>%
  # Sort by total_count in descending order
  arrange(desc(total_count)) %>%
  # Create a factor with levels ordered by count (for plotting)
  mutate(category = factor(category, levels = rev(category)))

plot_prot_KOs_metaG <- ggplot(plot_data, aes(x = category, y = total_count)) +
  geom_col(fill = "#80CDC1") +
  coord_flip() +
  labs(
    title = "KEGG Ortholog Distribution by Functional Category",
    x = "Functional category",
    y = "Times found in metagenome"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank()
  )
plot_prot_KOs_metaG

ggsave("KEGG_Ortholog_Distribution_by_Category_EVO.pdf", plot_prot_KOs_metaG, width = 8, height = 6, dpi = 300)

##########################
### Functional network ###
##########################

# Create abundance dataset
abundance_long <- KO_proteomics %>%
  pivot_longer(cols = -name, names_to = "sample", values_to = "abundance") %>%
  rename(KO = name)

# Add BRITE categories
abundance_with_categories <- abundance_long %>%
  left_join(proteomics_pathways %>% select(KO, category), by = "KO")

# Remove human-related categories
abundance_with_categories_filtered <- abundance_with_categories %>% 
  filter(!category %in% human_categories)

###############################################
### User defined filters - Adjust as needed ###
###############################################
# Minimum number of KOs per category
min_KOs_per_category <- 5

# Minimum Jaccard similarity to keep an edge (link)
min_jaccard <- 0.10

# Minimum number of shared KOs
min_shared_KOs <- 2

# Node size options:
#   "total_abundance"
#   "mean_abundance"
#   "n_KOs"

node_size_metric <- "n_KOs"

### Processing category mapping to KOs and sizes ###
category_ko <- proteomics_pathways_filtered %>%
  select(category, KO) %>%
  distinct()

# filter tiny categories
category_sizes <- category_ko %>%
  count(category, name = "n_KOs") %>%
  filter(n_KOs >= min_KOs_per_category)

# Check dropped categories
message("Dropped categories: ", 
            paste(setdiff(unique(proteomics_pathways_filtered$category), 
                  category_sizes$category), collapse = ", "))

category_ko <- category_ko %>%
  filter(category %in% category_sizes$category)

# Presence-absence matrix
presence_absence <- category_ko %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = KO, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "category") %>%
  as.matrix()

# Quick sanity check (should match number of KOs per category)
rowSums(presence_absence)

# Jaccard similarity
jaccard_sim <- 1 - as.matrix(vegdist(presence_absence, method = "jaccard", binary = TRUE))

# Table of edges (links) filtered (removed self-links and user-defined low similarity threshold)
edges <- as.data.frame(as.table(jaccard_sim)) %>%
  rename(from = Var1, to = Var2, jaccard = Freq) %>%
  filter(from != to, jaccard >= min_jaccard) %>%
  mutate(pair = map2_chr(from, to, ~paste(sort(c(.x, .y)), 
                                          collapse = "__to__"))) %>% 
  distinct(pair, .keep_all=TRUE) %>%
  mutate(overlap_category = case_when(
    jaccard < 0.3 ~ "Low KO overlap",
    jaccard < 0.5 ~ "Moderate KO overlap",
    TRUE ~ "High KO overlap"))

# Table of filtered nodes (and node size metric) 
nodes <- abundance_with_categories_filtered %>%
  filter(category %in% category_sizes$category) %>%
  filter(category %in% unique(c(edges$from, edges$to))) %>%
  group_by(category) %>%
  summarise(
    total_imputed_normalised_abundance = sum(abundance),
    mean_imputed_normalised_abundance = mean(abundance),
    stdev_imputed_normalised_abundance = sd(abundance),
    n_KOs = n_distinct(KO), .groups="drop")

# Sort nodes by average abundance for plotting
nodes <- nodes %>%
  arrange(desc(mean_imputed_normalised_abundance)) %>%
  mutate(node_ID = row_number())

############################
### Plotting the network ###
############################
graph_network <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Colour breaks for edges
edge_colour <- paletteer_d("khroma::highcontrast")

# Node_ID to Category
node_caption <- nodes %>%
  arrange(node_ID) %>%
  mutate(label = paste0(node_ID, " = ", category)) %>%
  pull(label) %>%
  paste(collapse = "   |   ")

plot_network_prot_categories <- ggraph(graph_network, 
                                       layout = "linear", circular = TRUE) + 
  geom_edge_arc(aes(alpha = overlap_category, colour = overlap_category),
                linewidth = 1.1) +
  geom_node_point(aes(size = .data[[node_size_metric]], 
                      colour = mean_imputed_normalised_abundance)) + 
  geom_node_text(aes(x = x * 1.1, y = y * 1.1, label = node_ID)) +
  # Edge colours
  scale_edge_colour_manual(values = edge_colour, name = "Jaccard similarity") +
  scale_edge_alpha_manual(values = c(
    "Low KO overlap" = 0.25,
    "Moderate KO overlap" = 0.6, 
    "High KO overlap" = 1),
    name = "Jaccard similarity") +
  # Node size and colours
  scale_size_continuous(
    name = "Num. of KOs",
    breaks = c(5, 10, 20, 30, 50),  # Adjust these to match your data range
    range = c(2, 8)  # Min and max point sizes
  ) +
    scale_colour_viridis_c(option = "turbo", name = "Mean KO abundance") +
  theme_void()
plot_network_prot_categories

ggsave("Functional_Network_Proteomics_EVO.pdf", plot_network_prot_categories, width = 10, height = 8, dpi = 300)


##############################################################
### Update first plot with node_id and highlighted colours ###
##############################################################
# Plot data with node IDs
plot_data <- proteomics_pathways_filtered %>%
  group_by(category) %>%
  summarise(total_count = sum(num_of_times_found),
            num_of_prot_KOs = n_distinct(KO), .groups = "drop") %>%
  # Join with node IDs from network
  left_join(nodes %>% select(category, node_ID), by = "category") %>%
  # Create label with node ID and bold formatting for connected categories
  mutate(category_label = case_when(
    !is.na(node_ID) ~ paste0("[", node_ID, "] ", category), TRUE ~ category)) %>%
  # Sort by total_count in descending order
  arrange(desc(total_count)) %>%
  # Create a factor with levels ordered by count (for plotting)
  mutate(category_label = factor(category_label, levels = rev(category_label)))

plot_prot_KOs_metaG_updated <- ggplot(plot_data, aes(x = category_label, y = total_count)) +
  geom_col(aes(fill = !is.na(node_ID))) +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#80CDC1", "FALSE" = "#D3D3D3"),
    labels = c("TRUE" = "Connected", "FALSE" = "Isolated"),
    name = "Functional overlap") +
  labs(
    title = "KEGG Ortholog Distribution by Functional Category",
    x = "Functional category",
    y = "Times found in metagenome") +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 9),  # Use element_markdown instead of element_text
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank()
  )
plot_prot_KOs_metaG_updated

ggsave("KEGG_Ortholog_Distribution_by_Category_EVO_Updated.pdf", plot_prot_KOs_metaG_updated, width = 8, height = 6, dpi = 300)

####################
### Final figure ###
####################

Figure_5 <- ggarrange(plot_network_prot_categories, plot_prot_KOs_metaG_updated, ncol = 1, heights = c(1.5, 1), labels = "AUTO")
Figure_5

ggsave("Figure_5_Proteomics_EVO.pdf", Figure_5, width = 12, height = 14, dpi = 300)
