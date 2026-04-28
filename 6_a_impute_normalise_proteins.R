devtools::install_github("mildpiggy/DEP2")

library(tidyverse)
library(DEP2)
library(compositions)
library(dplyr)
library(stringr)
library(purrr)
library(MASS)
library(broom)
library(rstatix)
setwd('~/Documents/MACE/EVO_study/warming_plaine_morte')

# Load data and add the Potential contaminant column based on gene name (Uniprot ID or megahit-prodigal contig_gene)
data = read.csv('data/proteomics/combined_protein.tsv', sep='\t') 
data$Potential.contaminant = '-'
data$Potential.contaminant[startsWith(data$Protein.ID, 'CON_')] = '+'
data <- filter(data, Reverse != "+", Potential.contaminant != "+")

samples_to_numbers = read.csv('data/proteomics/samples_to_numbers.csv', header = F, col.names = c('Sample', 'Number'))
samples_to_numbers$Name = paste0(samples_to_numbers$Number, '.mzML')
samples_to_numbers$Name[samples_to_numbers$Name %in% c('17.mzML', '23.mzML', '29.mzML', '35.mzML', '41.mzML')] = paste0('EVO_250604_0835_Ascend_A_S0', samples_to_numbers$Name[samples_to_numbers$Name %in% c('17.mzML', '23.mzML', '29.mzML', '35.mzML', '41.mzML')])
samples_to_numbers$Name[samples_to_numbers$Name %in% c('18.mzML', '24.mzML', '30.mzML', '36.mzML', '42.mzML')] = paste0('EVO_250604_0835_Ascend_B_S0', samples_to_numbers$Name[samples_to_numbers$Name %in% c('18.mzML', '24.mzML', '30.mzML', '36.mzML', '42.mzML')])
samples_to_numbers$Name[samples_to_numbers$Name %in% c('19.mzML', '25.mzML', '31.mzML', '37.mzML', '43.mzML')] = paste0('EVO_250604_0835_Ascend_C_S0', samples_to_numbers$Name[samples_to_numbers$Name %in% c('19.mzML', '25.mzML', '31.mzML', '37.mzML', '43.mzML')])
samples_to_numbers$Name[samples_to_numbers$Name %in% c('20.mzML', '26.mzML', '32.mzML', '38.mzML', '44.mzML')] = paste0('EVO_250604_0835_Ascend_D_S0', samples_to_numbers$Name[samples_to_numbers$Name %in% c('20.mzML', '26.mzML', '32.mzML', '38.mzML', '44.mzML')])
samples_to_numbers$Name[samples_to_numbers$Name %in% c('21.mzML', '27.mzML', '33.mzML', '39.mzML', '45.mzML')] = paste0('EVO_250604_0835_Ascend_E_S0', samples_to_numbers$Name[samples_to_numbers$Name %in% c('21.mzML', '27.mzML', '33.mzML', '39.mzML', '45.mzML')])
samples_to_numbers$Name[samples_to_numbers$Name %in% c('22.mzML', '28.mzML', '34.mzML', '40.mzML', '46.mzML')] = paste0('EVO_250604_0835_Ascend_F_S0', samples_to_numbers$Name[samples_to_numbers$Name %in% c('22.mzML', '28.mzML', '34.mzML', '40.mzML', '46.mzML')])
samples_to_numbers = samples_to_numbers[samples_to_numbers$Name != '47.mzML',]
samples_to_numbers$Name = paste0(samples_to_numbers$Name, '.MaxLFQ.Intensity')

for (name in samples_to_numbers$Name){
  colnames_to_change = colnames(data)[endsWith(colnames(data), name)]
  for (colname in colnames_to_change){
    colnames(data)[colnames(data) == colname] = gsub(name, samples_to_numbers$Sample[samples_to_numbers$Name == name], colname)
  }
}

dim(data)
metadata = data.frame(label = samples_to_numbers$Sample)
metadata$condition = 'Control'
metadata$condition[startsWith(metadata$label, 'L_D')] = 'Warming'
metadata$condition[startsWith(metadata$label, 'L_E')] = 'Warming'
metadata$condition[startsWith(metadata$label, 'L_F')] = 'Warming'
metadata$replicate = metadata$label
rownames(metadata) = metadata$label

# Load annotation and add PFAMs as gene names
annot = read.csv('data/proteomics/annotations_mags_EDU_2025.emapper.annotations', skip = 4, header = T, sep='\t')
data$Gene.names = map_chr(data$Entry.Name, function(x) ifelse(length(annot$KEGG_ko[annot$X.query %in% x]) > 0, paste0(annot$KEGG_ko[annot$X.query %in% x], collapse = ';'), ''))
data$Gene.names = gsub(',',';',data$Gene.names)
data$Gene.names = gsub('-','',data$Gene.names)
data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# additional: 4 contaminants were not removed...
data = data %>% filter(Protein.ID != '')

# Identify LFQ columns
LFQ_columns <- which(startsWith(colnames(data), "L_"))
LFQ_names <- colnames(data)[LFQ_columns]

# Keep relevant columns
data_ko <- data %>%
  select(Protein.ID, Gene.names, all_of(LFQ_names)) %>%
  filter(Gene.names != "" & !is.na(Gene.names))

# Split multiple KOs
data_ko <- data_ko %>%
  separate_rows(Gene.names, sep = ";")

# Sum intensities per KO
data_ko_summed <- data_ko %>%
  group_by(Gene.names) %>%
  summarise(
    across(
      all_of(LFQ_names),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Rename for DEP2
colnames(data_ko_summed)[1] <- "Gene.names"
data_ko_summed$Protein.ID <- data_ko_summed$Gene.names

# Replace original data
data <- data_ko_summed


# Make gene names unique, adds the "name" column
data_unique <- make_unique(data, "Gene.names", "Protein.ID", delim = ';')
data_unique$name %>% duplicated() %>% any()

# Load with the parsing function, specifying the LFQ intensity columns
LFQ_columns <- which(startsWith(colnames(data_unique), "L_"))
data_se <- make_se(data_unique, LFQ_columns, metadata)
plot_frequency(data_se)
plot_numbers(data_se)
plot_coverage(data_se)

# Filter proteins based on occurrence, max number of missing values is set to 50% of samples
data_filt <- filter_se(data_se, fraction = 0.5)
plot_frequency(data_filt)
plot_numbers(data_filt)
plot_coverage(data_filt)

# Normalisation
data_norm <- normalize_vsn(data_filt)
meanSdPlot(data_norm)
plot_normalization(data_filt, data_norm)

plot_missval(data_norm)
plot_detect(data_norm)

# We use KNN for imputation as based on the missval plot, it seems that missing values are missing at random (MAR)
# c.f. https://pmc.ncbi.nlm.nih.gov/articles/PMC10949645/#S10
data_imp <- DEP2::impute(data_norm, fun = "knn", rowmax = 0.95)
plot_imputation(data_norm, data_imp)

out_data = as.data.frame(data_imp@assays@data@listData)
out_data$name = rownames(out_data)
write.table(out_data, file = 'data/proteomics/metaP_norm_imp.tsv', quote = F, row.names = F, sep = '\t')

