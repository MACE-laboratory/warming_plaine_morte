library(tidyverse)
library(vegan)
library(ANCOMBC)
library(reshape2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(mgcv)
library(compositions)
library(ggpubr)
library(ggrepel)
library(KEGGREST)
library(performance)
library(gratia)
library(ggstance)
library(phyloseq)
library(microbiome)
library(cowplot)
library(PCAtest)
library(EcolUtils)

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
  "Citrate_scale" = 'Citrate')

load_mags_data <- function(){
  ab_tab = read.csv('data/mags/mags_abundances.tsv', sep='\t')
  tx_tab = read.csv('data/mags/gtdbtk.bac120.summary.tsv', sep='\t')
  ql_tab = read.csv('data/mags/quality_report.tsv', sep='\t')
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

load_meta_data <- function() {
  # Check for precomputed file
  if (file.exists("data/ko_table.tsv")) {
    message("Loading existing KO abundance table from data/ko_table.tsv")
    ko_tab = read_tsv("data/ko_table.tsv", show_col_types = FALSE)
    pfam_tab = read_tsv("data/pfam_table.tsv", show_col_types = FALSE)
    cazy_tab = read_tsv("data/cazy_table.tsv", show_col_types = FALSE)
    
    raw_ko_tab   <- ko_tab %>% column_to_rownames("KEGG_ko")
    raw_pfam_tab <- pfam_tab %>% column_to_rownames("PFAMs")
    raw_cazy_tab <- cazy_tab %>% column_to_rownames("CAZy")
    
    relab_ko_tab   <-   raw_ko_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
    relab_pfam_tab <- raw_pfam_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
    relab_cazy_tab <- raw_cazy_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
    
    relab_ko_tab_with_scalar = relab_ko_tab
    relab_ko_tab_with_scalar[relab_ko_tab_with_scalar == 0] <- min(relab_ko_tab_with_scalar[relab_ko_tab_with_scalar > 0])
    ko_clr_tab = clr(relab_ko_tab_with_scalar)
    
    relab_pfam_tab_with_scalar = relab_pfam_tab
    relab_pfam_tab_with_scalar[relab_pfam_tab_with_scalar == 0] <- min(relab_pfam_tab_with_scalar[relab_pfam_tab_with_scalar > 0])
    pfam_clr_tab = clr(relab_pfam_tab_with_scalar)
    
    relab_cazy_tab_with_scalar = relab_cazy_tab
    relab_cazy_tab_with_scalar[relab_cazy_tab_with_scalar == 0] <- min(relab_cazy_tab_with_scalar[relab_cazy_tab_with_scalar > 0])
    cazy_clr_tab = clr(relab_cazy_tab_with_scalar)
    
    return(list(ko_raw_table = raw_ko_tab, pfam_raw_table = raw_pfam_tab, cazy_raw_table = raw_cazy_tab, 
                ko_relab_table = relab_ko_tab, pfam_relab_table = relab_pfam_tab, cazy_relab_table = relab_cazy_tab, 
                ko_clr_table = ko_clr_tab, pfam_clr_table = pfam_clr_tab, cazy_clr_table = cazy_clr_tab))
  }
  
  message("Creating KO abundance table from input files...")
  
  # Load abundance
  ab_tab <- read_tsv("data/functional/all_contigs_trimmed_mean.tsv", show_col_types = FALSE)
  ab_tab[is.na(ab_tab)] <- 0
  colnames(ab_tab)[1] <- "contig"
  colnames(ab_tab)[-1] <- map_chr(colnames(ab_tab)[-1], ~strsplit(.x, "\\.\\.\\.")[[1]][1])
  colnames(ab_tab)[-1] <- gsub("^final\\.contigs\\.renamed\\.filtered\\.fa\\.", "", colnames(ab_tab)[-1])
  
  # Rename sample columns
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
  
  print('Samples sums in the abundance file:')
  print(quantile(colSums(ab_tab[,colnames(ab_tab) != 'contig'])))
  
  print('Samples sums for assembly 0:')
  print(quantile(colSums(ab_tab[startsWith(ab_tab$contig, 'sub_assembly_0'),colnames(ab_tab) != 'contig'])))
  
  print('Samples sums for assembly 1:')
  print(quantile(colSums(ab_tab[startsWith(ab_tab$contig, 'sub_assembly_1'),colnames(ab_tab) != 'contig'])))
  
  print('Samples sums for assembly 2:')
  print(quantile(colSums(ab_tab[startsWith(ab_tab$contig, 'sub_assembly_2'),colnames(ab_tab) != 'contig'])))
  
  print('Samples sums for assembly 3:')
  print(quantile(colSums(ab_tab[startsWith(ab_tab$contig, 'sub_assembly_3'),colnames(ab_tab) != 'contig'])))
  
  print('Range of contigs abundances (please check carefully there are no zeros):')
  print(quantile(rowSums(ab_tab[,colnames(ab_tab) != 'contig'])))
  
  # Load annotations
  ann_files <- list.files("data/functional", pattern = "assembly_\\d\\.emapper.annotations$", full.names = TRUE)
  ann_list <- lapply(ann_files, function(file) {
    prefix <- str_extract(basename(file), "assembly_\\d+")
    prefix <- paste0('sub_', prefix)
    lines <- read_lines(file)
    lines <- lines[!startsWith(lines, "##")]
    
    tmp <- tempfile()
    write_lines(lines, tmp)
    
    df <- read_tsv(tmp, show_col_types = FALSE)
    df %>%
      rename(query = `#query`) %>%
      mutate(query = paste0(prefix, "_", query)) %>%
      mutate(
        KEGG_ko = gsub("ko:", "", KEGG_ko),
        contig = str_remove(query, "_\\d+$")
      )
  })
  annotations <- bind_rows(ann_list)
  
  ab_tab_filtered = ab_tab[ab_tab$contig %in% annotations$contig,]
  rowsums_before = rowSums(ab_tab_filtered[,colnames(ab_tab_filtered) != 'contig'])
  ab_tab_filtered = ab_tab_filtered %>% column_to_rownames('contig')
  ab_tab_filtered = sapply(unique(colnames(ab_tab_filtered)), function(x) rowSums(ab_tab_filtered[,grepl(x, colnames(ab_tab_filtered))]))
  ab_tab_filtered = as.data.frame(ab_tab_filtered) %>% rownames_to_column(var='contig')
  rowsums_after = rowSums(ab_tab_filtered[,colnames(ab_tab_filtered) != 'contig'])
  print(mean(rowsums_after == rowsums_before)) # check that no counts were increased/reduced
  
  ko_annot = annotations %>% select(KEGG_ko, contig) %>% separate_longer_delim(KEGG_ko, ',') %>% filter(KEGG_ko != '-')
  pfam_annot = annotations %>% select(PFAMs, contig) %>% separate_longer_delim(PFAMs, ',') %>% filter(PFAMs != '-')
  cazy_annot = annotations %>% select(CAZy, contig) %>% separate_longer_delim(CAZy, ',') %>% filter(CAZy != '-')
  
  ko_merged = full_join(ab_tab_filtered, ko_annot, by='contig') %>% filter(!(is.na(KEGG_ko)))
  pfam_merged = full_join(ab_tab_filtered, pfam_annot, by='contig') %>% filter(!(is.na(PFAMs)))
  cazy_merged = full_join(ab_tab_filtered, cazy_annot, by='contig') %>% filter(!(is.na(CAZy)))
  
  ko_tab = ko_merged %>% select(-contig) %>% group_by(KEGG_ko) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))
  
  pfam_tab = pfam_merged %>% select(-contig) %>% group_by(PFAMs) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))
  
  cazy_tab = cazy_merged %>% select(-contig) %>% group_by(CAZy) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))
  
  write_tsv(ko_tab, "data/ko_table.tsv")
  write_tsv(pfam_tab, "data/pfam_table.tsv")
  write_tsv(cazy_tab, "data/cazy_table.tsv")
  
  raw_ko_tab   <- ko_tab %>% column_to_rownames("KEGG_ko")
  raw_pfam_tab <- pfam_tab %>% column_to_rownames("PFAMs")
  raw_cazy_tab <- cazy_tab %>% column_to_rownames("CAZy")
  
  relab_ko_tab   <-   raw_ko_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
  relab_pfam_tab <- raw_pfam_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
  relab_cazy_tab <- raw_cazy_tab %>% mutate(across(everything(), ~ .x / sum(.x)))
  
  relab_ko_tab_with_scalar = relab_ko_tab
  relab_ko_tab_with_scalar[relab_ko_tab_with_scalar == 0] <- min(relab_ko_tab_with_scalar[relab_ko_tab_with_scalar > 0])
  ko_clr_tab = clr(relab_ko_tab_with_scalar)
  
  relab_pfam_tab_with_scalar = relab_pfam_tab
  relab_pfam_tab_with_scalar[relab_pfam_tab_with_scalar == 0] <- min(relab_pfam_tab_with_scalar[relab_pfam_tab_with_scalar > 0])
  pfam_clr_tab = clr(relab_pfam_tab_with_scalar)
  
  relab_cazy_tab_with_scalar = relab_cazy_tab
  relab_cazy_tab_with_scalar[relab_cazy_tab_with_scalar == 0] <- min(relab_cazy_tab_with_scalar[relab_cazy_tab_with_scalar > 0])
  cazy_clr_tab = clr(relab_cazy_tab_with_scalar)
  
  return(list(ko_raw_table = raw_ko_tab, pfam_raw_table = raw_pfam_tab, cazy_raw_table = raw_cazy_tab, 
              ko_relab_table = relab_ko_tab, pfam_relab_table = relab_pfam_tab, cazy_relab_table = relab_cazy_tab, 
              ko_clr_table = ko_clr_tab, pfam_clr_table = pfam_clr_tab, cazy_clr_table = cazy_clr_tab))
}

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
  
  otu_table = otu_table[rowSums(otu_table[,c('negative_control', 'negative_control_24')]) == 0,]
  otu_table = otu_table %>% dplyr::select(-negative_control, -negative_control_24)
  otu_table = otu_table[rowSums(otu_table) > 1,]
  tax_table = tax_table[tax_table$OTU %in% rownames(otu_table),]
  
  relab_tab = sweep(otu_table, 2, colSums(otu_table), FUN = "/")
  
  relab_tab_with_scalar = otu_table
  relab_tab_with_scalar[relab_tab_with_scalar == 0] <- min(relab_tab_with_scalar[relab_tab_with_scalar > 0])
  clr_tab = clr(relab_tab_with_scalar)
  
  # rarefaction
  rarefy_depth <- 8000
  
  rar_table <- apply(otu_table, 2, function(x){
    if(sum(x) < rarefy_depth){
      # return NA or keep original (choose behavior)
      return(rep(NA, length(x)))
    } else {
      probs <- x / sum(x)
      as.vector(rmultinom(1, size = rarefy_depth, prob = probs))
    }
  })
  
  rownames(rar_table) <- rownames(otu_table)
  colnames(rar_table) <- colnames(otu_table)
  
  to_keep = colnames(rar_table)[endsWith(colnames(rar_table), suffix = '_24')]
  otu_table = otu_table[,to_keep]
  relab_tab = relab_tab[,to_keep]
  clr_tab = clr_tab[,to_keep]
  rar_table = rar_table[,to_keep]
  
  return(list(raw_table=otu_table, relab_table=relab_tab, clr_table=clr_tab, rar_table=rar_table, tax_table=tax_table)) 
}

load_metadata <- function(){
  metadata = read.csv('data/SiteL23_24_env_variables_final.csv')
  metadata$Sample = paste0('L_', metadata$Sample)
  metadata$Sample = gsub('\\.','_',metadata$Sample)
  return(metadata)
}

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

a16s_to_genus <- function(loaded_16sr){
  tab = as.data.frame(loaded_16sr$raw_table)
  genus_info = loaded_16sr$tax_table
  genus_info$genus_full = apply(genus_info[, c('Phylum','Class','Order','Family', 'Genus')], 1, paste, collapse = ";")
  genus_info$genus_full[endsWith(genus_info$genus_full, '; ')] = 'Others'
  tab$genus_full = map_chr(rownames(tab), function(x) genus_info$genus_full[genus_info$OTU == x])
  tab_sum <- tab %>%
    group_by(genus_full) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  tab_sum = tab_sum %>% column_to_rownames(var = 'genus_full')
  return(tab_sum)
}

mags_to_genus <- function(loaded_mags){
  tab = as.data.frame(loaded_mags$relab_table)
  genus_info = loaded_mags$mag_info
  genus_info$genus_full = apply(genus_info[, c('Phylum','Class','Order','Family','Genus')], 1, paste, collapse = ";")
  genus_info$genus_full[endsWith(genus_info$genus_full, '; ')] = 'Others'
  tab$genus_full = map_chr(rownames(tab), function(x) genus_info$genus_full[genus_info$MAG == x])
  tab_sum <- tab %>%
    group_by(genus_full) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  tab_sum = tab_sum %>% column_to_rownames(var = 'genus_full')
  return(tab_sum)
}

mags_to_order <- function(loaded_mags){
  tab = as.data.frame(loaded_mags$relab_table)
  class_info = loaded_mags$mag_info
  class_info$order_full = apply(class_info[, c('Phylum','Class','Order')], 1, paste, collapse = "; ")
  class_info$order_full[endsWith(class_info$order_full, '; ')] = 'Others'
  tab$order_full = map_chr(rownames(tab), function(x) class_info$order_full[class_info$MAG == x])
  tab_sum <- tab %>%
    group_by(order_full) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  tab_sum = tab_sum %>% column_to_rownames(var = 'order_full')
  return(tab_sum)
}

mags_to_class <- function(loaded_mags){
  tab = as.data.frame(loaded_mags$relab_table)
  class_info = loaded_mags$mag_info
  class_info$class_full = apply(class_info[, c('Phylum','Class')], 1, paste, collapse = "; ")
  class_info$class_full[endsWith(class_info$class_full, '; ')] = 'Others'
  tab$class_full = map_chr(rownames(tab), function(x) class_info$class_full[class_info$MAG == x])
  tab_sum <- tab %>%
    group_by(class_full) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  tab_sum = tab_sum %>% column_to_rownames(var = 'class_full')
  return(tab_sum)
}

a16s_to_class <- function(loaded_16sr){
  tab = as.data.frame(loaded_16sr$raw_table)
  class_info = loaded_16sr$tax_table
  class_info$class_full = apply(class_info[, c('Phylum','Class')], 1, paste, collapse = "; ")
  class_info$class_full[endsWith(class_info$class_full, '; ')] = 'Others'
  tab$class_full = map_chr(rownames(tab), function(x) class_info$class_full[class_info$OTU == x])
  tab_sum <- tab %>%
    group_by(class_full) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  tab_sum = tab_sum %>% column_to_rownames(var = 'class_full')
  return(tab_sum)
}

add_shannon <- function(metadata, loaded_mags, loaded_meta, loaded_16sr){
  metadata$Shannon_MAGs = diversity(loaded_mags$relab_table, index = "shannon")$shannon
  metadata$Shannon_KEGG = diversity(loaded_meta$ko_relab_table, index = "shannon")$shannon
  metadata$Shannon_PFAM = diversity(loaded_meta$pfam_relab_table, index = "shannon")$shannon
  metadata$Shannon_CAZy = diversity(loaded_meta$cazy_relab_table, index = "shannon")$shannon
  metadata$Shannon_16Sr = diversity(loaded_16sr$relab_table, index = "shannon")$shannon
  metadata$rec_a_copies = as.vector(unlist(loaded_meta$raw_table['K03553',]))
  return(metadata)
}

log_half_min <- function(x){
  half_min = min(x[x > 0], na.rm = T) / 2
  return(log(x + half_min))
}

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

