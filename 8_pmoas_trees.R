setwd('~/Documents/MACE/warming_plaine_morte')
library(purrr)
library(tidyverse)
library(seqinr)

pmoas_genes = read.csv('data/pmoas/subset_pmoA_incontigs.csv')
pmoas_genes$ID = map_chr(1:nrow(pmoas_genes), function(i) paste0(pmoas_genes$assembly[i], '_', pmoas_genes$gene[i]))

annot0 = read.fasta('data/functional/assembly_0.emapper.genepred.fasta')
names(annot0) = map_chr(names(annot0), function(x) paste0('sub_assembly_0_', x))
                        
annot1 = read.fasta('data/functional/assembly_1.emapper.genepred.fasta')
names(annot1) = map_chr(names(annot1), function(x) paste0('sub_assembly_1_', x))
                        
annot2 = read.fasta('data/functional/assembly_2.emapper.genepred.fasta')
names(annot2) = map_chr(names(annot2), function(x) paste0('sub_assembly_2_', x))
                        
annot3 = read.fasta('data/functional/assembly_3.emapper.genepred.fasta')
names(annot3) = map_chr(names(annot3), function(x) paste0('sub_assembly_3_', x))
                        
all_annot = c(annot0, annot1, annot2, annot3)

subset_annot <- all_annot[names(all_annot) %in% pmoas_genes$ID]
write.fasta(subset_annot, names(subset_annot), file.out = 'data/pmoas/pmoas_in_metag.fasta')

nuc <- readDNAStringSet('data/pmoas/pmoas_in_metag.fasta')
prot <- translate(nuc)
writeXStringSet(prot, 'data/pmoas/pmoas_in_metag_proteins.fasta')

#################################################################################
# First we do database sequences comparison
# pmoas are: https://www.kegg.jp/kegg-bin/uniprot_list?ko=K10944
# amoas are: https://www.kegg.jp/kegg-bin/uniprot_list?ko=K28504
# archaeal amoas as outgroups are: Thaumarchaeote_AmoA_KC735600.1 and Nitrososphaera_AmoA_KT351403.1
library(ape)
library(ggplot2)

# Read tree
tree <- read.tree("data/pmoas/pmoa_liebner_comparison/liebner_aligned_to_ref.nwk")

# Pairwise patristic distances
dist_mat <- cophenetic.phylo(tree)

tips <- tree$tip.label

# Reference sequences
outgroups <- grep("_archaeal", tips, value = TRUE)
amoa_refs <- grep("^amoa_", tips, value = TRUE)
amoa_refs = amoa_refs[!(amoa_refs %in% outgroups)]
pmoa_refs <- grep("^pmoa_", tips, value = TRUE)

# Keep all tips as focal sequences
focal_tips <- tips

# Mean distance to amoA references
phylo_dist_amoa <- sapply(focal_tips, function(x) {
  refs <- setdiff(amoa_refs, x)   # avoid self-distance
  if (length(refs) == 0) return(NA_real_)
  min(dist_mat[x, refs], na.rm = TRUE)
})

# Mean distance to pmoA references
phylo_dist_pmoa <- sapply(focal_tips, function(x) {
  refs <- setdiff(pmoa_refs, x)
  if (length(refs) == 0) return(NA_real_)
  min(dist_mat[x, refs], na.rm = TRUE)
})

# Output table
phylo_dist <- data.frame(
  tip = focal_tips,
  phylo_dist_amoa = phylo_dist_amoa,
  phylo_dist_pmoa = phylo_dist_pmoa,
  Group = ifelse(focal_tips %in% outgroups, 'Archaeal amoA', 'Liebner et al. 2016'),
  row.names = NULL
)

phylo_dist

phylo_dist = phylo_dist %>% dplyr::filter(!(tip %in% amoa_refs) | !(tip %in% pmoa_refs))

ggplot(phylo_dist, aes(x=phylo_dist_amoa, y=phylo_dist_pmoa, colour=Group)) + 
  geom_abline(intercept = 0, slope = 1, colour='dimgrey', linetype = 'dashed') + 
  geom_point() + xlab('Phylogenetic distance to closest reference amoA') +
  ylab('Phylogenetic distance to closest reference pmoA')



group <- dplyr::case_when(
  tips %in% outgroups ~ "Archaeal amoA",
  tips %in% amoa_refs ~ "amoA KEGG proteins",
  tips %in% pmoa_refs ~ "pmoA KEGG proteins",
  TRUE ~ "Liebner et al. 2016"
)

groups <- data.frame(
  tip = tips,
  group = group
)

library(ggtree)
tree_rooted <- root(
  tree,
  outgroup = outgroups,
  resolve.root = TRUE
)
groups$group <- factor(
  groups$group,
  levels = c(
    "Liebner et al. 2016",
    "Archaeal amoA",
    "amoA KEGG proteins",
    "pmoA KEGG proteins"
  )
)

p <- ggtree(tree_rooted) %<+% groups +
  geom_tree() +
  geom_tippoint(
    aes(color = group),
    size = 3, alpha=0.5
  ) +
  theme_tree2()
p

# All pairwise distances with group labels
dist_df <- as.data.frame(as.table(dist_mat)) %>%
  dplyr::rename(tip1 = Var1, tip2 = Var2, distance = Freq) %>%
  dplyr::left_join(groups, by = c("tip1" = "tip")) %>%
  dplyr::rename(group1 = group) %>%
  dplyr::left_join(groups, by = c("tip2" = "tip")) %>%
  dplyr::rename(group2 = group)

# Remove self comparisons
dist_df <- dist_df %>%
  dplyr::filter(tip1 != tip2)

# Mean distance between groups
mean_phylo_dist <- dist_df %>%
  dplyr::group_by(group1, group2) %>%
  dplyr::summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    sd_distance = sd(distance, na.rm = TRUE),
    median_distance = median(distance, na.rm = TRUE),
    iqr_distance = quantile(distance, na.rm = TRUE, probs = 0.75) - quantile(distance, na.rm = TRUE, probs = 0.25),
    n_comparisons = dplyr::n(),
    .groups = "drop"
  )

mean_phylo_dist




######### METAG ANALYIS ########################################################
# Read tree
tree <- read.tree("data/pmoas/pmoa_metag_comparison/metag_liebner_aligned_to_ref.nwk")

# Pairwise patristic distances
dist_mat <- cophenetic.phylo(tree)

tips <- tree$tip.label

# Reference sequences
outgroups <- grep("_archaeal", tips, value = TRUE)
amoa_refs <- grep("^amoa_", tips, value = TRUE)
amoa_refs = amoa_refs[!(amoa_refs %in% outgroups)]
pmoa_refs <- grep("^pmoa_", tips, value = TRUE)
metag <- grep("sub_assembly", tips, value = TRUE)


# Keep all tips as focal sequences
focal_tips <- tips

# Mean distance to amoA references
phylo_dist_amoa <- sapply(focal_tips, function(x) {
  refs <- setdiff(amoa_refs, x)   # avoid self-distance
  if (length(refs) == 0) return(NA_real_)
  min(dist_mat[x, refs], na.rm = TRUE)
})

# Mean distance to pmoA references
phylo_dist_pmoa <- sapply(focal_tips, function(x) {
  refs <- setdiff(pmoa_refs, x)
  if (length(refs) == 0) return(NA_real_)
  min(dist_mat[x, refs], na.rm = TRUE)
})

# Output table
phylo_dist <- data.frame(
  tip = focal_tips,
  phylo_dist_amoa = phylo_dist_amoa,
  phylo_dist_pmoa = phylo_dist_pmoa,
  row.names = NULL
)

phylo_dist

phylo_dist = phylo_dist %>% dplyr::filter(!(tip %in% amoa_refs) | !(tip %in% pmoa_refs))


group <- dplyr::case_when(
  tips %in% outgroups ~ "Archaeal amoA",
  tips %in% metag ~ "Metagenomic seqs.",
  tips %in% amoa_refs ~ "amoA KEGG proteins",
  tips %in% pmoa_refs ~ "pmoA KEGG proteins",
  TRUE ~ "Liebner et al. 2016"
)

groups <- data.frame(
  tip = tips,
  group = group
)


library(ggtree)
tree_rooted <- root(
  tree,
  outgroup = outgroups,
  resolve.root = TRUE
)
groups$group <- factor(
  groups$group,
  levels = c(
    "Liebner et al. 2016",
    "Archaeal amoA",
    "amoA KEGG proteins",
    "pmoA KEGG proteins",
    "Metagenomic seqs."
  )
)

p <- ggtree(tree_rooted) %<+% groups +
  geom_tree() +
  geom_tippoint(
    aes(color = group),
    size = 3, alpha=0.5
  ) +
  theme_tree2() + scale_colour_manual(values=c('grey','forestgreen','pink','tomato','steelblue'))
p













# Read fasta IDs
library(Biostrings)

ref <- readAAStringSet("data/pmoas/ref_pmoa_amoa.fasta")
meta <- readAAStringSet("data/pmoas/pmoas_in_metag.fasta")

ref_ids <- names(ref)
ref_ids = ref_ids[ref_ids %in% tree$tip.label]
ref_ids = ref_ids[startsWith(ref_ids, 'pmoa')]
meta_ids <- names(meta)
meta_ids = meta_ids[meta_ids %in% tree$tip.label]

# Tree distances
D <- cophenetic(tree)


#################################################
# Reference-reference distances
#################################################

ref_dist <- D[ref_ids, ref_ids]
ref_values <- ref_dist[upper.tri(ref_dist)]


#################################################
# Reference-metagenome distances
#################################################

meta_values <- as.vector(D[meta_ids, ref_ids])


#################################################
# Threshold range
#################################################

max_ref <- max(ref_values)

thresholds <- seq(
  1.5 * max_ref,
  0.1 * max_ref,
  length.out = 20
)


#################################################
# Run Wilcoxon tests
#################################################

set.seed(123)

results <- data.frame(
  threshold = thresholds,
  p_value = NA,
  median_ref = NA,
  median_meta = NA,
  iqr_ref = NA,
  iqr_meta = NA,
  n_ref = NA,
  n_meta = NA
)


for(i in seq_along(thresholds)) {
  
  t <- thresholds[i]
  
  # distances below threshold
  ref_sub <- ref_values[ref_values <= t]
  meta_sub <- meta_values[meta_values <= t]
  
  
  # require enough values
  if(length(ref_sub) >= 50 &
     length(meta_sub) >= 50) {
    
    # randomly sample equal numbers
    ref_sample <- sample(ref_sub, 50)
    meta_sample <- sample(meta_sub, 50)
    
    test <- wilcox.test(
      ref_sample,
      meta_sample
    )
    
    results$p_value[i] <- test$p.value
    results$median_ref[i] <- median(ref_sample)
    results$iqr_ref[i] <- quantile(ref_sample, probs = 0.75) - quantile(ref_sample, probs = 0.25)
    
    results$median_meta[i] <- median(meta_sample)
    results$iqr_meta[i] <- quantile(meta_sample, probs = 0.75) - quantile(meta_sample, probs = 0.25)
    
    results$n_ref[i] <- length(ref_sub)
    results$n_meta[i] <- length(meta_sub)
  }
}


results
results$padj = p.adjust(results$p_value, method='holm')
results %>% filter(padj > 0.05)

ggplot() + geom_vline(xintercept = 2, colour='black', linetype='dashed') +
  geom_point(results, mapping = aes(x=threshold, y=median_ref), colour='darkgrey') +
  geom_errorbar(results, mapping = aes(x=threshold, ymin=median_ref-iqr_ref, ymax=median_ref+iqr_ref), colour='darkgrey') +
  geom_line(results, mapping = aes(x=threshold, y=median_ref), colour='darkgrey') +
  geom_point(results, mapping = aes(x=threshold, y=median_meta), colour='tomato') +
  geom_errorbar(results, mapping = aes(x=threshold, ymin=median_meta-iqr_meta, ymax=median_meta+iqr_meta), colour='tomato') +
  geom_line(results, mapping = aes(x=threshold, y=median_meta), colour='tomato') + 
  ylab('Cophenetic distance') + xlab('Distance threshold')
ggsave('figures/SI_Figure_pmoa_threshold_phydist.pdf', width = 8, height = 7)

library(ape)
library(ggtree)
library(ggplot2)
library(Biostrings)
library(dplyr)

################################################
# Assign groups
################################################

tip_groups <- data.frame(
  label = tree$tip.label,
  group = "Other"
)
tip_groups$label = map_chr(tip_groups$label, function(x) gsub('_R_','',x))
amoa_ids = tree$tip.label[startsWith(tree$tip.label, 'amoa')]
tip_groups = tip_groups %>% dplyr::filter(label %in% c(meta_ids, ref_ids, amoa_ids))

# amoA outgroup
tip_groups$group[
  tip_groups$label %in% amoa_ids
] <- "amoA outgroups"

# Reference pmoA
tip_groups$group[
  tip_groups$label %in% ref_ids
] <- "ref_pmoa"


tip_groups$label[!(tip_groups$label %in% meta_ids) & (tip_groups$group == 'Other')]

# Metagenomes
meta_distance <- sapply(meta_ids, function(x){
  
  min(D[x, ref_ids])
  
})


# split metagenomes
threshold = 0.5
putative_meta <- names(meta_distance[
  meta_distance <= threshold
])

ambiguous_meta <- names(meta_distance[
  meta_distance > threshold
])


tip_groups$group[
  tip_groups$label %in% putative_meta
] <- "Metagenome putative pmoAs"


tip_groups$group[
  tip_groups$label %in% ambiguous_meta
] <- "Metagenome ambiguous"

table(tip_groups$group)



only_meta_tips = tip_groups %>% filter(label %in% meta_ids)

min_pmoa_distance <- sapply(only_meta_tips$label, function(x){
  
  min(D[x, ref_ids])
  
})

min_amoa_distance <- sapply(only_meta_tips$label, function(x){
  
  min(D[x, amoa_ids])
  
})
only_meta_tips$min_pmoa_dist = min_pmoa_distance
only_meta_tips$min_pmoa_dist = min_amoa_distance

only_meta_tips


amoa_amoa = D[amoa_ids, amoa_ids]
amoa_pmoa = sample(D[amoa_ids, ref_ids], 1000)
pmoa_pmoa = sample(D[ref_ids, ref_ids], 1000)
meta_pmoa = sample(D[ref_ids, meta_ids], 1000)
meta_amoa = sample(D[amoa_ids, meta_ids], 1000)
comp_df = data.frame(Distance=c(amoa_amoa, amoa_pmoa, pmoa_pmoa, meta_pmoa, meta_amoa),
                     Comparison=c(rep('amoa - amoa', length(amoa_amoa)),
                                  rep('amoa - pmoa', 1000),
                                  rep('within pmoa', 1000),
                                  rep('meta - pmoa', 1000),
                                  rep('meta - amoa', 1000)))
ggplot(comp_df, aes(x=Comparison,y=Distance)) + geom_boxplot()


library(ape)
library(Biostrings)
library(DECIPHER)
library(ggtree)

# alignment
aln <- readDNAStringSet("data/pmoas/amoas_pmoas_trim.fasta")

# 1. extract matching tips from tree
keep_tips <- tree$tip.label[tree$tip.label %in% ref_ids]

subtree <- keep.tip(tree, keep_tips)

# 2. extract matching sequences from alignment
seqs <- aln[names(aln) %in% keep_tips]

d <- DistanceMatrix(
  seqs,
  correction = "none"
)

# cluster at 95% identity
# 95% identity = 5% distance
seqs_nogap <- RemoveGaps(seqs)
clusters <- Clusterize(
  seqs_nogap,
  cutoff = 0.1,
  method = "longest"
)
length(unique(clusters$cluster))


cluster_id <- clusters[,1]   # first column containing cluster assignments
names(cluster_id) <- rownames(clusters)
rep_names <- sapply(split(names(cluster_id), cluster_id), function(x) {
  dm <- as.matrix(d)[x, x, drop = FALSE]
  x[which.min(rowMeans(dm))]
})

derep_seqs <- seqs[rep_names]
ref_to_keep = names(derep_seqs)

amoa_ids = names(seqs)[1:3]
new_align_ids = c(amoa_ids, ref_to_keep, tip_groups$label[tip_groups$group == 'Metagenome putative pmoAs'])
seqs <- readDNAStringSet("data/pmoas/amoas_pmoas.fasta")
subset_seqs <- seqs[names(seqs) %in% new_align_ids]
writeXStringSet(subset_seqs, "data/pmoas/pmoas_subset_final.fasta")

################################################
# Root tree by amoA
################################################

# Read tree
tree <- read.tree("data/pmoas/pmoas_subset_final_trim.fasta.treefile")

# Define outgroup
outgroup <- amoa_ids

# Reroot
tree_rooted <- root(
  tree,
  outgroup = outgroup,
  resolve.root = TRUE
)

# IDs to highlight
# ref_ids <- c(...)
# meta_ids <- c(...)

# Plot base tree
p <- ggtree(tree_rooted)

# Colour tips
p <- p +
  geom_tippoint(
    data = subset(p$data, label %in% ref_ids),
    colour = "grey50",
    size = 2
  ) +
  geom_tippoint(
    data = subset(p$data, label %in% meta_ids),
    colour = "blue",
    size = 2
  )

p