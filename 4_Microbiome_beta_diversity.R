setwd('~/Documents/MACE/EVO_study/warming_plaine_morte')
source('0_functions_and_packages.R')

# loading functions
loaded_mags = load_mags_data()
metadata = load_metadata()
metadata = check_normality_change(metadata, c("Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
                                              "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",              
                                              "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",          
                                              "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
                                              "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
                                              "Other_Cations","Grouped_OC"))
metadata$Aspect = 'North'
metadata$Aspect[metadata$Replicate %in% c('C','D','E')] = 'South'
metadata$OTC = 'No'
metadata$OTC[metadata$Replicate %in% c('A','B','C')] = 'Yes'

################################################################################
# Part 3. Microbiome - beta diversity
################################################################################

# On MAGs data
set.seed(23)
relab_table = loaded_mags$relab_table[,colnames(loaded_mags$relab_table) %in% metadata$Sample]
relab_table <- relab_table[, match(metadata$Sample, colnames(relab_table))]
all(metadata$Sample == colnames(relab_table))

nmds <- metaMDS(
  t(relab_table),
  distance = "bray",
  k = 2, trymax = 1000, try=500)

scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$Temp_soil <- metadata$Temp_soil
scores_df$Replicate <- metadata$Replicate
scores_df$Aspect <- metadata$Aspect
scores_df$OTC <- metadata$OTC

env_vars <- metadata[, c(
  "Temp_soil_scaled",
  "TC_scaled",
  "TOC_scaled",
  "TIN_scaled",
  "pH_scaled",
  "Water_Content_Lab_scaled"
)]

ef <- envfit(
  nmds,
  env_vars,
  permutations = 9999,
  strata = metadata$Replicate
)

print(ef)

ef_scores <- scores(ef, display = "vectors")
ef_r2 <- ef$vectors$r
ef_p  <- ef$vectors$pvals

ef_df <- data.frame(
  var = rownames(ef_scores),
  NMDS1 = ef_scores[,1],
  NMDS2 = ef_scores[,2],
  r = ef_r2,
  p = ef_p
)

ef_sig <- ef_df[ef_df$p < 0.1, ]
ef_sig$Significance = ifelse(ef_sig$p < 0.05, '< 0.05', '< 0.1')
ef_sig$var_renamed = ef_sig$var
ef_sig$var_renamed[ef_sig$var_renamed == 'TIN_scaled'] = 'TIN'
ef_sig$var_renamed[ef_sig$var_renamed == 'TOC_scaled'] = 'TOC'
ef_sig$var_renamed[ef_sig$var_renamed == 'Water_Content_Lab_scaled'] = 'Soil water content'
ef_sig$var_renamed[ef_sig$var_renamed == 'pH_scaled'] = 'pH'
ef_sig$var_renamed[ef_sig$var_renamed == 'Temp_soil_scaled'] = 'Soil temperature'

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
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      "Stress = ", round(nmds$stress, 3)
    ),
    hjust = 1.1, vjust = 1.1,
    size = 4
  ) + theme(legend.position = 'top')


# Compute dissimilarities
bray <- vegdist(t(relab_table), method = "bray")
jaccard <- vegdist(t(relab_table), method = "jaccard")
metadata$Season_year = map_chr(1:nrow(metadata), function(i) paste0(metadata$Season[i], metadata$Year[i], metadata$Replicate[i], collapse = '_'))

# Bray-Curtis
perm_bray <- adonis2(
  bray ~ Temp_soil_scaled,
  data = metadata,
  permutations = 999,
  strata = metadata$Season_year, by ='margin')
print(perm_bray)


# Jaccard
perm_jacc <- adonis2(
  jaccard ~ Temp_soil_scaled,
  data = metadata,
  permutations = 999,
  strata = metadata$Season_year, by ='margin')
print(perm_jacc)




library(tidyverse)

# Assume metadata has: Sample, Replicate, Year, Condition, Temp_soil_scaled

samples <- metadata$Sample

# Build pairwise combinations
pairs_bc <- expand_grid(Sample1 = samples, Sample2 = samples) %>%
  filter(Sample1 != Sample2)  %>%
  rowwise() %>%
  mutate(
    Bray = bray[Sample1, Sample2],
    TempDiff = abs(metadata$Temp_soil_scaled[metadata$Sample == Sample1] -
                     metadata$Temp_soil_scaled[metadata$Sample == Sample2]),
    Replicate1 = metadata$Replicate[metadata$Sample == Sample1],
    Replicate2 = metadata$Replicate[metadata$Sample == Sample2],
    Year1 = metadata$Year[metadata$Sample == Sample1],
    Year2 = metadata$Year[metadata$Sample == Sample2],
    Season1 = metadata$Season[metadata$Sample == Sample1],
    Season2 = metadata$Season[metadata$Sample == Sample2]
  ) %>%
  ungroup() %>%
  filter(Year1 == Year2, Season1 == Season2)

pairs_jc <- expand_grid(Sample1 = samples, Sample2 = samples) %>%
  filter(Sample1 != Sample2)  %>%
  rowwise() %>%
  mutate(
    Jaccard = jaccard[Sample1, Sample2],
    TempDiff = abs(metadata$Temp_soil_scaled[metadata$Sample == Sample1] -
                     metadata$Temp_soil_scaled[metadata$Sample == Sample2]),
    Replicate1 = metadata$Replicate[metadata$Sample == Sample1],
    Replicate2 = metadata$Replicate[metadata$Sample == Sample2],
    Year1 = metadata$Year[metadata$Sample == Sample1],
    Year2 = metadata$Year[metadata$Sample == Sample2],
    Season1 = metadata$Season[metadata$Sample == Sample1],
    Season2 = metadata$Season[metadata$Sample == Sample2]
  ) %>%
  ungroup() %>%
  filter(Year1 == Year2, Season1 == Season2)

# Plot
bc_temp = ggplot(pairs_bc, aes(x = TempDiff, y = Bray)) +
  geom_point(alpha = 0.2, colour='slategray4') +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  
  # Regression equation + R2
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = 0.95) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85) +
  theme_classic() +
  labs(
    x = "Difference in soil temperature (scaled)",
    y = "Bray–Curtis dissimilarity"
  ) + ylim(0.13, 0.8)

# Plot
jc_temp = ggplot(pairs_jc, aes(x = TempDiff, y = Jaccard)) +
  geom_point(alpha = 0.2, colour='slategray4') +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  
  # Regression equation + R2
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = 0.95) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85) +
  theme_classic() +
  labs(
    x = "Difference in soil temperature (scaled)",
    y = "Jaccard dissimilarity"
  ) + ylim(0.15, 0.95)


splot = ggarrange(bc_temp, jc_temp, labels = c('B', 'C'), nrow = 2, ncol = 1)
ggarrange(nmds_env, splot, ncol = 2, labels = c('A', ''), nrow = 1, widths = c(0.65, 0.35))
ggsave('figures/Figure3.pdf', width = 9, height = 7)







# On 16S data
set.seed(23)
metadata = metadata %>% filter(Year == 2024)
relab_table = loaded_16s$relab_table[,colnames(loaded_16s$relab_table) %in% metadata$Sample]
relab_table <- relab_table[, match(metadata$Sample, colnames(relab_table))]
all(metadata$Sample == colnames(relab_table))

nmds <- metaMDS(
  t(relab_table),
  distance = "bray",
  k = 2, trymax = 1000, try=500)

scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$Temp_soil <- metadata$Temp_soil
scores_df$Replicate <- metadata$Replicate
scores_df$Aspect <- metadata$Aspect
scores_df$OTC <- metadata$OTC

env_vars <- metadata[, c(
  "Temp_soil_scaled",
  "TC_scaled",
  "TOC_scaled",
  "TIN_scaled",
  "pH_scaled",
  "Water_Content_Lab_scaled"
)]

ef <- envfit(
  nmds,
  env_vars,
  permutations = 9999,
  strata = metadata$Replicate
)

print(ef)

ef_scores <- scores(ef, display = "vectors")
ef_r2 <- ef$vectors$r
ef_p  <- ef$vectors$pvals

ef_df <- data.frame(
  var = rownames(ef_scores),
  NMDS1 = ef_scores[,1],
  NMDS2 = ef_scores[,2],
  r = ef_r2,
  p = ef_p
)

ef_sig <- ef_df[ef_df$p < 0.1, ]
ef_sig$Significance = ifelse(ef_sig$p < 0.05, '< 0.05', '< 0.1')
ef_sig$var_renamed = ef_sig$var
ef_sig$var_renamed[ef_sig$var_renamed == 'TC_scaled'] = 'TC'
ef_sig$var_renamed[ef_sig$var_renamed == 'pH_scaled'] = 'pH'

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
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      "Stress = ", round(nmds$stress, 3)
    ),
    hjust = 1.1, vjust = 1.1,
    size = 4
  ) + theme(legend.position = 'top')


# Compute dissimilarities
bray <- vegdist(t(relab_table), method = "bray")
jaccard <- vegdist(t(relab_table), method = "jaccard")
metadata$Season_year = map_chr(1:nrow(metadata), function(i) paste0(metadata$Season[i], metadata$Year[i], metadata$Replicate[i], collapse = '_'))

# Bray-Curtis
perm_bray <- adonis2(
  bray ~ Temp_soil_scaled,
  data = metadata,
  permutations = 999,
  strata = metadata$Season_year, by ='margin')
print(perm_bray)


# Jaccard
perm_jacc <- adonis2(
  jaccard ~ Temp_soil_scaled,
  data = metadata,
  permutations = 999,
  strata = metadata$Season_year, by ='margin')
print(perm_jacc)




library(tidyverse)

# Assume metadata has: Sample, Replicate, Year, Condition, Temp_soil_scaled

samples <- metadata$Sample

# Build pairwise combinations
bray=as.matrix(bray)
pairs_bc <- expand_grid(Sample1 = samples, Sample2 = samples) %>%
  filter(Sample1 != Sample2)  %>%
  rowwise() %>%
  mutate(
    Bray = bray[Sample1, Sample2],
    TempDiff = abs(metadata$Temp_soil_scaled[metadata$Sample == Sample1] -
                     metadata$Temp_soil_scaled[metadata$Sample == Sample2]),
    Replicate1 = metadata$Replicate[metadata$Sample == Sample1],
    Replicate2 = metadata$Replicate[metadata$Sample == Sample2],
    Year1 = metadata$Year[metadata$Sample == Sample1],
    Year2 = metadata$Year[metadata$Sample == Sample2],
    Season1 = metadata$Season[metadata$Sample == Sample1],
    Season2 = metadata$Season[metadata$Sample == Sample2]
  ) %>%
  ungroup() %>%
  filter(Year1 == Year2, Season1 == Season2)

jaccard=as.matrix(jaccard)
pairs_jc <- expand_grid(Sample1 = samples, Sample2 = samples) %>%
  filter(Sample1 != Sample2)  %>%
  rowwise() %>%
  mutate(
    Jaccard = jaccard[Sample1, Sample2],
    TempDiff = abs(metadata$Temp_soil_scaled[metadata$Sample == Sample1] -
                     metadata$Temp_soil_scaled[metadata$Sample == Sample2]),
    Replicate1 = metadata$Replicate[metadata$Sample == Sample1],
    Replicate2 = metadata$Replicate[metadata$Sample == Sample2],
    Year1 = metadata$Year[metadata$Sample == Sample1],
    Year2 = metadata$Year[metadata$Sample == Sample2],
    Season1 = metadata$Season[metadata$Sample == Sample1],
    Season2 = metadata$Season[metadata$Sample == Sample2]
  ) %>%
  ungroup() %>%
  filter(Year1 == Year2, Season1 == Season2)

# Plot
bc_temp = ggplot(pairs_bc, aes(x = TempDiff, y = Bray)) +
  geom_point(alpha = 0.2, colour='slategray4') +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  
  # Regression equation + R2
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = 0.95) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85) +
  theme_classic() +
  labs(
    x = "Difference in soil temperature (scaled)",
    y = "Bray–Curtis dissimilarity"
  ) 

# Plot
jc_temp = ggplot(pairs_jc, aes(x = TempDiff, y = Jaccard)) +
  geom_point(alpha = 0.2, colour='slategray4') +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  
  # Regression equation + R2
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = 0.95) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85) +
  theme_classic() +
  labs(
    x = "Difference in soil temperature (scaled)",
    y = "Jaccard dissimilarity"
  ) 

splot = ggarrange(bc_temp, jc_temp, labels = c('B', 'C'), nrow = 2, ncol = 1)
ggarrange(nmds_env, splot, ncol = 2, labels = c('A', ''), nrow = 1, widths = c(0.65, 0.35))
ggsave('figures/Figure_S??_16s_betadiv.pdf', width = 9, height = 7)


pairs_jc$Condition1 = 'Control'
pairs_jc$Condition1[pairs_jc$Replicate1 %in% c('A', 'B', 'C')] = 'Warming'
pairs_jc$Condition2 = 'Control'
pairs_jc$Condition2[pairs_jc$Replicate2 %in% c('A', 'B', 'C')] = 'Warming'

pairs_jc$ConditionGroup = 'Across'
pairs_jc$ConditionGroup[(pairs_jc$Condition1 == pairs_jc$Condition2) & (pairs_jc$Condition1 == 'Control')] = 'Within Control'
pairs_jc$ConditionGroup[(pairs_jc$Condition1 == pairs_jc$Condition2) & (pairs_jc$Condition1 == 'Warming')] = 'Within Warming'

comparisons = list(c('Within Warming', 'Within Control'),
                   c('Across', 'Within Warming'),
                   c('Across', 'Within Control'))
sp1 = ggplot(pairs_jc, aes(x=ConditionGroup, fill=ConditionGroup, y=Jaccard)) + geom_boxplot() + stat_compare_means(label.y = 0.65) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = comparisons,
                     p.adjust.method = "BH") + theme_classic() +
  labs(
    x = "",
    y = "Jaccard dissimilarity"
  ) + scale_fill_manual(values = c('violet','#012169', '#C8102E'))



pairs_bc$Condition1 = 'Control'
pairs_bc$Condition1[pairs_bc$Replicate1 %in% c('A', 'B', 'C')] = 'Warming'
pairs_bc$Condition2 = 'Control'
pairs_bc$Condition2[pairs_bc$Replicate2 %in% c('A', 'B', 'C')] = 'Warming'

pairs_bc$ConditionGroup = 'Across'
pairs_bc$ConditionGroup[(pairs_bc$Condition1 == pairs_bc$Condition2) & (pairs_bc$Condition1 == 'Control')] = 'Within Control'
pairs_bc$ConditionGroup[(pairs_bc$Condition1 == pairs_bc$Condition2) & (pairs_bc$Condition1 == 'Warming')] = 'Within Warming'

sp2 = ggplot(pairs_bc, aes(x=ConditionGroup, fill=ConditionGroup, y=Bray)) + geom_boxplot() + 
  stat_compare_means(label.y = 0.45) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = comparisons,
                     p.adjust.method = "BH") + theme_classic() +
  labs(
    x = "",
    y = "Bray-Curtis dissimilarity"
  ) + scale_fill_manual(values = c('violet','#012169', '#C8102E'))

ggarrange(sp1 + theme(legend.position = 'none'), sp2 + theme(legend.position = 'none'), ncol = 1, nrow = 2, labels = c('A', 'B'))
ggsave('figures/Figure_S??_bray_jaccard_condition.pdf', width=5, height =7)


