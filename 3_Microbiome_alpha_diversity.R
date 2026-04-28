setwd('~/Documents/MACE/EVO_study/warming_plaine_morte')
source('0_functions_and_packages.R')

# loading functions
loaded_mags = load_mags_data()
loaded_16s = load_16sr_data()
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
# Part 2. Microbiome - alpha diversity
################################################################################
# 2.1 Taxonomic diversity
library(metacoder)
library(dplyr)

mag_data = loaded_mags$relab_table
mag_data$mag_id = rownames(mag_data)
mag_data$Taxonomy = map_chr(mag_data$mag_id, function(x) loaded_mags$mag_info$Taxonomy[loaded_mags$mag_info$MAG == x])
mag_data$Taxonomy = map_chr(mag_data$Taxonomy, function(x) strsplit(x, split = ';f__')[[1]][1])
mag_data$mean_rel_ab = rowMeans(loaded_mags$relab_table)

a16s_data = loaded_16s$relab_table
a16s_data$otu_id = rownames(a16s_data)
a16s_data$Taxonomy = map_chr(a16s_data$otu_id, function(x) loaded_16s$tax_table$Taxonomy[loaded_16s$tax_table$OTU == x])
a16s_data$Taxonomy = map_chr(a16s_data$Taxonomy, function(x) strsplit(x, split = ';f__')[[1]][1])
a16s_data$mean_rel_ab = rowMeans(loaded_16s$relab_table)

# MAGs stats section
write.table(loaded_mags$mag_info, file = 'stats/mags_checkm2_gtdb.csv', quote = F, row.names = F)

median_length = median(loaded_mags$mag_info$Length)
iqr_length = quantile(loaded_mags$mag_info$Length, probs = 0.75) - quantile(loaded_mags$mag_info$Length, probs = 0.25)

median_comp = median(loaded_mags$mag_info$Completeness)
iqr_comp = quantile(loaded_mags$mag_info$Completeness, probs = 0.75) - quantile(loaded_mags$mag_info$Completeness, probs = 0.25)

median_cont = median(loaded_mags$mag_info$Contamination)
iqr_cont = quantile(loaded_mags$mag_info$Contamination, probs = 0.75) - quantile(loaded_mags$mag_info$Contamination, probs = 0.25)

nhq = sum((loaded_mags$mag_info$Completeness > 90) & (loaded_mags$mag_info$Contamination < 5))
nmq = sum((loaded_mags$mag_info$Completeness > 50) & (loaded_mags$mag_info$Contamination < 10)) - nhq
nlq = nrow(loaded_mags$mag_info) - nhq - nmq
  
sink('stats/mags_checkm2_summary.txt')
print(paste0('Median (IQR) MAG length: ', median_length, 'bp (', iqr_length, 'bp)'))
print(paste0('Median (IQR) MAG completeness: ', median_comp, '% (', iqr_comp, '%)'))
print(paste0('Median (IQR) MAG contamination: ', median_cont, '% (', iqr_cont, '%)'))

print(paste0('Number of high quality MAGs (>90% completeness, <5% contamination): ', nhq))
print(paste0('Number of medium quality MAGs (>50% completeness, <10% contamination): ', nmq))
print(paste0('Number of low quality MAGs (>40% completeness): ', nlq))
sink()

# Taxonomic trees
taxmap_obj_mags <- parse_tax_data(
  mag_data,
  class_cols = "Taxonomy", # the column that contains taxonomic information
  class_sep = ";", # The character used to separate taxa in the classification
  class_regex = "^(.*)__(.*)$",
  class_key = c(tax_rank = "info",
                tax_name = "taxon_name"))

taxmap_obj_mags$data$tax_abund = calc_taxon_abund(taxmap_obj_mags, "tax_data", cols = c("mean_rel_ab"))
panel_left = taxmap_obj_mags %>% 
heat_tree(node_label = taxon_names,
          node_size = n_obs,
          node_color = taxmap_obj_mags$data$tax_abund$mean_rel_ab,
          node_size_axis_label = "MAG count",
          node_color_axis_label = "Mean relative abundance",
          layout = "da", # The primary layout algorithm
          initial_layout = "fr",  overlap_avoidance = 11,
          node_size_range = c(0.012, 0.04), aspect_ratio=0.75) # The layout algorithm that initializes node locations


# 16S near full-length
taxmap_obj_amplicon <- parse_tax_data(
  a16s_data,
  class_cols = "Taxonomy", # the column that contains taxonomic information
  class_sep = ";", # The character used to separate taxa in the classification
  class_regex = "^(.*)__(.*)$",
  class_key = c(tax_rank = "info",
                tax_name = "taxon_name"))

taxmap_obj_amplicon$data$tax_abund = calc_taxon_abund(taxmap_obj_amplicon, "tax_data", cols = c("mean_rel_ab"))
panel_left_16s = taxmap_obj_amplicon %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = taxmap_obj_amplicon$data$tax_abund$mean_rel_ab,
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Mean relative abundance",
            layout = "da", # The primary layout algorithm
            initial_layout = "fr",  overlap_avoidance = 11,
            node_size_range = c(0.012, 0.04), aspect_ratio=0.75) # The layout algorithm that initializes node locations


# 2.2 Temperature effect on alpha-diversity (MAGs)
min_non_zero = min(loaded_mags$raw_table[loaded_mags$raw_table > 0])
norm_table_cond = as.data.frame(lapply(loaded_mags$raw_table , function(x) as.integer(x / min_non_zero)))
norm_table_cond = norm_table_cond[,colnames(norm_table_cond) %in% metadata$Sample]
norm_table_cond <- norm_table_cond[, match(metadata$Sample, colnames(norm_table_cond))]
all(metadata$Sample == colnames(norm_table_cond))

shannon = vegan::diversity(t(norm_table_cond))
invsimp = vegan::diversity(t(norm_table_cond), index = 'invsimpson')
obs = colSums(norm_table_cond > 0)

metadata$Shannon_MAGs = shannon
metadata$Observed_MAGs = obs
metadata$InvSimp_MAGs = invsimp

metadata$Replicate = as.factor(metadata$Replicate)
mod_shannon = gam(formula = Shannon_MAGs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Aspect, 
                  data=metadata)
sink('stats/figure2_gams_shannon_soiltemp.txt')
summary(mod_shannon)
sink()

mod_observed = gam(formula = Observed_MAGs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Aspect, 
                   data=metadata)
sink('stats/figure2_gams_obs_rich_soiltemp.txt')
summary(mod_observed)
sink()

mod_invsimp = gam(formula = InvSimp_MAGs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Aspect, 
                   data=metadata)
sink('stats/figure2_gams_inv_simp_soiltemp.txt')
summary(mod_invsimp)
sink()

# Plotting
pred_grid <- expand.grid(
  Days_Since_Start = unique(metadata$Days_Since_Start),
  Replicate = unique(metadata$Replicate), 
  Aspect = unique(metadata$Aspect),
  Temp_soil = seq(-1, 25, length.out = 50)
)

pred_vals <- predict(mod_shannon, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit_shannon <- pred_vals$fit
pred_grid$se_shannon  <- pred_vals$se.fit

pred_vals <- predict(mod_observed, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit_observed <- pred_vals$fit
pred_grid$se_observed  <- pred_vals$se.fit

pred_vals <- predict(mod_invsimp, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit_invsimp <- pred_vals$fit
pred_grid$se_invsimp  <- pred_vals$se.fit

# upper/lower bounds for SE ribbon
pred_grid <- pred_grid %>%
  mutate(
    upper_shannon = fit_shannon + 1.96*se_shannon,
    lower_shannon = fit_shannon - 1.96*se_shannon
  )

pred_grid <- pred_grid %>%
  mutate(
    upper_observed = fit_observed + 1.96*se_observed,
    lower_observed = fit_observed - 1.96*se_observed
  )

pred_grid <- pred_grid %>%
  mutate(
    upper_invsimp = fit_invsimp + 1.96*se_invsimp,
    lower_invsimp = fit_invsimp - 1.96*se_invsimp
  )

pred_median <- pred_grid %>%
  group_by(Temp_soil) %>%
  summarise(
    fit_shannon    = mean(fit_shannon),
    upper_shannon  = mean(upper_shannon),
    lower_shannon  = mean(lower_shannon),
    fit_observed    = mean(fit_observed),
    upper_observed  = mean(upper_observed),
    lower_observed  = mean(lower_observed),
    fit_invsimp    = mean(fit_invsimp),
    upper_invsimp  = mean(upper_invsimp),
    lower_invsimp  = mean(lower_invsimp),
    .groups = "drop"
  )


a = ggplot() +
  # raw data points in the background
  geom_point(data = metadata, 
             aes(x = Temp_soil, y = Shannon_MAGs), 
             alpha = 0.5, colour='darkred') +
  # SE ribbon
  geom_ribbon(data = pred_median, 
              aes(x = Temp_soil, ymin = lower_shannon, ymax = upper_shannon), 
              alpha = 0.2, fill = 'dimgrey') +
  # median predicted line
  geom_line(data = pred_median, 
            aes(x = Temp_soil, y = fit_shannon), colour='black', size = 1) +
  labs(x = "Soil temperature [°C]", y = "Shannon Index") +
  theme_minimal()


b = ggplot() +
  # raw data points in the background
  geom_point(data = metadata, 
             aes(x = Temp_soil, y = Observed_MAGs), 
             alpha = 0.5, colour='royalblue') +
  # SE ribbon
  geom_ribbon(data = pred_median, 
              aes(x = Temp_soil, ymin = lower_observed, ymax = upper_observed), 
              alpha = 0.2, fill = 'dimgrey') +
  # median predicted line
  geom_line(data = pred_median, 
            aes(x = Temp_soil, y = fit_observed), colour='black', size = 1) +
  labs(x = "Soil temperature [°C]", y = "Observed MAGs") +
  theme_minimal()

c = ggplot() +
  # raw data points in the background
  geom_point(data = metadata, 
             aes(x = Temp_soil, y = InvSimp_MAGs), 
             alpha = 0.5, colour='forestgreen') +
  # SE ribbon
  geom_ribbon(data = pred_median, 
              aes(x = Temp_soil, ymin = lower_invsimp, ymax = upper_invsimp), 
              alpha = 0.2, fill = 'dimgrey') +
  # median predicted line
  geom_line(data = pred_median, 
            aes(x = Temp_soil, y = fit_invsimp), colour='black', size = 1) +
  labs(x = "Soil temperature [°C]", y = "Inverted Simpson") +
  theme_minimal()

panel_right = ggarrange(a, c, b, nrow = 3, ncol = 1, labels = c('B', 'C', 'D'))
both = ggarrange(panel_left, panel_right, widths = c(0.7, 0.3), ncol = 2, nrow = 1, labels = c('A', ''))
ggsave(filename = 'figures/Figure2.pdf', plot = both, width = 12, height = 8)



# 2.2 Temperature effect on alpha-diversity (16S amplicons)
norm_table_cond = loaded_16s$rar_table
metadata_2024 = metadata %>% filter(Year == 2024)
norm_table_cond = norm_table_cond[,colnames(norm_table_cond) %in% metadata_2024$Sample]
norm_table_cond <- norm_table_cond[, match(metadata_2024$Sample, colnames(norm_table_cond))]
all(metadata_2024$Sample == colnames(norm_table_cond))

shannon = vegan::diversity(t(norm_table_cond))
invsimp = vegan::diversity(t(norm_table_cond), index = 'invsimpson')
obs = colSums(norm_table_cond > 0)

metadata_2024$Shannon_OTUs = shannon
metadata_2024$Observed_OTUs = obs
metadata_2024$InvSimp_OTUs = invsimp

metadata_2024$Replicate = as.factor(metadata_2024$Replicate)
mod_shannon = gam(formula = Shannon_OTUS ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Aspect, 
                  data=metadata_2024)
sink('stats/16s_gams_shannon_soiltemp.txt')
summary(mod_shannon)
sink()

mod_observed = gam(formula = Observed_OTUs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Aspect, 
                   data=metadata_2024)
sink('stats/16s_gams_obs_rich_soiltemp.txt')
summary(mod_observed)
sink()

mod_invsimp = gam(formula = InvSimp_OTUs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Aspect, 
                  data=metadata_2024)
sink('stats/16s_gams_inv_simp_soiltemp.txt')
summary(mod_invsimp)
sink()

# Plotting
pred_grid <- expand.grid(
  Days_Since_Start = unique(metadata_2024$Days_Since_Start),
  Replicate = unique(metadata_2024$Replicate), 
  Aspect = unique(metadata_2024$Aspect),
  Temp_soil = seq(-1, 25, length.out = 50)
)

pred_vals <- predict(mod_shannon, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit_shannon <- pred_vals$fit
pred_grid$se_shannon  <- pred_vals$se.fit

pred_vals <- predict(mod_observed, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit_observed <- pred_vals$fit
pred_grid$se_observed  <- pred_vals$se.fit

pred_vals <- predict(mod_invsimp, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit_invsimp <- pred_vals$fit
pred_grid$se_invsimp  <- pred_vals$se.fit

# upper/lower bounds for SE ribbon
pred_grid <- pred_grid %>%
  mutate(
    upper_shannon = fit_shannon + 1.96*se_shannon,
    lower_shannon = fit_shannon - 1.96*se_shannon
  )

pred_grid <- pred_grid %>%
  mutate(
    upper_observed = fit_observed + 1.96*se_observed,
    lower_observed = fit_observed - 1.96*se_observed
  )

pred_grid <- pred_grid %>%
  mutate(
    upper_invsimp = fit_invsimp + 1.96*se_invsimp,
    lower_invsimp = fit_invsimp - 1.96*se_invsimp
  )

pred_median <- pred_grid %>%
  group_by(Temp_soil) %>%
  summarise(
    fit_shannon    = mean(fit_shannon),
    upper_shannon  = mean(upper_shannon),
    lower_shannon  = mean(lower_shannon),
    fit_observed    = mean(fit_observed),
    upper_observed  = mean(upper_observed),
    lower_observed  = mean(lower_observed),
    fit_invsimp    = mean(fit_invsimp),
    upper_invsimp  = mean(upper_invsimp),
    lower_invsimp  = mean(lower_invsimp),
    .groups = "drop"
  )


a = ggplot() +
  # raw data points in the background
  geom_point(data = metadata_2024, 
             aes(x = Temp_soil, y = Shannon_OTUs), 
             alpha = 0.5, colour='darkred') +
  # SE ribbon
  geom_ribbon(data = pred_median, 
              aes(x = Temp_soil, ymin = lower_shannon, ymax = upper_shannon), 
              alpha = 0.2, fill = 'dimgrey') +
  # median predicted line
  geom_line(data = pred_median, 
            aes(x = Temp_soil, y = fit_shannon), colour='black', size = 1) +
  labs(x = "Soil temperature [°C]", y = "Shannon Index") +
  theme_minimal()


b = ggplot() +
  # raw data points in the background
  geom_point(data = metadata_2024, 
             aes(x = Temp_soil, y = Observed_OTUs), 
             alpha = 0.5, colour='royalblue') +
  # SE ribbon
  geom_ribbon(data = pred_median, 
              aes(x = Temp_soil, ymin = lower_observed, ymax = upper_observed), 
              alpha = 0.2, fill = 'dimgrey') +
  # median predicted line
  geom_line(data = pred_median, 
            aes(x = Temp_soil, y = fit_observed), colour='black', size = 1) +
  labs(x = "Soil temperature [°C]", y = "Observed OTUs") +
  theme_minimal()

c = ggplot() +
  # raw data points in the background
  geom_point(data = metadata_2024, 
             aes(x = Temp_soil, y = InvSimp_OTUs), 
             alpha = 0.5, colour='forestgreen') +
  # SE ribbon
  geom_ribbon(data = pred_median, 
              aes(x = Temp_soil, ymin = lower_invsimp, ymax = upper_invsimp), 
              alpha = 0.2, fill = 'dimgrey') +
  # median predicted line
  geom_line(data = pred_median, 
            aes(x = Temp_soil, y = fit_invsimp), colour='black', size = 1) +
  labs(x = "Soil temperature [°C]", y = "Inverted Simpson") +
  theme_minimal()

panel_right = ggarrange(a, c, b, nrow = 3, ncol = 1, labels = c('B', 'C', 'D'))
both = ggarrange(panel_left_16s, panel_right, widths = c(0.7, 0.3), ncol = 2, nrow = 1, labels = c('A', ''))
ggsave(filename = 'figures/Figure_S??_alpha_16s.pdf', plot = both, width = 12, height = 8)





###### SPEC.GEN ANALYSIS
library(EcolUtils)
lev <- spec.gen(t(norm_table_cond), n=999, niche.width.method = 'levins')
table(lev$sign) # 26 Non Significant and 539 Specialists
lev$MAG = rownames(loaded_mags$relab_table)

spec_table = as.matrix(loaded_mags$relab_table) %>% melt() %>% filter(Var2 %in% colnames(norm_table_cond))
spec_table$Spec = map_chr(spec_table$Var1, function(x) as.character(lev$sign[lev$MAG == x]))
spec_table$Soil_temp = map_dbl(spec_table$Var2, function(x) as.numeric(metadata$Temp_soil[metadata$Sample == x]))
spec_table$Days_Since_Start = map_dbl(spec_table$Var2, function(x) metadata$Days_Since_Start[metadata$Sample == x])
spec_table$Replicate = map_chr(spec_table$Var2, function(x) as.character(metadata$Replicate[metadata$Sample == x]))
spec_table$Aspect = map_chr(spec_table$Var2, function(x) as.character(metadata$Aspect[metadata$Sample == x]))

spec_table_summary = spec_table %>% filter(Spec == 'SPECIALIST') %>% group_by(Var2, Replicate, Days_Since_Start, Soil_temp, Aspect) %>% summarise(rel_ab = sum(value))

spec_table_summary$Replicate = as.factor(spec_table_summary$Replicate)
mod_specialists = gam(formula = rel_ab ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Soil_temp, bs='ts', k=5), 
                  data=spec_table_summary)
summary(mod_specialists)


###### RAC ANALYSIS
relab_table = loaded_mags$relab_table[,colnames(loaded_mags$relab_table) %in% metadata$Sample]
relab_table <- relab_table[, match(metadata$Sample, colnames(relab_table))]
all(metadata$Sample == colnames(relab_table))

# Long format
rad_df <- relab_table %>%
  as.data.frame() %>%
  rownames_to_column(var = "MAG") %>%
  pivot_longer(-MAG, names_to = "Sample", values_to = "Abundance") %>%
  left_join(metadata, by = "Sample") %>%
  group_by(Sample) %>%
  filter(Abundance > 0) %>%
  mutate(
    rel_abund = Abundance / sum(Abundance),
    Rank = rank(-rel_abund, ties.method = "first"),
    log_rel_abund = log10(rel_abund + 1e-6)
  ) %>%
  ungroup()

library(dplyr)
library(broom)
boot_spearman <- function(df, n_boot = 1000) {
  replicate(n_boot, {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    cor(
      df$log_rel_abund[idx],
      df$Temp_soil[idx],
      method = "spearman"
    )
  })
}

n_samples <- n_distinct(rad_df$Sample)

rank_boot <- rad_df %>%
  group_by(Rank) %>%
  filter(n_distinct(Sample) >= 0.5 * n_samples) %>%
  group_modify(~ {
    rho_boot <- boot_spearman(.x, n_boot = 1000)
    
    tibble(
      rho = median(rho_boot, na.rm = TRUE),
      lower = quantile(rho_boot, 0.025, na.rm = TRUE),
      upper = quantile(rho_boot, 0.975, na.rm = TRUE)
    )
  }) %>%
  ungroup()

ggplot(rank_boot, aes(x = Rank, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "grey80",
    alpha = 0.6
  ) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Rank (within sample)",
    y = expression("Spearman's " * rho * " with soil temperature"),
    title = "Rank-wise temperature response"
  ) +
  theme_classic()
ggsave('figures/SuppFigure3_RAC.pdf')




























# Compute relative abundance and rank per sample
rad_df <- rad_df %>%
  group_by(Sample) %>% filter(Abundance > 0) %>%
  mutate(
    rel_abund = Abundance / sum(Abundance),
    Rank = rank(-rel_abund, ties.method = "first"),
    log_rel_abund = log10(rel_abund + 1e-6)   # log scale for plotting
  ) %>%
  ungroup()

rad_df$Condition = factor(rad_df$Condition, levels = c('control', 'warming'))
# GAM with smooth for Rank, separate by Condition, random effect for Sample
gam_model <- gam(
  log_rel_abund ~ s(Rank, bs='ts', k=5, by=Condition)+
    s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs = 're'),
  data = rad_df
)

summary(gam_model)

pred_out <- predict(
  gam_model,
  newdata = rad_df,
  se.fit = TRUE,
  exclude = "s(Days_Since_Start)"   # ignore replicate-by-day smooth
)

rad_df$pred <- pred_out$fit

rad_summary <- rad_df %>% 
  group_by(Condition, Rank) %>%
  summarise(
    median_pred = median(pred),
    lower = quantile(pred, 0.25),
    upper = quantile(pred, 0.75),
    .groups = "drop"
  )


pA <- ggplot(rad_summary, aes(x = Rank, y = median_pred, color = Condition, fill = Condition)) +
  geom_point(rad_df, mapping = aes(x = Rank, y = log_rel_abund, color = Condition), alpha=0.01) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    x = "Rank",
    y = "Predicted log10(Relative Abundance)",
  ) +
  theme_classic() + scale_color_manual(values = c('#012169','#C8102E')) +
  scale_fill_manual(values = c('#012169','#C8102E'))

# Assume two conditions: "Control" and "Warming"
diff <- difference_smooths(
  gam_model,
  select = "s(Rank)",
  n = 200
)
diff

# Plot Panel B with IQR ribbon
pB <- ggplot(diff, aes(x = Rank, y = -.diff)) +
  geom_line(size = 1.2, color = "purple") +
  geom_ribbon(aes(ymin = -.lower_ci, ymax = -.upper_ci), alpha = 0.2, fill = "violet") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Rank",
    y = "Difference (Warming - Control)",
  ) +
  theme_classic()

ggarrange(pA, pB, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = T)
ggsave(filename = 'figures/SuppFigure_RAC.pdf', width = 8, height = 6)

















# 2.2 Seasonal effect on alpha-diversity
min_non_zero = min(loaded_mags$raw_table[loaded_mags$raw_table > 0])
norm_table = as.data.frame(lapply(loaded_mags$raw_table , function(x) as.integer(x / min_non_zero)))
norm_table = norm_table[,colnames(norm_table) %in% metadata_timeseries$Sample]
norm_table <- norm_table[, match(metadata_timeseries$Sample, colnames(norm_table))]
all(metadata_timeseries$Sample == colnames(norm_table))

shannon = vegan::diversity(t(norm_table))
invsimp = vegan::diversity(t(norm_table), index = 'invsimpson')
obs = colSums(norm_table > 0)

metadata_timeseries$Shannon_MAGs = shannon
metadata_timeseries$Observed_MAGs = obs
metadata_timeseries$InvSimp_MAGs = invsimp

ggplot(metadata_timeseries, aes(x=Days_Since_Start, y=Shannon_MAGs)) + geom_point() + geom_smooth()
ggplot(metadata_timeseries, aes(x=Days_Since_Start, y=Observed_MAGs)) + geom_point() + geom_smooth()
ggplot(metadata_timeseries, aes(x=Days_Since_Start, y=InvSimp_MAGs)) + geom_point() + geom_smooth()

mod_shannon = gam(formula = Shannon_MAGs ~ s(DOY, k=5, bs = "ts") + s(Replicate, bs='re') + Year + Aspect, 
                  data=metadata_timeseries)
sink('stats/model_shannon_timeseries.txt')
summary(mod_shannon)
sink()

mod_observed = gam(formula = Observed_MAGs ~ s(DOY, k=5, bs = "ts") + s(Replicate, bs='re') + Year + Aspect, 
                   data=metadata_timeseries)
sink('stats/model_obs_rich_timeseries.txt')
summary(mod_observed)
sink()

mod_invsimp = gam(formula = InvSimp_MAGs ~ s(DOY, k=5, bs = "ts") + s(Replicate, bs='re') + Year + Aspect, 
                   data=metadata_timeseries)
sink('stats/model_inv_simp_timeseries.txt')
summary(mod_invsimp)
sink()

# PLOTTING 
# Unique values
replicates <- unique(metadata_timeseries$Replicate)
years <- unique(metadata_timeseries$Year)

# Sequence of DOY
doy_seq <- 190:310

# Create newdata for predictions
newdata <- expand.grid(
  DOY = doy_seq,
  Year = years,
  Replicate = replicates
)
get_median_pred <- function(model, newdata) {
  pred <- predict(model, newdata = newdata, se.fit = TRUE)
  
  newdata %>%
    mutate(
      fit = pred$fit,
      upr = fit + 2 * pred$se.fit,
      lwr = fit - 2 * pred$se.fit
    ) %>%
    group_by(DOY, Year) %>%
    summarise(
      fit_median = median(fit),
      upr_median = median(upr),
      lwr_median = median(lwr),
      .groups = "drop"
    )
}
pred_shannon <- get_median_pred(mod_shannon, newdata)
pred_observed <- get_median_pred(mod_observed, newdata)
pred_invsimp <- get_median_pred(mod_invsimp, newdata)
plot_median_spline <- function(pred_data, ylab_name) {
  ggplot(pred_data, aes(x = DOY, y = fit_median, color = as.factor(Year), fill = as.factor(Year))) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lwr_median, ymax = upr_median), alpha = 0.2, color = NA) +
    labs(x = "Day of Year", y = ylab_name, color = "Year", fill = "Year") +
    theme_minimal() +
    theme(text = element_text(size = 12)) + scale_colour_manual(values = c('coral4', 'gold')) +
    scale_fill_manual(values = c('coral4', 'gold'))
}
p_shannon <- plot_median_spline(pred_shannon, "Shannon Index")
p_observed <- plot_median_spline(pred_observed, "Observed Richness")
p_invsimp <- plot_median_spline(pred_invsimp, "Inverse Simpson")
ggarrange(
  p_shannon, p_observed, p_invsimp,
  ncol = 1, nrow = 3,  # stack vertically
  common.legend = TRUE, legend = "right"
)
