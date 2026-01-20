source('0_functions_and_packages.R')
setwd('~/Documents/MACE/EVO_study')

# loading functions
loaded_mags = load_mags_data()
loaded_meta = load_meta_data()
loaded_16sr = load_16sr_data()
metadata = load_metadata()
metadata = add_shannon(metadata, loaded_mags, loaded_meta, loaded_16sr)
metadata = check_normality_change(metadata, c("Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
                                              "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",              
                                              "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",          
                                              "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
                                              "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
                                              "Other_Cations","Grouped_OC"))

metadata_condition = metadata %>% filter(Days_Since_Start >= 90)
metadata_timeseries = metadata %>% filter(Condition == 'control')

scale_loadings <- function(pca, PCx, PCy){
  
  max_score <- max(abs(pca$x[, PCx]), abs(pca$x[, PCy]))
  
  max_loading <- max(abs(pca$rotation[, PCx]), abs(pca$rotation[, PCy]))
  
  scale_factor <- max_score / max_loading
  
  PCA_loadings_scaled <- as.data.frame(pca$rotation[, c(PCx, PCy)] * scale_factor)
  
  return(PCA_loadings_scaled)
  
}

################################################################################
# Part 2. Microbiome - alpha diversity
################################################################################
# 2.1 Condition effect on alpha-diversity
min_non_zero = min(loaded_mags$raw_table[loaded_mags$raw_table > 0])
norm_table_cond = as.data.frame(lapply(loaded_mags$raw_table , function(x) as.integer(x / min_non_zero)))
norm_table_cond = norm_table_cond[,colnames(norm_table_cond) %in% metadata_condition$Sample]
norm_table_cond <- norm_table_cond[, match(metadata_condition$Sample, colnames(norm_table_cond))]
all(metadata_condition$Sample == colnames(norm_table_cond))

shannon = vegan::diversity(t(norm_table_cond))
invsimp = vegan::diversity(t(norm_table_cond), index = 'invsimpson')
obs = colSums(norm_table_cond > 0)

metadata_condition$Shannon_MAGs = shannon
metadata_condition$Observed_MAGs = obs
metadata_condition$InvSimp_MAGs = invsimp

ggplot(metadata_condition, aes(x=Condition, y=Shannon_MAGs)) + geom_boxplot() + stat_compare_means()
ggplot(metadata_condition, aes(x=Condition, y=Observed_MAGs)) + geom_boxplot() + stat_compare_means()
ggplot(metadata_condition, aes(x=Condition, y=InvSimp_MAGs)) + geom_boxplot() + stat_compare_means()

mod_shannon = gam(formula = Shannon_MAGs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + Condition, 
                  data=metadata_condition)
sink('stats/model_shannon_condition.txt')
summary(mod_shannon)
sink()

mod_observed = gam(formula = Observed_MAGs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + Condition, 
                   data=metadata_condition)
sink('stats/model_obs_rich_condition.txt')
summary(mod_observed)
sink()

mod_invsimp = gam(formula = InvSimp_MAGs ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + Condition, 
                   data=metadata_condition)
sink('stats/model_inv_simp_condition.txt')
summary(mod_invsimp)
sink()

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

mod_shannon = gam(formula = Shannon_MAGs ~ s(DOY, k=5, bs = "ts") + s(Replicate, bs='re') + Year, 
                  data=metadata_timeseries)
sink('stats/model_shannon_timeseries.txt')
summary(mod_shannon)
sink()

mod_observed = gam(formula = Observed_MAGs ~ s(DOY, k=5, bs = "ts") + s(Replicate, bs='re') + Year, 
                   data=metadata_timeseries)
sink('stats/model_obs_rich_timeseries.txt')
summary(mod_observed)
sink()

mod_invsimp = gam(formula = InvSimp_MAGs ~ s(DOY, k=5, bs = "ts") + s(Replicate, bs='re') + Year, 
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
ggsave(filename = 'figures/alpha_div_timeseries.pdf', width = 6, height = 8)
