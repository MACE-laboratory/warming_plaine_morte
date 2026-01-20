source('0_functions_and_packages.R')
setwd('~/Documents/MACE/EVO_study')

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
# Part 1. PCA and environmental drivers 
################################################################################
variables_for_pca = c("Temp_soil_scaled","Temp_air_scaled","Luminosity_scaled","CH4_Flux_scaled",
                      "CO2_Flux_scaled","TC_scaled","TOC_scaled", "pH_scaled","Water_Content_Lab_scaled",
                      "F._scaled","NO3._scaled","Cl._scaled",
                      "NO2._scaled","PO4._scaled","SO4._scaled","Br._scaled","Li._scaled",
                      "NH4._scaled","Ca2._scaled","Sr2._scaled","Na._scaled","Mg2._scaled",
                      "K._scaled","Formate_scaled","Malate_scaled","Propionate_scaled",
                      "Lactate_scaled","Butyrate_scaled","Oxalate_scaled",
                      "Acetate_scaled")

# 1. PCA for condition
meta_cond <- metadata_condition %>%
  select(Sample, Condition, all_of(variables_for_pca))

pca_mat_cond <- meta_cond %>%
  select(-Sample, -Condition) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()
rownames(pca_mat_cond) <- meta_cond$Sample

pca_cond <- prcomp(pca_mat_cond, center = TRUE, scale. = TRUE)
variance_explained_cond <- ((pca_cond$sdev^2) / sum(pca_cond$sdev^2) * 100) [1:2]
pca_plot_df <- cbind(as.data.frame(pca_cond$x[, 1:2]), Sample = meta_cond$Sample, Condition = meta_cond$Condition)

pca_cond_test <- PCAtest(pca_mat_cond, 100, 100, 
                                          0.05, varcorr = TRUE, counter = FALSE, plot = TRUE)

sig_vars_PC1 = c(1, 2, 3, 7, 8, 12, 13, 14, 15, 17, 25, 26, 27, 28, 29)
sig_vars_PC2 = c(6, 10, 12, 15, 17, 18, 20, 21, 22, 23, 28)
sig_vars_cond = unique(c(sig_vars_PC1, sig_vars_PC2))

pca_cond_loadings <- scale_loadings(pca_cond, 1, 2)
pca_cond_loadings$name = map_chr(rownames(pca_cond_loadings), function(x) as.character(PCA_label_conversion[x]))
pca_cond_sig_var_names <- colnames(pca_mat_cond)[sig_vars_cond]
pca_cond_loadings = pca_cond_loadings[pca_cond_sig_var_names,]

pca_cond_plot <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_segment(
    data = pca_cond_loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "darkgrey", alpha = 0.75, inherit.aes = FALSE
  ) +
  stat_ellipse(level = 0.95, linewidth = 0.5, linetype = 'dashed') +
  geom_point(size = 3) +
  geom_text_repel(
    data = pca_cond_loadings,
    aes(x = PC1, y = PC2, label = name),
    size = 3, inherit.aes = FALSE, parse = TRUE
  ) +
  labs(
    x = paste0("PC1 (", round(variance_explained_cond[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained_cond[2], 2), "%)"),
    color = "Condition"
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_manual(values = c('#012169', '#C8102E'))
pca_cond_plot


# 2. PCA for timeseries
meta_time <- metadata_timeseries %>%
  select(Sample, Days_Since_Start, all_of(variables_for_pca)) %>% select(-CH4_Flux_scaled, -CO2_Flux_scaled)

pca_mat_time <- meta_time %>%
  select(-Sample, -Days_Since_Start) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()
rownames(pca_mat_time) <- meta_time$Sample

pca_time <- prcomp(pca_mat_time, center = TRUE, scale. = TRUE)
variance_explained_time <- ((pca_time$sdev^2) / sum(pca_time$sdev^2) * 100) [1:2]
pca_plot_df <- cbind(as.data.frame(pca_time$x[, 1:2]), Sample = meta_time$Sample, Days_Since_Start = meta_time$Days_Since_Start)

pca_time_test <- PCAtest(pca_mat_time, 100, 100, 
                         0.05, varcorr = TRUE, counter = FALSE, plot = TRUE)

sig_vars_PC1 = c(1, 2, 5, 6, 8, 12, 13, 14, 15, 16, 17, 18, 20, 21, 27)
sig_vars_PC2 = c(8, 9, 10, 11, 13, 16, 19, 21, 23, 24, 26, 27)
sig_vars_time = unique(c(sig_vars_PC1, sig_vars_PC2))

pca_time_loadings <- scale_loadings(pca_time, 1, 2)
pca_time_loadings$name = map_chr(rownames(pca_time_loadings), function(x) as.character(PCA_label_conversion[x]))
pca_time_sig_var_names <- colnames(pca_mat_time)[sig_vars_time]
pca_time_loadings = pca_cond_loadings[pca_time_sig_var_names,]

pca_time_plot <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Days_Since_Start)) +
  geom_segment(
    data = pca_time_loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "darkgrey", alpha = 0.75, inherit.aes = FALSE
  ) +
  geom_point(size = 3) +
  geom_text_repel(
    data = pca_time_loadings,
    aes(x = PC1, y = PC2, label = name),
    size = 3, inherit.aes = FALSE, parse = TRUE
  ) +
  labs(
    x = paste0("PC1 (", round(variance_explained_time[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained_time[2], 2), "%)"),
    color = "Days since the start of the experiment"
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_gradient(low = 'gold', high = 'violet')
pca_time_plot
ggsave(filename = 'figures/SuppFigure1.pdf', pca_time_plot, width = 10.5)

library(Hmisc)
# GAMs for env parameters, testing Condition and time
variables_for_gams = c("Temp_soil","TC","Water_Content_Lab","log(TIN)")

# A) Models Condition
metadata_condition <- metadata_condition %>%
  mutate(
    Condition = factor(Condition),
    Replicate = factor(Replicate)
  )

fit_gam_extract <- function(response, data) {
  
  formula_gam <- as.formula(
    paste0(
      response,
      "~ s(Days_Since_Start, bs='ts', k = 5) + s(Replicate, bs='re') + Condition"
    )
  )
  
  m <- gam(
    formula_gam,
    data = data,
    method = "REML"
  )
  
  disp <- sum(residuals(m, type = "pearson")^2) /
    df.residual(m)
  
  sm <- summary(m)

  # Extract Condition effect
  cond_table <- sm$p.table
  cond_rows <- grep("^Condition", rownames(cond_table))
  
  tibble(
    variable = response,
    dispersion = disp,
    deviance_explained = sm$dev.expl,
    condition_p_value = cond_table[cond_rows, "Pr(>|t|)"],
    condition_estimate = cond_table[cond_rows, "Estimate"],
    condition_se = cond_table[cond_rows, "Std. Error"]
  )
}

gam_results <- map_dfr(
  variables_for_gams,
  fit_gam_extract,
  data = metadata_condition
)

gam_results$condition_padj = p.adjust(gam_results$condition_p_value)
gam_results %>% filter(condition_padj < 0.05)
gam_results
write.csv(gam_results, file='stats/figure1_gams_condition.csv')


# B) Models year and season
metadata_timeseries <- metadata_timeseries %>%
  mutate(
    Condition = factor(Condition),
    Replicate = factor(Replicate)
  )

fit_gam_extract <- function(response, data, n_doy = 365) {
  
  formula_gam <- as.formula(
    paste0(
      response,
      " ~ s(DOY, bs = 'ts', k = 5) + ",
      "s(Replicate, bs = 're') + Year"
    )
  )
  
  m <- mgcv::gam(
    formula_gam,
    data = data,
    method = "REML"
  )
  
  disp <- sum(residuals(m, type = "pearson")^2) /
    df.residual(m)
  
  sm <- summary(m)
  
  ## ---- season (DOY spline): p-value ----
  season_row <- sm$s.table[grep("^s\\(DOY", rownames(sm$s.table)), , drop = FALSE]
  
  ## ---- year fixed effect ----
  year_row <- sm$p.table["Year", , drop = FALSE]
  
  ## ---- effect size: predicted seasonal range ----
  doy_seq <- seq(
    min(data$DOY, na.rm = TRUE),
    max(data$DOY, na.rm = TRUE),
    length.out = n_doy
  )
  
  newdat <- data.frame(
    DOY = doy_seq,
    Year = mean(data$Year, na.rm = TRUE),
    Replicate = data$Replicate[1]  # arbitrary; RE smooth is centred at 0
  )
  
  pred <- predict(m, newdata = newdat, type = "terms")
  
  season_range <- diff(range(pred[, "s(DOY)"]))
  
  tibble::tibble(
    variable = response,
    
    ## model fit
    deviance_explained = sm$dev.expl,
    
    dispersion = disp,

    ## season effect
    season_effect_range = season_range,
    season_p_value      = season_row[, "p-value"],
    
    ## year effect
    year_estimate = year_row[, "Estimate"],
    year_se       = year_row[, "Std. Error"],
    year_p_value  = year_row[, "Pr(>|t|)"]
  )
}


gam_results <- map_dfr(
  variables_for_gams,
  fit_gam_extract,
  data = metadata_timeseries
)

gam_results$season_effect_padj = p.adjust(gam_results$season_p_value)
gam_results$year_effect_padj = p.adjust(gam_results$year_p_value)

gam_results %>% filter((year_effect_padj < 0.05) | (season_effect_padj < 0.05))
gam_results
write.csv(gam_results, file='stats/figure1_gams_timeseries.csv')



# Figure 1: PCA + difference in temperatures
mod_temp_time <- mgcv::gam(
  Temp_soil ~ s(Days_Since_Start, bs = 'ts', k = 5) + s(Replicate, bs = 're') + Condition,
  data = metadata_condition,
  method = "REML"
)
# Expand grid over DOY, Year, and Replicate
newdata <- expand.grid(
  Days_Since_Start = 90:480,
  Condition = unique(metadata_condition$Condition),
  Replicate = unique(metadata_condition$Replicate)
)

pred <- predict(mod_temp_time, newdata = newdata, se.fit = TRUE)
newdata <- newdata %>%
  mutate(
    fit = pred$fit,
    upr = fit + 2 * pred$se.fit,
    lwr = fit - 2 * pred$se.fit
  )
median_pred <- newdata %>%
  group_by(Days_Since_Start, Condition) %>%
  summarise(
    fit_median = median(fit),
    upr_median = median(upr),
    lwr_median = median(lwr),
    .groups = "drop"
  )

b = ggplot(median_pred, aes(x = Days_Since_Start, y = fit_median, color = as.factor(Condition), fill = as.factor(Condition))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr_median, ymax = upr_median), alpha = 0.2, color = NA) +
  labs(
    x = "Day since the start of the experiment",
    y = "Predicted Soil Temperature",
    color = "Condition",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12)
  ) + scale_colour_manual(values = c('#012169', '#C8102E')) +
  scale_fill_manual(values = c('#012169', '#C8102E'))

mod_temp_time <- mgcv::gam(
  Temp_soil ~ s(DOY, bs = 'ts', k = 5) + s(Replicate, bs = 're') + Year,
  data = metadata_timeseries,
  method = "REML"
)
# Expand grid over DOY, Year, and Replicate
newdata <- expand.grid(
  DOY = 190:310,
  Year = unique(metadata_timeseries$Year),
  Replicate = unique(metadata_timeseries$Replicate)
)

pred <- predict(mod_temp_time, newdata = newdata, se.fit = TRUE)
newdata <- newdata %>%
  mutate(
    fit = pred$fit,
    upr = fit + 2 * pred$se.fit,
    lwr = fit - 2 * pred$se.fit
  )
median_pred <- newdata %>%
  group_by(DOY, Year) %>%
  summarise(
    fit_median = median(fit),
    upr_median = median(upr),
    lwr_median = median(lwr),
    .groups = "drop"
  )

c = ggplot(median_pred, aes(x = DOY, y = fit_median, color = as.factor(Year), fill = as.factor(Year))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr_median, ymax = upr_median), alpha = 0.2, color = NA) +
  labs(
    x = "Day of Year",
    y = "Predicted Soil Temperature",
    color = "Year",
    fill = "Year"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12)
  ) + scale_colour_manual(values = c('coral4', 'gold')) +
  scale_fill_manual(values = c('coral4', 'gold'))


bc = ggarrange(b, c, nrow = 2, ncol = 1, labels = c('B', 'C'))
fig = ggarrange(pca_cond_plot + theme(legend.position = 'none'), bc, ncol = 2, nrow = 1, labels = c('A',''))
fig
ggsave(filename = 'figures/Figure1.pdf', fig, width = 10.5, height = 6.5)