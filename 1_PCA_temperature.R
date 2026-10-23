################################################################################
# 1_PCA_temperature.R
#
# Environmental drivers of the experiment (Figure 1, SFig 1-3):
#   Part 1  PCA of the environmental variables (warming vs control, time series,
#           North vs South aspect) + PCAtest to choose which loadings to draw
#   Part 2  GAMs of key soil variables: warming effect, and year / season effect
#   Part 3  Soil-temperature trajectories and CH4 / CO2 flux comparisons
#           (van Elteren test, stratified by sampling date) -> Figure 1
#
# Inputs : data/SiteL23_24_env_variables_final.csv
# Outputs: figures/Fig_1, SFig_1, SFig_2, SFig_3, SFig_22-24   (.pdf + .png)
#          stats/01_*.csv|txt
################################################################################
setwd('~/Desktop/warming_plaine_morte')
source('0_functions_and_packages.R')   # run from the repository root

# Extra packages used only here
library(rstatix)    # identify_outliers()
library(ggsignif)
library(coin)       # van Elteren test (coin::wilcox_test with strata)

# Metadata with transformed/scaled environmental variables (+ Aspect, OTC) and
# the two analysis subsets:
#   condition  : samples >= 90 days after the start (warming vs control)
#   timeseries : control plots only (seasonal / yearly variation)
metadata <- prepare_metadata()
subsets <- subset_metadata(metadata)
metadata_condition  <- subsets$condition
metadata_timeseries <- subsets$timeseries

# Rotate loadings so that they span the same range as the sample scores
# (standard "scaled biplot" trick for plotting arrows on top of the scores)
scale_loadings <- function(pca, PCx, PCy){

  max_score <- max(abs(pca$x[, PCx]), abs(pca$x[, PCy]))

  max_loading <- max(abs(pca$rotation[, PCx]), abs(pca$rotation[, PCy]))

  scale_factor <- max_score / max_loading

  PCA_loadings_scaled <- as.data.frame(pca$rotation[, c(PCx, PCy)] * scale_factor)

  return(PCA_loadings_scaled)

}

# Axis label with the variance explained, e.g. "PC1 (22.3%)" (computed, not typed)
pc_label <- function(pca, pc) {
  var_expl <- summary(pca)$importance["Proportion of Variance", pc] * 100
  paste0("PC", pc, " (", formatC(var_expl, format = "f", digits = 1), "%)")
}

# Write PCA outputs (variance explained + full loadings) to stats/
save_pca_stats <- function(pca, name) {
  save_stat_table(as.data.frame(t(summary(pca)$importance)) %>% rownames_to_column("PC"),
                  paste0("01_pca_", name, "_variance_explained"))
  save_stat_table(as.data.frame(pca$rotation) %>% rownames_to_column("variable"),
                  paste0("01_pca_", name, "_loadings"))
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

pca_cond = prcomp(pca_mat_cond)
save_pca_stats(pca_cond, "condition")
pca_plot_df <- cbind(as.data.frame(pca_cond$x[, 1:2]), Sample = meta_cond$Sample, Condition = meta_cond$Condition, Aspect=metadata_condition$Aspect)

# PCAtest: permutation test of the PCA (100 bootstraps, 100 permutations, alpha = 0.05).
# It draws diagnostic plots, saved here instead of going to Rplots.pdf, and the
# console report is saved to stats/.
# The loadings drawn on the biplot are the variables PCAtest flags as significant
# on PC1 / PC2. The indices below (positions in `variables_for_pca`) were read
# from that report: RE-CHECK THEM against stats/01_pcatest_condition.txt if the
# data or the PCAtest version change.
set.seed(GLOBAL_SEED)   # PCAtest is permutation-based
pcatest_out <- utils::capture.output(
  save_base_fig({ pca_cond_test <- PCAtest(pca_mat_cond, 100, 100, 0.05, varcorr = TRUE,
                                           counter = FALSE, plot = TRUE) },
                'sfig_pcatest_condition', width = 10, height = 10))
save_stat_lines("01_pcatest_condition", pcatest_out)

sig_vars_PC1 = c(1, 2, 7, 8, 12, 13, 14, 15, 17, 25, 26, 27, 28, 29)
sig_vars_PC2 = c(6, 10, 12, 15, 17, 18, 21, 22, 23, 28)
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
    x = pc_label(pca_cond, 1),
    y = pc_label(pca_cond, 2),
    color = "Condition" 
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_manual(values = c('#012169', '#C8102E'))
pca_cond_plot    # panel A of Figure 1


# 2. PCA for timeseries
meta_time <- metadata_timeseries %>%
  select(Sample, Days_Since_Start, all_of(variables_for_pca)) %>% select(-CH4_Flux_scaled, -CO2_Flux_scaled)

pca_mat_time <- meta_time %>%
  select(-Sample, -Days_Since_Start) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()
rownames(pca_mat_time) <- meta_time$Sample

pca_time = prcomp(pca_mat_time)
save_pca_stats(pca_time, "timeseries")
pca_plot_df <- cbind(as.data.frame(pca_time$x[, 1:2]), Sample = meta_time$Sample, Days_Since_Start = meta_time$Days_Since_Start)

set.seed(GLOBAL_SEED)   # PCAtest is permutation-based
pcatest_out <- utils::capture.output(
  save_base_fig({ pca_time_test <- PCAtest(pca_mat_time, 100, 100, 0.05, varcorr = TRUE,
                                           counter = FALSE, plot = TRUE) },
                'sfig_pcatest_time', width = 10, height = 10))
save_stat_lines("01_pcatest_timeseries", pcatest_out)

# indices read from the PCAtest report (see note above); flux variables are excluded here
sig_vars_PC1 = c(1, 2, 5, 6, 8, 12, 14, 15, 16, 17, 18, 20, 21, 27)
sig_vars_PC2 = c(8, 9, 10, 11, 13, 16, 17, 19, 21, 23, 24, 26, 27)
sig_vars_time = unique(c(sig_vars_PC1, sig_vars_PC2))

pca_time_loadings <- scale_loadings(pca_time, 1, 2)
pca_time_loadings$name = map_chr(rownames(pca_time_loadings), function(x) PCA_label_conversion[x])
pca_time_sig_var_names <- colnames(pca_mat_time)[sig_vars_time]
pca_time_loadings = pca_time_loadings[pca_time_sig_var_names,]

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
    x = pc_label(pca_time, 1),
    y = pc_label(pca_time, 2),
    color = "Days since the start of the experiment"
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_gradient(low = 'gold', high = 'violet') + xlim(-7.5, 6.5) + ylim(-6.5, 6.5)
pca_time_plot
save_fig(pca_time_plot, 'sfig_pca_time', width = 10.5, height = 7)


# 3. PCA for aspect
meta_aspect <- metadata %>%
  select(Sample, Days_Since_Start, all_of(variables_for_pca), Aspect) %>% select(-CH4_Flux_scaled, -CO2_Flux_scaled)

pca_mat_aspect <- meta_aspect %>%
  select(-Sample, -Days_Since_Start, -Aspect) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()
rownames(pca_mat_aspect) <- meta_aspect$Sample

pca_aspect = prcomp(pca_mat_aspect)
save_pca_stats(pca_aspect, "aspect")
pca_plot_df <- cbind(as.data.frame(pca_aspect$x[, 1:2]), Sample = meta_aspect$Sample, 
                     Days_Since_Start = meta_aspect$Days_Since_Start, Aspect = meta_aspect$Aspect)

# NOTE: the original script re-ran PCAtest on `pca_mat_time` here (copy-paste);
# it now tests the matrix that is actually plotted (`pca_mat_aspect`).
# VERIFY that the indices below still match stats/01_pcatest_aspect.txt.
set.seed(GLOBAL_SEED)   # PCAtest is permutation-based
pcatest_out <- utils::capture.output(
  save_base_fig({ pca_aspect_test <- PCAtest(pca_mat_aspect, 100, 100, 0.05, varcorr = TRUE,
                                             counter = FALSE, plot = TRUE) },
                'sfig_pcatest_aspect', width = 10, height = 10))
save_stat_lines("01_pcatest_aspect", pcatest_out)

sig_vars_PC1 = c(1, 2, 5, 6, 12, 13, 14, 15, 16, 17, 18, 20, 21, 27)
sig_vars_PC2 = c(8, 9, 10, 11, 13, 16, 17, 19, 21, 23, 24, 26, 27)
sig_vars_aspect = unique(c(sig_vars_PC1, sig_vars_PC2))

pca_aspect_loadings <- scale_loadings(pca_aspect, 1, 2)
pca_aspect_loadings$name = map_chr(rownames(pca_aspect_loadings), function(x) as.character(PCA_label_conversion[x]))
pca_aspect_sig_var_names <- colnames(pca_mat_aspect)[sig_vars_aspect]
pca_aspect_loadings = pca_aspect_loadings[pca_aspect_sig_var_names,]

pca_aspect_plot <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Aspect)) +
  geom_segment(
    data = pca_aspect_loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "darkgrey", alpha = 0.75, inherit.aes = FALSE
  ) +
  geom_point(size = 3) + stat_ellipse() +
  geom_text_repel(
    data = pca_aspect_loadings,
    aes(x = PC1, y = PC2, label = name),
    size = 3, inherit.aes = FALSE, parse = TRUE
  ) +
  labs(
    x = pc_label(pca_aspect, 1),
    y = pc_label(pca_aspect, 2),
    color = "Aspect"
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_manual(values = c('turquoise','violet'))
pca_aspect_plot
save_fig(pca_aspect_plot, 'sfig_pca_aspect', width = 10.5, height = 7)




################################################################################
# Part 2. GAMs of soil variables
################################################################################
# GAMs for env parameters, testing Condition and time
variables_for_gams = c("Temp_soil","TC","TOC","pH","Water_Content_Lab")

# A) Warming effect: smooth of time + random replicate + fixed Condition and Aspect
metadata_condition <- metadata_condition %>%
  mutate(
    Condition = factor(Condition),
    Replicate = factor(Replicate),
    Aspect = factor(Aspect)
  )

fit_gam_extract_condition <- function(response, data) {
  
  formula_gam <- as.formula(
    paste0(
      response,
      "~ s(Days_Since_Start, bs='ts', k = 5) + s(Replicate, bs='re') + Condition + Aspect"
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
  fit_gam_extract_condition,
  data = metadata_condition
)

gam_results$condition_padj = p.adjust(gam_results$condition_p_value, method='holm')
gam_results %>% filter(condition_padj < 0.05)
gam_results
save_stat_table(gam_results, '01_env_gams_condition')


# B) Year and season effects in the control plots:
#    smooth of day-of-year + random replicate + fixed Year and Aspect
metadata_timeseries <- metadata_timeseries %>%
  mutate(
    Condition = factor(Condition),
    Replicate = factor(Replicate)
  )

fit_gam_extract_timeseries <- function(response, data, n_doy = 365) {
  
  formula_gam <- as.formula(
    paste0(
      response,
      " ~ s(DOY, bs = 'ts', k = 5) + ",
      "s(Replicate, bs = 're') + Year + Aspect"
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
  
  tibble::tibble(
    variable = response,
    
    ## model fit
    deviance_explained = sm$dev.expl,
    
    dispersion = disp,

    ## season effect
    season_p_value      = season_row[, "p-value"],
    
    ## year effect
    year_estimate = year_row[, "Estimate"],
    year_se       = year_row[, "Std. Error"],
    year_p_value  = year_row[, "Pr(>|t|)"]
  )
}


gam_results <- map_dfr(
  variables_for_gams,
  fit_gam_extract_timeseries,
  data = metadata_timeseries
)

# (p.adjust default = Holm, as for the condition models above)
gam_results$season_effect_padj = p.adjust(gam_results$season_p_value)
gam_results$year_effect_padj = p.adjust(gam_results$year_p_value)

gam_results %>% filter((year_effect_padj < 0.05) | (season_effect_padj < 0.05))
gam_results
save_stat_table(gam_results, '01_env_gams_timeseries')





################################################################################
# Part 3. Figure 1 - soil temperature trajectories + gas fluxes
################################################################################
# Panel B: GAM-predicted soil temperature, control vs warming. Predictions are
# made for every replicate / aspect combination and summarised by their median.
mod_temp_time <- mgcv::gam(
  Temp_soil ~ s(Days_Since_Start, bs = 'ts', k = 5) + s(Replicate, bs = 're') + Condition + Aspect,
  data = metadata_condition,
  method = "REML"
)
# Expand grid over DOY, Year, and Replicate
newdata <- expand.grid(
  Days_Since_Start = 90:480,
  Condition = unique(metadata_condition$Condition),
  Replicate = unique(metadata_condition$Replicate),
  Aspect = unique(metadata_condition$Aspect)
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

b = ggplot() +
  geom_point(metadata_condition, mapping = aes(x=Days_Since_Start, y = Temp_soil, color = as.factor(Condition)), alpha=0.2) +
  geom_line(median_pred, mapping = aes(x = Days_Since_Start, y = fit_median, color = as.factor(Condition)), linewidth = 1) +
  geom_ribbon(median_pred, mapping = aes(x = Days_Since_Start, ymin = lwr_median, ymax = upr_median, 
                                         fill = as.factor(Condition)), alpha = 0.2, color = NA) +
  labs(
    x = "Day since the start of the experiment",
    y = "Soil Temperature [°C]",
    color = "Condition",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12)
  ) + scale_colour_manual(values = c('#012169', '#C8102E')) +
  scale_fill_manual(values = c('#012169', '#C8102E'))

# SFig 3: seasonal soil temperature (control plots) in the two years
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

yearly = ggplot(median_pred, aes(x = DOY, y = fit_median, color = as.factor(Year), fill = as.factor(Year))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr_median, ymax = upr_median), alpha = 0.2, color = NA) +
  labs(
    x = "Day of Year",
    y = "Soil Temperature [°C]",
    color = "Year",
    fill = "Year"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12)
  ) + scale_colour_manual(values = c('coral4', 'gold')) +
  scale_fill_manual(values = c('coral4', 'gold'))
save_fig(yearly, 'sfig_yearly_temp', width = 10.5, height = 6.5)

# ---- Panels C and D: CH4 and CO2 fluxes, control vs warming --------------------
# Load the raw (untransformed) data
env_var <- read.csv("data/SiteL23_24_env_variables_final.csv", header = TRUE)
env_var <- rename(env_var, "Sample_ID" = "Sample")

custom_palette <- c("#0072B2", "#E69F00", "#009E73", "#CC3333", "#4EAEE8",
                    "#9A7DCC", "#D55E00", "#009CAD", "#FF7F6E", "#7E9A3C", "#999999")

# Subset data
gas_fluxes <- env_var[, c("Sample_ID", "CO2_Flux", "CH4_Flux", "Temp_soil", "Replicate", "Condition", "Season", "Year", "Date", "DOY", "Days_Since_Start")]

# Remove NAs
gas_fluxes <- gas_fluxes[complete.cases(gas_fluxes), ]

# Create new column Season_Year
gas_fluxes$Season_Year <- paste(gas_fluxes$Season, gas_fluxes$Year, sep = " ")

# Flag outliers (rstatix definition, within Condition x Season_Year). The flags are
# informational; only the two samples named below are removed from the analysis.
gas_fluxes <- gas_fluxes %>%
  group_by(Condition, Season_Year) %>%
  mutate(
    is_co2_outlier = CO2_Flux %in% identify_outliers(pick(everything()), variable = "CO2_Flux")$CO2_Flux,
    is_ch4_outlier = CH4_Flux %in% identify_outliers(pick(everything()), variable = "CH4_Flux")$CH4_Flux
  ) %>%
  ungroup()

# Manually identify and remove extreme outlier
outliers_filtered_CO2 <- gas_fluxes %>% filter(!Sample_ID == "B_t_25.07.24")
outliers_filtered_CH4 <- gas_fluxes %>% filter(!Sample_ID == "D_t_22.08.24")

# Fix order of Season_Year
levels_order <- c("Autumn 2023", "Summer 2024", "Autumn 2024")
outliers_filtered_CO2$Season_Year <- factor(outliers_filtered_CO2$Season_Year, levels = levels_order)
outliers_filtered_CH4$Season_Year <- factor(outliers_filtered_CH4$Season_Year, levels = levels_order)

##############################
### Connected Scatterplots ###
##############################

# Remove date missing warming Condition
# CO2
plot_data_CO2 <- outliers_filtered_CO2 %>%
  group_by(Date) %>% 
  filter(all(c("control", "warming") %in% Condition)) %>% ungroup()

# CH4
plot_data_CH4 <- outliers_filtered_CH4 %>%
  group_by(Date) %>%
  filter(all(c("control", "warming") %in% Condition)) %>% ungroup()

####################
### Effect sizes ###
####################

# Van Elteren test (Wilcoxon rank-sum stratified by sampling date), asymptotic.
# coin:: is explicit because rstatix also defines wilcox_test().
vE_co2 <- coin::wilcox_test(CO2_Flux ~ factor(Condition) | factor(Date),
                            data = plot_data_CO2,
                            distribution = "asymptotic")
vE_ch4 <- coin::wilcox_test(CH4_Flux ~ factor(Condition) | factor(Date),
                            data = plot_data_CH4,
                            distribution = "asymptotic")
save_stat_text("01_flux_vanElteren_co2", vE_co2)
save_stat_text("01_flux_vanElteren_ch4", vE_ch4)
save_stat_table(
  tibble(gas = c("CO2", "CH4"),
         n_samples = c(nrow(plot_data_CO2), nrow(plot_data_CH4)),
         n_dates = c(n_distinct(plot_data_CO2$Date), n_distinct(plot_data_CH4$Date)),
         Z = c(coin::statistic(vE_co2), coin::statistic(vE_ch4)),
         p_value = c(coin::pvalue(vE_co2), coin::pvalue(vE_ch4))),
  "01_flux_vanElteren_summary")


#######################
### Plotting fluxes ###
#######################

c <- plot_data_CH4 %>%
  ggplot(aes(x = Date, y = CH4_Flux, group = interaction(Date, Condition))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", show.legend = FALSE,
               aes(group = Condition, fill = Condition), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", aes(group = Condition, color = Condition),
               linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4,
               aes(group = Condition), linewidth = 1, position = position_dodge(width = 0.7)) +
  geom_boxplot(aes(fill = Condition), position = position_dodge2(padding = 0.2),
               coef = 0, outlier.shape = NA) +
  scale_y_continuous(breaks = extended_breaks(n = 10), expand = expansion(mult = c(0.1, 0.2))) +
  scale_colour_manual(values = custom_palette[c(1, 4)], aesthetics = c("colour", "fill")) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = -45, size = 10, hjust = 0)) +
  labs(x = NULL, y = expression(paste("CH"[4], " nmol · m"^-2, "· s"^-1))) + 
  # van Elteren test results
  annotate("text", x = 1, y = min(plot_data_CH4$CH4_Flux) * 1.15,
           label = paste0("van Elteren test: p = ", format.pval(coin::pvalue(vE_ch4), digits = 3)),
           hjust = 0, size = 3.5, fontface = "italic")



# Plot CO2
d <- plot_data_CO2 %>% 
  ggplot(aes(x = Date, y = CO2_Flux, group = interaction(Date, Condition))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", show.legend = FALSE,
               aes(group = Condition, fill = Condition), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", aes(group = Condition, color = Condition),
               linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4,
               aes(group = Condition), linewidth = 1, position = position_dodge(width = 0.7))  +
  geom_boxplot(aes(fill = Condition), position = position_dodge2(padding = 0.2),
               coef = 0, outlier.shape = NA) +
  scale_y_continuous(breaks = extended_breaks(n = 10), expand = expansion(mult=c(0.1,0.2))) +  
  scale_colour_manual(values = custom_palette[c(1, 4)], aesthetics = c("colour", "fill")) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = -45, size = 10,hjust = 0)) +
  labs(x = NULL, y = expression(paste("CO"[2], " μmol · m"^-2, " · s"^-1))) + 
  # van Elteren test results
  annotate("text", x = 1, y = min(plot_data_CO2$CO2_Flux) * 1.15,
           label = paste0("van Elteren test: p = ", format.pval(coin::pvalue(vE_co2), digits = 3)),
           hjust = 0, size = 3.5, fontface = "italic")






# ---- Assemble Figure 1: A = PCA, B = soil temperature, C = CH4, D = CO2 ---------
fig = ggarrange(pca_cond_plot, b, c, d, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'), common.legend = T)
fig
save_fig(fig, 'fig_pca_temp_fluxes', width = 11.5, height = 9.2)

save_session_info("01")

