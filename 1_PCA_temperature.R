setwd('~/Documents/MACE/EVO_study/warming_plaine_morte')
source('0_functions_and_packages.R')

dir.create('figures')
dir.create('stats')

metadata = load_metadata()
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

metadata$Aspect = 'North'
metadata$Aspect[metadata$Replicate %in% c('C','D','E')] = 'South'
metadata_condition$Aspect = 'North'
metadata_condition$Aspect[metadata_condition$Replicate %in% c('C','D','E')] = 'South'
metadata_timeseries$Aspect = 'North'
metadata_timeseries$Aspect[metadata_timeseries$Replicate %in% c('C','D','E')] = 'South'
metadata_condition$Aspect = as.factor(metadata_condition$Aspect)
metadata_timeseries$Aspect = as.factor(metadata_timeseries$Aspect)

metadata$OTC = 'No'
metadata$OTC[metadata$Replicate %in% c('A','B','C')] = 'Yes'
metadata_condition$OTC = 'No'
metadata_condition$OTC[metadata_condition$Replicate %in% c('A','B','C')] = 'Yes'
metadata_timeseries$OTC = 'No'
metadata_timeseries$OTC[metadata_timeseries$Replicate %in% c('A','B','C')] = 'Yes'
metadata_condition$OTC = as.factor(metadata_condition$OTC)
metadata_timeseries$OTC = as.factor(metadata_timeseries$OTC)

metadata_condition$Replicate = as.factor(metadata_condition$Replicate)
metadata_timeseries$Replicate = as.factor(metadata_timeseries$Replicate)

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
pca_plot_df <- cbind(as.data.frame(pca_cond$x[, 1:2]), Sample = meta_cond$Sample, Condition = meta_cond$Condition, Aspect=metadata_condition$Aspect)

pca_cond_test <- PCAtest(pca_mat_cond, 100, 100, 
                                          0.05, varcorr = TRUE, counter = FALSE, plot = TRUE)

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
    x = paste0("PC1 (22.3%)"),
    y = paste0("PC2 (17.9%)"),
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

pca_time = prcomp(pca_mat_time)
pca_plot_df <- cbind(as.data.frame(pca_time$x[, 1:2]), Sample = meta_time$Sample, Days_Since_Start = meta_time$Days_Since_Start)

pca_time_test <- PCAtest(pca_mat_time, 100, 100, 
                         0.05, varcorr = TRUE, counter = FALSE, plot = TRUE)

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
    x = paste0("PC1 (22.3% )"),
    y = paste0("PC2 (20%)"),
    color = "Days since the start of the experiment"
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_gradient(low = 'gold', high = 'violet') + xlim(-7.5, 6.5) + ylim(-6.5, 6.5)
pca_time_plot
ggsave(filename = 'figures/SuppFigure1_pca_time_plot.pdf', pca_time_plot, width = 10.5)


# 3. PCA for aspect
meta_aspect <- metadata %>%
  select(Sample, Days_Since_Start, all_of(variables_for_pca), Aspect) %>% select(-CH4_Flux_scaled, -CO2_Flux_scaled)

pca_mat_aspect <- meta_aspect %>%
  select(-Sample, -Days_Since_Start, -Aspect) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()
rownames(pca_mat_aspect) <- meta_aspect$Sample

pca_aspect = prcomp(pca_mat_aspect)
pca_plot_df <- cbind(as.data.frame(pca_aspect$x[, 1:2]), Sample = meta_aspect$Sample, 
                     Days_Since_Start = meta_aspect$Days_Since_Start, Aspect = meta_aspect$Aspect)

pca_time_test <- PCAtest(pca_mat_time, 100, 100, 
                         0.05, varcorr = TRUE, counter = FALSE, plot = TRUE)

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
    x = paste0("PC1 (22.3%)"),
    y = paste0("PC2 (20%)"),
    color = "Aspect"
  ) +
  theme_light() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  scale_colour_manual(values = c('turquoise','violet'))
pca_aspect_plot
ggsave(filename = 'figures/SuppFigure2_pca_aspect_plot.pdf', pca_aspect_plot, width = 10.5)




library(Hmisc)
# GAMs for env parameters, testing Condition and time
variables_for_gams = c("Temp_soil","TC","TOC","pH","Water_Content_Lab")

# A) Models Condition
metadata_condition <- metadata_condition %>%
  mutate(
    Condition = factor(Condition),
    Replicate = factor(Replicate),
    Aspect = factor(Aspect)
  )

fit_gam_extract <- function(response, data) {
  
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
  fit_gam_extract,
  data = metadata_condition
)

gam_results$condition_padj = p.adjust(gam_results$condition_p_value, method='holm')
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
  geom_line(median_pred, mapping = aes(x = Days_Since_Start, y = fit_median, color = as.factor(Condition)), size = 1) +
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
  geom_line(size = 1) +
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
ggsave(filename = 'figures/SuppFigure3_yearly_temp.pdf', yearly, width = 10.5, height = 6.5)

# panel C and D with gas fluxes
library(tidyverse)
library(rstatix)
library(paletteer)
library(ggpubr)
library(ggsignif)
library(cowplot)
library(scales)
library(rstatix)
library(coin)

# Load the data
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

# Identify outliers
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

# Van Elteren test stratified by Date
vE_co2 <- wilcox_test(CO2_Flux ~ factor(Condition) | factor(Date), 
                      data = plot_data_CO2,
                      distribution = "asymptotic")
summary(vE_co2)
pvalue(vE_co2)

vE_ch4 <- wilcox_test(CH4_Flux ~ factor(Condition) | factor(Date), 
                      data = plot_data_CH4,
                      distribution = "asymptotic")
summary(vE_ch4)
pvalue(vE_ch4)


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
  labs(x = NULL, y = expression(paste("CH"[4], " ppb/s"))) + 
  # van Elteren test results
  annotate("text", x = 1, y = min(plot_data_CH4$CH4_Flux) * 1.15,
           label = paste0("van Elteren test: p = ", format.pval(pvalue(vE_ch4), digits = 3)),
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
  labs(x = NULL, y = expression(paste("CO"[2], " ppm/s"))) + 
  # van Elteren test results
  annotate("text", x = 1, y = min(plot_data_CO2$CO2_Flux) * 1.15,
           label = paste0("van Elteren test: p = ", format.pval(pvalue(vE_co2), digits = 3)),
           hjust = 0, size = 3.5, fontface = "italic")






#bcd = ggarrange(b + theme(legend.position = 'none'), 
                #c + theme(legend.position = 'none'), 
                #d + theme(legend.position = 'none'), nrow = 3, ncol = 1, labels = c('B', 'C', 'D'))
#fig = ggarrange(pca_cond_plot + theme(legend.position = 'bottom'), bcd, ncol = 2, nrow = 1, labels = c('A',''))
fig = ggarrange(pca_cond_plot, b, c, d, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'), common.legend = T)
fig
ggsave(filename = 'figures/Figure1.pdf', fig, width = 11.5, height = 9.2)

