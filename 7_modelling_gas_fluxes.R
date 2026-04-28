setwd('~/Documents/MACE/EVO_study/warming_plaine_morte')

source('0_functions_and_packages.R')
library(performance)

loaded_mags = load_mags_data()
loaded_meta = load_meta_data()

metadata = load_metadata()
metadata = check_normality_change(metadata, c("Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
                                              "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",              
                                              "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",          
                                              "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
                                              "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
                                              "Other_Cations","Grouped_OC"))

run_gam_selection_backward <- function(metadata,
                              response_variable,
                              fixed_effects,
                              replicate_var = "Replicate") {
  
  # Subset only needed variables and remove NA
  vars_used <- c(response_variable, replicate_var, fixed_effects)
  dat <- metadata[, vars_used]
  dat <- na.omit(dat)
  
  dat[[replicate_var]] <- as.factor(dat[[replicate_var]])
  
  # Null model
  null_formula <- as.formula(
    paste0(
      response_variable,
      " ~ s(", replicate_var, ", bs='re')"
    )
  )
  
  null_model <- mgcv::gam(null_formula, data = dat, method = "REML")
  
  current_fixed <- fixed_effects
  models <- list()
  models[['null']] <- null_model
  
  results <- list()
  step <- 1
  last_removed <- NA
  
  while (TRUE) {
    
    rhs <- c(
      paste0("s(", replicate_var, ", bs='re')"),
      paste0("s(", current_fixed, ", bs='ts', k=5)")
    )
    
    formula_str <- paste(response_variable, "~", paste(rhs, collapse = " + "))
    f <- as.formula(formula_str)
    
    model_name <- paste0("model_", step)
    
    m <- mgcv::gam(f, data = dat, method = "REML")
    
    models[[model_name]] <- m
    
    results[[model_name]] <- data.frame(
      model = model_name,
      step = step,
      n_splines = length(current_fixed),
      removed = ifelse(step == 1, NA, last_removed),
      formula = formula_str,
      deviance_explained = summary(m)$dev.expl
    )
    
    if (length(current_fixed) <= 1) break
    
    sm <- summary(m)
    
    # Extract p-values from smooth terms
    smooth_table <- as.data.frame(sm$s.table)
    pvals <- smooth_table[, 'p-value']

    # Match smooth names back to variable names
    smooth_names <- rownames(smooth_table)
    var_map <- gsub("^s\\(|\\)$", "", smooth_names)
    
    pvals <- pvals[var_map %in% current_fixed]
    names(pvals) <- var_map[var_map %in% current_fixed]
    
    if (length(pvals) == 0) break
    
    worst_var <- names(which.max(pvals))
    last_removed <- worst_var
    current_fixed <- setdiff(current_fixed, worst_var)
    
    step <- step + 1
  }
  
  results_df <- dplyr::bind_rows(results)

  # Run BF comparison once
  bf <- performance::test_bf(models)
  
  bf_df <- as.data.frame(bf)
  bf_df$model <- names(models)
  print(bf_df)
  
  results_df <- results_df |>
    dplyr::left_join(bf_df[, c("model", "BF")], by = "model") |>
    dplyr::rename(BF_vs_null = BF) |>
    dplyr::arrange(desc(BF_vs_null))
  
  best <- results_df[1, ]
  
  cat("\n===== BEST MODEL =====\n")
  cat("Formula:\n", best$formula, "\n\n")
  cat("Bayes Factor vs null:", best$BF_vs_null, "\n")
  cat("Deviance explained:", best$deviance_explained, "\n")
  cat("Fixed effects:", best$n_fixed, "\n")
  
  cat("\n===== MODEL RANKING =====\n")
  print(results_df)
  
  return(list(
    best_model = best,
    results = results_df,
    models = models,
    null_model = null_model,
    bf_table = bf_df
  ))
}

run_gam_forward_selection <- function(metadata,
                                      response_variable,
                                      base_vars,
                                      try_effects,
                                      bf_threshold = 3) {
  
  vars_used <- unique(c(response_variable, try_effects, base_vars))
  dat <- metadata[, vars_used, drop = FALSE]
  dat <- na.omit(dat)
  
  # --- BASELINE MODEL (your new null) ---
  rhs_base <- c(
    paste0("s(", base_vars, ", bs='ts', k=5)")
  )
  
  base_formula_str <- paste(response_variable, "~", paste(rhs_base, collapse = " + "))
  current_model <- mgcv::gam(as.formula(base_formula_str), data = dat, method = "REML")
  
  models <- list(null_model = current_model)
  
  selected_vars <- c()
  remaining_vars <- setdiff(try_effects, base_vars)
  
  step <- 1
  improvement <- TRUE
  
  while (improvement && length(remaining_vars) > 0) {
    
    improvement <- FALSE
    best_var <- NULL
    best_model <- NULL
    best_bf <- -Inf
    
    for (var in remaining_vars) {
      
      print(var)
      
      feature_values <- dat[[var]]
      if (length(unique(feature_values[!is.na(feature_values)])) <= 1) next
      
      rhs_try <- c(
        paste0("s(", base_vars, ", bs='ts', k=5)"),
        if (length(selected_vars) > 0) paste0("s(", selected_vars, ", bs='ts', k=5)"),
        paste0("s(", var, ", bs='ts', k=5)")
      )
      
      formula_str <- paste(response_variable, "~", paste(rhs_try, collapse = " + "))
      f_try <- as.formula(formula_str)
      
      m_try <- mgcv::gam(f_try, data = dat, method = "REML")
      
      bf <- tryCatch({
        performance::test_bf(m_try, current_model)$BF[2]
      }, error = function(e) NA)
      
      print(bf)
      
      if (!is.na(bf) && bf > best_bf) {
        best_bf <- bf
        best_var <- var
        best_model <- m_try
      }
    }
    
    if (!is.null(best_var) && best_bf > bf_threshold) {
      selected_vars <- c(selected_vars, best_var)
      remaining_vars <- setdiff(remaining_vars, best_var)
      models[[paste0("step_", step, "_add_", best_var)]] <- best_model
      current_model <- best_model
      improvement <- TRUE
      step <- step + 1
    }
  }
  
  final_model <- current_model
  
  bf_final <- tryCatch({
    performance::test_bf(final_model, models[["null_model"]])$BF[2]
  }, error = function(e) NA)
  
  # If BF < 1 → keep baseline as best
  if (!is.na(bf_final) && bf_final < 1) {
    final_model <- models[["null_model"]]
    selected_vars <- c()
  }
  
  results_df <- data.frame(
    model = "final_model",
    formula = paste(deparse(formula(final_model)), collapse = ""),
    deviance_explained = summary(final_model)$dev.expl,
    BF_vs_null = bf_final,
    selected_vars = paste(selected_vars, collapse = ", "),
    stringsAsFactors = FALSE
  )
  
  cat("\n===== FINAL MODEL =====\n")
  cat("Formula:\n", results_df$formula, "\n\n")
  cat("Bayes Factor vs null:", results_df$BF_vs_null, "\n")
  cat("Deviance explained:", results_df$deviance_explained, "\n")
  cat("Selected splines:", if(nchar(results_df$selected_vars) > 0) results_df$selected_vars else "none", "\n")
  
  return(list(
    final_model = final_model,
    results = results_df,
    null_model = models[["null_model"]],
    selected_vars = selected_vars
  ))
}

abiotic_ch4_soilt = run_gam_forward_selection(metadata, response_variable = 'CH4_Flux',
                                        base_vars = c("Temp_soil_scaled", "Water_Content_Lab_scaled"),
                                        try_effects = c("TC_scaled","TOC_scaled","TIN_scaled",
                                                      "pH_scaled","NH4._scaled"))

abiotic_ch4_airt = run_gam_forward_selection(metadata, response_variable = 'CH4_Flux',
                                        base_vars = c("Temp_air_scaled", "Water_Content_Lab_scaled"),
                                        try_effects = c("TC_scaled","TOC_scaled","TIN_scaled",
                                                        "pH_scaled","NH4._scaled"))

abiotic_co2_soilt = run_gam_forward_selection(metadata, response_variable = 'CO2_Flux',
                                        base_vars = c("Temp_soil_scaled", "Water_Content_Lab_scaled"),
                                        try_effects = c("TC_scaled","TOC_scaled","TIN_scaled",
                                                        "pH_scaled","NH4._scaled"))

abiotic_co2_airt = run_gam_forward_selection(metadata, response_variable = 'CO2_Flux',
                                        base_vars = c("Temp_air_scaled", "Water_Content_Lab_scaled"),
                                        try_effects = c("TC_scaled","TOC_scaled","TIN_scaled",
                                                        "pH_scaled","NH4._scaled"))

summary(abiotic_ch4_airt$final_model)
summary(abiotic_co2_soilt$final_model)

mags_table = as.data.frame(loaded_mags$clr_table)
meta_table = as.data.frame(loaded_meta$ko_clr_table)
prot_table = read.table('data/proteomics/metaP_norm_imp.tsv', header = T)
prot_table = prot_table %>% column_to_rownames('name')
rownames(prot_table) = map_chr(rownames(prot_table), function(x) gsub('ko:','',x))
colnames(prot_table) = map_chr(colnames(prot_table), function(x) gsub('Control_','',x))
colnames(prot_table) = map_chr(colnames(prot_table), function(x) gsub('Warming_','',x))

correlate_features_with_gam <- function(feature_table, 
                                        metadata, 
                                        sample_col = "Sample",
                                        best_model_result,
                                        n_perm = 40) {
  
  best_abiotic <- best_model_result$final_model
  base_formula <- formula(best_abiotic)
  
  meta <- metadata
  
  samples_common <- intersect(colnames(feature_table), meta[[sample_col]])
  feature_table <- feature_table[, samples_common, drop = FALSE]
  meta <- meta[meta[[sample_col]] %in% samples_common, , drop = FALSE]
  
  # Order metadata to match feature table columns
  meta <- meta[match(colnames(feature_table), meta[[sample_col]]), , drop = FALSE]
  
  output_list <- vector("list", nrow(feature_table))
  
  for (i in seq_len(nrow(feature_table))) {
    cat("\rProgress:", i, "/", nrow(feature_table))
    flush.console()
    
    feature_name <- rownames(feature_table)[i]
    
    meta_tmp <- meta
    meta_tmp[[feature_name]] <- as.numeric(feature_table[i, ])
    
    meta_tmp <- meta_tmp[complete.cases(meta_tmp), ]
    
    base_refit <- mgcv::gam(base_formula, data = meta_tmp, method = "REML")
    
    new_formula <- update(
      base_formula,
      paste0(". ~ . + s(", feature_name, ", bs='ts', k=5)")
    )
    
    new_model <- mgcv::gam(new_formula, data = meta_tmp, method = "REML")
    
    bf <- performance::test_bf(base_refit, new_model)$BF[2]
    
    dev_base <- summary(base_refit)$dev.expl
    dev_new <- summary(new_model)$dev.expl
    
    # Extract spline edf and p-value
    smooth_table <- summary(new_model)$s.table
    row_id <- grep(paste0("^s\\(", feature_name, "\\)"), rownames(smooth_table))
    
    if (length(row_id) > 0) {
      edf_feature <- smooth_table[row_id, "edf"]
      p_val <- smooth_table[row_id, "p-value"]
    } else {
      edf_feature <- NA
      p_val <- NA
    }
    
    output_list[[i]] <- data.frame(
      feature = feature_name,
      edf_feature = edf_feature,
      pval_feature = p_val,
      BF_vs_base = bf,
      deviance_base = dev_base,
      deviance_with = dev_new,
      delta_deviance = dev_new - dev_base,
      logLik_base = as.numeric(logLik(base_refit)),
      logLik_new = as.numeric(logLik(new_model)),
      delta_logLik = as.numeric(logLik(new_model) - logLik(base_refit)),
      AIC_base = AIC(base_refit),
      AIC_new = AIC(new_model),
      delta_AIC = AIC(new_model) - AIC(base_refit),
      stringsAsFactors = FALSE
    )
  }
  
  cat("\n")
  output_df <- dplyr::bind_rows(output_list)
  output_df$padj_feature <- p.adjust(output_df$pval_feature, method = "fdr")
  
  return(output_df)
}

############################################
# CH4
# For MAGs table
mags_stats_ch4 <- correlate_features_with_gam(
  feature_table = mags_table,
  metadata = metadata,
  sample_col = "Sample",
  best_model_result = abiotic_ch4_airt
)

# For Metagenome table
meta_stats_ch4 <- correlate_features_with_gam(
  feature_table = meta_table,
  metadata = metadata,
  sample_col = "Sample",
  best_model_result = abiotic_ch4_airt
)

# For Proteins table
prot_stats_ch4 <- correlate_features_with_gam(
  feature_table = prot_table,
  metadata = metadata,
  sample_col = "Sample",
  best_model_result = abiotic_ch4_airt
)

mags_stats_ch4 %>% filter(BF_vs_base > 3, padj_feature < 0.05)
meta_stats_ch4 %>% filter(BF_vs_base > 3, padj_feature < 0.05)
prot_stats_ch4 %>% filter(BF_vs_base > 3, padj_feature < 0.05)

############################################
# C02

# For MAGs table
mags_stats_co2 <- correlate_features_with_gam(
  feature_table = mags_table,
  metadata = metadata,
  sample_col = "Sample",
  best_model_result = abiotic_co2_soilt
)

# For Metagenome table
meta_stats_co2 <- correlate_features_with_gam(
  feature_table = meta_table,
  metadata = metadata,
  sample_col = "Sample",
  best_model_result = abiotic_co2_soilt
)

# For Proteins table
prot_stats_co2 <- correlate_features_with_gam(
  feature_table = prot_table,
  metadata = metadata,
  sample_col = "Sample",
  best_model_result = abiotic_co2_soilt
)

mags_stats_co2 %>% filter(BF_vs_base > 3, padj_feature < 0.05)
meta_stats_co2 %>% filter(BF_vs_base > 3, padj_feature < 0.05)
prot_stats_co2 %>% filter(BF_vs_base > 3, padj_feature < 0.05)

mags_stats_ch4$Dataset = 'MAG'
meta_stats_ch4$Dataset = 'Meta KO'
prot_stats_ch4$Dataset = 'Prot KO'

mags_stats_co2$Dataset = 'MAG'
meta_stats_co2$Dataset = 'Meta KO'
prot_stats_co2$Dataset = 'Prot KO'

ch4_res = rbind(mags_stats_ch4, meta_stats_ch4, prot_stats_ch4)
co2_res = rbind(mags_stats_co2, meta_stats_co2, prot_stats_co2)

write.table(ch4_res, file = 'stats/modelling_ch4.csv', sep=',', quote = F, row.names = F)
write.table(co2_res, file = 'stats/modelling_co2.csv', sep=',', quote = F, row.names = F)





plot_gam_partial_effects <- function(gam_model,
                                     metadata,
                                     sample_col = "Sample",
                                     replicate_var = "Replicate",
                                     colour = "#2C7FB8",
                                     n_points = 200) {
  
  # Data used in the model fit
  data <- gam_model$model
  
  # Pretty variable names
  var_labels <- c(
    Temp_soil_scaled = "Scaled soil temperature",
    Water_Content_Lab_scaled = "Scaled soil water content",
    pH_scaled = "Scaled pH"
  )
  
  sm <- summary(gam_model)$s.table
  smooth_terms <- rownames(sm)
  
  # remove random effect
  smooth_terms <- smooth_terms[!grepl(replicate_var, smooth_terms)]
  
  vars <- gsub("^s\\(|,.*|\\)$", "", smooth_terms)
  
  res <- residuals(gam_model, type = "working")
  term_pred <- predict(gam_model, type = "terms", se.fit = TRUE)
  
  plots <- list()
  
  for (v in vars) {
    
    term_name <- grep(paste0("s\\(", v), colnames(term_pred$fit), value = TRUE)[1]
    
    smooth_contrib <- term_pred$fit[, term_name]
    partial_resid <- smooth_contrib + res
    
    df_points <- data.frame(
      x = data[[v]],
      partial_resid = partial_resid,
      replicate = data[[replicate_var]]
    )
    
    rng <- range(data[[v]], na.rm = TRUE)
    grid <- seq(rng[1], rng[2], length.out = n_points)
    
    newdata <- data[rep(1, n_points), ]
    newdata[[v]] <- grid
    
    pred_grid <- predict(gam_model, newdata = newdata, type = "terms", se.fit = TRUE)
    
    df_smooth <- data.frame(
      x = grid,
      fit = pred_grid$fit[, term_name],
      se = pred_grid$se.fit[, term_name]
    )
    
    df_smooth$upper <- df_smooth$fit + 2 * df_smooth$se
    df_smooth$lower <- df_smooth$fit - 2 * df_smooth$se
    
    # Apply pretty name if available
    x_label <- if (v %in% names(var_labels)) var_labels[[v]] else v
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = df_points,
        ggplot2::aes(x = x, y = partial_resid),
        colour = colour,
        alpha = 0.35,
        size = 1.5
      ) +
      ggplot2::geom_ribbon(
        data = df_smooth,
        ggplot2::aes(x = x, ymin = lower, ymax = upper),
        fill = colour,
        alpha = 0.2
      ) +
      ggplot2::geom_line(
        data = df_smooth,
        ggplot2::aes(x = x, y = fit),
        colour = colour,
        linewidth = 1.2
      ) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        x = x_label,
        y = "Partial effect"
      )
    
    plots[[v]] <- p
  }
  
  plots
}

subplots_ch4 = plot_gam_partial_effects(abiotic_ch4$final_model, metadata)
subplots_co2 = plot_gam_partial_effects(abiotic_co2$final_model, metadata, colour = 'gold')

ch4p = ggarrange(subplots_ch4$Temp_soil_scaled + labs(title = 'Methane'), subplots_ch4$Water_Content_Lab_scaled, subplots_ch4$pH_scaled,nrow=1, align = 'h')
co2p = ggarrange(subplots_co2$Temp_soil_scaled + labs(title = 'Carbon dioxide') + ylim(-0.95,0.5), subplots_co2$Water_Content_Lab_scaled + ylim(-0.95,0.5), ggplot(), nrow=1, align = 'h')
ggarrange(ch4p, 
          co2p, 
          nrow = 2, align = 'v', label.x = 0)
ggsave('figures/Figure6.pdf', width = 8, height = 5)

###############

# For the carbon dioxide-associated MAG
loaded_mags$mag_info[loaded_mags$mag_info$MAG == 'sub_assembly_3_bin_752717',] # MAG information, taxonomy

quantile(as.numeric(loaded_mags$relab_table['sub_assembly_3_bin_752717',])) # MAG relative abundance

data %>% filter(grepl('sub_assembly_3_bin_752717', Protein)) %>% pull(Gene.names) # Protein groups from that MAG in the metaproteomes

annot %>% filter(startsWith(X.query, 'sub_assembly_3_bin_752717')) %>% filter(grepl('ko00195', KEGG_Pathway)) # Photosynthesis pathway? no genes except ATP synthases

annot %>% filter(startsWith(X.query, 'sub_assembly_3_bin_752717')) %>% filter(grepl('map01120', KEGG_Pathway)) # Microbial metabolism in diverse environments

annot %>% filter(startsWith(X.query, 'sub_assembly_3_bin_752717')) %>% filter(grepl('map00720', KEGG_Pathway)) # Other carbon fixation pathways

mag_annot = annot %>% filter(startsWith(X.query, 'sub_assembly_3_bin_752717')) 
write.table(mag_annot, file = 'stats/sub_assembly_3_bin_752717_annotation.csv', sep=',', quote = F)


# Check which protein groups have the methane-associated KO assigned to them
data = read.csv('data/proteomics/combined_protein.tsv', sep='\t') 
data$Potential.contaminant = '-'
data$Potential.contaminant[startsWith(data$Protein.ID, 'CON_')] = '+'
data <- filter(data, Reverse != "+", Potential.contaminant != "+")
annot = read.csv('data/proteomics/annotations_mags_EDU_2025.emapper.annotations', skip = 4, header = T, sep='\t')
data$Gene.names = map_chr(data$Entry.Name, function(x) ifelse(length(annot$KEGG_ko[annot$X.query %in% x]) > 0, paste0(annot$KEGG_ko[annot$X.query %in% x], collapse = ';'), ''))
genes = data %>% filter(grepl("K03286", Gene.names)) %>% pull(Protein)
bins_names <- sub("_k127.*", "", genes)
methane_ko_mags = loaded_mags$mag_info[loaded_mags$mag_info$MAG %in% bins_names,]
