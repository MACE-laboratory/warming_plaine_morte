library(jsonlite)
library(ggplot2)
library(patchwork)

setwd('~/Documents/MACE/warming_plaine_morte')
source('0_functions_and_packages.R')
## ---- 1. Parameters -------------------------------------------------------------
category_level          <- "sub"   # "sub" (Recommended, ~40 BRITE functional sub-categories)
# or "top" (~7 broad categories: Metabolism, ...)
prevalence_threshold_meta <- 0.5   # Panel A + Panel B background: KO must be present
# (abundance > 0) in MORE THAN this fraction of
# meta_KO samples to count as "detected"
prot_detection_threshold  <- 0     # Panel B test set: min. fraction of prot_KO samples
# a KO must be present in to count as "detected" in
# the metaproteome. Default 0 = present in >=1 sample
# (proteomics is typically far sparser than
# metagenomics); set to 0.5 to use the same ">half
# the samples" rule as the metaG side for consistency.
padj_threshold           <- 0.05   # significance cutoff used for panel B's fading
panelA_top_n              <- 20    # panel A shows only the N most abundant categories
out_prefix                <- "figures/KO_category_panels"
## ---- 2. Load your data ----------------------------------------------------------
metadata = load_metadata()
metadata = check_normality_change(metadata, c("Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
                                              "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",
                                              "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",
                                              "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
                                              "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
                                              "Other_Cations","Grouped_OC"))
loaded_meta = load_meta_data()
meta_KO <- loaded_meta$ko_raw_table
prot_KO <- read.delim("data/proteomics/metaP_norm_imp.tsv", sep = "\t", header = TRUE)
rownames(prot_KO) = prot_KO$name
prot_KO$name = NULL
meta_KO <- as.matrix(meta_KO); storage.mode(meta_KO) <- "numeric"
prot_KO <- as.matrix(prot_KO); storage.mode(prot_KO) <- "numeric"
rownames(meta_KO) <- sub("^ko:", "", rownames(meta_KO))
rownames(prot_KO) <- sub("^ko:", "", rownames(prot_KO))
## ---- 3. KO -> KEGG BRITE category mapping (ko00001), via the KEGG REST API -----
## The ko00001 hierarchy is 4 levels deep: A = top category, B = sub-category,
## C = pathway, D = individual KO. We fetch it as JSON (much easier to parse
## reliably than the flat-text "htext" format) and walk it recursively, so a
## KO that sits under more than one branch produces one row per branch.
message("Fetching the KEGG BRITE ko00001 functional hierarchy ...")
brite <- tryCatch(
  jsonlite::fromJSON("https://rest.kegg.jp/get/br:ko00001/json", simplifyVector = FALSE),
  error = function(e) {
    stop(
      "Could not reach the KEGG REST API (https://rest.kegg.jp/get/br:ko00001/json).\n",
      "Check your internet connection / firewall. Original error: ", conditionMessage(e),
      call. = FALSE
    )
  }
)
extract_ko_categories <- function(node, ancestors = character(0)) {
  nm <- node$name
  if (grepl("^K[0-9]{5}\\b", nm)) {
    ko <- sub("^([K][0-9]{5}).*$", "\\1", nm)
    cat_top <- if (length(ancestors) >= 1) ancestors[1] else NA_character_
    cat_sub <- if (length(ancestors) >= 2) ancestors[2] else NA_character_
    return(list(data.frame(ko = ko, category_top = cat_top, category_sub = cat_sub,
                           stringsAsFactors = FALSE)))
  }
  kids <- node$children
  if (is.null(kids) || length(kids) == 0) return(list())
  new_ancestors <- c(ancestors, nm)
  do.call(c, lapply(kids, extract_ko_categories, ancestors = new_ancestors))
}
## Start from brite$children, not brite itself -- the root node's own name
## ("ko00001") must NOT become part of the ancestor chain, or every KO's
## category/sub-category would be picked up one level too high.
ko_cat_rows <- do.call(c, lapply(brite$children, extract_ko_categories, ancestors = character(0)))
ko_cat_map  <- unique(do.call(rbind, ko_cat_rows))
clean_label <- function(x) trimws(sub("^[0-9]+\\s+", "", x))
ko_cat_map$category_top_label <- clean_label(ko_cat_map$category_top)
ko_cat_map$category_sub_label <- clean_label(ko_cat_map$category_sub)
cat_col <- if (category_level == "top") "category_top_label" else "category_sub_label"
ko_cat_map <- unique(ko_cat_map[!is.na(ko_cat_map[[cat_col]]), c("ko", cat_col)])
colnames(ko_cat_map) <- c("ko", "category")

## Drop categories tied to human/animal physiology & disease -- KEGG's
## "Organismal Systems" and "Human Diseases" top-level branches, plus
## antimicrobial/antineoplastic drug resistance -- since they aren't
## meaningful for an environmental soil microbiome. A KO that ALSO maps to a
## non-human category (most do) keeps that mapping; only the human-specific
## rows are removed, not the KO itself.
human_categories <- c(
  # Organismal Systems
  "Immune system", "Endocrine system", "Circulatory system", "Digestive system",
  "Excretory system", "Nervous system", "Sensory system",
  "Development and regeneration", "Aging", "Environmental adaptation",
  # Human Diseases
  "Cancer: overview", "Cancer: specific types", "Infectious disease: bacterial",
  "Infectious disease: viral", "Infectious disease: parasitic",
  "Neurodegenerative disease", "Substance dependence", "Cardiovascular disease",
  "Endocrine and metabolic disease", "Immune disease",
  # Drug Development
  "Drug resistance: antimicrobial", "Drug resistance: antineoplastic"
)
n_cat_before <- length(unique(ko_cat_map$category))
ko_cat_map <- ko_cat_map[!ko_cat_map$category %in% human_categories, ]
message("Dropped ", n_cat_before - length(unique(ko_cat_map$category)),
        " human-related categories (of ", n_cat_before, ").")

message("Parsed ", length(unique(ko_cat_map$ko)), " KOs across ",
        length(unique(ko_cat_map$category)), " '", category_level, "'-level categories.")
## ---- 4. Panel A -- KOs per category, prevalent (>half samples) in the metagenome ----
detect_kos <- function(mat, min_prevalence) {
  prevalence <- rowSums(mat > 0, na.rm = TRUE) / ncol(mat)
  rownames(mat)[prevalence > min_prevalence]
}
meta_detected <- detect_kos(meta_KO, prevalence_threshold_meta)
message(length(meta_detected), " / ", nrow(meta_KO),
        " KOs detected in >", prevalence_threshold_meta * 100, "% of metagenome samples.")
meta_map <- unique(ko_cat_map[ko_cat_map$ko %in% meta_detected, ])
n_unmapped_A <- length(setdiff(meta_detected, meta_map$ko))
if (n_unmapped_A) {
  message(n_unmapped_A, " prevalent metaG KOs had no (non-human) KEGG category and were excluded from panel A.")
}
panelA_data <- as.data.frame(table(meta_map$category), stringsAsFactors = FALSE)
colnames(panelA_data) <- c("category", "n_KO")
panelA_data <- panelA_data[order(-panelA_data$n_KO), ]
panelA_data <- panelA_data[seq_len(min(panelA_top_n, nrow(panelA_data))), ]  # top N only
panelA_data$category <- factor(panelA_data$category, levels = rev(panelA_data$category))
BLUE  <- "#2a78d6"
RED   <- "#e34948"
GREY  <- "#c3c2b7"
panelA <- ggplot(panelA_data, aes(x = n_KO, y = category)) +
  geom_col(fill = BLUE, width = 0.75) +
  geom_text(aes(label = n_KO), hjust = -0.2, size = 3, color = "#52514e") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of KOs", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
## ---- 5. Panel B -- category enrichment: metaP (test) vs. metaG (background) --------
prot_detected <- detect_kos(prot_KO, prot_detection_threshold)
message(length(prot_detected), " / ", nrow(prot_KO),
        " KOs detected in metaproteome (>", prot_detection_threshold * 100, "% of samples).")
test_map <- unique(ko_cat_map[ko_cat_map$ko %in% prot_detected, ])
bg_map   <- meta_map   # same prevalence-filtered metaG set used in panel A (top-N slice
# above only affects the panelA_data summary, not bg_map/test_map)
n_test_total <- length(unique(test_map$ko))
n_bg_total   <- length(unique(bg_map$ko))
categories <- union(unique(test_map$category), unique(bg_map$category))
enrich_one <- function(cat) {
  a <- length(unique(test_map$ko[test_map$category == cat]))          # test, in category
  b <- n_test_total - a                                                 # test, not in category
  c <- length(unique(bg_map$ko[bg_map$category == cat]))               # background, in category
  d <- n_bg_total - c                                                   # background, not in category
  
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft  <- fisher.test(tab)
  
  ## Haldane-Anscombe correction (add 0.5) so log2FC is defined when a or c is 0
  test_frac <- (a + 0.5) / (n_test_total + 0.5)
  bg_frac   <- (c + 0.5) / (n_bg_total + 0.5)
  
  data.frame(
    category   = cat,
    n_test     = a, n_test_total = n_test_total,
    n_bg       = c, n_bg_total   = n_bg_total,
    odds_ratio = unname(ft$estimate),
    log2FC     = log2(test_frac / bg_frac),
    p_value    = ft$p.value,
    stringsAsFactors = FALSE
  )
}
panelB_data <- do.call(rbind, lapply(categories, enrich_one))
panelB_data$padj <- p.adjust(panelB_data$p_value, method = "BH")
panelB_data$significant <- panelB_data$padj < padj_threshold
panelB_data <- panelB_data[order(panelB_data$log2FC), ]
panelB_data$category <- factor(panelB_data$category, levels = panelB_data$category)
panelB <- ggplot(panelB_data, aes(x = log2FC, y = category, fill = log2FC, alpha = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#898781", linewidth = 0.4) +
  geom_col(width = 0.75) +
  scale_fill_gradient2(low = RED, mid = GREY, high = BLUE, midpoint = 0, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(x = expression(log[2]*" fold-enrichment (metaP / metaG)"), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
## ---- 6. Combine, draw, and export --------------------------------------------------
fig <- panelA / panelB +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
fig_height <- max(8, 0.22 * max(nrow(panelA_data), nrow(panelB_data)) + 4)
ggsave(paste0(out_prefix, ".pdf"), fig, width = 9, height = fig_height, limitsize = FALSE)
ggsave(paste0(out_prefix, ".png"), fig, width = 9, height = fig_height, dpi = 300, limitsize = FALSE)
message("Wrote ", out_prefix, ".pdf and ", out_prefix, ".png")
write.csv(panelA_data, paste0('stats/panelA_category_counts.csv'), row.names = FALSE)
write.csv(panelB_data, paste0('stats/_panelB_category_enrichment.csv'), row.names = FALSE)
message("Wrote ", out_prefix, "_panelA_category_counts.csv and ",
        out_prefix, "_panelB_category_enrichment.csv")