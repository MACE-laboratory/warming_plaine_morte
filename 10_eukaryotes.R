library(data.table)
library(taxonomizr)

setwd('~/Documents/MACE/warming_plaine_morte')
source('0_functions_and_packages.R')

## ---------------------------------------------------------------------
## 0. USER-EDITABLE PATHS
## ---------------------------------------------------------------------
tax_files <- c(
  sub_assembly_0 = "data/contigs_taxo/sub_assembly_0_mmseqs2_tax_lca.tsv",
  sub_assembly_1 = "data/contigs_taxo/sub_assembly_1_mmseqs2_tax_lca.tsv",
  sub_assembly_2 = "data/contigs_taxo/sub_assembly_2_mmseqs2_tax_lca.tsv",
  sub_assembly_3 = "data/contigs_taxo/sub_assembly_3_mmseqs2_tax_lca.tsv"
)
coverage_file <- "data/functional/all_contigs_trimmed_mean.tsv"

taxa_sql <- "data/accessionTaxa.sql"   # local lineage lookup db, built once (see below)
prepareDatabase(taxa_sql)

tax_cols <- c("contig", "taxid", "rank", "name",
              "frags_retained", "frags_assigned", "frags_agree", "support")

## ---------------------------------------------------------------------
## 1. Read + combine the 4 mmseqs2 LCA tables, rebuild the
##    "sub_assembly_N_<contig>" id used in the coverage table
## ---------------------------------------------------------------------
tax_list <- lapply(names(tax_files), function(sa) {
  dt <- fread(tax_files[[sa]], header = FALSE, col.names = tax_cols, sep = "\t")
  dt[, contig_full := paste0(sa, "_", contig)]
  dt
})
tax_dt <- rbindlist(tax_list)
cat("Total contigs with a taxonomic call:", nrow(tax_dt), "\n")

## ---------------------------------------------------------------------
## 2. Resolve full lineage for every unique taxid
## ---------------------------------------------------------------------
uid <- unique(tax_dt$taxid)
lineage <- getTaxonomy(uid, taxa_sql)              # matrix: superkingdom, phylum, class, ...
lineage_dt <- data.table(taxid = uid, as.data.table(lineage))

tax_dt <- merge(tax_dt, lineage_dt, by = "taxid", all.x = TRUE)

## ---------------------------------------------------------------------
## 3. Flag eukaryotes, label unresolved phylum/class within Eukaryota
## ---------------------------------------------------------------------
tax_dt[, is_eukaryote := !is.na(domain) &  domain == "Eukaryota"]

cat("Eukaryotic contigs:", sum(tax_dt$is_eukaryote), "out of", nrow(tax_dt), "\n")

tax_dt[is_eukaryote == TRUE & is.na(phylum), phylum := "Unclassified_Eukaryota"]
tax_dt[is_eukaryote == TRUE & is.na(class),  class  := "Unclassified_Eukaryota"]

## ---------------------------------------------------------------------
## 4. Read coverage table, reshape to long format
## ---------------------------------------------------------------------
cov_wide <- fread(coverage_file, header = TRUE, sep = "\t")
setnames(cov_wide, names(cov_wide)[1], "contig_full")
sample_cols <- setdiff(names(cov_wide), "contig_full")

cov_long <- melt(cov_wide, id.vars = "contig_full", measure.vars = sample_cols,
                 variable.name = "sample", value.name = "coverage")

## ---------------------------------------------------------------------
## 5. Join taxonomy calls onto coverage
##    (contigs with no mmseqs2 call, or with a call that isn't Eukaryota,
##    are kept for the total but excluded from the Euk tables)
## ---------------------------------------------------------------------
merged <- merge(cov_long,
                tax_dt[, .(contig_full, is_eukaryote, phylum, class)],
                by = "contig_full", all.x = TRUE)
merged[is.na(is_eukaryote), is_eukaryote := FALSE]

## ---------------------------------------------------------------------
## 6. Phylum-level and Class-level abundance tables (Eukaryota only)
##
##    NOTE on the abundance metric: this sums each contig's CoverM
##    trimmed-mean depth across all contigs assigned to the same
##    phylum/class, per sample. This is NOT length-weighted -- a taxon
##    spread across many short contigs and one on a few long contigs get
##    summed the same way. If you have a contig-length table, tell me and
##    I'll switch this to a length-weighted sum instead.
## ---------------------------------------------------------------------
euk <- merged[is_eukaryote == TRUE]
euk$coverage[is.na(euk$coverage)] = 0

phylum_abund <- dcast(euk, phylum ~ sample, value.var = "coverage",
                      fun.aggregate = sum, fill = 0)

setnames(phylum_abund,
         old = names(phylum_abund),
         new = gsub("^final\\.contigs\\.renamed\\.filtered\\.fa\\.(\\d{2})(\\d{2})(\\d{2})_L([A-Za-z])_([A-Za-z]+)$",
                    "L_\\4_\\5_\\3_\\2_\\1",
                    names(phylum_abund)))
colnames(phylum_abund)
rownames(phylum_abund) = phylum_abund$phylum
phylum_abund$phylum = NULL

fwrite(phylum_abund, "eukaryote_phylum_abundance.tsv", sep = "\t")

## ---------------------------------------------------------------------
## 7. Relative abundance of eukaryotes, per sample
##    (total eukaryotic coverage / total coverage across ALL contigs,
##    classified or not, eukaryotic or not)
## ---------------------------------------------------------------------
totals     <- merged[, .(total_coverage = sum(coverage, na.rm = TRUE)), by = sample]
euk_totals <- merged[is_eukaryote == TRUE,
                     .(euk_coverage = sum(coverage, na.rm = TRUE)), by = sample]

euk_rel <- merge(totals, euk_totals, by = "sample", all.x = TRUE)
euk_rel[is.na(euk_coverage), euk_coverage := 0]
euk_rel[, relative_abundance := euk_coverage / total_coverage]

median(euk_rel$relative_abundance) *100 # 0.14 %
100 * (quantile(euk_rel$relative_abundance, probs = 0.75) - quantile(euk_rel$relative_abundance, probs = 0.25)) # 0.09 %

euk_rel$sample = gsub("^final\\.contigs\\.renamed\\.filtered\\.fa\\.(\\d{2})(\\d{2})(\\d{2})_L([A-Za-z])_([A-Za-z]+)$",
                    "L_\\4_\\5_\\3_\\2_\\1", euk_rel$sample)

euk_rel$Temp_soil = map_dbl(euk_rel$sample, function(x) metadata$Temp_soil_scaled[metadata$Sample == x])
euk_rel$Days_Since_Start = map_dbl(euk_rel$sample, function(x) metadata$Days_Since_Start[metadata$Sample == x])
euk_rel$Replicate = map_chr(euk_rel$sample, function(x) metadata$Replicate[metadata$Sample == x])
euk_rel$Replicate = as.factor(euk_rel$Replicate)

mod_euk_relab = gam(formula = relative_abundance ~ s(Days_Since_Start, k=5, bs = "ts") + s(Replicate, bs='re') + s(Temp_soil, bs='ts', k=5) + Condition, 
                  data=euk_rel)
summary(mod_euk_relab)

fwrite(euk_rel, "eukaryote_relative_abundance.tsv", sep = "\t")

## ---------------------------------------------------------------------
## Done -- quick look
## ---------------------------------------------------------------------
cat("\n--- Class-level Eukaryota abundance (data.frame `class_abund`) ---\n")
print(class_abund)
cat("\n--- Phylum-level Eukaryota abundance (data.frame `phylum_abund`) ---\n")
print(phylum_abund)
cat("\n--- Relative abundance of Eukaryota per sample (`euk_rel`) ---\n")
print(euk_rel)




#####
metadata = load_metadata()
metadata = check_normality_change(metadata, c("Temp_soil","Temp_air","Luminosity","CH4_Flux","CO2_Flux","Water_Content_Field",
                                              "TC","TOC","pH","Water_Content_Lab","F.","NO3.","Cl.",              
                                              "NO2.","PO4.","SO4.","Br.","Li.","NH4.","Ca2.",          
                                              "Sr2.","Na.","Mg2.","K.","Formate","Malate","Propionate",
                                              "Citrate","Lactate","Butyrate","Oxalate","Acetate","TIN","Other_Anions",
                                              "Other_Cations","Grouped_OC"))

phylum_abund <- phylum_abund[, metadata$Sample, drop = FALSE]
all(colnames(phylum_abund) == metadata$Sample)
rownames(metadata) = metadata$Sample
metadata$year_season = map_chr(1:nrow(metadata), function(i) paste0(metadata$Year[i],
                                                                    metadata$Season[i], 
                                                                    collapse = '_'))

phylum_abund = phylum_abund[rowSums(phylum_abund) > 0,]

out_mags = ancombc2(data = as.data.frame(phylum_abund), 
                    meta_data = metadata,
                    fix_formula = "Temp_soil_scaled", rand_formula = "(1|Replicate) + (1 | Replicate/year_season)",
                    p_adj_method = 'fdr', struc_zero = FALSE, 
                    prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                    alpha = 0.05, n_cl = 5, verbose = TRUE)
table(out_mags$res$diff_Temp_soil_scaled)
table(out_mags$res$passed_ss_Temp_soil_scaled)
table(out_mags$res$passed_ss_Temp_soil_scaled, out_mags$res$diff_Temp_soil_scaled)
write.table(out_mags$res, file = 'stats/euk_diff_abundance_bysoiltemp_fdr.csv', 
            quote = F, row.names = T, sep = ',')

library(ggplot2)
library(scales)

sample_order <- colnames(phylum_abund)

rel_dt <- as.data.table(phylum_abund, keep.rownames = "phylum")

## relative abundance = each sample column normalized to sum to 1
rel_dt[, (sample_order) := lapply(.SD, function(x) x / sum(x)), .SDcols = sample_order]

## fold everything past the top 7 phyla (by mean relative abundance) into "Other"
mean_abund <- rel_dt[, rowMeans(.SD), .SDcols = sample_order]
names(mean_abund) <- rel_dt$phylum
top_phyla <- names(sort(mean_abund, decreasing = TRUE))[seq_len(min(9, length(mean_abund)))]

rel_dt[, phylum_grouped := fifelse(phylum %in% top_phyla, phylum, "Other")]

rel_long <- data.table::melt(rel_dt, id.vars = "phylum_grouped", measure.vars = sample_order,
                 variable.name = "sample", value.name = "relative_abundance")
rel_long <- rel_long[, .(relative_abundance = sum(relative_abundance)),
                     by = .(phylum_grouped, sample)]

phylum_levels <- c(top_phyla, "Other")
rel_long[, phylum_grouped := factor(phylum_grouped, levels = phylum_levels)]
rel_long[, sample := factor(sample, levels = sample_order)]

rel_long$Treatment = map_chr(rel_long$sample, function(x) metadata$Condition[metadata$Sample == x])

p_phylum <- ggplot(rel_long, aes(x = sample, y = relative_abundance, fill = phylum_grouped)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative abundance",
       title = NULL) +
  theme_minimal(base_size = 11) + facet_grid(~Treatment, scales ='free', space='free') +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) + scale_fill_paletteer_d("RColorBrewer::Set3", name = 'Phylum') 

p_phylum 
ggsave('figures/SFig_5.pdf', width = 8, height = 6)
