################################################################################
# 6_a1_extract_mag_genes.R                      [ENV: warming_plaine_morte]
#
# Builds the protein-search database and annotation table for the
# metaproteomics: keeps only the predicted genes that sit on contigs belonging
# to MAGs, renames them 'sub_assembly_N_<contig>_<gene>', and writes the matching
# eggNOG annotation subset.
#
# Inputs : data/all_mags/MAGs_to_contigs.txt
#          data/functional/assembly_{0..3}.emapper.genepred.fasta
#          data/functional/assembly_{0..3}.emapper.annotations
# Outputs: data/functional/EVO_mags_genes_082026.fasta
#          data/proteomics/annotations_mags_EDU_2026.emapper.annotations
#            (this file used to be written by 0_functions_and_packages.R, which
#             depended on `all_genes` created here - a circular dependency)
# Next   : run 6_a2_impute_normalise_proteins.R in the DEP2 env
################################################################################
source('0_functions_and_packages.R')   # run from the repository root
library(seqinr)

# MAG -> contig table; gene predictions of each sub-assembly
mags_contigs <- read.table('data/all_mags/MAGs_to_contigs.txt',
                           header = FALSE, stringsAsFactors = FALSE,
                           col.names = c('mag_id', 'contig_id'))

all_genes <- list()

for (i in 0:3) {
  prefix <- paste0('sub_assembly_', i)
  
  # contigs (from MAGs_to_contigs.txt) that belong to this sub-assembly
  sub_contigs <- mags_contigs$contig_id[startsWith(mags_contigs$mag_id, prefix)]
  
  # gene predictions for this sub-assembly
  genes <- read.fasta(paste0('data/functional/assembly_', i, '.emapper.genepred.fasta'))
  
  # prodigal-style gene IDs are "<contig_id>_<gene_number>" -> strip the suffix to recover the contig
  gene_contigs <- sub('_[0-9]+$', '', names(genes))
  
  keep <- genes[gene_contigs %in% sub_contigs]
  names(keep) = map_chr(names(keep), function(x) paste0(prefix, '_', x))
  all_genes <- c(all_genes, keep)
  
  cat(prefix, ':', length(keep), 'of', length(genes), 'genes kept (in MAGs)\n')
}

write.fasta(sequences = all_genes, names = names(all_genes),
            file.out = 'data/functional/EVO_mags_genes_082026.fasta')

ann_files <- list.files("data/functional", pattern = "assembly_\\d\\.emapper.annotations$", full.names = TRUE)
ann_list <- lapply(ann_files, function(file) {
  prefix <- str_extract(basename(file), "assembly_\\d+")
  prefix <- paste0('sub_', prefix)
  lines <- read_lines(file)
  lines <- lines[!startsWith(lines, "##")]
  
  tmp <- tempfile()
  write_lines(lines, tmp)
  
  df <- read_tsv(tmp, show_col_types = FALSE)
  df %>%
    rename(query = `#query`) %>%
    mutate(query = paste0(prefix, "_", query)) %>%
    mutate(
      KEGG_ko = gsub("ko:", "", KEGG_ko),
      contig = str_remove(query, "_\\d+$")
    )
})
annotations <- bind_rows(ann_list)

cat('Total genes written to EVO_mags_genes_082026.fasta:', length(all_genes), '\n')

# eggNOG annotations restricted to the MAG genes
mags_annotations <- annotations[annotations$query %in% names(all_genes), ]
write.table(mags_annotations, file = 'data/proteomics/annotations_mags_EDU_2026.emapper.annotations',
            sep = '\t', quote = FALSE)
save_session_info("06a1")

