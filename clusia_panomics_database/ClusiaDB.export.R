EXPORT_DIR <- paste(RESULTS_DIR, "tables", sep = '/')

# *** Export dataframes ***

data %>%
  write_tsv(paste(EXPORT_DIR, "data.tsv", sep = '/'), na = "")

data.diploidization %>%
  write_tsv(paste(EXPORT_DIR, "data.diploidization.tsv", sep = '/'), na = "")
data.diploidization.plot %>%
  write_tsv(paste(EXPORT_DIR, "data.diploidization.plot.tsv", sep = '/'), na = "")

data.transcriptomics %>%
  write_tsv(paste(EXPORT_DIR, "data.transcriptomics.tsv", sep = '/'), na = "")
data.transcriptomics.plot %>%
  write_tsv(paste(EXPORT_DIR, "data.transcriptomics.plot.tsv", sep = '/'), na = "")

data.proteomics %>%
  write_tsv(paste(EXPORT_DIR, "data.proteomics.tsv", sep = '/'), na = "")
data.proteomics.plot %>%
  write_tsv(paste(EXPORT_DIR, "data.proteomics.plot.tsv", sep = '/'), na = "")

# Export features for Supplemental Table 3
feature.pseudogenes %>%
  write_tsv(paste(EXPORT_DIR, "feature.pseudogenes.tsv", sep = '/'), na = "")

# TPM matrices for NCBI GEO submission
feature.expression %>%
  filter(str_starts(sample, "F")) %>%
  pivot_wider(id_cols = gene_id, names_from = sample, values_from = TPM) %>%
  arrange(str_sub(gene_id, 1, 5), nchar(gene_id), gene_id)%>%
  write_tsv(paste(EXPORT_DIR, "feature.expression.multiflora.tsv", sep = '/'), na = "")

feature.expression %>%
  filter(str_starts(sample, "R")) %>%
  pivot_wider(id_cols = gene_id, names_from = sample, values_from = TPM) %>%
  arrange(nchar(gene_id), gene_id) %>%
  write_tsv(paste(EXPORT_DIR, "feature.expression.rosea.tsv", sep = '/'), na = "")