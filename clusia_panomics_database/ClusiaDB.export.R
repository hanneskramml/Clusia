EXPORT_DIR <- paste(RESULTS_DIR, "tables", sep = '/')

# *** Export dataframes ***

data %>%
  mutate(across(c(Group.ref, Gene, Transcript), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  mutate(Chr = str_replace(Chr, pattern = "CMU(..)",  replacement ="CMA\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "data.tsv", sep = '/'), na = "")


data.diploidization %>%
  mutate(across(c(Group.ref, Gene, Transcript, Pseudo.id, Pseudo.parent), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  mutate(Chr = str_replace(Chr, pattern = "CMU(..)",  replacement ="CMA\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "data.diploidization.tsv", sep = '/'), na = "")

data.diploidization.plot %>%
  write_tsv(paste(EXPORT_DIR, "data.diploidization.plot.tsv", sep = '/'), na = "")


data.transcriptomics %>%
  mutate(across(c(Group.ref, Gene), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  write_tsv(paste(EXPORT_DIR, "data.transcriptomics.tsv", sep = '/'), na = "")

data.transcriptomics.plot %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "data.transcriptomics.plot.tsv", sep = '/'), na = "")


data.proteomics %>%
  mutate(across(c(Group.ref, Gene), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  write_tsv(paste(EXPORT_DIR, "data.proteomics.tsv", sep = '/'), na = "")

data.proteomics.plot %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "data.proteomics.plot.tsv", sep = '/'), na = "")


# Export features for Supplemental Table 3
feature.pseudogenes %>%
  mutate(chr = str_replace(chr, pattern = "CMU(..)",  replacement ="CMA\\1")) %>%
  mutate(across(c(pid, parent, overlap), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  write_tsv(paste(EXPORT_DIR, "feature.pseudogenes.tsv", sep = '/'), na = "")


# TPM matrices for NCBI GEO submission
feature.expression %>%
  filter(str_starts(sample, "F")) %>%
  pivot_wider(id_cols = gene_id, names_from = sample, values_from = TPM) %>%
  mutate(gene_id = str_replace(gene_id, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  arrange(str_sub(gene_id, 1, 5), nchar(gene_id), gene_id) %>%
  write_tsv(paste(EXPORT_DIR, "feature.expression.major.tsv", sep = '/'), na = "")

feature.expression %>%
  filter(str_starts(sample, "R")) %>%
  pivot_wider(id_cols = gene_id, names_from = sample, values_from = TPM) %>%
  arrange(nchar(gene_id), gene_id) %>%
  write_tsv(paste(EXPORT_DIR, "feature.expression.rosea.tsv", sep = '/'), na = "")


# TPM matrices for PRIDE
feature.proteins %>%
  filter(str_starts(sample, "F")) %>%
  pivot_wider(id_cols = transcript, names_from = sample, values_from = regulation) %>%
  mutate(transcript = str_replace(transcript, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  arrange(str_sub(transcript, 1, 5), nchar(transcript), transcript) %>%
  write_tsv(paste(EXPORT_DIR, "feature.protein_abundance.major.tsv", sep = '/'), na = "")

feature.proteins %>%
  filter(str_starts(sample, "M")) %>%
  pivot_wider(id_cols = transcript, names_from = sample, values_from = regulation) %>%
  arrange(str_sub(transcript, 1, 5), nchar(transcript), transcript) %>%
  write_tsv(paste(EXPORT_DIR, "feature.protein_abundance.minor.tsv", sep = '/'), na = "")

feature.proteins %>%
  filter(str_starts(sample, "R")) %>%
  pivot_wider(id_cols = transcript, names_from = sample, values_from = regulation) %>%
  arrange(str_sub(transcript, 1, 5), nchar(transcript), transcript) %>%
  write_tsv(paste(EXPORT_DIR, "feature.protein_abundance.rosea.tsv", sep = '/'), na = "")
