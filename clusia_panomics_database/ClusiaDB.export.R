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

# Export features for supplement
feature.pseudogenes %>%
  write_tsv(paste(EXPORT_DIR, "feature.pseudogenes.tsv", sep = '/'), na = "")
