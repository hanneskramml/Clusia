EXPORT_DIR <- paste(RESULTS_DIR, "tables", sep = '/')

# *** Export dataframes (Source Data, NCBI GEO, EMBL PRIDE) ***

data %>%
  mutate(across(c(Group.ref, Gene, Transcript), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  mutate(Chr = str_replace(Chr, pattern = "CMU(..)",  replacement ="CMA\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "data.tsv", sep = '/'), na = "")


# Figure 1
phenotyping.gasexchange %>%
  write_tsv(paste(EXPORT_DIR, "Fig1b.GasExchange.tsv", sep = '/'), na = "")
phenotyping.acidity.delta %>%
  write_tsv(paste(EXPORT_DIR, "Fig1c.DeltaTA.tsv", sep = '/'), na = "")
phenotyping.acidity.delta.aov %>%
  write_tsv(paste(EXPORT_DIR, "Fig1c.ANOVA.tsv", sep = '/'), na = "")

# Figure 3
feature.repeats.landscape %>%
  write_tsv(paste(EXPORT_DIR, "Fig3b.RepeatLandscape.tsv", sep = '/'), na = "")
feature.pseudogenes %>%
  mutate(chr = str_replace(chr, pattern = "CMU(..)",  replacement ="CMA\\1")) %>%
  mutate(across(c(pid, parent, overlap), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  write_tsv(paste(EXPORT_DIR, "Fig3c.ListOfPseudogenes.tsv", sep = '/'), na = "")
feature.counts.og.vitis %>%
  rename("Clusia_major_H1" = "Clusia_multiflora_H1", "Clusia_major_H2" = "Clusia_multiflora_H2") %>%
  write_tsv(paste(EXPORT_DIR, "Fig3d.GeneFamilies.tsv", sep = '/'), na = "")

# Figure 4
data.diploidization.plot %>%
  write_tsv(paste(EXPORT_DIR, "Fig4a.SummarizedDiploidization.tsv", sep = '/'), na = "")
data.diploidization %>%
  mutate(across(c(Group.ref, Gene, Transcript, Pseudo.id, Pseudo.parent), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  mutate(Chr = str_replace(Chr, pattern = "CMU(..)",  replacement ="CMA\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Fig4c.GenewiseDiploidization.tsv", sep = '/'), na = "")

# Figure 5
opengreenhouse.arduino.soil %>%
  write_tsv(paste(EXPORT_DIR, "Fig5a.SoilWater.tsv", sep = '/'), na = "")
opengreenhouse.photosynq.par %>%
  write_tsv(paste(EXPORT_DIR, "Fig5bc.PAR.tsv", sep = '/'), na = "")
opengreenhouse.photosynq %>%
  write_tsv(paste(EXPORT_DIR, "Fig5d.PhotosynQ.tsv", sep = '/'), na = "")
data.proteomics.plot %>%
  filter(Pathway == "Carboxylation") %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Fig5g.CarboxylatingProteins.tsv", sep = '/'), na = "")
data.transcriptomics %>%
  filter(Function == "PEPC-kinase") %>%
  mutate(across(c(Group.ref, Gene), ~str_replace(., pattern = "Cmu(..)",  replacement ="Cma\\1") )) %>%
  write_tsv(paste(EXPORT_DIR, "Fig5h.PPCK.tsv", sep = '/'), na = "")
data.proteomics.plot %>%
  filter(Pathway == "Decarboxylation") %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Fig5i.DecarboxylatingProteins.tsv", sep = '/'), na = "")

# Figure 6
data.transcriptomics.plot %>%
  filter(GeneFamily == "OG0005653", Homoeolog == "BAM3") %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Fig6a.BAM3.tsv", sep = '/'), na = "")
data.proteomics.plot %>%
  filter(GeneFamily == "OG0004814", Homoeolog == "PHS1") %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Fig6b.PHS1.tsv", sep = '/'), na = "")
data.transcriptomics.plot %>%
  filter(GeneFamily == "OG0008062", Homoeolog == "PGMP") %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Fig6c.PGMP.tsv", sep = '/'), na = "")


# Supplementary Figure 1
phenotyping.env.water %>%
  filter(DateTime >= ymd("2025-11-01"), DateTime < ymd("2025-12-05")) %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig1a.SoilWater.tsv", sep = '/'), na = "")
phenotyping.env.par %>%
  filter(DateTime >= ymd("2025-11-01"), DateTime < ymd("2025-12-05")) %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig1b.PAR.tsv", sep = '/'), na = "")
phenotyping.env.fytotron %>%
  filter(DateTime >= ymd("2025-11-01"), DateTime < ymd("2025-12-05")) %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig1b.ENV.tsv", sep = '/'), na = "")
phenotyping.acidity %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig1c.TA.tsv", sep = '/'), na = "")
phenotyping.acidity.delta %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig1d.DeltaTA.tsv", sep = '/'), na = "")

# Supplementary Figure 6
feature.goterms.diploidization %>%
  as_tibble() %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig6.GenicDiploidization.tsv", sep = '/'), na = "")
feature.goterms.repeats %>%
  as_tibble() %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig6.Repeats.tsv", sep = '/'), na = "")
feature.goterms.conserved %>%
  as_tibble() %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig6.ConservedGenes.tsv", sep = '/'), na = "")

# Supplementary Figure 7
data.diploidization.plot %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig7.SummarizedDiploidization.tsv", sep = '/'), na = "")

# Supplementary Figure 9 and 10
data.transcriptomics.plot %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig9_10.GeneExpression.tsv", sep = '/'), na = "")

# Supplementary Figure 11 and 12
data.proteomics.plot %>%
  mutate(Genes = str_replace_all(Genes, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  write_tsv(paste(EXPORT_DIR, "Supplementary_Fig11_12.ProteinAbundance.tsv", sep = '/'), na = "")


# TPM matrices for NCBI GEO submission
feature.expression %>%
  filter(str_starts(sample, "F")) %>%
  pivot_wider(id_cols = gene_id, names_from = sample, values_from = TPM) %>%
  mutate(gene_id = str_replace(gene_id, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  arrange(str_sub(gene_id, 1, 5), nchar(gene_id), gene_id) %>%
  write_tsv(paste(EXPORT_DIR, "Cmajor.TPM.tsv", sep = '/'), na = "")

feature.expression %>%
  filter(str_starts(sample, "R")) %>%
  pivot_wider(id_cols = gene_id, names_from = sample, values_from = TPM) %>%
  arrange(nchar(gene_id), gene_id) %>%
  write_tsv(paste(EXPORT_DIR, "Crosea.TPM.tsv", sep = '/'), na = "")


# Protein matrices for PRIDE
feature.proteins %>%
  filter(str_starts(sample, "F")) %>%
  pivot_wider(id_cols = transcript, names_from = sample, values_from = regulation) %>%
  mutate(transcript = str_replace(transcript, pattern = "Cmu(..)",  replacement ="Cma\\1")) %>%
  arrange(str_sub(transcript, 1, 5), nchar(transcript), transcript) %>%
  write_tsv(paste(EXPORT_DIR, "Cmajor.PROT.tsv", sep = '/'), na = "")

feature.proteins %>%
  filter(str_starts(sample, "M")) %>%
  pivot_wider(id_cols = transcript, names_from = sample, values_from = regulation) %>%
  arrange(str_sub(transcript, 1, 5), nchar(transcript), transcript) %>%
  write_tsv(paste(EXPORT_DIR, "Cminor.PROT.tsv", sep = '/'), na = "")

feature.proteins %>%
  filter(str_starts(sample, "R")) %>%
  pivot_wider(id_cols = transcript, names_from = sample, values_from = regulation) %>%
  arrange(str_sub(transcript, 1, 5), nchar(transcript), transcript) %>%
  write_tsv(paste(EXPORT_DIR, "Crosea.PROT.tsv", sep = '/'), na = "")
