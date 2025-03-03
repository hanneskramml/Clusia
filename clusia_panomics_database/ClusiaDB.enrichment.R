if(!requireNamespace('clusterProfiler', quietly = TRUE))
  BiocManager::install('clusterProfiler')
if(!requireNamespace('org.At.tair.db', quietly = TRUE))
  BiocManager::install('org.At.tair.db')
if(!requireNamespace('enrichplot', quietly = TRUE))
  BiocManager::install('enrichplot')

library(enrichplot)


# Genic diploidization
tmp.goterms.select <-
  data.diploidization %>%
    filter(!is.na(Pathway)) %>%
    filter(Pseudo.evidence >= 1 | CDS.cov <= 50) %>%
    distinct(Orthogroup) %>%
    pull(Orthogroup)

tmp.goterms.genes <- data %>%
  filter(!is.na(Pathway), Species == "Arabidopsis_thaliana") %>%
  filter(Orthogroup %in% tmp.goterms.select) %>%
  select(Pathway, Type, Function, GeneFamily, Orthogroup, Gene, Location)

feature.goterms.diploidization <- clusterProfiler::enrichGO(gene = tmp.goterms.genes %>% pull(Gene), OrgDb = org.At.tair.db::org.At.tair.db,
                                                       keyType = "TAIR", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

pdf(paste(RESULTS_DIR, "figs", "goterms.diploidization.pdf", sep = '/'), height = 10, width = 6)
barplot(feature.goterms.diploidization, showCategory=20)
dev.off()


# Intron/repeat lengths
tmp.goterms.select <-
  data.diploidization %>%
    filter(!is.na(Pathway)) %>%
    filter(Repeat.length.z >= 1 | Intron.length.z >= 1) %>%
    distinct(Orthogroup) %>%
    pull(Orthogroup)

tmp.goterms.genes <- data %>%
  filter(!is.na(Pathway), Species == "Arabidopsis_thaliana") %>%
  filter(Orthogroup %in% tmp.goterms.select) %>%
  select(Pathway, Type, Function, GeneFamily, Orthogroup, Gene, Location)

feature.goterms.repeats <- clusterProfiler::enrichGO(gene = tmp.goterms.genes %>% pull(Gene), OrgDb = org.At.tair.db::org.At.tair.db,
                                                       keyType = "TAIR", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

pdf(paste(RESULTS_DIR, "figs", "goterms.repeats.pdf", sep = '/'), height = 10, width = 6)
barplot(feature.goterms.repeats, showCategory=20)
dev.off()


# Conserved genes
tmp.goterms.select <-
  data.diploidization %>%
    filter(!is.na(Pathway)) %>%
    filter(Pseudo.evidence < 1 & CDS.cov > 50 & Repeat.length.z < 1 & Intron.length.z < 1) %>%
    distinct(Orthogroup) %>%
    pull(Orthogroup)

tmp.goterms.genes <- data %>%
  filter(!is.na(Pathway), Species == "Arabidopsis_thaliana") %>%
  filter(Orthogroup %in% tmp.goterms.select) %>%
  select(Pathway, Type, Function, GeneFamily, Orthogroup, Gene, Location)

feature.goterms.conserved <- clusterProfiler::enrichGO(gene = tmp.goterms.genes %>% pull(Gene), OrgDb = org.At.tair.db::org.At.tair.db,
                                 keyType = "TAIR", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

pdf(paste(RESULTS_DIR, "figs", "goterms.conserved.pdf", sep = '/'), height = 10, width = 6)
barplot(feature.goterms.conserved, showCategory=20)
dev.off()



AnnotationDbi::keytypes(org.At.tair.db::org.At.tair.db)
AnnotationDbi::select(org.At.tair.db::org.At.tair.db, keys = tmp.goterms.genes %>% pull(Gene), keytype = "TAIR", column = c("ENTREZID", "SYMBOL")) %>%
  group_by(TAIR, ENTREZID) %>%
  summarise(SYMBOLS = paste0(SYMBOL, collapse = "/"))
