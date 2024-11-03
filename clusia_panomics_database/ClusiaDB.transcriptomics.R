library(tidyverse)
library(ggtext)
library(data.table)   # TODO: replace rbindlist() with tidyverse::bind_rows(..., .id = NULL)


# *** Transcript expression analysis ***

# Import transcript expression data (TPM), DEseq transformation (vsd) not applicable
load_expression <- function(path) {
  file <- basename(path)
  sample <- substr(file, 1, 6)
  replicate <- substr(file, 1, 3)
  condition <- substr(file, 2, 2)
  timepoint <- substr(file, 4, 6)

  expr <- cbind(
    sample = sample,
    replicate = replicate,
    condition = ifelse(condition == 'D', "C", ifelse(condition == 'H', "T", NA)),
    timepoint = timepoint,
    read_tsv(path),
    stringsAsFactors = TRUE)
}

feature.expression <-
  list.files(path = paste(DATA_ROOT, "Transcriptomics", sep = '/'), pattern = "*.genes.results", full.names = TRUE) %>%
    lapply(load_expression) %>%
    rbindlist()

# Create data matrix
data.transcriptomics <- data %>%
  filter(!Homoeolog.og %in% (meta.homoeologs.blacklist %>% pull(SynOG))) %>%
  left_join(feature.expression, by = join_by(Gene == gene_id)) %>%
  filter(!is.na(TPM)) %>%
  mutate(rep = case_when(
    replicate %in% c("FD1", "FH2", "RD6", "RH1") ~ "rep1",
    replicate %in% c("FD2", "FH3", "RD7", "RH5") ~ "rep2",
    replicate %in% c("FD3", "FH5", "RD8", "RH9") ~ "rep3",
    .default = NA)) %>%
  mutate(GeneName = if_else(!is.na(AltName) & AltName != Name, paste(Name, AltName, sep="/"), Name)) %>%
  mutate(TPM.mean = mean(TPM), TPM.sd = sd(TPM), .by = c(Gene, condition, timepoint)) %>%   # average over replicates
  pivot_wider(id_cols = c(1:13, GeneName, condition, timepoint, TPM.mean, TPM.sd), names_from = rep, values_from = TPM) %>%
  relocate(TPM.rep1 = rep1, TPM.rep2 = rep2, TPM.rep3 = rep3, .before = TPM.mean) %>%
  rename(Condition = condition, Timepoint = timepoint)

data.transcriptomics.plot <-
  data.transcriptomics %>%
  filter(!is.na(Pathway)) %>%
  group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group, Species, Condition, Timepoint) %>%
  summarise(Genes = paste0(unique(Gene), collapse = ", "), Genes.n = n(), TPM.sum = sum(TPM.mean)) %>%
  mutate(TPM.z = (TPM.sum-mean(TPM.sum))/sd(TPM.sum)) %>%
  relocate(Condition, Timepoint, .before = TPM.sum) %>%
  ungroup()


# Supplemental Figure 6

pdf(paste(RESULTS_DIR, "figs", "expression.cam.control.pdf", sep = '/'), width = 9, height = 15)
data.transcriptomics.plot %>%
  filter(Condition == "C") %>%
  #filter(Condition == "T") %>%
  mutate(Order = max(TPM.sum), .by = c(Pathway, Type, Function, GeneFamily, Homoeolog)) %>%
  mutate(Strip = paste0(Group, " *C.", str_split_i(Species, "_", 2), "*")) %>%
  mutate(y = paste0("**", Homoeolog, "** (", GeneFamily, ")")) %>%
  arrange(Order) %>%
  ggplot(aes(x = Timepoint, y = factor(y, levels = unique(y)), size = TPM.sum, color = TPM.z)) +
  geom_point() +
  facet_grid(rows = vars(Pathway), cols = vars(Strip), scales = "free", space = "free", switch = "y") +
  scale_size_area(name = "**Gene expression**<br>summed TPM", limits = c(0, 500), oob = scales::squish) +
  scale_color_gradientn(name = "**Circadian expression**<br>z-score within group", colours = viridis::viridis(20), limits = c(-1,1.5), oob = scales::squish) +
  guides(shape = guide_legend(order = 1), size = guide_legend(order = 2)) +
  xlab("Timepoint") +
  ylab("CAM-related genes (GeneFamily)") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", angle=45, vjust=1, hjust=1), strip.text = element_markdown(face = "bold", size = rel(0.6)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
  theme(axis.title.x = element_markdown(), axis.text.y = element_markdown(), legend.title = element_markdown())
dev.off()

pdf(paste(RESULTS_DIR, "figs", "expression.cam.treatment.pdf", sep = '/'), width = 9, height = 15)
data.transcriptomics.plot %>%
  #filter(Condition == "C") %>%
  filter(Condition == "T") %>%
  mutate(Order = max(TPM.sum), .by = c(Pathway, Type, Function, GeneFamily, Homoeolog)) %>%
  mutate(Strip = paste0(Group, " *C.", str_split_i(Species, "_", 2), "*")) %>%
  mutate(y = paste0("**", Homoeolog, "** (", GeneFamily, ")")) %>%
  arrange(Order) %>%
  ggplot(aes(x = Timepoint, y = factor(y, levels = unique(y)), size = TPM.sum, color = TPM.z)) +
  geom_point() +
  facet_grid(rows = vars(Pathway), cols = vars(Strip), scales = "free", space = "free", switch = "y") +
  scale_size_area(name = "**Gene expression**<br>summed TPM", limits = c(0, 500), oob = scales::squish) +
  scale_color_gradientn(name = "**Circadian expression**<br>z-score within group", colours = viridis::viridis(20), limits = c(-1,1.5), oob = scales::squish) +
  guides(shape = guide_legend(order = 1), size = guide_legend(order = 2)) +
  xlab("Timepoint") +
  ylab("CAM-related genes (GeneFamily)") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", angle=45, vjust=1, hjust=1), strip.text = element_markdown(face = "bold", size = rel(0.6)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
  theme(axis.title.x = element_markdown(), axis.text.y = element_markdown(), legend.title = element_markdown())
dev.off()


# Figure 5a
pdf(paste(RESULTS_DIR, "figs", "expression.bam3.pdf", sep = '/'), width = 9, height = 6)
data.transcriptomics.plot %>%
  filter(GeneFamily == "OG0005653", Homoeolog == "BAM3") %>%
  #filter(GeneFamily == "OG0004814", Homoeolog == "PHS1") %>%
  #filter(GeneFamily == "OG0008062", Homoeolog == "PGMP") %>%
  ggplot(aes(x = Timepoint, y = TPM.sum, color = Condition, group = Condition)) +
  geom_point() +
  geom_line() +
  #geom_errorbar(aes(ymin = TPM.sum-TPM.sd, ymax=TPM.sum+TPM.sd, width = .1)) +
  facet_grid(Group ~ Species, scales = "free") +
  theme_classic()
dev.off()


# Figure 4h
pdf(paste(RESULTS_DIR, "figs", "expression.pepc-kinase.pdf", sep = '/'))
data.transcriptomics %>%
  filter(Function == "PEPC-kinase") %>%
  filter(Species == "Clusia_multiflora") %>%
  #filter(Species == "Clusia_rosea") %>%
  mutate(Strip = paste(Homoeolog.og, Group)) %>%
  ggplot(aes(x=Timepoint, y=Gene, fill=TPM.mean)) +
  geom_tile(colour="black", size=0.25) +
  #coord_fixed() +
  facet_grid(Strip ~ Condition, scale = "free_y") +
  scale_fill_viridis_b() +
  theme_classic()
dev.off()