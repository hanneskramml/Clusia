library(tidyverse)
library(ggtext)
library(data.table)   # TODO: replace rbindlist() with tidyverse::bind_rows(..., .id = NULL)


# *** Protein regulation analysis ***

# Import protein data
feature.proteins <- list.files(path = paste(DATA_ROOT, "Proteomics", sep = '/'), pattern = "C.*.regulation.xls", full.names = TRUE) %>%
  lapply(read_excel) %>%
  rbindlist(fill = TRUE) %>%
  drop_na(Accession) # remove contaminants

feature.proteins %<>%
  pivot_longer(cols = where(is.numeric), names_to = "SID", values_to = "regulation", values_drop_na = TRUE) %>%
  mutate(SID = as.integer(SID)) %>%
  left_join(meta.samples, by = join_by(SID == SID)) %>%
  mutate(replicate = substr(Code, 1, 3),
         condition = ifelse(substr(Code, 2, 2) == 'D', "C", ifelse(substr(Code, 2, 2) == 'H', "T", NA)),
         timepoint = substr(Code, 4, 6)) %>%
  select(sample = Code, replicate, condition, timepoint, transcript = Accession, regulation) %>%
  arrange(sample, replicate, condition, timepoint, transcript)

# Create data matrix
data.proteomics <- data %>%
  filter(!Homoeolog.og %in% (meta.homoeologs.blacklist %>% pull(SynOG))) %>%
  left_join(feature.proteins, by = join_by(Transcript == transcript)) %>%
  filter(!is.na(regulation)) %>%
  mutate(rep = case_when(
    replicate %in% c("FD1", "FH2", "RD6", "RH1", "MD1", "MH2") ~ "rep1",
    replicate %in% c("FD2", "FH3", "RD7", "RH5", "MD4", "MH3") ~ "rep2",
    replicate %in% c("FD3", "FH5", "RD8", "RH9", "MD6", "MH6") ~ "rep3",
    .default = NA)) %>%
  mutate(GeneName = if_else(!is.na(AltName) & AltName != Name, paste(Name, AltName, sep="/"), Name)) %>%
  mutate(REG.mean = mean(regulation), REG.sd = sd(regulation), .by = c(Gene, condition, timepoint)) %>%   # average over replicates
  pivot_wider(id_cols = c(1:13, GeneName, condition, timepoint, REG.mean, REG.sd), names_from = rep, values_from = regulation) %>%
  relocate(REG.rep1 = rep1, REG.rep2 = rep2, REG.rep3 = rep3, .before = REG.mean) %>%
  rename(Condition = condition, Timepoint = timepoint)

data.proteomics.plot <-
  data.proteomics %>%
    filter(!is.na(Pathway)) %>%
    group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group, Species, Condition, Timepoint) %>%
    summarise(Genes = paste0(unique(Gene), collapse = ", "), Genes.n = n(), REG.sum = sum(REG.mean)) %>%
    mutate(REG.z = (REG.sum-mean(REG.sum))/sd(REG.sum)) %>%
    relocate(Condition, Timepoint, .before = REG.sum) %>%
    ungroup()


# Supplemental Figure 7

pdf(paste(RESULTS_DIR, "figs", "proteins.cam.control.pdf", sep = '/'), width = 11, height = 8)
data.proteomics.plot %>%
  filter(Condition == "C") %>%
  #filter(Condition == "T") %>%
  mutate(Order = max(REG.sum), .by = c(Pathway, Type, Function, GeneFamily, Homoeolog)) %>%
  mutate(Strip = paste0(Group, " *C.", str_split_i(Species, "_", 2), "*")) %>%
  mutate(y = paste0("**", Homoeolog, "** (", GeneFamily, ")")) %>%
  arrange(Order) %>%
  ggplot(aes(x = Timepoint, y = factor(y, levels = unique(y)), size = REG.sum, color = REG.z)) +
  geom_point() +
  facet_grid(rows = vars(Pathway), cols = vars(Strip), scales = "free", space = "free", switch = "y") +
  scale_size_area(name = "**Protein abundance**<br>summed relative regulation", limits = c(0, 0.01), oob = scales::squish) +
  scale_color_gradientn(name = "**Circadian abundance**<br>z-score within group", colours = viridis::viridis(20), limits = c(-1,1.5), oob = scales::squish) +
  guides(shape = guide_legend(order = 1), size = guide_legend(order = 2)) +
  xlab("Timepoint") +
  ylab("CAM-related genes (GeneFamily)") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", angle=45, vjust=1, hjust=1), strip.text = element_markdown(face = "bold", size = rel(0.6)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
  theme(legend.title = element_markdown(), axis.title.x = element_markdown(), axis.text.y = element_markdown())
dev.off()

pdf(paste(RESULTS_DIR, "figs", "proteins.cam.treatment.pdf", sep = '/'), width = 11, height = 8)
data.proteomics.plot %>%
  #filter(Condition == "C") %>%
  filter(Condition == "T") %>%
  mutate(Order = max(REG.sum), .by = c(Pathway, Type, Function, GeneFamily, Homoeolog)) %>%
  mutate(Strip = paste0(Group, " *C.", str_split_i(Species, "_", 2), "*")) %>%
  mutate(y = paste0("**", Homoeolog, "** (", GeneFamily, ")")) %>%
  arrange(Order) %>%
  ggplot(aes(x = Timepoint, y = factor(y, levels = unique(y)), size = REG.sum, color = REG.z)) +
  geom_point() +
  facet_grid(rows = vars(Pathway), cols = vars(Strip), scales = "free", space = "free", switch = "y") +
  scale_size_area(name = "**Protein abundance**<br>summed relative regulation", limits = c(0, 0.01), oob = scales::squish) +
  scale_color_gradientn(name = "**Circadian abundance**<br>z-score within group", colours = viridis::viridis(20), limits = c(-1,1.5), oob = scales::squish) +
  guides(shape = guide_legend(order = 1), size = guide_legend(order = 2)) +
  xlab("Timepoint") +
  ylab("CAM-related genes (GeneFamily)") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", angle=45, vjust=1, hjust=1), strip.text = element_markdown(face = "bold", size = rel(0.6)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
  theme(legend.title = element_markdown(), axis.title.x = element_markdown(), axis.text.y = element_markdown())
dev.off()


# Figure 5b
pdf(paste(RESULTS_DIR, "figs", "proteins.phs1.pdf", sep = '/'), width = 10, height = 4)
data.proteomics.plot %>%
  group_by(Species) %>%
  filter(GeneFamily == "OG0004814", Homoeolog == "PHS1") %>%
  #filter(Homoeolog == "PGMP") %>%
  ggplot(aes(x = Timepoint, y = REG.sum, color = Species, group = Species)) +
  geom_point() +
  geom_line() +
  facet_grid(~ Condition) +
  theme_classic()
dev.off()


# Figure 4gi
pdf(paste(RESULTS_DIR, "figs", "proteins.carboxylation.heatmap.pdf", sep = '/'))
data.proteomics.plot %>%
  filter(Pathway == "Carboxylation") %>%
  mutate(y = paste(Homoeolog, Group, sep = '.')) %>%
  ggplot(aes(x=Timepoint, y=y, fill=REG.z)) +
  geom_tile(colour="black", size=0.25) +
  coord_fixed() +
  facet_grid(Species ~ Condition) +
  scale_fill_viridis_b() +
  theme_classic() +
  ylab("Carboxylating proteins")
dev.off()

pdf(paste(RESULTS_DIR, "figs", "proteins.decarboxylation.heatmap.pdf", sep = '/'))
data.proteomics.plot %>%
  filter(Pathway == "Decarboxylation") %>%
  mutate(y = paste(Homoeolog, Group, sep = '.')) %>%
  ggplot(aes(x=Timepoint, y=y, fill=REG.z)) +
  geom_tile(colour="black", size=0.25) +
  coord_fixed() +
  facet_grid(Species ~ Condition) +
  scale_fill_viridis_b() +
  theme_classic() +
  ylab("Decarboxylating proteins")
dev.off()


pdf(paste(RESULTS_DIR, "figs", "proteins.carboxylation.abundance.pdf", sep = '/'))
data.proteomics.plot %>%
  filter(Pathway == "Carboxylation") %>%
  group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group, Species, Condition) %>%
  summarise(x = sum(REG.sum)) %>%
  mutate(y = paste(Homoeolog, Group, sep = '.')) %>%
  ggplot(aes(x=x, y=y)) +
  facet_grid(Species ~ Condition) +
  geom_bar(stat='identity') +
  theme_classic() +
  xlab("Relative protein abundance") +
  ylab("Carboxylating proteins")
dev.off()

pdf(paste(RESULTS_DIR, "figs", "proteins.decarboxylation.abundance.pdf", sep = '/'))
data.proteomics.plot %>%
  filter(Pathway == "Decarboxylation") %>%
  group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group, Species, Condition) %>%
  summarise(x = sum(REG.sum)) %>%
  mutate(y = paste(Homoeolog, Group, sep = '.')) %>%
  ggplot(aes(x=x, y=y)) +
  facet_grid(Species ~ Condition) +
  geom_bar(stat='identity') +
  theme_classic() +
  xlab("Relative protein abundance") +
  ylab("Decarboxylating proteins")
dev.off()