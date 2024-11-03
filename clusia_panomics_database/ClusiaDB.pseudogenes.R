library(tidyverse)


src.pseudogenes <- read_tsv(paste(DATA_ROOT, "Pseudogenes", "Cmultiflora.pseudogenes.tsv", sep = '/'), col_type = "ciiccdiiiiddicc") %>%
  mutate(length = as.integer(end - start + 1)) %>%
  select(chr, start, end, pid = id, parent = query,strand, length, 6:14)

src.overlaps <- read_tsv(paste(DATA_ROOT, "Pseudogenes", "Cmultiflora.pseudogenes.overlaps.bed", sep = '/'), col_type = "ciicciicicii", col_names = FALSE, na = c("", "NA", ".")) %>%
  select(chr = X1, start = X2, end = X3, pid = X4, eid = X8, estart = X6, eend = X7, estrand = X10, elength = X11, eoverlap = X12) %>%
  mutate(start = as.integer(start+1)) %>%
  group_by(chr, start, end, pid) %>%
  summarise(
    overlap = dplyr::first(str_replace(eid, pattern = "(.*)\\.exon.*",  replacement ="\\1")),
    ostart = min(estart), oend = max(eend), ostrand = dplyr::first(estrand),
    exons = n(), elength = sum(elength), olength = sum(eoverlap)) %>%
  ungroup() %>%
  mutate(overlap = case_when(
    pid == "Cmu25.p9651" ~ "Cmu25.g1112.t1",  #GPT2 => Genebrowser: Pseudogene wrongly assigned to upstream gene
    pid == "Cmu20.p7305" ~ "Cmu20.g54.t1",    #NHD => Genebrowser: Pseudogene wrongly mapped to another upstream NHD fragment
    pid == "Cmu20.p7307" ~ "Cmu20.g61.t1",    #NHD => similiar
    .default = overlap
  ))

feature.pseudogenes <- src.pseudogenes %>%
  left_join(src.overlaps, by = join_by(chr, start, end, pid))


# quality control
src.pseudogenes %>% nrow()
src.pseudogenes %>% filter(expect <= 1e-5) %>% filter(ident >= 0.2) %>% filter(frac >= 0.05) %>% filter(length >= 90) %>% nrow()

# overlap with multiple genes?
src.overlaps %>% nrow()
src.overlaps %>% distinct(pid) %>% nrow()

# Number of pseudogenes overlapping gene models
feature.pseudogenes %>%
  mutate(model = !is.na(overlap)) %>%
  count(model)


# Figure 2c
pdf(paste(RESULTS_DIR, "figs", "pseudogenes.pdf", sep = '/'), width = 5)

feature.pseudogenes %>%
  mutate(model = if_else(!is.na(overlap), "Yes", "No")) %>%
  ggplot(aes(type, fill = model)) +
  geom_bar(alpha=0.8) +
  theme_classic() +
  scale_fill_manual(name="Gene models", breaks=c("No", "Yes"), labels=c("Not annotated", "Overlapping"), values=c("grey80", "#CC6666")) +
  xlab("Classification of pseudogenes") +
  ylab("Number of pseudogenes (#)")

dev.off()
