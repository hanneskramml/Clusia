library(Rsamtools)
library(tidyverse)
# Script requires ClusiaDB.init.R


# gather list of pseudogenes
src.pseudogenes <- read_tsv(paste(DATA_ROOT, "Pseudogenes", "Cmajor.pseudogenes.tsv", sep = '/'), col_type = "ciiccdiiiiddicc") %>%
  mutate(length = as.integer(end - start + 1)) %>%
  select(chr, start, end, pid = id, parent = query,strand, length, 6:14)

# handle overlaps to annotated gene models
src.overlaps <- read_tsv(paste(DATA_ROOT, "Pseudogenes", "Cmajor.pseudogenes.overlaps.bed", sep = '/'), col_type = "ciicciicicii", col_names = FALSE, na = c("", "NA", ".")) %>%
  select(chr = X1, start = X2, end = X3, pid = X4, eid = X8, estart = X6, eend = X7, estrand = X10, elength = X11, eoverlap = X12) %>%
  mutate(start = as.integer(start+1)) %>%
  group_by(chr, start, end, pid) %>%
  summarise(
    overlap = dplyr::first(str_replace(eid, pattern = "(.*)\\.exon.*",  replacement ="\\1")),
    ostart = min(estart), oend = max(eend), ostrand = dplyr::first(estrand),
    exons = n(), elength = sum(elength), olength = sum(eoverlap)) %>%
  ungroup() %>%
  mutate(overlap = case_when(
    pid == "Cmu25.p9651" ~ "Cmu25.g1114.t1",  #GPT2 => Genebrowser: Pseudogene wrongly assigned to upstream gene VUP1 Cmu25.g1110
    pid == "Cmu20.p7305" ~ "Cmu20.g54.t1",    #NHD => Genebrowser: Pseudogene wrongly mapped to another upstream NHD fragment
    pid == "Cmu20.p7307" ~ "Cmu20.g61.t1",    #NHD => similiar
    .default = overlap
  ))

# Generate pileup from CANU alignment to filter for haplotype confusions
src.pileup <-
  pileup(paste(DATA_ROOT, "Alignment", "Cmajor.scaffolds.alignment.poly.bam", sep = '/'),
    scanBamParam = ScanBamParam(
      which = GRanges(
        src.pseudogenes %>%
          mutate(position = paste0(chr, ':', start, '-', end)) %>%
          pull(position)),
      flag = scanBamFlag(
        isSecondaryAlignment = FALSE,
        isUnmappedQuery = FALSE,
        isDuplicate = FALSE)),
    pileupParam = PileupParam(
      min_mapq = 20,
      include_insertions = TRUE,
      distinguish_strands = FALSE))

src.pileup %<>%
  pivot_wider(id_cols = c(which_label, pos), names_from = nucleotide, values_from = count) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(match = if_else(rowSums(pick(A:G)) > 0, 1, 0)) %>%
  relocate(ins = '+', del = '-', .after = match) %>%
  mutate(ins = if_else(match == 0 & ins > 0, 1, 0), del = if_else(match == 0 & del > 0, 1, 0)) %>%
  group_by(position = which_label) %>%
  summarise(pmatch = sum(match), pins = sum(ins), pdel = sum(del)) %>%
  mutate(position = str_extract(position, ".*(?=\\.)")) %>%
  distinct()


# create feature matrix and mark putative haplotype confusions
feature.pseudogenes <- src.pseudogenes %>%
  left_join(src.overlaps, by = join_by(chr, start, end, pid)) %>%
  mutate(position = paste0(chr, ':', start, '-', end)) %>%
  left_join(src.pileup, by = join_by(position)) %>%
  mutate(pfilter = if_else(pmatch > pdel, TRUE, FALSE)) %>%   # filter for large deletions while SNPs are still allowed
  replace_na(list(pfilter = FALSE)) %>%    # filter for lacking alignments
  mutate(mode = if_else(!pfilter, "filtered", if_else(!is.na(overlap), "annotated", "unannotated")))


# quality control
src.pseudogenes %>% nrow()
src.pseudogenes %>% filter(expect <= 1e-5) %>% filter(ident >= 0.2) %>% filter(frac >= 0.05) %>% filter(length >= 90) %>% nrow()

# overlap with multiple genes?
src.overlaps %>% nrow()
src.overlaps %>% distinct(pid) %>% nrow()

# Number of total, overlapping, and filtered pseudogenes
feature.pseudogenes %>%
  count(mode) %>%
  mutate(perc = n/sum(n)*100)


# Figure 3c
pdf(paste(RESULTS_DIR, "figs", "pseudogenes.pdf", sep = '/'), width = 5)
feature.pseudogenes %>%
  ggplot(aes(type, fill = mode)) +
  geom_bar(alpha=0.8) +
  theme_classic() +
  scale_fill_manual(name="Pseudogenes", breaks=c("annotated", "unannotated", "filtered"), labels=c("Pseudogenized genes", "Genic fragments", "Haplotype confusion"), values=c("#BE1E2D", "#133F66", "grey80")) +
  xlab("Classification of pseudogenes") +
  ylab("Number of pseudogenes (#)")
dev.off()


# Number of hapotigs per gene => 75% encounter two pseudo-haplotypes
# Pseudogenes follow same curve => haplotype collapsing/duplication unlikely for 3/4 of all genes/pseudogenes
pdf(paste(RESULTS_DIR, "figs", "counts.haplotigs.genes.pdf", sep = '/'), width = 4)
tmp.alignment %>%
  filter(Contig == str_sub(qname, end = 10)) %>%
  count(Gene) %>%
  mutate(haplotigs = if_else(n>=4, ">=4", as.character(n))) %>%
  mutate(haplotigs = fct_relevel(haplotigs, c("1","2","3",">=4"))) %>%
  ggplot(aes(x = haplotigs)) +
  geom_bar() +
  xlab("Number of haplotigs (#)") +
  ylab("Gene counts (#)") +
  theme_minimal()
dev.off()
