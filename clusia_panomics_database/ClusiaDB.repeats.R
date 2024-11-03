library(tidyverse)


# Load and reclassify repeats based on Wicker et al. (2007)
src.repeats <- read_tsv(paste(DATA_ROOT, "Repeats", "Cmultiflora_v2.scaffolds.repeats.bed.gz", sep = '/'), col_types = "ciicicdddcccici", col_names = FALSE, na = c("", "NA", ".")) %>%
  mutate(order = str_split_i(X11,"/", 1), superfamily = str_split_i(X11,"/", 2), rstart = as.integer(if_else(X6 == "+", X12, X14))) %>%
  select(chr = X1, start = X2, end = X3, repeats = X4, order, superfamily, rstart, rend = X13, rstrand = X6, score = X5, div = X7, del = X8, ins = X9, id = X15) %>%
  mutate(
    order = case_when(
      order %in% c("DNAnona", "DNAauto") ~ "DNA",
      order == "Retroelement" ~ "LTR",  #Class1 (LTR/non-LTR)
      order == "TRIM" | superfamily == "TRIM" ~ "LTR/TRIM",
      order == "LARD" ~ "LTR/LARD",
      order == "non-LTR(SINE)" ~ "SINE",
      order == "LINE?" ~ "LINE",
      order %in% c("Unclassified", "Unknown", "nonLTR", "non", "Other", "RathE1_cons", "RathE2_cons", "RathE3_cons", "DIRS|TIR", "Helitron|LARD", "Helitron|TRIM", "LARD|TRIM", "LTR|DIRS", "LTR|TIR", "SINE|LARD", "SINE|TRIM") ~ "Unknown_Order",
      is.na(order) ~ "Unknown_Order",
      .default = order),
    superfamily = case_when(
      order == "hAT" | superfamily %in% c("hAT-Ac", "hAT-Tip100", "hAT-Tag1", "HAT") ~ "hAT",
      superfamily %in% c("Harbinger", "PIF-Harbinge") ~ "PIF-Harbinger",
      superfamily %in% c("TcMar", "Mariner", "Tc1") ~ "Tc1-Mariner",
      order == "Mutator" | superfamily %in% c("MuDR", "MULE", "MULE-MuDR", "MULEtir") ~ "Mutator",
      superfamily %in% c("CACTA", "En-Spm", "EnSpm", "CMC-EnSpm") ~ "CACTA/EnSpm",
      superfamily == "Stow" ~ "Stowaway",
      superfamily == "ERTBV-A" ~ "ERTBV",
      order == "Maverick" ~ "Maverick",
      ((is.na(superfamily) & ! order %in% c("Simple_repeat", "Low_complexity", "ARTEFACT", "subtelomere")) |  superfamily %in% c("Unknown", "Ukn", "rice", "Simple")) ~ "Unknown_Superfam",
      .default = superfamily)) %>%
  mutate(
    order = case_when(
      order %in% c("RC", "Helitron") | superfamily == "Helitron" ~ "DNA/RC",
      order == "TIR" | superfamily %in% c("PIF-Harbinger", "Tc1-Mariner", "Mutator", "hAT", "CACTA/EnSpm") ~ "DNA/TIR",
      order == "MITE" | superfamily %in% c("Stowaway", "Tourist", "Mite") ~ "DNA/MITE",
      order == "Maverick" ~ "DNA",
      .default = order),
    superfamily = case_when(
      order == "DNA/RC" ~ "Helitron",
      superfamily == "Unknown_Superfam" & order == "DNA" ~ "Unknown_DNA_Superfam",
      superfamily %in% c("Unknown_Superfam", "Mite") & order == "DNA/MITE" ~ "Unknown_MITE_Superfam",
      superfamily == "Unknown_Superfam" & order == "DNA/TIR" ~ "Unknown_TIR_Superfam",
      superfamily == "Unknown_Superfam" & order == "LTR" ~ "Unknown_LTR_Superfam",
      superfamily %in% c("Unknown_Superfam", "TRIM") & order == "LTR/TRIM" ~ "Unknown_TRIM_Superfam",
      superfamily == "Unknown_Superfam" & order == "LTR/LARD" ~ "Unknown_LARD_Superfam",
      superfamily == "Unknown_Superfam" & order == "LINE" ~ "Unknown_LINE_Superfam",
      superfamily == "Unknown_Superfam" & order == "SINE" ~ "Unknown_SINE_Superfam",
      superfamily == "Unknown_Superfam" & order == "DIRS" ~ "Unknown_DIRS_Superfam",
      .default = superfamily)) %>%
  arrange(chr, start)



# *** Kimura repeat landscapes ***

gs_colors <- function(n = 10){
  cols <- c("#C4645C", "#F5915C", "#FFC765", "#FCF8D5", "#BEF3F9", "#66B8FF", "#6666FF", "#9C63E1", "#F4BDFF")
  pal <- colorRampPalette(cols)
  return(pal(n))
}

src.repeats.kimura <- read_tsv(paste(DATA_ROOT, "Repeats", "Cmultiflora_v2.repeatlandscape.tsv.gz", sep = '/'), col_types = "ccii") %>%
  filter(!repeats == "combined")

feature.repeats.landscape <- src.repeats %>%
  distinct(repeats, order) %>%
  inner_join(
    src.repeats.kimura,
    by = join_by(repeats),
    relationship = "one-to-many") %>%
  left_join(
    read_tsv(paste(DATA_ROOT, "Assembly", "Cmultiflora_v2.scaffolds.size.bed", sep = '/'), col_types = "cii", col_names = FALSE) %>%
      select(chr = X1, size = X3),
    by = join_by(chr)) %>%
  inner_join(
    src.repeats.kimura %>%
      distinct(chr, repeats) %>%
      count(repeats),
    by = join_by(repeats)) %>%
  mutate(perc = bp/size/n)


pdf(paste(RESULTS_DIR, "figs", "repeat.landscape.pdf", sep = '/'), height = 4)
feature.repeats.landscape %>%
  mutate(order = str_split_i(order,"/", 1)) %>%
  mutate(order = if_else(order %in% c("Crypton", "DIRS", "Evirus", "MobileElement", "PLE", "rRNA", "subtelomere", "Satellite", "SINE"), "Other", order)) %>%
  group_by(order, div) %>%
  summarise(bp = sum(bp), perc = sum(perc)) %>%
  #filter(order == "LINE") %>%
  ggplot(aes(fill=order, y=perc, x=div)) +
  geom_bar(position="stack", stat="identity",color="black", alpha=0.8, linewidth = 0) +
  scale_fill_manual(values = gs_colors(12)) +
  #scale_fill_manual(values = "#E7845C") +
  #viridis::scale_fill_viridis(discrete = T, option = "E") +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") +
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12))
dev.off()

pdf(paste(RESULTS_DIR, "figs", "repeat.landscape.line.pdf", sep = '/'), height = 4)
feature.repeats.landscape %>%
  mutate(order = str_split_i(order,"/", 1)) %>%
  mutate(order = if_else(order %in% c("Crypton", "DIRS", "Evirus", "MobileElement", "PLE", "rRNA", "subtelomere", "Satellite", "SINE"), "Other", order)) %>%
  group_by(order, div) %>%
  summarise(bp = sum(bp), perc = sum(perc)) %>%
  filter(order == "LINE") %>%
  ggplot(aes(fill=order, y=perc, x=div)) +
  geom_bar(position="stack", stat="identity",color="black", alpha=0.8, linewidth = 0) +
  #scale_fill_manual(values = gs_colors(12)) +
  scale_fill_manual(values = "#E7845C") +
  #viridis::scale_fill_viridis(discrete = T, option = "E") +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") +
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12))
dev.off()



# *** Analysis of intronic/surrounding repeats ***

# Load features (*.bed files)
load_introns <- function(file) {
  file %>%
    read_tsv(col_types = "ciic-ci", col_names = FALSE) %>%
    rename_all(~ c("chr", "start", "end", "parent", "strand", "length")) %>%
    mutate(num = 1, id = paste(parent, ".i", cumsum(num), sep = ""), num = NULL, .by = "parent", .before = "parent")
}

load_intronic_repeats <- function(file) {
  file %>%
    read_tsv(col_types = "ciiccciciicicdddcccccii", col_names = FALSE, na = c("", "NA", ".")) %>%
    select(chr = X1, start = X2, end = X3, parent = X4, strand = X6, length = X7, rid = X22, rstart = X9, rend = X10, repeats = X11, roverlap = X23) %>%
    drop_na(repeats)
}

src.introns <- list.files(path = paste(DATA_ROOT, "Features", sep = '/'), pattern = "C*.introns.bed", full.names = TRUE) %>%
  lapply(load_introns) %>%
  bind_rows()

src.introns.repeats <-
  list.files(path = paste(DATA_ROOT, "Features", sep = '/'), pattern = "C*.introns.repeats.bed", full.names = TRUE) %>%
    lapply(load_intronic_repeats) %>%
    bind_rows() %>%
    arrange(chr, start)

# TODO: Create feature matrix of intronic repeats used for export and diploidization
