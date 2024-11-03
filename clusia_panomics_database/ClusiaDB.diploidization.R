if(!requireNamespace('viridis', quietly = TRUE))
  install.packages('viridis')
if (!requireNamespace("ggnewscale", quietly = TRUE))
  install.packages("ggnewscale")
if (!requireNamespace("ggtext", quietly = TRUE))
  install.packages("ggtext")

library(tidyverse)
library(ggnewscale)
library(ggtext)


# Script requires ClusiaDB.R, ClusiaDB.pseudogenes.R, ClusiaDB.repeats.R


# Create feature matrix and calculate genome-wide z-scores
data.diploidization <- data %>%
  filter(Haplotype == "PREDOM") %>%
  filter(!Homoeolog.og %in% (meta.homoeologs.blacklist %>% pull(SynOG))) %>%
  mutate(GeneName = if_else(!is.na(AltName) & AltName != Name, paste(Name, AltName, sep="/"), Name)) %>%
  select(1:15, last_col()) %>%
  left_join(
    src.introns %>%
      group_by(parent) %>%
      summarise(Intron.n = n(), Intron.max = max(length), Intron.length = sum(length)),
    by = join_by(Transcript == parent)) %>%
  left_join(
    src.introns.repeats %>%
      group_by(parent) %>%
      summarise(Repeat.n = n(), Repeat.length = sum(roverlap)),
    by = join_by(Transcript == parent)) %>%
  left_join(
    feature.pseudogenes %>%
      group_by(Gene = str_replace(overlap, pattern = "(^[Cmu|Cmi|Cro].*)\\.t\\d+",  replacement ="\\1")) %>%
      summarise(Pseudo.id = paste0(pid, collapse = ", "), Pseudo.parent = paste0(parent, collapse = ", "), Pseudo.frac = mean(frac), Pseudo.ins = sum(ins), Pseudo.del = sum(del), Pseudo.shift = sum(shift), Pseudo.stop = sum(stop), Pseudo.polya = sum(polya), Pseudo.ident = mean(ident), Pseudo.type = paste0(type, collapse = ", ")),
    by = join_by(Gene)) %>%
  replace_na(list(Intron.n = 0, Intron.length = 0, Intron.max = 0, Repeat.n = 0, Repeat.length = 0, Pseudo.frac = 0, Pseudo.ident = 0)) %>%
  mutate(Intron.length.z = (Intron.length-mean(Intron.length))/sd(Intron.length), .after = Intron.length) %>%
  mutate(Repeat.length.z = (Repeat.length-mean(Repeat.length))/sd(Repeat.length), .after = Repeat.length) %>%
  rowwise() %>%
  mutate(Pseudo.evidence = sum(c_across(c(Pseudo.shift, Pseudo.stop, Pseudo.polya)), na.rm = TRUE), .after = Pseudo.polya) %>%
  ungroup() %>%
  mutate(Pseudo.evidence.z = (Pseudo.evidence-mean(Pseudo.evidence))/sd(Pseudo.evidence), .after = Pseudo.evidence) %>%
  mutate(XFunct.non = if_else(!is.na(Pseudo.id), TRUE, FALSE)) %>%
  mutate(XFunct.neo = FALSE) %>%
  mutate(XFunct.sub = FALSE) %>%
  arrange(Pathway, Type, Function, GeneFamily, Homoeolog, Homoeolog.og, Species, Haplotype, Gene, Transcript)


# Number of pseudgenized gene models
data.diploidization %>%
  filter(!is.na(Pseudo.id)) %>%
  count(Group)

# Genome-wide intron and intronic repeat lengths
data.diploidization %>% arrange(desc(Intron.length))   #GWD3/PWD => Top1
data.diploidization %>% arrange(desc(Repeat.length))   #GWD3/PWD => Top10 (total intronic repeat length)

# Figure 3c
pdf(paste(RESULTS_DIR, "figs", "repeat.lengths.pdf", sep = '/'))
data.diploidization %>%
  ggplot(aes(x=Repeat.length, y=Intron.length)) +
  geom_point(alpha=0.1) +
  geom_point(data=. %>% drop_na(Pathway), aes(color=Pathway), na.rm=TRUE) +
  theme_classic() +
  xlab("Intronic TE length (bp)") +
  ylab("Intron length (bp)")
dev.off()


# Plot figures for genic diploidization
data.diploidization.plot <-
  data.diploidization %>%
  filter(!is.na(Pathway)) %>%
  pivot_wider(id_cols = c(Pathway, Type, Function, GeneFamily, Homoeolog), names_from = Group, values_from = Pseudo.frac, values_fn = mean) %>%
  replace_na(list(H1 = 1, H2 = 1)) %>%
  mutate(H1 = if_else(H1 > 0 & H1 < 1, 1-H1, H1)) %>%
  mutate(H2 = if_else(H2 > 0 & H2 < 1, 1-H2, H2)) %>%
  pivot_longer(cols = c(H1, H2), names_to = "Group", values_to = "Frac") %>%
  left_join(
    data.diploidization %>%
      filter(!is.na(Pathway)) %>%
      group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group) %>%
      summarise(Genes = n(), Pseudo = max(Pseudo.evidence), Intron = max(Intron.length.z), Repeat = max(Repeat.length.z), XFunct.non = sum(XFunct.non), XFunct.neo = sum(XFunct.neo), XFunct.sub = sum(XFunct.sub)),
    by = join_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group),
    relationship = "one-to-one") %>%
  replace_na(list(Genes = 0,  Pseudo = 0, Intron = 0, Repeat = 0, XFunct.non = 0, XFunct.neo = 0, XFunct.sub = 0)) %>%
  mutate(Order = max(Intron), .by = c(Pathway, Type, Function, GeneFamily, Homoeolog)) %>%
  mutate(Strip = max(Repeat) >= 1 | max(Pseudo) > 0, .by = c(Pathway, Type, Function, GeneFamily, Homoeolog)) %>%
  relocate(Genes, .after = Group)


# Create plot variants
pdf(paste(RESULTS_DIR, "figs", "diploidization.main.pdf", sep = '/'), width = 6, height = 8)

data.diploidization.plot %>%
  filter(Type != "TF") %>%
  pivot_longer(cols = c(Frac, Pseudo, Intron, Repeat, XFunct.non, XFunct.neo, XFunct.sub), names_to = "type", values_to = "value") %>%
  mutate(size = if_else(type %in% c("Intron", "Frac"), value, NA)) %>%
  mutate(fill = if_else(type %in% c("Repeat", "Pseudo"), value, NA)) %>%
  mutate(type = if_else(type %in% c("Intron", "Repeat"), "IntronicTEs", type)) %>%
  mutate(type = if_else(type == "Pseudo", "Frac", type)) %>%
  mutate(shape = if_else(str_starts(type, "XFunct"), value, NA)) %>%
  group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group, Genes, Order, Strip, type) %>%
  summarise(size = sum(size, na.rm = TRUE), fill = sum(fill, na.rm = TRUE), shape = if_else(sum(shape, na.rm = TRUE) > 0, TRUE, FALSE)) %>%
  mutate(type = fct_relevel(type, c("IntronicTEs", "Frac", "XFunct.non", "XFunct.neo", "XFunct.sub"))) %>%
  filter(str_starts(type, "XFunct") & shape | !str_starts(type, "XFunct")) %>%
  mutate(Group = if_else(str_starts(type, "XFunct"), "XF", Group)) %>%
  mutate(Homoeolog = if_else(Function == "AGPase", paste0(Function, "/", Homoeolog), Homoeolog)) %>%
  mutate(y = paste0("**", Homoeolog, "** (", GeneFamily, ")", sep = "")) %>%
  arrange(Order) %>%
  ggplot(mapping = aes(Group, factor(y, levels = unique(y)))) +
  geom_point(data = . %>% filter(type == "IntronicTEs"), aes(size = size, fill = fill, shape = type)) +
  scale_size_area(name = "**Intron size**<br>max. z-scored length (bp)", limits = c(0,4), breaks = c(0, 2, 4), oob = scales::squish) +
  scale_fill_gradientn(name = "**Intronic TEs**<br>max. z-scored length (bp)", colours = viridis::viridis(20), limits = c(-2,4), breaks = c(-2, 0, 2, 4), oob = scales::squish) +
  guides(size = guide_legend(order = 1)) +
  new_scale_fill() +
  new_scale("size") +
  geom_point(data = . %>% filter(type == "Frac"), aes(size = size, fill = fill, shape = type, alpha = if_else(fill > 0, 1, 0.5))) +
  scale_size_area(name = "Homoeolog<br>**Fractionation** (%)", limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_fill_gradientn(name = "Evidence for<br>**pseudogenization** (#)", colours = viridis::viridis(20, begin = 0.5), limits = c(0,4), breaks = c(0, 2, 4), oob = scales::squish) +
  scale_alpha_continuous(range = c(0.5, 1), guide = NULL) +
  guides(size = guide_legend(order = 1)) +
  new_scale_fill() +
  new_scale("size") +
  geom_point(data = . %>% filter(str_starts(type, "XFunct")), aes(size = 1, stroke = 1, shape = type), color = "#fde725", show.legend = F) +
  scale_shape_manual(name="Type of<br>**Diploidization**", breaks=c("IntronicTEs", "Frac", "XFunct.non", "XFunct.neo", "XFunct.sub"), labels=c("TEs within Introns", "Fractionation", "Nonfunctionalization", "Neofunctionalization<br>(not implemented)", "Subfunctionalization<br>(not implemented)"), values=c(21, 23, 4, 3, 8)) +
  guides(shape = guide_legend(order = 1)) +
  #facet_grid(Pathway ~ type, scales = "free", space = "free", switch = "y") +
  facet_grid(factor(Strip, levels = c(TRUE, FALSE), labels = c("Highly affected","Below threshold")) ~ type, scales = "free", space = "free", switch = "y") +
  xlab("Signals of genic diploidization") +
  ylab("CAM-related genes (GeneFamily)") +
  theme_classic() +
  theme(legend.title = element_markdown(), axis.title.x = element_markdown(), axis.text.y = element_markdown()) +
  theme(axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = rel(1)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
  theme(axis.text.x = element_text(face = "bold"), strip.text.x = element_blank(), strip.background.x = element_blank(), strip.text.y = element_text(face = "bold", size = rel(1)), strip.background.y = element_rect(fill = "gray95", colour = "black", linewidth = 0))
  #theme(axis.text.x = element_text(face = "bold"), strip.text.x = element_blank(), strip.background.x = element_blank(), strip.text.y = element_text(face = "bold", size = rel(0.6)), strip.background.y = element_rect(fill = "gray95", colour = "black", linewidth = 0))

# TODO: gene counts next to every homology group incl. Cmi & Cro => clearly shows single-copy gene families (disomic inheritance) and gene copies in agreement with polyploidy levels for all clusia species
dev.off()


pdf(paste(RESULTS_DIR, "figs", "diploidization.pathways.pdf", sep = '/'), width = 6, height = 8.5)

data.diploidization.plot %>%
  filter(Type != "TF") %>%
  pivot_longer(cols = c(Frac, Pseudo, Intron, Repeat, XFunct.non, XFunct.neo, XFunct.sub), names_to = "type", values_to = "value") %>%
  mutate(size = if_else(type %in% c("Intron", "Frac"), value, NA)) %>%
  mutate(fill = if_else(type %in% c("Repeat", "Pseudo"), value, NA)) %>%
  mutate(type = if_else(type %in% c("Intron", "Repeat"), "IntronicTEs", type)) %>%
  mutate(type = if_else(type == "Pseudo", "Frac", type)) %>%
  mutate(shape = if_else(str_starts(type, "XFunct"), value, NA)) %>%
  group_by(Pathway, Type, Function, GeneFamily, Homoeolog, Group, Genes, Order, Strip, type) %>%
  summarise(size = sum(size, na.rm = TRUE), fill = sum(fill, na.rm = TRUE), shape = if_else(sum(shape, na.rm = TRUE) > 0, TRUE, FALSE)) %>%
  mutate(type = fct_relevel(type, c("IntronicTEs", "Frac", "XFunct.non", "XFunct.neo", "XFunct.sub"))) %>%
  filter(str_starts(type, "XFunct") & shape | !str_starts(type, "XFunct")) %>%
  mutate(Group = if_else(str_starts(type, "XFunct"), "XF", Group)) %>%
  mutate(Homoeolog = if_else(Function == "AGPase", paste0(Function, "/", Homoeolog), Homoeolog)) %>%
  mutate(y = paste0("**", Homoeolog, "** (", GeneFamily, ")", sep = "")) %>%
  arrange(Order) %>%
  ggplot(mapping = aes(Group, factor(y, levels = unique(y)))) +
  geom_point(data = . %>% filter(type == "IntronicTEs"), aes(size = size, fill = fill, shape = type)) +
  scale_size_area(name = "**Intron size**<br>max. z-scored length (bp)", limits = c(0,4), breaks = c(0, 2, 4), oob = scales::squish) +
  scale_fill_gradientn(name = "**Intronic TEs**<br>max. z-scored length (bp)", colours = viridis::viridis(20), limits = c(-2,4), breaks = c(-2, 0, 2, 4), oob = scales::squish) +
  guides(size = guide_legend(order = 1)) +
  new_scale_fill() +
  new_scale("size") +
  geom_point(data = . %>% filter(type == "Frac"), aes(size = size, fill = fill, shape = type, alpha = if_else(fill > 0, 1, 0.5))) +
  scale_size_area(name = "Homoeolog<br>**Fractionation** (%)", limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_fill_gradientn(name = "Evidence for<br>**pseudogenization** (#)", colours = viridis::viridis(20, begin = 0.5), limits = c(0,4), breaks = c(0, 2, 4), oob = scales::squish) +
  scale_alpha_continuous(range = c(0.5, 1), guide = NULL) +
  guides(size = guide_legend(order = 1)) +
  new_scale_fill() +
  new_scale("size") +
  geom_point(data = . %>% filter(str_starts(type, "XFunct")), aes(size = 1, stroke = 1, shape = type), color = "#fde725", show.legend = F) +
  scale_shape_manual(name="Type of<br>**Diploidization**", breaks=c("IntronicTEs", "Frac", "XFunct.non", "XFunct.neo", "XFunct.sub"), labels=c("TEs within Introns", "Fractionation", "Nonfunctionalization", "Neofunctionalization<br>(not implemented)", "Subfunctionalization<br>(not implemented)"), values=c(21, 23, 4, 3, 8)) +
  guides(shape = guide_legend(order = 1)) +
  facet_grid(Pathway ~ type, scales = "free", space = "free", switch = "y") +
  #facet_grid(factor(Strip, levels = c(TRUE, FALSE), labels = c("Highly affected","Below threshold")) ~ type, scales = "free", space = "free", switch = "y") +
  xlab("Signals of genic diploidization") +
  ylab("CAM-related genes (GeneFamily)") +
  theme_classic() +
  theme(legend.title = element_markdown(), axis.title.x = element_markdown(), axis.text.y = element_markdown()) +
  theme(axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = rel(1)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
  #theme(axis.text.x = element_text(face = "bold"), strip.text.x = element_blank(), strip.background.x = element_blank(), strip.text.y = element_text(face = "bold", size = rel(1)), strip.background.y = element_rect(fill = "gray95", colour = "black", linewidth = 0))
  theme(axis.text.x = element_text(face = "bold"), strip.text.x = element_blank(), strip.background.x = element_blank(), strip.text.y = element_text(face = "bold", size = rel(0.6)), strip.background.y = element_rect(fill = "gray95", colour = "black", linewidth = 0))

dev.off()
