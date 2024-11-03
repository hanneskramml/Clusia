library(tidyverse)


# Orthogroup counts based on diploid outgroup
# Select for orthogroups having at least half of all eudicot species present (remove lineage/species-specific likely paralogous orthogroups)
# Filter for overdispersed orthogroups having not more than 8*ploidy genes per homoeolog/species
# Outgroup having at least one gene present, TODO: Include others like highly diploidized Ath
feature.counts.og.vitis <- feature.counts.og %>%
  filter(Species > 25) %>%
  filter(Clusia_multiflora_H1 <= 8, Clusia_multiflora_H2 <= 8, Clusia_minor <= 24, Clusia_rosea <= 32, between(Vitis_vinifera, 1, 8)) %>%
  #filter((Clusia_multiflora_H1 > 0 | Clusia_multiflora_H2 > 0) & (Clusia_minor > 0 | Clusia_rosea > 0)) %>%
  mutate(Clusia_multiflora_H1 = 1.0 - (Vitis_vinifera - Clusia_multiflora_H1) / Vitis_vinifera) %>%
  mutate(Clusia_multiflora_H2 = 1.0 - (Vitis_vinifera - Clusia_multiflora_H2) / Vitis_vinifera) %>%
  mutate(Clusia_minor = 1.0 - (Vitis_vinifera - Clusia_minor) / Vitis_vinifera) %>%
  mutate(Clusia_rosea = 1.0 - (Vitis_vinifera - Clusia_rosea) / Vitis_vinifera) %>%
  select(Orthogroup, Clusia_multiflora_H1, Clusia_multiflora_H2, Clusia_minor, Clusia_rosea)

# Summary of genes included
feature.counts.og %>%
  filter(Species > 25) %>%
  filter(Clusia_multiflora_H1 <= 8, Clusia_multiflora_H2 <= 8, Clusia_minor <= 24, Clusia_rosea <= 32, between(Vitis_vinifera, 1, 8)) %>%
  select(Orthogroup, Clusia_multiflora_H1, Clusia_multiflora_H2, Clusia_minor, Clusia_rosea, Vitis_vinifera, Arabidopsis_thaliana) %>%
  pivot_longer(cols = !Orthogroup, names_to = "Species", values_to = "n") %>%
  group_by(Species) %>%
  summarise(counts = sum(n))

# Plot gene copy densities per outgroup
pdf(paste(RESULTS_DIR, "figs", "counts.og.vitis.density.pdf", sep = '/'))
feature.counts.og.vitis %>%
  pivot_longer(cols = !Orthogroup, names_to = "Species", values_to = "ratio") %>%
  mutate(Species = fct_relevel(Species, c("Clusia_multiflora_H1", "Clusia_multiflora_H2", "Clusia_minor", "Clusia_rosea"))) %>%
  ggplot(aes(ratio, stat(count))) +
  geom_density(alpha=0.8) +
  facet_grid(vars(Species), scales = "free_y") +
  scale_x_continuous(name="Relative gene copies per outgroup", breaks=0:10, limits=c(0, 10))
dev.off()

# Plot gene familiy expansion/contraction, Figure 2d
pdf(paste(RESULTS_DIR, "figs", "counts.og.vitis.variation.pdf", sep = '/'), height = 4, width = 8)
feature.counts.og.vitis %>%
  pivot_longer(cols = !Orthogroup, names_to = "Species", values_to = "ratio") %>%
  mutate(Species = fct_relevel(Species, rev(c("Clusia_multiflora_H1", "Clusia_multiflora_H2", "Clusia_minor", "Clusia_rosea")))) %>%
  mutate(type = if_else(ratio < 1, "lost", if_else(ratio == 1, "retained", "gained"))) %>%
  count(Species, type) %>%
  mutate(perc = paste(round(n/sum(n)*100, 1), "%"), .by = Species) %>%
  ggplot(aes(n, type, fill = Species)) +
  geom_col(position="dodge", alpha=0.8) +
  geom_text(aes(label = perc), hjust = 1, nudge_x = -1) +
  scale_fill_manual(values = gs_colors(12), limits = c("Clusia_multiflora_H1", "Clusia_multiflora_H2", "Clusia_minor", "Clusia_rosea")) +
  theme_minimal() +
  scale_y_discrete(name="Gene copy variation relative to outgroup") +
  scale_x_continuous(name="Number of gene families")
dev.off()
