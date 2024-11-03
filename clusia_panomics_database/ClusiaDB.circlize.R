if(!requireNamespace('circlize', quietly = TRUE))
  install.packages("circlize")

library(tidyverse)
library(circlize)


# Visualization => opposite pairs of chromosomes
ticks <- c(0,20,40,60,80)
chrom_order <- c("CMU01", "CMU02", "CMU03", "CMU04", "CMU06", "CMU07", "CMU10", "CMU11", "CMU12", "CMU14", "CMU16", "CMU19", "CMU20", "CMU25",
                 "CMU26", "CMU23", "CMU29", "CMU27", "CMU18", "CMU13", "CMU17", "CMU21", "CMU24", "CMU28", "CMU08", "CMU15", "CMU05", "CMU22", "CMU09")
color <- c('#C4645C', '#E27F5C', '#F79D5E', '#FDBE63', '#FDDD98', '#F7F7D7', '#D1F4ED', '#A2E0FA', '#6CBCFE', '#668BFF', '#6E65FA', '#8F63E7', '#BD85EC', '#F4BDFF',
           '#F4BDFF', '#BD85EC', '#8F63E7', '#6E65FA', '#668BFF', '#6CBCFE', '#A2E0FA', '#D1F4ED', '#F7F7D7', '#FDDD98', '#FDBE63', '#F79D5E', '#E27F5C', '#E27F5C', '#C4645C',
           '#999999', '#999999', '#999999'
)

genome <- read.table(paste(DATA_ROOT, "Assembly", "Cmultiflora_v2.scaffolds.size.bed", sep = '/'), col.names = c("chr","start","end")) %>%
  mutate(chr = fct_relevel(chr, chrom_order),
         chr = fct_relevel(chr, c("CMU00.1", "CMU00.2", "CMU00.3"), after = Inf)) %>%
  arrange(chr)

synteny <- read.csv(paste(DATA_ROOT, "Pangenes/clusia.outgroup/results/syntenicBlock_coordinates.csv", sep = '/')) %>%
  filter(genome1 == "Clusia_multiflora_H1", genome2 == "Clusia_multiflora_H2", chr1 != chr2) %>%
  mutate(color = case_when(
    chr1 %in% c("CMU01", "CMU09") & chr2 %in% c("CMU01", "CMU09") ~ "#C4645C",
    chr1 %in% c("CMU02", "CMU05", "CMU22") & chr2 %in% c("CMU02", "CMU05", "CMU22") ~ "#E27F5C",
    chr1 %in% c("CMU03", "CMU15") & chr2 %in% c("CMU03", "CMU15") ~ "#F79D5E",
    chr1 %in% c("CMU04", "CMU08") & chr2 %in% c("CMU04", "CMU08") ~ "#FDBE63",
    chr1 %in% c("CMU06", "CMU28") & chr2 %in% c("CMU06", "CMU28") ~ "#FDDD98",
    chr1 %in% c("CMU07", "CMU24") & chr2 %in% c("CMU07", "CMU24") ~ "#F7F7D7",
    chr1 %in% c("CMU10", "CMU21") & chr2 %in% c("CMU10", "CMU21") ~ "#D1F4ED",
    chr1 %in% c("CMU11", "CMU17") & chr2 %in% c("CMU11", "CMU17") ~ "#A2E0FA",
    chr1 %in% c("CMU12", "CMU13") & chr2 %in% c("CMU12", "CMU13") ~ "#6CBCFE",
    chr1 %in% c("CMU14", "CMU18") & chr2 %in% c("CMU14", "CMU18") ~ "#668BFF",
    chr1 %in% c("CMU16", "CMU27") & chr2 %in% c("CMU16", "CMU27") ~ "#6E65FA",
    chr1 %in% c("CMU19", "CMU29") & chr2 %in% c("CMU19", "CMU29") ~ "#8F63E7",
    chr1 %in% c("CMU20", "CMU23") & chr2 %in% c("CMU20", "CMU23") ~ "#BD85EC",
    chr1 %in% c("CMU25", "CMU26") & chr2 %in% c("CMU25", "CMU26") ~ "#F4BDFF",
    .default = "#00000010"
  ))


genes <- genomicDensity(read.table(paste(DATA_ROOT, "Annotation/Cmultiflora_v2.2.annotation.longestIsoform.bed.gz", sep = '/')), window.size = 1e6)
pseudogenes <- genomicDensity(read.table(paste(DATA_ROOT, "Pseudogenes/Cmultiflora.pseudogenes.bed.gz", sep = '/')), window.size = 1e6)
repeats <- genomicDensity(read.table(paste(DATA_ROOT, "Repeats/Cmultiflora_v2.scaffolds.repeats.composite_masked.bed.gz", sep = '/')), window.size = 1e6)



pdf(paste(RESULTS_DIR, "figs", "circusplot.pdf", sep = '/'))

circos.clear()
circos.par(
  start.degree = 90,
  track.height = 0.1,
  cell.padding = c(0,0,0,0),
  gap.degree = c(rep(1, 13), 5, rep(1, 14), 2, 1, 1, 2)
)
circos.genomicInitialize(genome, plotType = NULL)

# Axes & Labels
circos.trackPlotRegion(ylim = c(0, 1), bg.border = color, bg.col = color, track.height = 0.05, panel.fun = function(x, y) {
  chr <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.genomicAxis(h = "top", sector.index = chr, major.at = ticks*10^6, labels = ticks, minor.ticks = 3, tickLabelsStartFromZero = TRUE, labels.cex = 0.4)
  circos.text(mean(xlim), mean(ylim), labels = if_else(!startsWith(chr, "CMU00"), chr, NA), facing = "bending", niceFacing = TRUE, cex = 0.5, font = 2, col = "white")
})

# Gene Density
circos.genomicTrackPlotRegion(genes, ylim = c(0, 0.4), track.height = 0.08, bg.border = NA, panel.fun = function(region, value, ...) {
  i <- get.cell.meta.data("sector.numeric.index")
  circos.genomicLines(region, value, area = TRUE, border = NA, baseline = 0, col = color[i])
})

# Putative pseudogenes that are fragmented and/or containing permature stop codons or framshifts
circos.genomicTrackPlotRegion(pseudogenes, bg.border = NA, panel.fun = function(region, value, ...) {
  i <- get.cell.meta.data("sector.numeric.index")
  circos.genomicLines(region, value, area = TRUE, border = NA, baseline = 0, col = color[i])
})

# Repeat Density
circos.genomicTrackPlotRegion(repeats, bg.border = NA, panel.fun = function(region, value, ...) {
  i <- get.cell.meta.data("sector.numeric.index")
  circos.genomicLines(region, value, area = TRUE, border = NA, baseline = "bottom", col = color[i])
})

# Rearrangements (syntenic links)
circos.genomicLink(
  synteny %>% select(chr1, startBp1, endBp1),
  synteny %>% select(chr2, startBp2, endBp2),
  col = synteny %>% pull(color),
  #col = scales::alpha(synteny %>% pull(color), alpha=0.4),
  border = NA
)

# Label groups (putative subgenomes)
text(0.9, 0.8, "H1", cex = 1.5)
text(-0.9, 0.8, "H2", cex = 1.5)

dev.off()
