library(tidyverse)
library(GENESPACE)

setwd(paste(DATA_ROOT, "Pangenes", sep = '/'))
load("clusia.subgenomes/results/gsParams.rda")

# Riparian plots of homoeologous chromosomes in C. multiflora (Fig. 2a, Supplemental Fig. 4)
pdf(paste(RESULTS_DIR, "figs", "riparianplot.homoeologs.generank.pdf", sep = '/'), width = 10, height = 3)
ripd <- plot_riparian(gsParam = gsParam, useRegions = FALSE, useOrder = TRUE, forceRecalcBlocks = TRUE,
                      refGenome = "Clusia_multiflora_H1",
                      genomeIDs = c("Clusia_multiflora_H2", "Clusia_multiflora_H1"),
                      chrFill = "grey90", addThemes = theme_classic(), chrLabFontSize = 8,
                      invertTheseChrs = data.frame(genome = "Clusia_multiflora_H2", chr = c("CMU08","CMU28","CMU24","CMU21","CMU18","CMU23","CMU26")),
                      chrLabFun = function(x) gsub("^0", "", gsub("chr|cmu|cmi.*|cro.*", "", tolower(x))))
dev.off()

pdf(paste(RESULTS_DIR, "figs", "riparianplot.homoeologs.bp.pdf", sep = '/'), width = 10, height = 3)
ripd <- plot_riparian(gsParam = gsParam, useRegions = FALSE, useOrder = FALSE, forceRecalcBlocks = TRUE,
                      refGenome = "Clusia_multiflora_H1",
                      genomeIDs = c("Clusia_multiflora_H2", "Clusia_multiflora_H1"),
                      chrFill = "grey90", addThemes = theme_classic(), chrLabFontSize = 8,
                      invertTheseChrs = data.frame(genome = "Clusia_multiflora_H2", chr = c("CMU08","CMU28","CMU24","CMU21","CMU18","CMU23","CMU26")),
                      chrLabFun = function(x) gsub("^0", "", gsub("chr|cmu|cmi.*|cro.*", "", tolower(x))))
dev.off()


# Riparian plots of Clusia spp.
pdf(paste(RESULTS_DIR, "figs", "riparianplot.clusia.pdf", sep = '/'), width = 10, height = 3)
ripd <- plot_riparian(gsParam = gsParam, useRegions = FALSE, useOrder = TRUE, forceRecalcBlocks = FALSE,
                      refGenome = "Clusia_multiflora_H1", syntenyWeight = 1,
                      chrFill = "grey90", addThemes = theme_classic(), chrLabFontSize = 8,
                      chrLabFun = function(x) gsub("^0", "", gsub("chr|cmu|cmi.*|cro.*", "", tolower(x))))
dev.off()


# Dotplot of homoeologous chromosomes (Supplemental Fig. 4)
load("clusia.outgroup/results/gsParams.rda")

pdf(paste(RESULTS_DIR, "figs", "dotplot.H1.pdf", sep = '/'))
hits <- read_allBlast(filepath = file.path(gsParam$paths$syntenicHits, "Clusia_multiflora_H1_vs_Clusia_multiflora_H1.allBlast.txt.gz"))
ggdotplot(hits = hits, type = "raw", verbose = FALSE)
dev.off()

pdf(paste(RESULTS_DIR, "figs", "dotplot.H2.pdf", sep = '/'))
hits <- read_allBlast(filepath = file.path(gsParam$paths$syntenicHits, "Clusia_multiflora_H2_vs_Clusia_multiflora_H2.allBlast.txt.gz"))
ggdotplot(hits = hits, type = "raw", verbose = FALSE)
dev.off()


setwd(CODE_DIR)
