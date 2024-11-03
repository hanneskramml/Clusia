# Plot MSAs using msavisr or ggmsa package

library(tidyverse)
library(seqinr)
library(msa)
library(ggmsa)
library(stringi)

# define function for saving msa as fasta file
alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

# load motifs data
feature.motifs <- read_tsv(paste(DATA_ROOT, "Genomics/cis_regulatory_elements/TSS_run/Clusia.motifs.tsv", sep='/')) %>%
  select(Motif = motif_name, Gene = gene, Start = start, Stop = stop, pvalue.bonferroni = p_value_bonf, MatchedSequence = matched_sequence)

data.motifs <- clusia.data.cam %>%
  left_join(feature.motifs, by = join_by(Gene == Gene)) 
data.motifs <- data.motifs %>% filter(Motif != "CCA1" | is.na(Motif)) %>% 
  mutate(Product = str_replace_all(Product, "%2C", ",")) %>% 
  mutate(Motif = str_replace_all(Motif, "RVE1", "EE")) %>%
  mutate(Species = str_replace_all(Species, c("C3" = "C", "CAM" = "A")))
data.motifs$Repeat <- str_detect(data.motifs$MatchedSequence, regex("[:lower:]"))

# load promoter sequences 
Clusia.promoters <- readDNAStringSet("Clusia.cam_promoters.fna")

# select genes of interest for alignment
motifs.select <- data.motifs %>% filter(og.clusia %in% og) %>% arrange(Gene) 
genes <- motifs.select$Gene # list of GOIs
promoters.select <- Clusia.promoters[names(Clusia.promoters) %in% genes]

# create msa with ClustalOmega algorithm
alignment <- msa(promoters.select, method = "ClustalOmega") 

## msavisr ------------------------------------------------------------------------------------------------------------------------------------------------
# save msa as fasta (required for msavisr package)
filename <- # define filename
alignment2Fasta(alignment, filename)

# get regions of interest and colors for highlighting them
seq_fnc = Vectorize(function(x,y) seq(x,y)) # function for creating sequence of numbers

motifs.select <- motifs.select %>% mutate(ROI = seq_fnc(Start, Stop)) # add column containing positions of roi

# create list of vectors for rois
rois <- motifs.select %>% filter(!is.na(Motif)) %>% select(Gene, ROI, Motif) %>%
  split(1:nrow(motifs.select)) %>% lapply(transpose) %>% lapply(deframe) %>% lapply(unname) %>% lapply(unlist)

# create list of colors for rois (LHY, EE, LUX)
motifs.select <- motifs.select %>% mutate(Color = case_when((Motif == "LUX") ~ "#E69F00",
                                                            (Motif == "LHY") ~ "#339999",
                                                            (Motif == "EE") ~ "#0072B2")) 
roicolors <- motifs.select$Color %>% na.omit()

# # define rois and colors manually - example
# rois <- list(c("Cmi0.g42655", 1245:1254, "LHY"), c("Cmi0.g50089", 1224:1233, "LHY"),
#              c("Cmu00.g14109", 1223:1232, "LHY"), c("Cmu08.g106", 79:92, "LUX"),
#              c("Cro0.g33738", 1081:1090, "LHY"), c("Cro0.g47321", 1247:1256, "LHY"))
# roicolors <- c("#339999", "#339999", "#339999","#E69F00", "#339999", "#339999")

# plot msa
myref <- "" # define reference sequence
filename <- "" # define filename
aln <- msavisr(mymsa = filename, myref = myref,
               myroi = rois, roicolors = roicolors, basecolors = c("gray", "black", "white")) +
  scale_x_reverse()  + 
  ggplot2::theme(legend.position = "none",
                 axis.line.y = ggplot2::element_blank(), 
                 axis.line.x = ggplot2::element_blank(), 
                 axis.text.y = ggplot2::element_text(face = "bold"), 
                 axis.ticks.y = ggplot2::element_blank())
ggsave(filename, plot = aln, width = 5.0)

### NOTE!!!
# Reference sequence can only be at top or bottom. For the refseq to be at any position in the plot, change the source code of msavisr() function
# using trace(msavisr, edit = TRUE). 
# Remove lines 48 - 57 with the following code:
# olvls <- fasdf %>% dplyr::filter(curhead != myref) %>% dplyr::distinct(curhead)
# if (isTRUE(refontop)) {
#   olvls <- rbind(olvls, myref)
# }
# else {
#   olvls <- rbind(myref, olvls)
# }
# olvls <- olvls$curhead
# fasdf$curhead <- factor(fasdf$curhead, levels = olvls)
# rm(olvls)


## ggmsa ---------------------------------------------------------------------------------------------------------------------------------------------
# change class of alignment from MSADNAMultipleAlignment to DNAMultipleAlignment
class(alignment) <- "DNAMultipleAlignment"

highlights <- 79:92 # define region to be highlighted
NT_aln <- ggmsa(alignment, start, end, seq_name = TRUE, position_highlight = highlights, color = "Chemistry_NT")



