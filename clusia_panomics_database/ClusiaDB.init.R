# Clusia Panomics/feature database

# *** Package/library requirements ***
if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if(!requireNamespace('tidyverse', quietly = TRUE))
  install.packages('tidyverse')

if(!requireNamespace('cogeqc', quietly = TRUE))
  BiocManager::install('cogeqc')
if(!requireNamespace('Rsamtools', quietly = TRUE))
  BiocManager::install('Rsamtools')

if(!requireNamespace('GENESPACE', quietly = TRUE))
  devtools::install_github('jtlovell/GENESPACE')

#devtools::install_github("jtlovell/GENESPACE@v1.2.3", upgrade = F)
#devtools::load_all("~/git/Clusia/Submodules/GENESPACE")  # local dev => overwrite threshold for rerunning OF in syntenic orthogroups (Cmu), base: v1.2.3

library(cogeqc)
library(GENESPACE)
library(Rsamtools)
library(tidyverse)
library(readxl)


# *** Global parameters/settings ***
DATA_ROOT <- "~/git/Clusia/data"
CODE_DIR <- "~/git/Clusia/clusia_panomics_database"
RESULTS_DIR <- "~/git/Clusia/clusia_panomics_database"

SPECIES_CLUSIA <- c("Clusia_multiflora", "Clusia_minor", "Clusia_rosea")
SPECIES_OUTGROUP <- c(SPECIES_CLUSIA, "Vitis_vinifera", "Arabidopsis_thaliana")

options(dplyr.summarise.inform = FALSE)


# *** Import data sources ***

# Orthogroups / gene families (50+ eudicots sampling)
setwd(paste(DATA_ROOT, "Orthogroups/eudicots.haplotypes", sep = '/'))

src.orthogroups <- read_orthogroups("Orthogroups.tsv") %>%
  mutate(Transcript = Gene, Gene = str_replace(Gene, pattern = "(^[Cmu|Cmi|Cro].*)\\.t\\d+",  replacement ="\\1")) %>%     # longestIsoform => primaryTranscript
  mutate(Species = fct_relevel(Species, c("Clusia_multiflora_H1", "Clusia_multiflora_H2", "Clusia_multiflora_ALLELIC", "Clusia_minor", "Clusia_rosea"))) %>%    # reorder clusia species
  arrange(Orthogroup, Species, Gene, Transcript) %>%
  as_tibble()


# Pangenes / syntenic orthogroups / homoeologs
setwd(paste(DATA_ROOT, "Pangenes", sep = '/'))

load_pangenes <- function (gsParam, refGenome) {
  query_pangenes(gsParam = gsParam, refGenome = refGenome, transform = FALSE, showNSOrtho = FALSE, showUnPlacedPgs = FALSE) %>%
    mutate(repGene = id[pgRepID == ofID][1], og = og[flag == "PASS"][1], .by = "pgID") %>%
    mutate(genome = fct_relevel(genome, c("Clusia_multiflora_H1", "Clusia_multiflora_H2", "Clusia_minor", "Clusia_rosea"))) %>%
    left_join(
      read_combBed(file.path(gsParam$paths$results, "combBed.txt")) %>%
        select(pgRepID = ofID, repGenome = genome, repChr = chr),
      by = join_by(pgRepID)) %>%
    select(pgID, interpChr, interpOrd, og, repGene, repGenome, repChr, flag, id, genome, chr, ord, start, end) %>%
    as_tibble()
}

load("clusia.outgroup/results/gsParams.rda")
src.pangenes.H1 <- load_pangenes(gsParam, "Clusia_multiflora_H1")
src.pangenes.H2 <- load_pangenes(gsParam, "Clusia_multiflora_H2")


# Clusia annotation (*.tsv)
setwd(paste(DATA_ROOT, "Annotation", sep = '/'))

src.annotation <-
  c("Cmultiflora_v2.2.annotation.tsv", "Cminor_v1.2.annotation.tsv", "Crosea_v1.2.annotation.tsv") %>%
  read_tsv(col_types = "ccciiciicccccccccccddd", na = c("", "NA", "."))


# Join all data sources and create columns to be subsequently modified
data <- src.orthogroups %>%
  filter(str_starts(Species, "Clusia_multiflora") | Species %in% SPECIES_OUTGROUP) %>%
  left_join(
    src.pangenes.H1 %>%
      dplyr::union(src.pangenes.H2) %>%
      filter(interpChr == repChr) %>%
      distinct(id, og, flag),
    by = join_by(Transcript == id),
    relationship = "one-to-one") %>%
  left_join(
    src.annotation,
    by = join_by(Gene == Gene, Transcript == Transcript)) %>%
  relocate(Homoeolog.og = og, Homoeolog.flag = flag, .after = Orthogroup) %>%
  mutate(Homoeolog.flag = str_to_upper(Homoeolog.flag)) %>%
  mutate(Haplotype = if_else(str_starts(Species, "Clusia_multiflora"), str_replace_all(Species, c("Clusia_multiflora_H1" = "PREDOM", "Clusia_multiflora_H2" = "PREDOM", "Clusia_multiflora_ALLELIC" = "ALLELIC")), NA), .after = Species) %>%
  mutate(Group = if_else(Haplotype == "PREDOM", str_replace(Species, pattern = "Clusia_multiflora_(..)",  replacement ="\\1"), NA), .after = Homoeolog.flag) %>%
  mutate(Group.ref = if_else(Haplotype == "PREDOM", Gene, NA), .after = Group) %>%
  mutate(Species = if_else(str_starts(Species, "Clusia_multiflora"), "Clusia_multiflora", Species)) %>%
  mutate(Species = fct_relevel(Species, SPECIES_OUTGROUP)) %>%
  mutate(Haplotype = fct_relevel(Haplotype, c("PREDOM", "ALLELIC"))) %>%
  arrange(Orthogroup, Homoeolog.og, Species, Haplotype, Gene, Transcript)


# Load data mappings/metadata (map.xlsx)
setwd(CODE_DIR)

meta.genes <- read_excel("metadata.xlsx", sheet = "Genes", range = cell_cols("A:L"))
meta.families.revision <- read_excel("metadata.xlsx", sheet = "Families", range = cell_cols("A:C"), col_names = TRUE)
meta.families.blacklist <- read_excel("metadata.xlsx", sheet = "Families", range = cell_cols("E:F"), col_names = TRUE) %>% select(Orthogroup = Blacklist, Note)
meta.families.gene <- read_excel("metadata.xlsx", sheet = "Families", range = cell_cols("H:J"), col_names = TRUE)
meta.homoeologs.name <- read_excel("metadata.xlsx", sheet = "Homoeologs", range = cell_cols("A:B"))
meta.homoeologs.blacklist <- read_excel("metadata.xlsx", sheet = "Homoeologs", range = cell_cols("D:E")) %>% select(SynOG = Blacklist, Note)
meta.samples <- read_excel("metadata.xlsx", sheet = "Samples")



# *** Homoeolog resolution ***

# Get gene counts per orthogroup for eudicots sampling
# and further create gene families from orthogroups having at least half of all species present
# (removes lineage/species-specific likely paralogous orthogroups containing neofunctionalized or pseudogenized/fragmented genes)
if (!exists("feature.counts.og") || is.null(feature.counts.og))
  feature.counts.og <- src.orthogroups %>%
    select(-Transcript) %>%
    pivot_wider(names_from = Species, values_from = Gene, values_fn = length) %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    rowwise(Orthogroup) %>%
    mutate(Species = sum(c_across(1:53) > 0, na.rm = TRUE), Genes = sum(c_across(1:53), na.rm = TRUE)) %>%
    relocate(Species, Genes, .after = Orthogroup) %>%
    ungroup()

data %<>%
  left_join(
    feature.counts.og %>%
      select(Orthogroup, n_species = Species),
    by = join_by(Orthogroup)) %>%
  mutate(GeneFamily = if_else(n_species > 25, Orthogroup, NA), .before = Orthogroup)

# Incorporate manually curated gene families
data <- data %>%
  rows_update(meta.families.revision %>% select(Orthogroup, GeneFamily), by = "Orthogroup")
data <- data %>%
  rows_update(meta.families.gene %>% select(Gene, GeneFamily), by = "Gene")


# Assign underlying contigs to chromosomal genes (incorporate Cmu alleles based on unzipping/phasing information)
feature.contigs <- list.files(path = paste(DATA_ROOT, "Features", sep = '/'), pattern = "Cmultiflora.genes.contigs.bed", full.names = TRUE) %>%
  read_tsv(col_select = c(4, 1, 11, 15), col_names = FALSE, na = c("", "NA", ".")) %>%
  rename_all(~c("gene", "chr", "contig", "length")) %>%
  drop_na(contig) %>%
  group_by(gene) %>%
  filter(length == max(length)) %>%
  summarise(chr, contig)

data %<>%
  left_join(feature.contigs, by = join_by(Gene == gene, Chr == chr)) %>%
  mutate(Contig = if_else(!is.na(contig), contig, Chr), contig = NULL) %>%
  relocate(Contig, .after = Chr)


# Process clusia sequence alignments for every C. multiflora gene resolved at chrom-level (takes some time...)
load_alignment <- function (file) {
  cols <- c('qname', 'flag', 'strand', 'pos', 'qwidth', 'mapq', 'cigar')
  ranges <- data %>%
    filter(Haplotype == "PREDOM") %>%
    mutate(Position = paste0(Chr, ':', Start, '-', End)) %>%
    pull(Position)

  param <- ScanBamParam(
    what = cols,
    which = GRanges(ranges),
    flag = scanBamFlag(isSecondaryAlignment = FALSE))

  cat("Processing file: ", basename(file), "\n")

  scanBam(file, param = param) %>%
    enframe(name = "Position", value = "value") %>%
    unnest_wider(value) %>%
    unnest(cols = all_of(cols))
}

if (!exists("feature.bam") || is.null(feature.bam))
  feature.bam <- list.files(path = paste(DATA_ROOT, "Alignment", sep = '/'), pattern = "Cmultiflora.scaffolds.alignment.[contigs|Cminor|Crosea]*.bam$", full.names = TRUE) %>%
    lapply(load_alignment) %>%
    bind_rows()


# Resolve syntenic orthogroup hierarchy (Clusia-specific paralogs resulting from neofunctionalization or most likely artefacts (genic fragments))
propagate_synteny <- function (data, from, to) {
  data %>%
    distinct(pick(all_of(c(from, to)))) %>%
    drop_na() %>%
    group_by(pick(all_of(from))) %>%
    filter(n() == 1) %>%
    select(all_of(c(from, to)))
}


# Resolve alignment-based homologs of Cmu_ALLELIC, C. minor & C. rosea (handles "non-syntenic"/unresolved genes due to short sequences, gene artefacts, etc.)
# Propagate syntenic relationship to next higher level (group & gene family)
# Iterate until all unambigous homologous relationships are resolved (or after max. 5 runs)
for (i in 1:5) {

  tmp.alignment <- data %>%
    filter(Haplotype == "PREDOM") %>%
    mutate(Position = paste0(Chr, ':', Start, '-', End)) %>%
    select(GeneFamily, Homoeolog.og, Group, Gene, Position) %>%
    left_join(feature.bam, by = join_by(Position)) %>%
    arrange(Gene, qname)

  tmp.update <- data %>%
    filter(Species %in% SPECIES_CLUSIA, is.na(Homoeolog.og) | is.na(Group)) %>%
    mutate(id = str_sub(Contig, end = 10)) %>%
    inner_join(
      tmp.alignment %>%
        filter(!is.na(GeneFamily), !is.na(Homoeolog.og)) %>%
        group_by(GeneFamily, Homoeolog.og, Group, Gene, qname) %>%
        summarise(breaks = n()-1, width = sum(qwidth)) %>%
        mutate(id = str_sub(qname, end = 10)) %>%
        rename_with(~ paste0("aln.", .), -id & -GeneFamily),
      by = join_by(id, GeneFamily),
      relationship = "many-to-many") %>%
    group_by(GeneFamily, Gene) %>%
    filter(id == aln.qname) %>%     # map genes to predominant chrom-level haplotype
    filter(Gene != aln.Gene) %>%    # remove self matching alignments
    filter(n() == 1) %>%            # filter for unambigous homologues, TODO: handle tandem gene duplications located on a single contig
    mutate(Homoeolog.og = if_else(is.na(Homoeolog.og), aln.Homoeolog.og, Homoeolog.og)) %>%
    mutate(Homoeolog.flag = if_else(is.na(Homoeolog.flag), "ALN", Homoeolog.flag)) %>%
    mutate(Group = if_else(is.na(Group), aln.Group, Group)) %>%
    mutate(Group.ref = if_else(is.na(Group.ref), aln.Gene, Group.ref)) %>%
    select(GeneFamily, Gene, Homoeolog.og, Homoeolog.flag, Group, Group.ref)

  if (nrow(tmp.update) > 0) {
    cat("Iteration ", i, ": Resolving ", nrow(tmp.update), " gene(s)...\n")

    data <- data %>%
      rows_update(tmp.update, by = c("GeneFamily", "Gene"))

    data <- data %>%
      rows_update(propagate_synteny(.,"Homoeolog.og", "Group"), by = "Homoeolog.og")

    data <- data %>%
      rows_update(propagate_synteny(.,"Homoeolog.og", "GeneFamily"), by = "Homoeolog.og")

  } else { break }
}

# Define homoeologs / homoeologous groups using the most frequent gene names and manually curated pairs
data %<>%
  left_join(
    data %>%
      drop_na(Homoeolog.og) %>%
      select(Homoeolog.og, Name, AltName) %>%
      pivot_longer(cols = c(Name, AltName), names_to = NULL, values_to = "Homoeolog", values_drop_na = TRUE) %>%
      mutate(Homoeolog = str_to_upper(Homoeolog)) %>%
      count(Homoeolog.og, Homoeolog) %>%
      group_by(Homoeolog.og) %>%
      filter(n == max(n)) %>%
      #summarise(SynOG.name = paste0(SynOG.name, collapse = "/")),
      summarise(Homoeolog = dplyr::first(Homoeolog)),
    by = join_by(Homoeolog.og)) %>%
  relocate(Homoeolog, .before = Homoeolog.og)

# Incorporate manually curated homoeolog names (=> investigate mapping rules in case of error)
data <- data %>%
  rows_update(meta.homoeologs.name %>% select(Homoeolog.og = SynOG, Homoeolog = Name), by = "Homoeolog.og")

# TODO: Check feasibility to additionally propagate from Homoeolog to gene family or reintegrate single-copy genespace run (clusia.subgenomes)


# Map Ath-based prediction of subcellular location
feature.cellular_location <- read_excel(paste(DATA_ROOT, "Features", "TAIR10.subcellular_predictions.xlsx", sep = '/'))

feature.cellular_location %<>%
  add_column(Transcript = str_replace(feature.cellular_location$Gene, pattern = "^>(.*)",  replacement ="\\1"), .after = "Gene") %>%
  mutate(Gene = str_replace(feature.cellular_location$Gene, pattern = "^>(.*)\\.\\d+",  replacement ="\\1")) %>%
  rename(Prediction = "Final Prediction") %>%
  arrange(Gene, Transcript)

data %<>%
  left_join(
    feature.cellular_location %>%
      distinct(Gene, Prediction) %>%
      group_by(Gene) %>%
      summarise(Location = if_else(n() <= 1, paste0(Prediction, collapse = "/"), "Various")),
    by = join_by(Gene)) %>%
  relocate(Location, .after = Contig)

# TODO: TAIR10 prediction not accurate (e.g. PGMP) => run TargetP for our selected CAM genes



# *** Selection of CAM relevant gene families ***


# Select based on multiple identifiers: EC number (enzymes only), Arabidopsis locus & Uniprot ID
# in case of relationship error => Mapping problem due to constraint violation => investigate metadata.xlsx
select.ec <- meta.genes %>%
  filter(!is.na(ec_number)) %>%
  distinct(pathway, type, name, ec_number) %>%
  left_join(data, by = join_by(ec_number == EC), relationship = "one-to-many")

select.uniprot <- meta.genes %>%
  filter(!is.na(uniprot_id)) %>%
  distinct(pathway, type, name, gene, protein_description, cellular_location, uniprot_id) %>%
  left_join(data, by = join_by(uniprot_id == UniprotID), relationship = "one-to-many")

select.ath_locus <- meta.genes %>%
  filter(!is.na(ath_locus)) %>%
  distinct(pathway, type, name, gene, protein_description, cellular_location, ath_locus) %>%
  left_join(data, by = join_by(ath_locus == Gene), relationship = "one-to-many", keep = TRUE)


# Associate tibble/dataframe with CAM genes
# Propagate selection to all genes within selected orthogroups that are not blacklisted (e.g. transposase)
data %<>%
  left_join(
    select.ec %>% distinct(pathway, type, name, Orthogroup) %>%
      dplyr::union(select.ath_locus %>% distinct(pathway, type, name, Orthogroup)) %>%
      dplyr::union(select.uniprot %>% distinct(pathway, type, name, Orthogroup)) %>%
      select(Pathway = pathway, Type = type, Function = name, Orthogroup) %>%
      filter(!Orthogroup %in% (meta.families.blacklist %>% pull(Orthogroup))) %>%
      drop_na(Orthogroup),
    by = join_by(Orthogroup)) %>%
  mutate(Pathway = factor(Pathway, levels = unique(meta.genes$pathway)), Type = factor(Type, levels = c("enzyme", "transporter", "TF"))) %>%
  relocate(Pathway, Type, Function, .before = GeneFamily) %>%
  arrange(Pathway, Type, Function, GeneFamily, Homoeolog, Homoeolog.og, Group, Group.ref, Species, Haplotype, Gene, Transcript)


# *** Statistiscs ***

# Show number of selected genes per pathway/type/function
data %>% filter(Species %in% SPECIES_CLUSIA, !is.na(Pathway)) %>% nrow()
data %>% filter(Species %in% SPECIES_CLUSIA, !is.na(Pathway)) %>% distinct(GeneFamily) %>% drop_na() %>% nrow()
data %>% filter(Species %in% SPECIES_CLUSIA) %>% count(Pathway)
data %>% filter(Species %in% SPECIES_CLUSIA) %>% count(Type)
data %>% filter(Species %in% SPECIES_CLUSIA) %>% count(Pathway, Type)
data %>% filter(Species %in% SPECIES_CLUSIA) %>% count(Pathway, Function)

# CAM genes having an unresolved gene family
data %>% filter(is.na(GeneFamily), !is.na(Pathway)) %>% count(Species)

# Show list of ignored genes
select.ec %>%
  filter(Orthogroup %in% (meta.families.blacklist %>% pull(Orthogroup)), !is.na(Gene))

# not found, TODO: combine all selection types
select.ec %>% filter(is.na(Gene))
select.uniprot %>% filter(is.na(Gene))


# Total number of gene models (excluding likely genic fragments that were not assigned into an orthogroup)
data %>%
  count(Species, Haplotype)

# Gene models without functional annotation in either of the data sources
data %>%
  filter(Species %in% SPECIES_CLUSIA) %>%
  filter(Product == "hypothetical protein" & is.na(em_target)) %>%
  count(Species, Haplotype)

# Number of resolved gene models in Clusia spp.
data %>%
  filter(Species %in% SPECIES_CLUSIA) %>%
  count(Species, Haplotype, Group, Homoeolog.flag)
