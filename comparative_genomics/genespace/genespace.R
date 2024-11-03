#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")

#detach("package:GENESPACE", unload = TRUE)
#devtools::install_github("jtlovell/GENESPACE@v1.1.8", upgrade = F)
#devtools::install_github("jtlovell/GENESPACE@v1.2.3", upgrade = F)

#library(GENESPACE)
devtools::load_all("~/git/GENESPACE")  # overwrite threshold for rerunning OF in syntenic orthogroups (Cmu), base: v1.2.3

source("parse_peptides.R")
#source("syntenic_orthogroups.R")

args <- commandArgs(trailingOnly = TRUE)
if(!is.null(args) && !is.na(args[1])) {
  wd <- args[1]
} else wd <- "clusia.vitis"


genomes_phytozome <- c(
  "Aquilegia_coerulea",
  "Arabidopsis_thaliana",
  "Manihot_esculenta",
  "Populus_trichocarpa",
  "Vitis_vinifera",
  "Ricinus_communis",
  "Linum_usitatissimum",
  "Salix_purpurea",
  "Lactuca_sativa",
  "Daucus_carota",
  "Carica_papaya",
  "Beta_vulgaris",
  "Chenopodium_quinoa",
  "Portulaca_amilis",
  "Arachis_hypogaea",
  "Glycine_max",
  "Medicago_truncatula",
  "Gossypium_hirsutum",
  "Theobroma_cacao",
  "Eucalyptus_grandis",
  "Kalanchoe_fedtschenkoi",
  "Kalanchoe_laxiflora",
  "Solanum_lycopersicum"
)

genomes_ncbi <- c(
  "Nelumbo_nucifera",
  "Artemisia_annua",
  "Brassica_napus",
  "Cucumis_sativus",
  "Coffea_arabica",
  "Camellia_sinensis",
  "Quercus_suber",
  "Juglans_regia",
  "Erythranthe_guttata",
  "Hevea_brasiliensis",
  "Jatropha_curcas",
  "Morus_notabilis",
  "Ziziphus_jujuba",
  "Rosa_chinensis",
  "Citrus_sinensis",
  "Ipomoea_nil"
)

genomes_custom <- c(
  "Fraxinus_excelsior",
  "Basella_alba",
  "Mesembryanthemum_cristallinum",
  "Sedum_album",
  "Garcinia_oblongifolia",
  "Hypericum_perforatum",
  "Ochna_serrulata",
  "Elaeocarpus_sylvestris",
  "Oxalis_sp",
  "Clusia_multiflora",
  "Clusia_minor",
  "Clusia_rosea",
  "Clusia_multiflora_PREDOM",
  "Clusia_multiflora_ALLELIC",
  "Clusia_multiflora_H1",
  "Clusia_multiflora_H2"
)

# NEW
clusia.subgenomes <- data.frame (
  genomes = genomes_custom[c(11,12,15,16)],
  ploidy = c(3, 4, 1, 1),
  outgroup = NA
)

clusia.outgroup <- data.frame (
  genomes = c(genomes_phytozome[5], genomes_custom[c(11,12,15,16)]),
  ploidy = c(1, 3, 4, 1, 1),
  outgroup = c(1, replicate(4, 0))
)

clusia.vitis <- data.frame (
  genomes = c(genomes_phytozome[5], genomes_custom[c(11,12,15,16)]),
  ploidy = c(1, 3, 4, 1, 1),
  #outgroup = replicate(5, 0)
  outgroup = NA
)

repo <- "repo"
path2mcscanx <- "~/git/MCScanX/"
path2orthofinder <- "orthofinder"
path2diamond <- "diamond"
#path2orthofinder <- "/home/user/kramml/git/OrthoFinder/orthofinder"
#path2diamond <- "/home/apps/conda/miniconda3/envs/diamond-2.1.3/bin/diamond"
dotplots <- "check"
nCores <- 10

data <- get(wd)
print(data)


# prepare input files for orthofinder and genespace
if(!dir.exists(paste(wd, "peptide", sep = "/"))) {

  # Sources: clusia, sedum (CoGe), mesembryanthemum (TVIR), Fraxinus (ensembl), 1KP
  if(any(genomes_custom %in% data[!data$outgroup | is.na(data$outgroup), 1])) {
    parsedPaths <- parse_annotations(
      genomeDirs = genomes_custom[genomes_custom %in% data[!data$outgroup | is.na(data$outgroup), 1]],
      faString = "faa$|faa\\.gz$",
      gffString = "gff3$|gff3\\.gz$",
      headerEntryIndex = 1,
      gffIdColumn = "ID",
      rawGenomeRepo = repo,
      genespaceWd = wd)
  }
  if(any(genomes_custom %in% data[!!data$outgroup, ]$genomes)) {
    parsedPaths <- parse_peptides(
      genomeDirs = genomes_custom[genomes_custom %in% data[!!data$outgroup, ]$genomes],
      faString = "faa$|faa\\.gz$",
      headerEntryIndex = 1,
      rawGenomeRepo = repo,
      genespaceWd = wd)
  }

  # phytozome
  if(any(genomes_phytozome %in% data[!data$outgroup, ]$genomes)) {
    parsedPaths <- parse_annotations(
      genomeDirs = genomes_phytozome[genomes_phytozome %in% data[!data$outgroup, ]$genomes],
      presets = "phytozome",
      faString = "\\.protein_primaryTranscriptOnly\\.fa",
      gffString = "\\.gene\\.gff3",
      rawGenomeRepo = repo,
      genespaceWd = wd)
  }
  if(any(genomes_phytozome %in% data[!!data$outgroup, ]$genomes)) {
    parsedPaths <- parse_peptides(
      genomeDirs = genomes_phytozome[genomes_phytozome %in% data[!!data$outgroup, ]$genomes],
      presets = "phytozome",
      faString = "\\.protein_primaryTranscriptOnly\\.fa",
      rawGenomeRepo = repo,
      genespaceWd = wd)
  }

  # ncbi
  if(any(genomes_ncbi %in% data[!data$outgroup, ]$genomes)) {
    parsedPaths <- parse_annotations(
      genomeDirs = genomes_ncbi[genomes_ncbi %in% data[!data$outgroup, ]$genomes],
      faString = ".*protein\\.faa$|.*protein\\.faa\\.gz$",
      headerEntryIndex = 1,
      gffIdColumn = "Name",
      rawGenomeRepo = repo,
      genespaceWd = wd)
  }
  if(any(genomes_ncbi %in% data[!!data$outgroup, ]$genomes)) {
    parsedPaths <- parse_peptides(
      genomeDirs = genomes_ncbi[genomes_ncbi %in% data[!!data$outgroup, ]$genomes],
      faString = ".*protein\\.faa$|.*protein\\.faa\\.gz$",
      headerEntryIndex = 1,
      rawGenomeRepo = repo,
      genespaceWd = wd)
  }
}

# -- initalize the run and QC the inputs
gsParam <- init_genespace(
  wd = wd,
  genomeIDs = data$genomes,
  ploidy = data$ploidy,
  ignoreTheseGenomes = if(all(is.na(data$outgroup))) NA else data[!!data$outgroup, ]$genomes,
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder,
  path2diamond = path2diamond,
  useHOGs = TRUE,               # use HOGs only in case of user-supplied species tree (Cminor/Crosea speciation event)
  #onlyOgAnchors = FALSE,       # default: TRUE,
  #onlyOgAnchorsSelf = FALSE,   # should only hits in orthogroups be considered for anchors in self-hits (particularly polyploids), default: TRUE
  onlyOgAnchorsSecond = TRUE,  # only orthogroup hits for homeolog block anchors in case of any ploidy, default: FALSE
  #nSecondaryHits = 1,           # search for secondary hits in case no outgroup defined, default: 0
  #nGapsSecond = 10,            # default
  dotplots = dotplots,
  nCores = nCores
)

# -- accomplish the run
out <- run_genespace(gsParam)
