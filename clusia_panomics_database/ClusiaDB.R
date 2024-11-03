DATA_ROOT <- "~/git/Clusia/data"
CODE_DIR <- "~/git/Clusia/clusia_panomics_database"
RESULTS_DIR <- "~/git/Clusia/clusia_panomics_database"

cat("Initializing panomics database...\n")
source("ClusiaDB.init.R", echo = TRUE)

cat("Incorporating pseudogenes...\n")
source("ClusiaDB.pseudogenes.R", echo = TRUE)

cat("Incorporating repeats...\n")
source("ClusiaDB.repeats.R", echo = TRUE)

cat("Creating dataframes and plots for diploidization...\n")
source("ClusiaDB.diploidization.R", echo = TRUE)

cat("Processing RNA-Seq...\n")
source("ClusiaDB.transcriptomics.R", echo = TRUE)

cat("Processing proteins...\n")
source("ClusiaDB.proteomics.R", echo = TRUE)

cat("Making circos plot...\n")
source("ClusiaDB.circlize.R", echo = TRUE)

cat("Making riparian plots...\n")
source("ClusiaDB.synteny.R", echo = TRUE)

cat("Making gene copy plot...\n")
source("ClusiaDB.counts.R", echo = TRUE)

cat("Exporting selected dataframes...\n")
source("ClusiaDB.export.R", echo = TRUE)

cat("DONE!\n")
