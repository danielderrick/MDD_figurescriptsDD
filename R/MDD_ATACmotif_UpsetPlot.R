# The purpose of this script is to generate an UpSet plot 
# of the ATAC-seq motif features associated with each ligand
# treatment
#
library(eulerr)
library(tidyverse)
library(UpSetR)
library(xlsx)

colFile           <- "../misc/MDD_color_codes.csv"
sharedFeatureFile <- "../misc/Sup Table 3 MDD_multiomics_shared_features.xlsx"
uniqueFeatureFile <- "../misc/Sup Table 4 MDD_multiomics_unique_features.xlsx"
outDirPlots       <- "../plots/MDD_manuscript_figures"

col_raw <- read.csv(colFile,
                    stringsAsFactors = FALSE)

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

###############################################################################

col <- list(
  ligand = dplyr::slice(col_raw, 1:8),
  experimentalTimePoint = dplyr::slice(col_raw, 10:15),
  secondLigand = dplyr::slice(col_raw, 17:18),
  collection = dplyr::slice(col_raw, 26:28)
)

col <-
  lapply(col, function(x) {
    x <- setNames(x[, 2], x[, 1])
  })

###############################################################################

feats <- read.xlsx(featureFile, sheetIndex = 1) %>% 
  filter(Type == "motifs")

sharedFeatures <- 
  read.xlsx(sharedFeatureFile, sheetIndex = 1, 
            colClasses = c("character", "character",
                           "numeric", "character"),
            stringsAsFactors = FALSE) %>% 
  filter(Type == "motifs") %>% 
  filter(Log2FoldChange > 0)

uniqueFeatures <- 
  read.xlsx(uniqueFeatureFile, sheetIndex = 1,
            colClasses = c("character", "character",
                           "numeric", "character", "character", "character"),
            stringsAsFactors = FALSE) %>% 
  filter(Type == "motifs") %>% 
  filter(Direction == "Positive") %>%
  dplyr::select(-contains("Direction"))

motifFeatures <- 
  bind_rows(sharedFeatures,
            uniqueFeatures) %>% 
  split(.$ligand, .$feature)

###############################################################################

if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots, recursive = TRUE)
}

pdf(sprintf("%s/MDD_F4C_DC.pdf", outDirPlots),
    height = 4.5, width = 4.5)
lapply(motifFeatures, pull, feature) %>% 
  fromList() %>% 
  upset(nintersects = NA, nsets = 100, 
        sets.bar.color = col$ligand[rev(c("BMP2", "HGF", "EGF",
                                          "OSM", "TGFB", "IFNG"))])

lapply(motifFeatures, pull, feature) %>% 
  fromList() %>% 
  upset(nintersects = NA, nsets = 100, 
        order.by = "freq", 
        sets.bar.color = col$ligand[rev(c("BMP2", "HGF", "EGF",
                                          "OSM", "TGFB", "IFNG"))])
dev.off()  
