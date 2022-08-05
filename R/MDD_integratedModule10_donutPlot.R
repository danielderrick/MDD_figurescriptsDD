# The purpose of this script is to generate a donut plot showing feature type
# annotations for integrated module 10.
# This script generates Fig 6A.

###############################################################################
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

importRNAscript <- "MDD_import_RNAseq_ensembl.R"
geneListFile    <- "../misc/MDD_multiomics14Modules_TFannotations_module10_2.csv"
outDirPlots     <- "../plots/MDD_manuscript_figures/"

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}
# source(importRNAscript)

###############################################################################

module10 <- read_csv(geneListFile)
assay_cols <- RColorBrewer::brewer.pal(6, "Paired")

donutPlot <-
  module10 %>% 
  group_by(Category) %>% 
  summarize(n = n()) %>% 
  ggdonutchart("n", "n", 
                  fill = "Category") +
  ggtitle("Module 10 Features") +
  scale_fill_manual(values = c("RNA" = brewer.pal(6, "Paired")[3],
                               "lncRNA"       = brewer.pal(4, "Pastel2")[2],
                               "Kinase" = brewer.pal(4, "Pastel2")[3],
                               "Epigenetic" = brewer.pal(6, "Pastel2")[6],
                               "Transcription Factor"     = brewer.pal(4, "Pastel2")[4],
                               "RPPA Features" = brewer.pal(6, "Paired")[2]))

if (!dir.exists(outDirPlots)) {dir.create(outDirPlots, recursive = TRUE)}

pdf(sprintf("%s/MDD_F6A.pdf", outDirPlots), height = 6, width = 8)
donutPlot
dev.off()
