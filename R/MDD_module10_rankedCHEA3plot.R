# The purpose of this script is to generate a donut plot showing feature
# This script generates Fig 6A.

###############################################################################
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(cowplot)
theme_set(theme_cowplot())

outDir <- "../plots/MDD_manuscript_figures"
TF_enrichmentsFile <- "../misc/CHEA3_ReMap_TFEnrichments/MDD_CHEA3_ReMap.csv"

###############################################################################
if(!grepl("R$", getwd())) {setwd("R")}


CHEA3_in <- 
  read_csv(TF_enrichmentsFile) %>% 
  filter(`Query Name` == "module_10")

###############################################################################

if (!dir.exists(outDir)) {dir.create(outDir, recursive = TRUE)}

pdf(sprintf("%s/MDD_F6C.pdf", outDir), height = 5.5, width = 5.5, 
    useDingbats = FALSE)
CHEA3_in %>% 
  filter(Rank <= 15) %>% 
  ggplot(aes(Rank, -log10(FDR))) +
  geom_point() +
  geom_text_repel(aes(label = TF)) +
  ggtitle("Module 10\nTop 15 Enriched ReMap Transcription Factors",
          subtitle = "Rank vs -log10(FDR)")
dev.off()

