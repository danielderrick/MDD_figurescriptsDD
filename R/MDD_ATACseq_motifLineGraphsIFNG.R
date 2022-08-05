# The purpose of this script is to generate a line graph of
# IRF and STAT family motif enrichment scores for EGF and IFNG
# treated samples.
library(tidyverse)
library(ggforce)
library(cowplot)
theme_set(theme_cowplot())

colFile     <- "../misc/MDD_color_codes.csv"
saRNAFile   <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"
motifFile   <- "../ATACseq/Data/motif/MDD_ATACseq_MotifFamilyScores.csv"
outDir      <- "../plots/MDD_manuscript_figures"
families    <- c("STAT factors", "Interferon-regulatory factors")

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

col <- read.csv(colFile, stringsAsFactors = FALSE)
col <- list(
  ligand = dplyr::slice(col, 1:8),
  experimentalTimePoint = dplyr::slice(col, 10:15),
  secondLigand = dplyr::slice(col, 17:18),
  replicate = dplyr::slice(col, 19:24),
  collection = dplyr::slice(col, 26:28)
)

col <-
  lapply(col, function(x) {
    x <- setNames(x[, 2], x[, 1])
  })

## Functions
renameTime0 <- function(data, sampleAnno, newLigand) {
  x <- 
    data %>% 
    left_join(sampleAnno) %>% 
    filter(ligand == "ctrl") %>% 
    dplyr::select(1:8) %>% 
    dplyr::select(-contains("specimen")) %>% 
    mutate(ligand = newLigand,
           experimentalCondition = paste0(newLigand, "_0"))
  x
}

joinAndFilterFeatures <- function(data, annoTable, featureName, myFilter) {
  quo_featureName <- enquo(featureName)
  x <- 
    data %>% 
    left_join(annoTable) %>% 
    filter(!!quo_featureName %in% myFilter) %>% 
    dplyr::select(!!quo_featureName, contains("sid"))
  
  x <- 
    x %>% 
    pivot_longer(cols = -!!quo_featureName,
                 names_to = "specimenID",
                 values_to = "expression")
  x
}

joinBindGroupNoSummary <- function(data, sampleAnno, ligands, featureName, bind_with1, bind_with2) {
  x <-
    data %>% 
    left_join(sampleAnno) %>% 
    filter(ligand %in% ligands) %>% 
    dplyr::select(1:8) %>% 
    dplyr::select(-contains("specimen")) %>% 
    bind_rows(bind_with1) %>% 
    bind_rows(bind_with2) %>% 
    group_by(experimentalCondition, family, ligand, experimentalTimePoint) %>%
    arrange(experimentalTimePoint) %>%
    data.frame() %>%
    arrange(experimentalTimePoint) %>%
    mutate(clear = case_when(ligand == "ctrl" ~ 1,
                             ligand != "ctrl" ~ 0))
  x
}



joinBindGroupAndSummarize <- function(data, sampleAnno, ligands, featureName, bind_with1, bind_with2) {
  x <-
    data %>% 
    left_join(sampleAnno) %>% 
    filter(ligand %in% ligands) %>% 
    dplyr::select(1:8) %>% 
    dplyr::select(-contains("specimen")) %>% 
    bind_rows(bind_with1) %>% 
    bind_rows(bind_with2) %>% 
    group_by(experimentalCondition, family, ligand, experimentalTimePoint) %>%
    arrange(experimentalTimePoint) %>%
    dplyr::summarize(meanExpr = mean(expression),
              stdevExpr = var(expression)^(1/2),
              stderrExpr = (var(expression)^(1/2))/(length(expression)^(1/2))) %>%
    data.frame() %>%
    arrange(experimentalTimePoint) %>%
    mutate(clear = case_when(ligand == "ctrl" ~ 1,
                             ligand != "ctrl" ~ 0))
  x
}

###############################################################################
# importing ATACseq motifs

if(!grepl("R$", getwd())) {setwd("R")}

sa.RNA    <- read.csv(saRNAFile, 
                      header = TRUE, stringsAsFactors = FALSE)

motifs <- read.csv(motifFile, 
                   header = TRUE, stringsAsFactors = FALSE)


motifs.L3.EGF0 <-
  motifs %>%
  filter(family %in% families) %>% 
  pivot_longer(cols = -family, 
               names_to = "specimenID",
               values_to = "expression") %>% 
  renameTime0(sa.RNA, "EGF")


motifs.L3.IFNG0 <-
  motifs %>%
  filter(family %in% families) %>% 
  pivot_longer(cols = -family, 
               names_to = "specimenID",
               values_to = "expression") %>% 
  renameTime0(sa.RNA, "IFNG")

indMotifScores <- lapply(c(1:2), function(x) {
  motifs %>% 
    filter(family %in% families) %>% 
    pivot_longer(cols = -family, 
                 names_to = "specimenID",
                 values_to = "expression") %>% 
    f2(sa.RNA,
                              c("EGF", "IFNG", "ctrl"),
                              family,
                              motifs.L3.EGF0,
                              motifs.L3.IFNG0) 
})


meanMotifScores <- lapply(c(1:2), function(x) {
  motifs %>% 
    filter(family %in% families) %>% 
    pivot_longer(cols = -family, 
                 names_to = "specimenID",
                 values_to = "expression") %>% 
    joinBindGroupAndSummarize(sa.RNA,
       c("EGF", "IFNG", "ctrl"),
       family,
       motifs.L3.EGF0,
       motifs.L3.IFNG0) 
})

pdf(sprintf("%s/MDD_SF3E.pdf", outDir), height = 5, width = 6)
a[[1]] %>% left_join(b[[1]]) %>%
  ggplot(aes(experimentalTimePoint,
             meanExpr,
             color = ligand)) +
  geom_line(aes(group = ligand)) +
  geom_point(size = .5,aes(experimentalTimePoint, expression), color = "black") +
  geom_errorbar(aes(min = meanExpr - stdevExpr,
                    max = meanExpr + stdevExpr),
                width = 2) +
  # geom_point(aes(alpha = clear)) +
  facet_wrap_paginate(~family, scales = "free", nrow = 1, ncol = 1, page = 1) +
  scale_color_manual(values = col$ligand) +
  ylab("Motif Enrichment") +
  ggtitle("ATACseq Motif Enrichment") +
  scale_x_continuous("time", breaks = c(0, 24, 48), labels = c("0", "24", "48"))
a[[1]] %>% left_join(b[[1]]) %>%
  # filter(family == "STAT factors") %>% 
  ggplot(aes(experimentalTimePoint,
             meanExpr,
             color = ligand)) +
  geom_line(aes(group = ligand)) +
  geom_point(aes(experimentalTimePoint, expression), color = "black") +
  geom_errorbar(aes(min = meanExpr - stdevExpr,
                    max = meanExpr + stdevExpr),
                width = 2) +
  facet_wrap_paginate(~family, scales = "free", nrow = 1, ncol = 1, page = 2) +
  scale_color_manual(values = col$ligand) +
  ylab("Motif Enrichment") +
  ggtitle("ATACseq Motif Enrichment") +
  scale_x_continuous("time", breaks = c(0, 24, 48), labels = c("0", "24", "48"))
dev.off()
