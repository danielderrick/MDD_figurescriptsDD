library(ComplexHeatmap)
library(tidyverse)
library(DESeq2)
library(plotly)
library(cowplot)
library(DT)
library(biomaRt)
library(circlize)

mSfile <- "../ATACseq/Data/motif/MDD_ATACseq_MotifScores.csv"
mAfile <- "../ATACseq/Data/motif/MDD_ATACseq_MotifAnno.csv"
annoFile  <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
ssGSEAFile  <- "../misc/MDD_RNAseq_ssGSEAScores.csv"

RNANormfile <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"

colFile <- "../misc/MDD_color_codes.csv"

col_raw <- read.csv(colFile, stringsAsFactors = FALSE)
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

sampleAnno <- read.csv(annoFile, stringsAsFactors = FALSE) 
sampleAnno

motifMat   <- read.csv(mSfile, stringsAsFactors = FALSE)

motifAnno  <- read.csv(mAfile, stringsAsFactors = FALSE) %>% 
  mutate(family = str_remove_all(family, "[ ]*$"))

ssGSEA <- 
  read.csv(ssGSEAFile, header = TRUE, 
           stringsAsFactors = FALSE, row.names = 1) %>% 
  t() %>% 
  scale(scale = FALSE) %>%
  t() %>% 
  data.frame %>%
  rownames_to_column("hgnc_symbol")

motifAnno <-
  motifAnno %>% 
  mutate(ssGSEA_exists = hgnc_symbol %in% ssGSEA$hgnc_symbol)

mmTP <- column_to_rownames(motifMat, "motif")

colnames(mmTP) <-
  sampleAnno %>% 
  filter(specimenID %in% colnames(mmTP)) %>% 
  pull(specimenName)

Heatmap(mmTP, 
        show_row_names = FALSE,
        cluster_columns = FALSE)

motifVardf <- 
  motifMat %>% 
  column_to_rownames("motif") %>% 
  apply(., 1, function(x) {x <- var(x)^(1/2)}) %>% 
  data.frame(motif = names(.),
             stdev = .,
             stringsAsFactors = FALSE)

ggplot(motifVardf, aes(stdev)) +
  geom_histogram(binwidth = .5) +
  ggtitle("Distribution of motif standard deviations")

mostVarMotifs <-
  motifVardf %>% 
  left_join(motifAnno) %>% 
  group_by(family) %>% 
  arrange(desc(stdev)) %>% 
  filter(ssGSEA_exists) %>% 
  distinct(family, .keep_all = TRUE) %>% 
  ungroup %>% 
  mutate(rankStDev = order(stdev, decreasing = TRUE)) %>% 
  filter(rankStDev <= 10)


mvmAll <-
  mostVarMotifs %>% 
  distinct(motif, .keep_all =TRUE) %>% 
  left_join(motifMat)

mvmAllToPlot <-
  mvmAll %>% 
  dplyr::select(motif, contains("sid")) %>% 
  column_to_rownames("motif")

colnames(mvmAllToPlot) <-
  sampleAnno %>% 
  filter(ATACseq_QCpass) %>% 
  pull(specimenName)


mAll <-
  motifVardf %>% 
  left_join(motifAnno) %>% 
  distinct(motif, .keep_all =TRUE) %>% 
  left_join(motifMat)

mAllToPlot <-
  mvmAll %>% 
  dplyr::select(motif, contains("sid")) %>% 
  column_to_rownames("motif")

colnames(mAllToPlot) <-
  sampleAnno %>% 
  filter(ATACseq_QCpass) %>% 
  pull(specimenName)

mvmByFamily <-
  mostVarMotifs %>% 
  filter(ssGSEA_exists) %>% 
  group_by(family) %>% 
  filter(stdev == max(stdev)) %>% 
  ungroup %>% 
  distinct(motif, .keep_all =TRUE) %>% 
  left_join(motifMat)

mvmToPlot <-
  mvmByFamily %>% 
  dplyr::select(hgnc_symbol, contains("sid")) %>% 
  column_to_rownames("hgnc_symbol")

colnames(mvmToPlot) <-
  sampleAnno %>% 
  filter(ATACseq_QCpass) %>% 
  pull(specimenName)

mAllByFamily <-
  motifVardf %>% 
  left_join(motifAnno) %>%   
  filter(ssGSEA_exists) %>% 
  group_by(family) %>% 
  filter(stdev == max(stdev)) %>% 
  ungroup %>% 
  distinct(motif, .keep_all =TRUE) %>% 
  left_join(motifMat)

mAllToPlot <-
  mAllByFamily %>% 
  dplyr::select(hgnc_symbol, contains("sid")) %>% 
  column_to_rownames("hgnc_symbol")

colnames(mAllToPlot) <-
  sampleAnno %>% 
  filter(ATACseq_QCpass) %>% 
  pull(specimenName)

mAllPeriod <-
  motifVardf %>% 
  left_join(motifAnno) %>%   
  filter(ssGSEA_exists) %>% 
  arrange(desc(stdev)) %>% 
  distinct(hgnc_symbol, .keep_all =TRUE) %>% 
  left_join(motifMat)

mAllPeriodToPlot <-
  mAllPeriod %>% 
  dplyr::select(hgnc_symbol, contains("sid")) %>% 
  column_to_rownames("hgnc_symbol")

colnames(mAllPeriodToPlot) <-
  sampleAnno %>% 
  filter(ATACseq_QCpass) %>% 
  pull(specimenName)

## Scatterplots: most variant motifs vs. ssGSEA scores.

The most variant motifs, filtered by family, are joined with ssGSEA scores and plotted
in the scatterplots below.

mvmByFamily[mvmByFamily$hgnc_symbol == "MAF", ]$hgnc_symbol <- "NFE2"

mvmToJoin <-
  mvmByFamily %>% 
  dplyr::select(hgnc_symbol, family, contains("sid")) %>% 
  pivot_longer(-c(hgnc_symbol, family), names_to = "specimenID", values_to = "motifEnrichment")

ssGSEAtoJoin <-
  ssGSEA %>% 
  pivot_longer(-hgnc_symbol, names_to = "specimenID", values_to = "ssGSEA")

joined <-
  inner_join(mvmToJoin, ssGSEAtoJoin) %>% 
  left_join(sampleAnno)

pdf("../plots/MDD_manuscript_figures/MDD_SF4F_a.pdf", height = 12, width = 8)
ggplot(joined, aes(motifEnrichment, ssGSEA, color = ligand, shape = as.factor(experimentalTimePoint))) +
  geom_point(size = 4) +
  facet_wrap(~hgnc_symbol, scale = "free", ncol = 2) +
  scale_shape("Time") +
  xlab("Motif Enrichment Score\n(ATACseq)") +
  ylab("ssGSEA Enrichment Score\n(RNAseq)") +
  scale_color_manual("Ligand", values = col$ligand)
dev.off()

pdf("../plots/MDD_manuscript_figures/MDD_SF4F_b.pdf", height = 6, width = 15)
ggplot(joined, aes(motifEnrichment, ssGSEA, color = ligand, shape = as.factor(experimentalTimePoint))) +
  geom_point(size = 3.5) +
  facet_wrap(~hgnc_symbol, scale = "free", ncol = 5) +
  scale_shape("Time") +
  xlab("Motif Enrichment Score\n(ATACseq)") +
  ylab("ssGSEA Enrichment Score\n(RNAseq)") +
  scale_color_manual("Ligand", values = col$ligand)
dev.off()
