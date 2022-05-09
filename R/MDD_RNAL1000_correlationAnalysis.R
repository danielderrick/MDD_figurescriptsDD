library(tidyverse)
library(ComplexHeatmap)
library(cmapR)
library(eulerr)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(circlize)

colScript            <- "MDD_importColors_pretty.R"
RNAseqFile           <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"
RNAseqAnnoFile       <- "../RNAseq/Metadata/MDD_RNAseq_geneAnnotations.csv"
RNAseqSampleAnnoFile <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"

L1000File           <- "../L1000/Data/MDD_L1000_Level3.csv"
L1000AnnoFile       <- "../L1000/Metadata/MDD_L1000_probeAnnotations.csv"
L1000SampleAnnoFile <- "../L1000/Metadata/MDD_L1000_sampleMetadata.csv"

# Setup
## Importing data

RNAseq            <- read.csv(RNAseqFile)
RNAseqAnno        <- read.csv(RNAseqAnnoFile)
RNAseqSampleAnno  <- read.csv(RNAseqSampleAnnoFile)

# Importing L1000 data matrix and metadata:
L1000            <- read.csv(L1000File)
L1000Anno        <- read.csv(L1000AnnoFile)
L1000SampleAnno  <- read.csv(L1000SampleAnnoFile)

## Importing colors
source(colScript)

## Converting data to HGNC symbols

### RNAseq
RNAseq_symbol <-
  RNAseq %>% 
  left_join(RNAseqAnno) %>% 
  # dplyr::select(hgnc_symbol, everything()) %>%
  filter(!duplicated(hgnc_symbol)) %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  column_to_rownames("hgnc_symbol")

### L1000
L1000_symbol <- 
  L1000 %>% 
  left_join(L1000Anno) %>% 
  dplyr::select(pr_gene_symbol, contains("sid")) %>% 
  arrange(pr_gene_symbol) %>% 
  dplyr::rename(hgnc_symbol = pr_gene_symbol) %>% 
  column_to_rownames("hgnc_symbol")

## Filtering data to overlapping samples and intersecting genes
### Selecting overlapping samples

which_samples <-
  L1000SampleAnno %>% 
  filter(!outlierSample) %>% 
  pull(specimenID) %>% 
  as.character() %>% 
  intersect(., colnames(RNAseq))

which_samples

### Selecting intersecting genes
filt <- apply(RNAseq_symbol, 1, anyNA)
RNAseq_symbol_f1 <- RNAseq_symbol[!filt, ]

intersect_genes <- base::intersect(rownames(L1000_symbol),
                                   rownames(RNAseq_symbol_f1))

intersect_genes %>% summary()

### Subsetting objects to overlapping genes/samples
L1000_symbol_filt  <- L1000_symbol[intersect_genes, which_samples]
RNAseq_symbol_filt <- RNAseq_symbol_f1[intersect_genes, which_samples]

## Z-scoring data
L1000_symbol_Z <- 
  L1000_symbol_filt %>% 
  t() %>% 
  scale() %>% 
  t()

RNAseq_symbol_Z <- 
  RNAseq_symbol_filt %>% 
  t() %>% 
  scale() %>% 
  t()

filt <- apply(RNAseq_symbol_Z, 1, anyNA)

RNAseq_symbol_Z_filt <- RNAseq_symbol_Z[!filt, ]
L1000_symbol_Z_filt  <- L1000_symbol_Z[!filt, ]

## Computing pearson's correlation between samples

corZ <- 
  cor((L1000_symbol_Z_filt),
      (RNAseq_symbol_Z_filt), 
      method = "pearson")

# Heatmaps
## Heatmap - between-sample correlations, unclustered {.tabset}

L1000SA_Filt <- 
  L1000SampleAnno %>% 
  filter(specimenID %in% rownames(corZ))

colnames(corZ) == L1000SA_Filt$specimenID

rA_sampleAnno1 <-
  L1000SA_Filt %>% 
  dplyr::rename(Time = experimentalTimePoint,
                Ligand = ligand) %>% 
  dplyr::select(Time, Ligand) %>% 
  mutate(Ligand = fct_recode(Ligand, 
                             "CTRL" = "ctrl",
                             "BMP2+EGF" = "BMP2",
                             "IFNG+EGF" = "IFNG",
                             "TGFB+EGF" = "TGFB")) %>% 
  rowAnnotation(df = ., col = col_MDD)

rA_sampleAnno2 <-
  L1000SA_Filt %>% 
  dplyr::rename(Time = experimentalTimePoint,
                Ligand = ligand) %>% 
  dplyr::select(Time, Ligand) %>% 
  mutate(Ligand = fct_recode(Ligand, 
                             "CTRL" = "ctrl",
                             "BMP2+EGF" = "BMP2",
                             "IFNG+EGF" = "IFNG",
                             "TGFB+EGF" = "TGFB")) %>% 
  HeatmapAnnotation(df = ., col = col_MDD, show_legend = FALSE)


outDirPlots <- "../plots/MDD_manuscript_figures"

if (!dir.exists(outDirPlots)) {dir.create(outDirPlots)}

# The plots below compare samples with matched/unmatched experimentalConditions, ligand treatments, and sampleIDs.
## Matched conditions


```{r, echo = FALSE, message = FALSE}
corMatchedTrt <- 
  corZ %>% 
  data.frame() %>% 
  rownames_to_column("sample1") %>%
  pivot_longer(-sample1, names_to = "specimenID") %>% 
  left_join(dplyr::select(RNAseqSampleAnno, specimenID, experimentalCondition)) %>% 
  dplyr::rename(sample2 = specimenID,
                trt2 = experimentalCondition) %>% 
  dplyr::rename(specimenID = sample1) %>% 
  left_join(dplyr::select(RNAseqSampleAnno, specimenID, experimentalCondition)) %>% 
  dplyr::rename(sample1 = specimenID,
                trt1 = experimentalCondition) %>% 
  mutate(matched = (trt1 == trt2))

corMatchedTrt %>% 
  ggplot(aes(value, fill = fct_rev(as.factor(matched)))) +
  # geom_density(alpha = .5, adjust = 2) +
  geom_density(alpha = .5) +
  xlab("Pearson's r") +
  scale_fill_manual("Matched\nConditions",
                    values = c("TRUE" = "steelblue2",
                               "FALSE" = "gray80")) +
  xlim(c(-1, 1)) +
  ggtitle("Pearson's R Distribution", 
          sprintf("matched (n = %s) & unmatched (n = %s) pairs",
                  nrow(filter(corMatchedTrt, matched)),
                  nrow(filter(corMatchedTrt, !matched))))

pdf(sprintf("%s/MDD_SF4D.pdf", outDirPlots), height = 6, width = 7.5)
corMatchedTrt %>% 
  ggplot(aes(value, fill = fct_rev(as.factor(matched)))) +
  # geom_density(alpha = .5, adjust = 2) +
  geom_density(alpha = .5) +
  xlab("Pearson's r") +
  scale_fill_manual("Matched\nConditions",
                    values = c("TRUE" = "steelblue2",
                               "FALSE" = "gray80")) +
  xlim(c(-1, 1)) +
  ggtitle("Pearson's R Distribution", 
          sprintf("matched (n = %s) & unmatched (n = %s) pairs",
                  nrow(filter(corMatchedTrt, matched)),
                  nrow(filter(corMatchedTrt, !matched))))
dev.off()

# Significance testing for matched conditions

# Testing for a difference in Pearson's R between matched and unmatched samples using
# a Wilcoxon rank sum test.
sig <- wilcox.test(value ~ matched, data = corMatchedTrt)
round(sig$p.value)


corMatchedTrt %>% 
  ggplot(aes(value, fill = fct_rev(as.factor(matched)))) +
  # geom_density(alpha = .5, adjust = 2) +
  geom_density(alpha = .5) +
  xlab("Pearson's r") +
  scale_fill_manual("Matched\nConditions",
                    values = c("TRUE" = "steelblue2",
                               "FALSE" = "gray80")) +
  xlim(c(-1, 1)) +
  ggtitle("Pearson's R Distribution", 
          sprintf("matched (n = %s) & unmatched (n = %s) pairs\np = %s (Mann-Whitney U Test)",
                  nrow(filter(corMatchedTrt, matched)),
                  nrow(filter(corMatchedTrt, !matched)),
                  signif(sig$p.value, 2)))
