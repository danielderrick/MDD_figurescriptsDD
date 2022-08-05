# The purpose of this script is to examine the expression of genes from
# integrated MDD module 10 in a dataset of CDKi-treated cell lines
# from Hafner, et al.

library(tidyverse)
library(biomaRt)
library(cowplot)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(Hmisc)
library(RColorBrewer)

grmetricsFile    <- "../misc/TableS6_GRmetrics_panel.tsv"
geneListFile     <- "../../MDD/misc/MDD_multiomics14Module_allFeatures.csv"
sorgerFile       <- "../misc/GSE99116_RPKM_acute.csv"

# mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                 dataset = "hsapiens_gene_ensembl")
#                 # host = "uswest.ensembl.org")
# 
# 
# annoTable <- getBM(attributes = c("ensembl_gene_id",
#                                   "hgnc_symbol"),
#                    mart = mart)

# Reading in BiomaRt annotation table from 03/04/2021

annoTable <- read_csv(
  c("../misc/symbol_id_table_03042021.csv")
  )

colClassifications <- list("Classification.3class" = c(
  "Luminal" = "steelblue",
  "Basal"   = "darkgreen",
  "Claudin-low" = "firebrick",
  "Non-malignant" = "purple",
  "none" = "gray70"
  ),
  "dt.median" = colorRamp2(c(400, 100, 0), c("white", "yellow", "red")))

outDir <- "../plots/MDD_manuscript_figures"

if (!dir.exists(outDir)) {dir.create(outDir, recursive = TRUE)}

###############################################################################

## Importing Sorger RNA-seq data

# Importing Sorger RNAseq data, in units of log2(rpkm + 1), and
# creating sample metadata based on column names:

sorgerRPKM <- read.csv(sorgerFile, 
                       stringsAsFactors = FALSE) %>% 
  dplyr::rename(ensembl_gene_id = X) %>% 
  left_join(annoTable, by = "ensembl_gene_id") %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  filter(hgnc_symbol != "") %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  dplyr::select(hgnc_symbol, everything())

sampleAnno_sorger <-
  data.frame(sample = colnames(sorgerRPKM)[-1],
             stringsAsFactors = FALSE) %>% 
  mutate(CellLine  = str_split_fixed(sample, "_", 4)[, 1],
         Treatment = str_split_fixed(sample, "_", 4)[, 2],
         Dose      = str_split_fixed(sample, "_", 4)[, 3],
         Time      = str_split_fixed(sample, "_", 5)[, 4],
  ) %>% 
  mutate(Treatment = fct_relevel(Treatment, "Control"),
         Time      = as.numeric(Time),
         Dose      = as.numeric(Dose)) %>% 
  arrange(CellLine, Time, Treatment) %>% 
  mutate(sample = fct_inorder(as.factor(sample)))

sorgerColors <- list("Treatment" = setNames(brewer.pal(4, "Set2"),
                                            unique(sampleAnno_sorger$Treatment)),
                     "Dose"      = setNames(brewer.pal(4, "Reds"),
                                            sort(as.numeric(unique(sampleAnno_sorger$Dose)))),
                     "Time"      = setNames(brewer.pal(3, "Blues")[-2],
                                            sort(as.numeric(unique(sampleAnno_sorger$Time)))),
                     "CellLine"  = setNames(brewer.pal(7, "Set1"),
                                            unique(sampleAnno_sorger$CellLine))
                     )

# The 50 most variant genes are plotted below.

sorger_geneExpression_metadata <-
  data.frame(hgnc_symbol = sorgerRPKM$hgnc_symbol,
             expressed = TRUE,
             stringsAsFactors = FALSE)

mostVar50 <-
  sorgerRPKM %>% 
  column_to_rownames("hgnc_symbol") %>% 
  apply(1, function(x) {
    X <- var(x)^(1/2)
  }) %>% 
  sort(decreasing = TRUE) %>% 
  head(50) %>% 
  names

sorgerRPKM_TP <-
  column_to_rownames(sorgerRPKM, "hgnc_symbol") %>% 
  t() %>% 
  scale(scale = FALSE) %>% 
  t() %>% 
  .[, as.character(sampleAnno_sorger$sample)]

ha.Sorger <- HeatmapAnnotation(df = dplyr::select(sampleAnno_sorger, -c(sample)),
                               col = sorgerColors)


Heatmap(sorgerRPKM_TP[mostVar50, ], 
        bottom_annotation = ha.Sorger,
        name = "mean-centered\nlog2(rpkm + 1)",
        column_names_gp = gpar(fontsize = 9),
        show_column_names = FALSE,
        cluster_columns = FALSE,
        column_title = "sorger RNAseq Data\n50 Most Variant Genes by variance")

# **Sample HS578T_Control_0_24 is an outlier sample and is removed.**

sampleAnno_sorger <-
  sampleAnno_sorger %>% 
  # filter(Dose == 0.0 |
  # Dose == 3.0) %>% 
  filter(Time == 24) %>%
  filter(CellLine != "HS578T")

samples_order <-
  sampleAnno_sorger %>% 
  mutate(Time = as.numeric(Time),
         Dose = as.numeric(Dose)) %>% 
  pull(sample) %>% 
  as.character

### Removing outlier sample
sorgerRPKM <- sorgerRPKM[, c("hgnc_symbol", samples_order)]

sorgerRPKM_TP <-
  column_to_rownames(sorgerRPKM, "hgnc_symbol") %>% 
  t() %>% 
  scale(scale = FALSE) %>% 
  t()

## Importing MDD module RNAseq features
modulesDF <- 
  read.csv(geneListFile, stringsAsFactors = FALSE) %>% 
  filter(Type == "RNAseq") %>% 
  dplyr::select(Cluster, feature) %>% 
  dplyr::rename(hgnc_symbol = feature) %>% 
  mutate(Cluster = as.factor(Cluster)) %>% 
  mutate(Cluster = fct_inorder(Cluster)) %>% 
  left_join(sorger_geneExpression_metadata, by = "hgnc_symbol")  %>% 
  mutate(expressed = replace_na(expressed, FALSE))

modulesDF %>% 
  group_by(Cluster, expressed) %>% 
  dplyr::summarize(n_genes = n()) %>% 
  ggplot(aes(Cluster, n_genes, fill = expressed)) +
  geom_col() +
  ylab("Module") +
  xlab("# Genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("RNAseq features in MDD Modules") +
  scale_fill_manual("Expressed\nin sorger",
                    values = c("TRUE" = "gray80",
                               "FALSE" = "black"))
modulesDF %>% 
  group_by(Cluster, expressed) %>% 
  dplyr::summarize(n_genes = n()) %>% 
  arrange(Cluster, fct_rev(as.factor(expressed))) %>% 
  pivot_wider(names_from = expressed,
              values_from = n_genes) %>% 
  mutate(pct_expressed = 100*(`TRUE`/(`TRUE` + `FALSE`))) %>% 
  mutate(pct_expressed = replace_na(pct_expressed, 100)) %>% 
  mutate(pct_expressed = round(pct_expressed, 1))

# The majority of module genes are found in the Sorger data. 

mDF_means <-
  modulesDF %>% 
  filter(expressed) %>% 
  dplyr::select(-expressed) %>% 
  left_join(sorgerRPKM, by = "hgnc_symbol") %>% 
  pivot_longer(-c(Cluster, hgnc_symbol),
               names_to = "sample",
               values_to = "expr") %>% 
  group_by(Cluster, sample) %>% 
  dplyr::summarize(meanExpr = mean(expr, na.rm = TRUE))

grMetricsDF <- read_tsv(grmetricsFile)

## Boxplots - mean module 10 expression
### Ordered by Abemaciclib GR50
module10_meansBoxplot_abemaOrder <-
  mDF_means %>% 
  filter(Cluster == "module_10") %>% 
  left_join(sampleAnno_sorger) %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(Treatment = fct_relevel(Treatment, "Control")) %>% 
  left_join(grMetricsDF, by = c("CellLine" = "Cell line")) %>%
  arrange(`Abemaciclib GR50`) %>%
  mutate(CellLine = fct_inorder(as.factor(CellLine))) %>%
  ggplot(aes(CellLine, meanExpr, fill = Treatment)) +
  geom_boxplot() +
  ggtitle("Mean Module 10 Gene Expression",
          subtitle = "ordered by Abemaciclib GR50") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1))

module10_meansBoxplot_abemaOrder

mDF_means %>% pull(sample) %>% unique %>% sort

pdf(sprintf("%s/MDD_F6G.pdf", outDir), 
    height = 7, width = 7)
module10_meansBoxplot_abemaOrder
dev.off()

# Mann-Whitney U test
# Mean Module 10 expression in CDKi vs ctrl, 24H only
# Does mean module 10 gene expression differ between control and cdki?

tempMeans <-
  mDF_means %>% 
  filter(Cluster == "module_10") %>% 
  left_join(sampleAnno_sorger) %>% 
  filter(Time == 24)

meansControl <- 
  tempMeans %>% 
  filter(Treatment == "Control") %>% 
  pull(meanExpr)

meansCDKi <- 
  tempMeans %>% 
  filter(Treatment != "Control") %>% 
  pull(meanExpr)

wilcox.test(meansControl,
            meansCDKi)

# Yes; a Welch's t-test (unpaired, unequal variances) shows a significant 
# difference between the module 10 gene expression means.

