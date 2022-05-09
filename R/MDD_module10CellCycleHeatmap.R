# The purpose of this script is to generate a heatmap of the genes
# in the intersection of the Module 10 RNA-seq features and the Seurat
# cell cycle gene list.
# 
# The EdU positive proportion (from IF) is used as an annotation 
# and to order samples.

###############################################################################
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

seuratGeneList    <- "../RNAseq/misc/seurat_cellCycle_geneList.txt"
moduleGeneFile    <- "../misc/MDD_multiomics14Module_allFeatures.csv"
RNA_resultsFile   <- "../RNAseq/Data/MDD_RNAseq_shrunkenResults_df.Rdata"
# combinedTableFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"
combinedTableFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable_2.csv"
dataIF            <- "../IF/Data/MDD_IF_Level4.csv"


heatmapScriptFile    <- "MDD_RNAseq_heatmapFunctions.R"
importRNAScriptFile  <- "MDD_import_RNAseq_ensembl.R"

outDir            <- "../plots/MDD_manuscript_figures"

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

source(heatmapScriptFile)
source(importRNAScriptFile)

colnames(RNAseqL4)[1] <- "CTRL_0"

if (!dir.exists(outDir)) {dir.create(outDir)}
###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

expressedGenes <- read.csv(combinedTableFile, stringsAsFactors = FALSE)

seuratList <- read.table(seuratGeneList,
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::rename(hgnc_symbol = "Gene") %>% 
  left_join(distinct(at.RNA, hgnc_symbol, ensembl_gene_id)) %>% 
  left_join(exprFilt) %>% 
  filter(expressedInRNA)

sa.RNA.L3$experimentalCondition[sa.RNA.L3$experimentalCondition == "ctrl_0"] <- "CTRL_0"
sa.RNA.L3$Ligand[sa.RNA.L3$experimentalCondition == "ctrl_0"] <- "CTRL"

# Mean-centering RNAseq Level4 data and changing identifiers to gene symbols
RNAseqL4C_s <-
  RNAseqL4[, setdiff(colnames(RNAseqL4), "CTRL_0")] %>% 
  t %>% 
  scale %>% 
  t %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(at.RNA) %>% 
  filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  dplyr::select(-ensembl_gene_id, -gene_biotype) %>% 
  column_to_rownames("hgnc_symbol")

eduDataC2 <-
  read.csv(dataIF, row.names = 1, stringsAsFactors = FALSE) %>%
  t() %>%
  data.frame %>%
  rownames_to_column("specimenName") %>%
  filter(collection == "C2") %>%
  mutate(experimentalCondition = sprintf("%s_%s",
                                         ligand, experimentalTimePoint)) %>% 
  dplyr::select(experimentalCondition, contains("EdU"))

eduDataC2$experimentalCondition[eduDataC2$experimentalCondition == "ctrl_ 0"] <- "CTRL_0"

sa.RNA.joined <-
  sa.RNA.L3 %>% 
  left_join(eduDataC2) %>% 
  arrange((Edu_Positive_Proportion))

ordL4 <- sa.RNA.joined %>% 
  filter(Ligand != "CTRL") %>%
  pull(experimentalCondition) %>% 
  unique


whichModule10 <- read_csv(moduleGeneFile) %>% 
  filter(Cluster == "module_10") %>% 
  filter(Type == "RNAseq") %>% 
  pull(feature)

whichToPlot <- seuratList %>% 
  filter(Class ==  "G1/S") %>% 
  pull(hgnc_symbol) %>% 
  intersect(whichModule10)

ordL4 <- sa.RNA.joined %>% 
  filter(Ligand != "CTRL") %>%
  pull(experimentalCondition) %>% 
  unique

col_fun <- colorRamp2(c(0, 25), c("white", "red"))

ha.RNA.L4 <-
  HeatmapAnnotation(df = sa.RNA.joined %>% 
                      distinct(Ligand, Time, .keep_all = TRUE) %>% 
                      filter(Ligand != "CTRL") %>%
                      dplyr::select(Ligand, Time, Edu_Positive_Proportion) %>% 
                      dplyr::rename("EdU Positive Proportion" = Edu_Positive_Proportion) %>% 
                      dplyr::mutate(`EdU Positive Proportion` = as.numeric(`EdU Positive Proportion`)),
                    col = c(col_MDD, list("EdU Positive Proportion" = col_fun)))

if (!dir.exists(outDir)) {dir.create(outDir, recursive = TRUE)}

pdf(sprintf("%s/MDD_F6F.pdf", outDir),
    height = 7, width = 9.5)
HeatmapL4C(RNAseqL4C_s[whichToPlot, ordL4], 
           featureName = "Mean-centered\nlog2(fpkm + 1)",
           showRow = TRUE, ca = ha.RNA.L4,
           column_title = "G1/S Genes in Module 10")
dev.off()

