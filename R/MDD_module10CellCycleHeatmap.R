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

RNAseqMetadataFile <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"
geneMetadataFile <- "../RNAseq/Metadata/MDD_RNAseq_geneAnnotations.csv"
RNAseqL3DataFile   <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"
RNAseqL4DataFile   <- "../RNAseq/Data/MDD_RNAseq_Level4.csv"

seuratGeneList    <- "../RNAseq/misc/seurat_cellCycle_geneList.txt"
moduleGeneFile    <- "../misc/MDD_multiomics14Modules_TFannotations_module10.csv"
dataIF            <- "../IF/Data/MDD_IF_Level4.csv"
colFile           <- "../misc/MDD_color_codes.csv"

outDir            <- "../plots/MDD_manuscript_figures"

###############################################################################

HeatmapL4C <- function(mat, ca = ha.RNA.L4, featureName = "Z-score", bks = c(-3, 0, 3), showRow = FALSE, clust = FALSE, ...) {
  Heatmap(mat, 
          top_annotation = ca,
          cluster_columns = clust,
          show_row_names = showRow,
          cluster_row_slices = FALSE,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          name = featureName,
          col = colorRamp2(bks, c("blue", "white", "red")),
          ...)
}
###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

###############################################################################

col_raw <- read.csv(colFile,
                    stringsAsFactors = FALSE)
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

col

###############################################################################

moduleGenes <- read.csv(moduleGeneFile, stringsAsFactors = FALSE)

seuratList <- read.table(seuratGeneList,
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::rename(hgnc_symbol = "Gene") %>% 
  filter(Class == "G1/S")

seuratG1S_module10 <- base::intersect(moduleGenes$feature,
                                      seuratList$hgnc_symbol)

###############################################################################

at.RNA <- read_csv(geneMetadataFile)

sa.RNA.L3 <-
  read.csv(RNAseqMetadataFile, 
           stringsAsFactors = FALSE) %>% 
  mutate(ligand = fct_inorder(ligand),
         experimentalTimePoint = fct_inorder(as.factor(experimentalTimePoint)),
         experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>% 
  filter(RNAseq_QCpass) %>% 
  dplyr::select(-contains("RNA")) %>% 
  dplyr::select(-contains("sequencing")) %>% 
  dplyr::select(-contains("File")) %>% 
  dplyr::rename(Time = experimentalTimePoint,
                Ligand = ligand) %>% 
  arrange(experimentalCondition) %>% 
  mutate(experimentalCondition = as.character(experimentalCondition)) %>% 
  mutate(Ligand = fct_recode(Ligand, 
                             "CTRL" = "ctrl",
                             "BMP2+EGF" = "BMP2",
                             "IFNG+EGF" = "IFNG",
                             "TGFB+EGF" = "TGFB"),
         experimentalCondition = fct_recode(experimentalCondition,
                                            "CTRL_0" = "ctrl_0")) %>% 
  mutate(Ligand = fct_relevel(Ligand, "CTRL", "PBS",
                              "HGF", "OSM", "EGF",
                              "BMP2+EGF", "IFNG+EGF", "TGFB+EGF")) %>% 
  arrange(Ligand, Time) %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition))

sa.RNA.L3$experimentalCondition[sa.RNA.L3$experimentalCondition == "ctrl_0"] <- "CTRL_0"
sa.RNA.L3$Ligand[sa.RNA.L3$experimentalCondition == "ctrl_0"] <- "CTRL"

RNAseqL4 <- 
  read.csv(RNAseqL4DataFile,
           header = TRUE, row.names = 1) %>% 
  data.matrix
colnames(RNAseqL4)[1] <- "CTRL_0"

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
  dplyr::select(-ensembl_gene_id) %>% 
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

names(col)[1:2] <- c("Ligand", "Time")
names(col$ligand) <- levels(sa.RNA.joined$Ligand)

col_fun <- colorRamp2(c(0, 25), c("white", "red"))

ha.RNA.L4 <-
  HeatmapAnnotation(df = sa.RNA.joined %>% 
                      distinct(Ligand, Time, .keep_all = TRUE) %>% 
                      filter(Ligand != "CTRL") %>%
                      dplyr::select(Ligand, Time, Edu_Positive_Proportion) %>% 
                      dplyr::rename("EdU Positive Proportion" = Edu_Positive_Proportion) %>% 
                      dplyr::mutate(`EdU Positive Proportion` = as.numeric(`EdU Positive Proportion`)),
                    col = c(col, list("EdU Positive Proportion" = col_fun)))

if (!dir.exists(outDir)) {dir.create(outDir, recursive = TRUE)}

pdf(sprintf("%s/MDD_F6F.pdf", outDir),
    height = 7, width = 9.5)
HeatmapL4C(RNAseqL4C_s[seuratG1S_module10, ordL4], 
           featureName = "Mean-centered\nlog2(fpkm + 1)",
           showRow = TRUE, ca = ha.RNA.L4,
           column_title = "G1/S Genes in Module 10")
dev.off()
