# The purpose of this script is to examine the expression of genes from 
# the integrated MDD modules in the BCCL RNA-seq data set.
library(tidyverse)
library(biomaRt)
library(cowplot)
library(Hmisc)
library(RColorBrewer)

dtFile         <- "../misc/JWGray_doubling_time_all_v2.txt"
bcclFile       <- "../misc/JWGray_BCCL_rnaseq_matrix_v3_tatlowVST.csv"
bcclAnnotFile  <- "../misc/JWGray_BCCL_classifications_v5.txt"

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

annoTable <- getBM(attributes = c("ensembl_gene_id",
                                  "hgnc_symbol"),
                   mart = mart)

colClassifications <- list("Classification" = c(
                                                  "Luminal" = "steelblue",
                                                  "Basal"   = "darkgreen",
                                                  "Claudin-low" = "firebrick",
                                                  "Non-malignant" = "purple",
                                                  "none" = "gray70"),
                          "dt.median" = colorRamp2(c(400, 100, 0), 
                                                   c("white", "yellow", "red")))

outDir <- "../plots/MDD_manuscript_figures"

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

# Importing BCCL RNAseq data in counts 
# These data have been VST-normalizing using DESeq2.
# The 50 most variant genes (by CV) are plotted below.

bcclVST <- read.csv(bcclFile, 
                    stringsAsFactors = FALSE) %>% 
  left_join(annoTable, by = "ensembl_gene_id") %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  filter(hgnc_symbol != "") %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  dplyr::select(hgnc_symbol, everything())

bccl_geneExpression_metadata <-
  data.frame(hgnc_symbol = bcclVST$hgnc_symbol,
             expressed = TRUE,
             stringsAsFactors = FALSE)

# Importing BCCL doubling times, and annotating with classifications. 
# Some cell lines lack classifications.

dt <- 
  read_tsv(dtFile)

dtAnno <- read_tsv(bcclAnnotFile) %>% 
  dplyr::rename(CellLineName = `Cell Line`)

dt <- 
  dt %>%
  left_join(dtAnno, by = "CellLineName") %>% 
  filter(!grepl("1%", CellLineName)) %>% 
  mutate(startsWithNumber = (grepl("^[0-9]", CellLineName))) %>% 
  mutate(CellLineName = case_when(startsWithNumber ~ sprintf("X_%s", CellLineName),
                                  !startsWithNumber ~ CellLineName)) %>% 
  dplyr::select(-startsWithNumber) %>% 
  mutate(CellLineName = toupper(CellLineName))

# Fixing differences between cell line names in datasets and comparing
# cell lines.

dt <- 
  dt %>%
  mutate(CellLineName = str_remove(CellLineName, "1%")) %>% 
  mutate(startsWithNumber = (grepl("^[0-9]", CellLineName))) %>% 
  mutate(CellLineName = case_when(startsWithNumber ~ sprintf("X_%s", CellLineName),
                                  !startsWithNumber ~ CellLineName)) %>% 
  dplyr::select(-startsWithNumber) %>% 
  mutate(CellLineName = toupper(CellLineName))

rnaNames <- 
  bcclVST %>% 
  colnames %>% 
  .[-1] %>% 
  str_split_fixed("_", 3) %>% 
  .[, 1] %>% 
  str_replace("^[X]", "X_") %>% 
  toupper()

# These cell line names are are present in the BCCL cell line data, 
# but are absent in the doubling time data:
setdiff(rnaNames, dt$CellLineName)

# Now looking at names in the DT data that are missing from RNA
setdiff(dt$CellLineName,
        rnaNames)

# Changing cell line name `SUM1315` in doubling time data
# to `SUM1315MO2`, to match RNAseq data.
# Other  cell lines appear to be missing from the 
# doubling time data and will be removed: 

if(rnaNames[55] == "SUM1315") {rnaNames[55] <- "SUM1315MO2"}

colnames(bcclVST)[-1] <- rnaNames

# Filtering both RNAseq data and doubling time data
# to select only overlapping cell lines.

overlapCellLines <- base::intersect(colnames(bcclVST), 
                                    dt$CellLineName) %>% 
  sort

bcclVST <- bcclVST[, c("hgnc_symbol", overlapCellLines)]
dt <- dt[match(overlapCellLines, dt$CellLineName), ]

dt.filt <- dt %>%
  dplyr::rename(Classification = `Classification-3class`)

bcclVST <- bcclVST[, c("hgnc_symbol", dt.filt$CellLineName)]

###############################################################################
### Importing MDD module RNAseq features
# Importing RNAseq features from 14 module analysis and filtering to
# genes that are found in the BCCL RNA-seq data.

geneListFile    <- "../misc/MDD_multiomics14Module_allFeatures.csv"

modulesDF <- 
  read.csv(geneListFile, stringsAsFactors = FALSE) %>% 
  filter(Type == "RNAseq") %>% 
  dplyr::select(Cluster, feature) %>% 
  dplyr::rename(hgnc_symbol = feature) %>% 
  mutate(Cluster = as.factor(Cluster)) %>% 
  mutate(Cluster = fct_inorder(Cluster)) %>% 
  left_join(bccl_geneExpression_metadata, by = "hgnc_symbol")  %>% 
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
  scale_fill_manual("Expressed\nin BCCL",
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

# The majority of genes are found in the BCCL data. 

mDF_means <-
  modulesDF %>% 
  filter(expressed) %>% 
  dplyr::select(-expressed) %>% 
  left_join(bcclVST, by = "hgnc_symbol") %>% 
  pivot_longer(-c(Cluster, hgnc_symbol),
               names_to = "CellLineName",
               values_to = "expr") %>% 
  group_by(Cluster, CellLineName) %>% 
  dplyr::summarize(meanExpr = mean(expr, na.rm = TRUE))

dt_subset <- 
  dt.filt %>% 
  dplyr::select(dt.median, Classification) %>% 
  data.frame

rA <-
  rowAnnotation(df = dt_subset,
                col = colClassifications)

heatmap_meanOrdered <-
  mDF_means %>% 
  pivot_wider(names_from = Cluster, values_from = meanExpr) %>% 
  column_to_rownames("CellLineName") %>% 
  scale(scale = FALSE) %>% 
  Heatmap(cluster_columns = FALSE,
          name = "mean\nexpression",
          column_title = "Mean Expression of Genes in MDD Modules\nOrdered",
          left_annotation = rA)

heatmap_meanOrdered


if (!dir.exists(outDir)) {dir.create(outDir, recursive = TRUE)}


module10Expr <- 
  modulesDF %>% 
  filter(expressed) %>% 
  filter(Cluster == "module_10") %>% 
  dplyr::select(-expressed, -Cluster) %>% 
  left_join(bcclVST, by = "hgnc_symbol") %>% 
  column_to_rownames("hgnc_symbol") %>% 
  data.matrix

moduleMeansWide <-
  mDF_means %>% 
  pivot_wider(names_from = "Cluster", values_from = "meanExpr")

if(all(dt.filt$CellLineName == moduleMeansWide$CellLineName)) {
   moduleMeansCorrS <- rcorr((dt.filt$dt.median), data.matrix(moduleMeansWide[, -1]),
                                 type = "spearman")
}

mMeanCorrS <-
  data.frame(
    Cluster = names(moduleMeansCorrS$r[1, -1]),
    r = moduleMeansCorrS$r[1, -1],
    p = p.adjust(moduleMeansCorrS$P[1, -1], method = "fdr")
  )

mMeanCorrS %>% filter(Cluster == "module_10")

theme_set(theme_cowplot())

meanScatterplot10Log10 <-
    mDF_means %>% 
    left_join(dt.filt) %>% 
    left_join(mMeanCorrS) %>% 
    ungroup %>% 
    filter(Cluster == "module_10") %>% 
    mutate(Cluster = str_replace(Cluster, "^m", "M")) %>% 
    mutate(Cluster = str_replace(Cluster, "_", " ")) %>% 
    mutate(Cluster = sprintf("%s\nPearson's r = %s, q = %s", Cluster, round(r, 3), round(p, 5))) %>% 
    mutate(Cluster = fct_inorder(as.factor(Cluster))) %>% 
    ggplot(aes(meanExpr, dt.median)) +
    # geom_text(aes(label = r)) +
    geom_smooth(method = "lm", linetype = 3, color = "black") +
    geom_point(size = 2.75, alpha = .75, aes(color = Classification)) +
    # facet_wrap(~Cluster, scales = "free_x", ncol = 2) +
    # ggtitle("BCCL doubling time vs mean MDD module expression") +
    scale_color_manual(values = colClassifications$Classification, 
                       na.value = "gray50") +
    # xlab("mean gene expression\n(VST)") +
    # ylab("median doubling time\n(hours)")
    xlab("Mean Gene Expression\n(VST)") +
    ylab("Median Doubling Time\n(hours)") +
    scale_y_log10()

pdf(sprintf("%s/MDD_F6H.pdf", outDir), height = 5, width = 7)
meanScatterplot10Log10
dev.off()
