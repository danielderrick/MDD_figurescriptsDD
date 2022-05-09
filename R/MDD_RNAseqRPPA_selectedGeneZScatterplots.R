library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(DT)
library(ggforce)
library(biomaRt)

saFile   <- "../Metadata/MDD_sample_annotations.csv"
atFile   <- "../RPPA/Metadata/MDD_RPPA_antibodyAnnotations.csv"
color_annoFile <- "../misc/MDD_color_codes.csv"

# Importing RPPA normLog2 data
RPPA_normLog2 <- read.csv("../RPPA/Data/MDD_RPPA_Level3.csv",
                          header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE) %>% 
  data.matrix()

sampleAnno <-  read.csv(saFile, 
                        stringsAsFactors = FALSE) %>% 
  filter(specimenID %in% colnames(RPPA_normLog2)) %>% 
  dplyr::select(specimenID, specimenName, experimentalCondition,
                ligand, experimentalTimePoint, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>% 
  dplyr::rename(Time = "experimentalTimePoint",
                Ligand = "ligand") %>% 
  mutate(Ligand = str_replace(Ligand, "ctrl", "CTRL"))

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "uswest.ensembl.org")

aTrna <-
  getBM(attributes = c("ensembl_gene_id",
                       "hgnc_symbol"),
        mart = mart)

col <- read.csv(color_annoFile, stringsAsFactors = FALSE)
col <- list(
  Ligand = dplyr::slice(col, 1:8),
  Time = dplyr::slice(col, 10:15),
  secondLigand = dplyr::slice(col, 17:18),
  collection = dplyr::slice(col, 26:28)
)

col$Ligand[1, 1] <- "CTRL"
col <-
  lapply(col, function(x) {
    x <- setNames(x[, 2], x[, 1])
  })

## Setup
### Importing RNAseq data

RNAseqL3 <- read.csv("../RNAseq/Data/MDD_RNAseq_Level3.csv", 
                     row.names = 1, header = TRUE) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(aTrna, by = "ensembl_gene_id") %>% 
  filter(hgnc_symbol != "",
         !is.na(hgnc_symbol)) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  column_to_rownames("hgnc_symbol")

### Importing and filtering RPPA data
aTrppa <- 
  read.csv(atFile, header = TRUE, stringsAsFactors = FALSE) %>% 
  mutate(antibody = str_remove(MDD, "-[:alnum:]-[:alnum:]$")) %>% 
  dplyr::rename(hgnc_symbol = Symbols) %>% 
  dplyr::select(antibody, hgnc_symbol, Sites)

# First, need to remove phospho-antibodies:

aTrppa <- 
  aTrppa %>% 
  filter(is.na(Sites))

# This filtered the list of 295 down to 231. 

# Antibodies with multiple target proteins also need to be removed (e.g. Akt antibody, which targets AKT1, AKT2, AKT3)

aTrppa <-
  aTrppa %>% 
  mutate(symbolFilt = !str_detect(hgnc_symbol, " ")) %>% 
  filter(symbolFilt)


any(duplicated(aTrppa$hgnc_symbol))

# The list of 231 is now reduced to 225, and the hgnc_symbols in the list are now unique

### Identifying out unexpressed genes/proteins

The aTrppa data frame is now used to select features for comparison between RNAseq and RPPA datasets.

```{r filterDatasets}
RNAseqL3filt <- RNAseqL3[aTrppa$hgnc_symbol, ]
RPPAfilt <- RPPA_normLog2[aTrppa$antibody, ]
rownames(RPPAfilt) <- aTrppa$hgnc_symbol

### Filtering datasets to selected samples
selectedIDs <- base::intersect(colnames(RNAseqL3), colnames(RPPA_normLog2))
RNAseqL3filt <- RNAseqL3filt[, selectedIDs]
RPPAfilt <- RPPAfilt[, selectedIDs]

Heatmap(RPPAfilt, name = "RPPA\nz-score", show_row_names = FALSE, cluster_columns = FALSE) + 
  Heatmap(t(scale(t(RNAseqL3filt), scale = FALSE)), 
          name = "RNAseq\nz-score", 
          cluster_columns = FALSE, show_row_names = FALSE)

# Both datasets are z-scored (independently).
RNA_Z  <- t(scale(t(RNAseqL3filt))) %>% data.frame 
RPPA_Z <- t(scale(t(RPPAfilt))) %>% data.frame

### Reshaping data and merging
# The RNAseq and RPPA data are reshaped to a "long" format, and the two are merged.

RNA_toJoin <- 
  RNA_Z %>% 
  rownames_to_column("hgnc_symbol") %>% 
  pivot_longer(-hgnc_symbol, names_to = "specimenID", values_to = "RNA_Z")

RPPA_toJoin <- 
  RPPA_Z %>% 
  rownames_to_column("hgnc_symbol") %>% 
  pivot_longer(-hgnc_symbol, names_to = "specimenID", values_to = "RPPA_Z")

joined <- 
  RNA_toJoin %>% 
  left_join(RPPA_toJoin) %>% 
  left_join(sampleAnno) %>% 
  mutate(Ligand = fct_inorder(as_factor(Ligand)))

### Scatterplots: RNAseq Z-score vs RPPA Z-score
# RNA-seq Z-score is plotted against the RPPA Z-score.

ggplot(joined, aes(RNA_Z, RPPA_Z)) + 
  geom_point(size = 1, alpha = .25) +
  geom_smooth(method = "lm") +
  facet_wrap(~experimentalCondition)

symbols  <- c("EGFR", "CCNB1", "JAK2", "DUSP4", "SOD2", "MYC", "HES1", "WWTR1", "CD274", "RPS6", "IRS1", "RB1")
library(cowplot)
theme_set(theme_cowplot())

outDir <- "../plots/MDD_manuscript_figures"

if (!dir.exists(outDir)) {
  dir.create(outDir)
}

plots <- list()
for (i in 1:length(symbols)) {
  plots[[i]] <-
    joined %>% 
    mutate(ligand = fct_inorder(as_factor(Ligand))) %>% 
    filter(hgnc_symbol %in% symbols) %>% 
    ggplot(aes(x = RNA_Z, y = RPPA_Z, color = Ligand, shape = as.factor(Time))) + 
    geom_point(size = 4, alpha = .9) +
    scale_color_manual("Ligand", values = col$Ligand) +
    scale_shape("Time") +
    facet_wrap_paginate(~hgnc_symbol, nrow = 1, ncol = 1, page = i, scale = "free_y") +
    xlab("RNA-seq (z-score)") +
    ylab("RPPA (z-score)") +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
}

pdf(sprintf("%s/MDD_S4A.pdf", dir), height = 5, width = 5.5)
plots
dev.off()
