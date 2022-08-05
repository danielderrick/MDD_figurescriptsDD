library(jsonlite)
library(tidyverse)
library(GSVA)
library(cmapR)
library(GSEABase)

geneSets <- getGmt("../misc/CHEA3_ReMap_ChIP-seq.gmt")
RNAseqFile           <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"
RNAseqAnnoFile           <- "../RNAseq/Metadata/MDD_RNAseq_geneAnnotations.csv"
RNAseqMetadataFile       <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"



## Calculating gene set enrichment for ReMap ChIP-seq TF gene sets

# In this document, transcription factor enrichment scores are calculated from
# MCF10A RNAseq values using the R package `GSVA` and transcription factor gene
# sets derived from ReMap ChIP-seq data. Enrichment scores will be calculated using
# multiple parameter options. The performance of the analysis will be evaluated
# A) by examining enrichment of known TFs (e.g. STAT1, IRF1, etc), and B) by 
# evaluating replicate correlation.

## Importing RNAseq data

at.RNA <- read_csv(RNAseqAnnoFile)

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


RNAseq.L3 <- read.csv("../RNAseq/Data/MDD_RNAseq_Level3.csv", 
                      header = TRUE, stringsAsFactors = FALSE) %>% 
  left_join(at.RNA) %>% 
  filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  dplyr::select(-ensembl_gene_id) %>% 
  column_to_rownames("hgnc_symbol") %>% 
  data.matrix()

filt <- apply(RNAseq.L3, 1, function(x) x <- sum(x >= .5) >= 3)
RNAseq.L3 <- RNAseq.L3[filt, ]
RNAseq.L3 <- data.matrix(RNAseq.L3)


## Running GSVA

GSVA_results <- list()

GSVA_results$L3.ss <- gsva(expr = RNAseq.L3, gset.idx.list = geneSets, method = "ssgsea")

outDir <- "../misc"

if (!dir.exists(dir)) {
  dir.create(dir)
}

ssGSEA_L4 <-
  GSVA_results$L3.ss %>% 
  data.frame %>% 
  rownames_to_column("hgnc_symbol") %>% 
  pivot_longer(-hgnc_symbol, values_to = "enrichment", names_to = "specimenID") %>% 
  left_join(sa.RNA.L3) %>% 
  group_by(hgnc_symbol, experimentalCondition) %>% 
  summarize(medianEnrichment = median(enrichment)) %>% 
  pivot_wider(names_from = "experimentalCondition", values_from = "medianEnrichment") %>% 
  column_to_rownames("hgnc_symbol")

write_csv(
  rownames_to_column(
    data.frame(GSVA_results$L3.ss),
    "hgnc_symbol"), 
  path = "../misc/MDD_RNAseq_ssGSEAScores.csv")

