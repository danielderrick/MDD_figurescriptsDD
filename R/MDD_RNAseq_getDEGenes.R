# The purpose of this script is to run differential expression analyses
# for each condition vs ctrl_0. For each condition, the shrunken log2FoldChange
# estimates will then be used to rank genes. The top n (25-50) for each condition 
# will be selected. These lists will be combined (and any duplicates removed)
# to produce a list of genes for use in the "mega matrix" heatmap.

library(tidyverse)
library(DESeq2)

RNAannoFile <- "../RNAseq/Metadata/MDD_RNAseq_geneAnnotations.csv"

###############################################################################

if( str_extract(getwd(), "[:alnum:]+$") != "R" ) {
  setwd("R")
}

load("../rnaSeq/MDD_RNAseq_shrunkenResults.Rdata")

annoTable <- read_csv(RNAannoFile)

res <- list()
for (i in 2:15) {
  res[[i - 1]] <- lfcShrink(dds, coef = i, type = "apeglm")
}

names(res) <- str_extract(resultsNames(dds)[2:15], "[:alnum:]+_[248]{2}")

getSigGenes <- function(res, lfc = 1.5) {
  res <- data.frame(res)
  res <- res %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    left_join(annoTable) %>% 
    filter(hgnc_symbol != "") %>% 
    filter(!duplicated(hgnc_symbol)) %>% 
    dplyr::select(-ensembl_gene_id) %>% 
    filter(padj < 0.01) %>% 
    filter(abs(log2FoldChange) > lfc) %>% 
    pull(hgnc_symbol)
}

getUniqueList <- function(resList, lfc = 1.5) {
  mostVarList <- lapply(resList, 
                        getSigGenes, 
                        lfc = lfc)
  mostVar <- unlist(mostVarList)
  mostVarDF <-
    data.frame(hgnc_symbol = mostVar,
               condition = str_extract(names(mostVar), "[:alnum:]+_[248]{2}"),
               stringsAsFactors = FALSE) %>% 
    distinct(hgnc_symbol, .keep_all = TRUE)
  mostVarDF
}


mvL <- lapply(res_RNA, 
              getSigGenes, 
              lfc = 1.5)

mv <- unlist(mvL)

uniq <- getUniqueList(res_RNA, 1.5)

write_csv(uniq %>% dplyr::select(hgnc_symbol) %>% unique, path = "../misc/MDD_geneList_lfc15.csv")
