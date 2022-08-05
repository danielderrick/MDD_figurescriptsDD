# The purpose of this script is to do the following annotations on the MDD
# integrated module features:
# RNA-seq features are annotated if they represent a transcription factor,
# kinase/epigenetic modifier, or a lncRNA.
#
# TF gene symbols are taken from JASPAR 2018 database and 
# from GO annotation "DNA-binding transcription factor activity"
#
# Epigenetic modifier annotations from dbEM (Nanda, 2016)
#
# kinase annotations from uniProt
#
library(tidyverse)
library(biomaRt)
library(cmapR)
library(eulerr)
library(ComplexHeatmap)
library(circlize)

geneListFile      <- "../misc/MDD_multiomics14Module_allFeatures.csv"
RNA_resultsFile   <- "../RNAseq/Data/MDD_RNAseq_shrunkenResults_df.Rdata"
combinedTableFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"

epigeneticModifiersFile <- "../misc/exp_cnv_final.csv"
kinaseListFile          <- "../misc/kinases-uniProt.txt"

outDir <- "../misc/"
###########################
jaspar <- function (...) 
  # This function is edited from the getJasparMotifs function in 
  # the R package chromVAR.
  # It has been changed to retrieve JASPAR 2018 motifs, rather than
  # the JASPAR 2016 motifs retrieved by the function in chromVAR
{
  opts <- list()
  opts["species"] <- "Homo sapiens"
  opts["collection"] <- "CORE"
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), 
                        sep = "_")
  return(out)
}

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

###############################################################################
# TF ANNOTATIONS
# Importing module table
modules <- read.csv(geneListFile, stringsAsFactors = FALSE)

# Importing genes with GO annotation `DNA-binding transcription factor activity` 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "useast.ensembl.org")

at.RNA <- getBM(attributes = c("ensembl_gene_id",
                               "hgnc_symbol",
                               "gene_biotype"), 
                mart = mart)

GO_genes <- getBM(attributes = "hgnc_symbol",
                  filters = "go",
                  values = "GO:0003700",
                  mart = mart)

# Importing and formatting JASPAR 2018 TFs
myJ <- jaspar()
jTFs <- lapply(myJ, function(x) {x <- x@name}) %>%
  unlist %>% 
  as.character %>% 
  set_names(., NULL) %>% 
  str_split_fixed("[:(]", 2) %>% 
  .[, 1]

# Euler diagram: JASPAR and GO transcription factors
myListForEuler <- list(
  JASPAR = unique(jTFs),
  GO = unique(GO_genes$hgnc_symbol))

euler(myListForEuler) %>% 
  plot(quantities = TRUE)

# Manually add in TCF19 (Transcription Factor 19), a TF that 
# is missing from these lists
allTFs <- c(unique(c(myListForEuler$JASPAR, myListForEuler$GO)), "TCF19")

modules <- 
  modules %>% 
  mutate(TFany = feature %in% allTFs) %>% 
  mutate(TF = feature %in% allTFs &
           Type == "RNAseq") %>% 
  dplyr::select(feature_type, Cluster, Set, feature, Type, TF, everything())

# separating module 10
module10 <-
  modules %>% 
  filter(Cluster == "module_10") %>% 
  dplyr::select(feature_type, Cluster, Set, feature, Type, TF, everything()) %>% 
  arrange(desc(TF))

module10_TFs <-
  modules %>% 
  filter(Cluster == "module_10") %>% 
  filter(TF) %>% 
  dplyr::select(feature_type, Cluster, Set, feature, Type, TF, everything())

###############################################################################
# EPIGENETIC MODIFIER ANNOTATIONS
epigeneticModifiers <-
  read.table(epigeneticModifiersFile, sep = "\t")$V1 %>% unique

module10 <- 
  module10 %>% 
  mutate(epigeneticModifer = feature %in% epigeneticModifiers)

###############################################################################
# KINASE ANNOTATIONS
at.RNA.kinases_uniProt <- 
  read_table(kinaseListFile, col_names = FALSE) %>% 
  dplyr::mutate(hgnc_symbol1 = X1,
                hgnc_symbol2 = str_extract(X2, "[:alnum:]+_")) %>% 
  mutate(hgnc_symbol2 = str_remove(hgnc_symbol2, "_")) %>% 
  dplyr::select(contains("hgnc_symbol")) %>% 
  unlist %>% 
  unique

module10 <- 
  module10 %>% 
  mutate(kinase = feature %in% at.RNA.kinases_uniProt)

module10 <-
  module10 %>% 
  left_join(at.RNA, by = c("feature" = "hgnc_symbol")) %>% 
  mutate(gene_biotype = case_when(Type == "RPPA" ~ "RPPA Features",
                                  TF             ~ "Transcription Factor",
                                  kinase             ~ "Kinase",
                                  epigeneticModifer         ~ "Epigenetic",
                                  Type != "RPPA" &
                                    !kinase &
                                    !epigeneticModifer &
                                    !TF ~ gene_biotype)) %>% 
  mutate(gene_biotype = replace_na(gene_biotype, "intronic transcript")) %>%
  mutate(gene_biotype = as.factor(gene_biotype)) %>% 
  mutate(gene_biotype = fct_collapse(gene_biotype,
                                     "RNA" = c("intronic transcript",
                                               "processed_pseudogene",
                                               "unprocessed_pseudogene",
                                               "protein_coding",
                                               "transcribed_processed_pseudogene",
                                               "transcribed_unprocessed_pseudogene"))) %>% 
  mutate(gene_biotype = fct_recode(gene_biotype, 
                                   "lncRNA" = "lncRNA") )%>% 
  dplyr::rename(Category = "gene_biotype") %>% 
  dplyr::select(1:6, "kinase", everything()) %>% 
  distinct(feature_type, .keep_all = TRUE)
  
module10Summary <- 
  module10 %>% 
  group_by(Category) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>% 
  pull(n) %>% sum()

if (!dir.exists(outDir)) {dir.create(outDir)}

# writing to file
write_csv(modules, sprintf("%s/MDD_multiomics14Modules_TFannotations_allModules.csv", outDir))
write_csv(module10, sprintf("%s/MDD_multiomics14Modules_TFannotations_module10.csv", outDir))
