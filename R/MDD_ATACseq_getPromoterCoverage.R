# The purpose of this script is to get the promoter regions of 
# hg38 transcripts
library(biomaRt)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(genomation)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2019.archive.ensembl.org")

# aT <-
#   getBM(attributes = c("ensembl_gene_id", 
#                        "ensembl_transcript_id",
#                        "entrezgene"),
#         mart = mart) %>% 
#   mutate(entrezgene = as.character(entrezgene))

aT <-
  getBM(attributes = c("ensembl_gene_id", 
                       "ensembl_transcript_id"),
        mart = mart)

sampleAnno <- read.csv("../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv",
                       header = TRUE) %>% 
  mutate(specimenName = fct_inorder(specimenName))

##############################################################################

getProm <- function(txdb, annoTable, genesOnly = FALSE) {
  txdb <- keepStandardChromosomes(txdb)
  txdb <- dropSeqlevels(txdb, c("chrY", "chrM"))
  
  if (genesOnly) {
    txdb <- genes(txdb)
    prom <- promoters(txdb, upstream = 500, downstream = 500)
    prom$ensembl_gene_id <- aT[match(prom$gene_id, aT$entrezgene), 1]
  } else {
    prom <- promoters(txdb, upstream = 500, downstream = 500)
    prom$tx_name         <- str_remove(prom$tx_name, "[.][0-9]+")
    prom$ensembl_gene_id <- aT[match(prom$tx_name, aT$ensembl_transcript_id), 1]
  }
  
  filt <- is.na(prom$ensembl_gene_id)
  prom <- prom[!filt]
  
  return(prom)
  }


writeBed <- function(granges, outFile) {
  granges <- data.frame(granges) %>% 
    mutate(score = "") %>% 
    dplyr::select(seqnames, start, end, tx_name, score, strand)
  
  write.table(granges, file = outFile, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

##############################################################################
prom <- getProm(txdb, aT)

# writeBed(prom, "promoters.bed")
# Accessibility of the promoter regions was quantified using the multiBamCov
# function of bedtools-2.26.0.
# Quantification was performed in groups for faster processing. 
# Samples were split by sample replicate (e.g. 01/D, 02/A, 03/B, 04/C) into groups
# The output files ("coverages_0[1-4].bed") are re-combined into a single object
# Sample order was recorded in sampleOrder_promoterCoverage.txt

promCov <- readBed("../misc/promoterCoverage/ATACseq_promoterCov_01.bed")
promCov$score <- NULL

C1 <- mcols(promCov)
coverage_02 <- read.table("../misc/promoterCoverage/ATACseq_promoterCov_02.tsv")
coverage_03 <- read.table("../misc/promoterCoverage/ATACseq_promoterCov_03.tsv")
coverage_04 <- read.table("../misc/promoterCoverage/ATACseq_promoterCov_04.tsv")
order <- str_remove(readLines("../misc/promoterCoverage/sampleOrder_promoterCoverage.txt"), "../")

mcolsAll <- bind_cols(data.frame(C1), coverage_02, coverage_03, coverage_04)

# adding specimenIDs
colnames(mcolsAll) <- c("ensembl_transcript_id", 
                         sampleAnno %>%
                           arrange(match(bamReads, order)) %>% 
                           pull(specimenID)
                         )

# normalizing for total 
countSums <- apply(mcolsAll[, -1], 2, sum)
mcolsAll[, -1] <- 1000000*(mcolsAll[, -1]/countSums)

mcAfilt <- mcolsAll[, c("ensembl_transcript_id",
                   sampleAnno %>% 
                     filter(ATACseq_QCpass) %>%
                     pull(specimenID))] %>% 
  left_join(aT)

sumAcc <- apply(mcAfilt[, -c(1, ncol(mcAfilt))],
                1, sum)

mostAccessibleTranscripts <- data.frame(
  ensembl_gene_id = mcAfilt$ensembl_gene_id,
  ensembl_transcript_id = mcAfilt$ensembl_transcript_id,
  totAcc = sumAcc) %>% 
  arrange(desc(totAcc)) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  dplyr::select(ensembl_transcript_id,
                ensembl_gene_id)

write_csv(mostAccessibleTranscripts,
          file = "../misc/MDD_ATACseq_mostAccessibleTranscripts.csv")
