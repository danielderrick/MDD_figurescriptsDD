library(DESeq2)
library(tidyverse)
library(plotly)
library(biomaRt)
library(genomation)
library(cowplot)

RNAseqFile           <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"
RNAseqAnnoFile           <- "../RNAseq/Metadata/MDD_RNAseq_geneAnnotations.csv"
RNAseqMetadataFile       <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"

transcriptsFile <- "../misc/MDD_ATACseq_mostAccessibleTranscripts.csv"
HLAGeneListFile     <- "../misc/MDD_ATACseq_HLAGeneList.csv"
combinedTableFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"

saFile     <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
promCovDir <- "../misc/promoterCoverage"
orderFile  <- "../misc/promoterCoverage/sampleOrder_promoterCoverage.txt"


outDir <- "../plots/MDD_manuscript_figures/"
theme_set(theme_cowplot())

###############################################################################

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'https://apr2019.archive.ensembl.org')

at.RNA <- getBM(attributes = c("ensembl_gene_id",
                               "gene_biotype"), 
                mart = mart)

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

###############################################################################

sa.ATAC.L3 <- read.csv(saFile,
                       header = TRUE,
                       stringsAsFactors = FALSE) %>% 
  mutate(specimenName = fct_inorder(as.factor(specimenName)))

###############################################################################

promCov <- readBed(sprintf("%s/ATACseq_promoterCov_01.bed", promCovDir))
promCov$score <- NULL
coverage_01 <- mcols(promCov)
coverage_02 <- read.table(sprintf("%s/ATACseq_promoterCov_02.tsv", promCovDir))
coverage_03 <- read.table(sprintf("%s/ATACseq_promoterCov_03.tsv", promCovDir))
coverage_04 <- read.table(sprintf("%s/ATACseq_promoterCov_04.tsv", promCovDir))

myOrder <- str_replace(readLines(orderFile), "../", "ATACseq/")
mcolsAll <- bind_cols(data.frame(coverage_01), coverage_02, coverage_03, coverage_04)

# adding specimenIDs
colnames(mcolsAll) <- c("ensembl_transcript_id", 
                        sa.ATAC.L3 %>%
                          arrange(match(bamReads, myOrder)) %>% 
                          pull(specimenID)
                        )

# adding specimenIDs
MDD_ATACseq_promAcc <- mcolsAll[, c("ensembl_transcript_id",
                                    sa.ATAC.L3 %>% 
                                      filter(ATACseq_QCpass) %>%
                                      pull(specimenID))] %>% 
  column_to_rownames("ensembl_transcript_id") %>% 
  data.matrix


# table of the most-accessible promoter for every ensembl_gene_id
bestTranscripts <- read_csv(transcriptsFile) %>% 
  filter(!is.na(ensembl_gene_id))

# Selecting most acc. promoter and renaming with gene id
MDD_ATACseq_promAcc <- MDD_ATACseq_promAcc[as.character(bestTranscripts$ensembl_transcript_id), ]
rownames(MDD_ATACseq_promAcc) <- bestTranscripts$ensembl_gene_id

# Creating DESeq2 dataset from promoter accessibility matrix
# the purpose is to create a matrix of promoter accessibility 
# after DESeq2 variance-stabilized transformation

dds <- DESeqDataSetFromMatrix(MDD_ATACseq_promAcc,
                              sa.ATAC.L3 %>% filter(ATACseq_QCpass),
                              ~experimentalCondition)
dds <- DESeq(dds)

ATACseqPromotersL3 <- 
  vst(dds) %>% 
  assay() %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id")

ATACseqL3Long <- 
  ATACseqPromotersL3 %>% 
  pivot_longer(-ensembl_gene_id, names_to = "specimenID", values_to = "accessibility") %>% 
  left_join(dplyr::select(sa.RNA.L3, specimenID, experimentalCondition, Ligand, Time))

###############################################################################

RNAseqL3Long <-
  read_csv(RNAseqFile) %>% 
  data.frame %>% 
  pivot_longer(-ensembl_gene_id, names_to = "specimenID", values_to = "expression") %>% 
  left_join(dplyr::select(sa.RNA.L3, specimenID, experimentalCondition, Ligand, Time))

###############################################################################

L3JoinedLong <-
  inner_join(ATACseqL3Long, RNAseqL3Long) 

L3JFilt <-
  L3JoinedLong %>% 
  filter(experimentalCondition == "EGF_48")

###############################################################################

at.RNA <-
  at.RNA %>% 
  left_join(read_csv(HLAGeneListFile)) %>% 
  mutate(HLAregion = replace_na(HLAregion, FALSE))

L3JFiltAnno <- 
  L3JFilt %>% 
  left_join(dplyr::select(at.RNA, ensembl_gene_id, gene_biotype, HLAregion))

L3JFiltAnnoSummarized <-
  L3JFiltAnno %>% 
  group_by(ensembl_gene_id, gene_biotype, HLAregion) %>% 
  dplyr::summarize(medianExpr = median(expression),
            medianAcc = median(accessibility)) %>% 
  ungroup

###############################################################################
if(!dir.exists(outDir)) {dir.create(outDir)}

pdf(sprintf("%s/MDD_SF4E.pdf", outDir))
L3JFiltAnnoSummarized %>% 
  filter(gene_biotype == "protein_coding") %>% 
  filter(!HLAregion) %>% 
  ggplot(aes(x = medianAcc, y = medianExpr)) +
  geom_point(alpha = .25, size = .75, shape = 16) +
  xlab("promoter accessibility\n(vst counts)") +
  ylab("gene expression\n(log2(fpkm + 1))") +
  ggtitle("TSS accessibililty vs gene expression", subtitle = "EGF_48, protein-coding genes") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
dev.off()
