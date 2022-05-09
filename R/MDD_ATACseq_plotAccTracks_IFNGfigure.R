# 

library(tidyverse)
library(biomaRt)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Gviz)
library(motifmatchr)

annoFile     <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"

###############################################################################
# Plotting functions
makeDBSampleSheet <- function(sampleAnno) {
  # renames the columns of a sample annotation dataframe
  # to create a new dataframe that can be used as input to Diffbind.
  ss <- sampleAnno %>% 
    dplyr::rename("SampleID"   = specimenID,
                  "Collection" = collection,
                  "Time"       = experimentalTimePoint) %>% 
    mutate(Replicate           = as.factor(replicate),
           Factor              = fct_inorder(as.factor(experimentalCondition)),
           PeakCaller          = "narrow",
           Peaks               = paste0("../", Peaks),
           bamReads            = paste0("../", bamReads)
    )
  ss
}

bufferCoords <- function(START, END, BUFFER) {
  START <- as.numeric(START)
  END   <- as.numeric(END)
  
  if (length(BUFFER) == 2) {
    START <- START - BUFFER[1]
    END   <- END   + BUFFER[2]
  } else {
    START <- START - BUFFER
    END   <- END   + BUFFER
  }
  return(c(START, END))
}

makeGranges <- function(CHR, START, END, GENOME = "hg38", BINWIDTH = NULL) {
  
  if (!is.null(BINWIDTH)) {
    END   <- START + BINWIDTH*round((END - START)/BINWIDTH)
    STARTS <- seq(START, END - BINWIDTH, by = BINWIDTH) + 1
    ENDS   <- seq(START + BINWIDTH, END, by = BINWIDTH)
    DF <- data.frame(chr = CHR, start = STARTS, end = ENDS)
  } else {
    DF <- data.frame(chr = CHR, start = START, end = END)
  }
  
  GR <- makeGRangesFromDataFrame(DF)
  genome(GR) <- GENOME
  return(GR)
}

makeGrangesAndCountSplit <- function(dob, chr, start, end, binwidth = 20, score = DBA_SCORE_RPKM) {
  gr_temp   <- makeGranges(CHR = chr, START = start, 
                           END = end, BINWIDTH = binwidth)
  
  which <- as.character(unique(dba.show(dob)$Factor))
  
  dobsByFactor <-
    lapply(which, function(x) {
      print(x)
      dobX <- dba(dob, mask = dob$masks[[x]])
      print(dba.show(dobX))
      return(dobX)
    })
  
  dobGranges     <- lapply(dobsByFactor, dba.count, minOverlap = 1,
                           peaks = gr_temp, score = score)
  grangesCounted <- lapply(dobGranges,
                           dba.peakset, bRetrieve = TRUE)
  
  names(grangesCounted) <- which
  return(grangesCounted)
}

makeGrangesFromGene <- function(symbol, buffer, basicDob, mart, aTrna = NULL) {
  
  if(is.null(aTrna)) {
    aTrna <- getBM(attributes = c("ensembl_gene_id",
                                  "hgnc_symbol",
                                  "chromosome_name",
                                  "start_position",
                                  "end_position"), 
                   mart = mart,
                   filters = "hgnc_symbol",
                   values = symbol) %>% 
      mutate(chromosome_name = paste0("chr", chromosome_name))
  }
  
  if(nrow(aTrna) > 1) {
    stop("biomaRt returned multiple matches to symbol")
  }
  
  range <- bufferCoords(aTrna$start_position, aTrna$end_position, buffer)
  chr   <- aTrna$chromosome_name
  start <- range[1]
  end   <- range[2]
  
  geneGranges <- makeGrangesAndCountSplit(dob = basicDob,
                                          chr = chr,
                                          start = start,
                                          end = end,
                                          binwidth = 20)
  return(geneGranges)
}

makeGeneTrackFromGene <- function(symbol, mart, grayCol = "gray70") {
  biomTrack <- BiomartGeneRegionTrack(genome="hg38",
                                      biomart = mart,
                                      filters=list(hgnc_symbol=symbol),
                                      name=symbol,
                                      collapseTranscripts = "longest",
                                      col = "black",
                                      col.line = "black",
                                      fill = grayCol,
                                      protein_coding = grayCol,
                                      "utr3" = grayCol,
                                      "composite" = grayCol,
                                      "utr5" = grayCol)
}

geneToStacks <- function(symbol, buffer, basicDob, mart, aTrna = NULL) {
  
  colorLookup <- c("ctrl_0"  = "#7A4A2A",
                   "PBS_24" = "#8dd3c7",
                   "PBS_48" = "#8dd3c7",
                   "HGF_24"  = "#80b1d3",
                   "HGF_48"  = "#80b1d3",
                   "OSM_24"  = "#fdb462",
                   "OSM_48"  = "#fdb462",
                   "EGF_24"  = "#FB8072",
                   "EGF_48"  = "#FB8072",
                   "BMP2_24"  = "#B3DE69",
                   "BMP2_48"  = "#B3DE69",
                   "IFNG_24" = "#BEBADA",
                   "IFNG_48" = "#BEBADA",
                   "TGFB_24" = "#FFD92F",
                   "TGFB_48" = "#FFD92F",
                   "ctrl_0"  = "#7A4A2A")
  
  biomTrack <- makeGeneTrackFromGene(symbol, mart = mart)
  grList <- makeGrangesFromGene(symbol, buffer, basicDob, mart, aTrna = aTrna)
  
  dTracks <- mapply(function(x, x_nm) {
    DTx <- DataTrack(x,
                     range             = x,
                     genome            = "hg38",
                     isPaired          = TRUE,
                     name = x_nm,
                     fill              = colorLookup[x_nm],
                     col.histogram     = colorLookup[x_nm]
    )
    return(DTx)
  } , x = grList, x_nm = names(grList)
  )
  
  return(c(dTracks, biomTrack))
}

geneToPlot <- function(symbol, buffer, basicDob, mart, sA, extraGenes = NULL) {
  aTrna <- getBM(attributes = c("ensembl_gene_id",
                                "hgnc_symbol",
                                "chromosome_name",
                                "start_position",
                                "end_position"), 
                 mart = mart,
                 filters = "hgnc_symbol",
                 values = symbol) %>% 
    mutate(chromosome_name = paste0("chr", chromosome_name))
  
  plotRange <- bufferCoords(aTrna$start_position, aTrna$end_position, 
                            buffer)
  
  stack <- geneToStacks(symbol, buffer, basicDob, mart, aTrna)
  
  axisTrack   <- GenomeAxisTrack()
  
  if(!is.null(extraGenes)) {
    extraGeneTracks <- lapply(extraGenes, 
                              makeGeneTrackFromGene, 
                              mart = mart)
    
    plotTracks(c(stack, extraGeneTracks, axisTrack),
               from = plotRange[1],
               to   = plotRange[2],
               type = "histogram",
               scale = .1)
  } else {
    plotTracks(c(stack, axisTrack),
               from = plotRange[1],
               to   = plotRange[2],
               type = "histogram",
               scale = .1)
  }
}


###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

###############################################################################
# for annotating peaks with nearest gene, etc.
###############################################################################
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

sampleAnno <- read.csv(annoFile, stringsAsFactors = FALSE) %>% 
  filter(ATACseq_QCpass) %>% 
  filter(ligand %in% c("ctrl", "EGF", "IFNG")) %>% 
  filter(experimentalTimePoint %in% c(0, 48))

###############################################################################
# Creating DiffBind samplesheet. Low-quality samples are removed.
sSheet    <- makeDBSampleSheet(sampleAnno) %>% 
  mutate(Factor = fct_recode(Factor, CTRL = "ctrl_0", EGF = "EGF_48", "IFNG+EGF" = "IFNG_48")) %>% 
  mutate(Factor = fct_relevel(Factor, "CTRL", "EGF", "IFNG+EGF")) %>% 
  arrange(Factor)
dob_temp  <- dba(sampleSheet = sSheet)

###############################################################################
# Genes of interest
outDir <- "../plots/MDD_manuscript_figures"

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

myGene <- "IRF1"
eg <- "IRF1-AS1"
pdf(sprintf("%s/MDD_S3F_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp, 
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()

myGene <- "STAT1"
eg <- c("GLS", "STAT4")
pdf(sprintf("%s/MDD_S3F_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 8)
geneToPlot(symbol = myGene, 
           buffer = c(50000),
           basicDob = dob_temp,
           mart = mart,
           sA = sampleAnno,
           extraGenes = eg)

dev.off()  


