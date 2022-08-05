# The purpose of this script is to submit the RNA-seq features
# of the MDD integrated modules to CHEA3 and write the
# TF enrichment results for the ReMap library to CSV
# and Excel files.
#
###############################################################################

# Libraries
library(tidyverse)
library(DT)
library(jsonlite)
library(httr)
library(readxl)   # for excel_sheets
library(openxlsx) # for write.xlsx
library(cowplot)

libaryName                <- "ReMap--ChIP-seq"

# Paths
integratedModuleExcelFile <- "../misc/MDD_multiomics_14.xlsx"
outDir                    <- "../misc/CHEA3_ReMap_TFEnrichments2"
outDirPlots               <- "../misc/CHEA3_ReMap_TFEnrichments2/plots"

# Functions
excelToReMap <- function(x) {
  # This function imports the RNA-seq features from an excel file of
  # integrated modules, submits the gene lists to CHEA3, and returns
  # results the ReMap ChIP-seq library
  
  x  <- importRNAFeatures(x)
  rL <- getCHEA3LibraryList(x, libraryName == "ReMap--ChIP-seq")
  
  return(rL)
}

importRNAFeatures <- function(excelFileName){
  # The purpose of this function is to import the RNA features
  # from an excel sheet of multiome features
  
  readExcelSheets <- function(eFN) {
    # The purpose of this function is to import all sheets of an Excel file
    # as a named list of data.frames.
    sheets <- excel_sheets(eFN)
    
    myT <- lapply(1:length(sheets), function(x) {
      X <- read_excel(eFN, sheet = x)
      X
    })
    names(myT) <- sheets
    
    return(myT)
  }
  
  getGeneSymbols <- function(iM){
    # The purpose of this script is to return the RNA features 
    # from a data frame of integrated MDD modules 
    x <- iM %>% 
      filter(Type == "RNAseq") %>% 
      distinct(feature) %>% 
      pull(feature) %>% 
      as.character()
    
    return(x)
  }
  
  eS <- readExcelSheets(excelFileName)
  geneSymbolList <- lapply(eS, getGeneSymbols)
  
  return(geneSymbolList)
}

getCHEA3LibraryList <- function(geneSetList, libraryName){
  # The purpose of this script is to submit each gene set in
  # a list of genesets to CHEA3 and return only the results
  # from the named library, as a named list
  
  getCHEA3Library <- function(geneList, queryName, libraryName) {
    # The purpose of this function is to submit a list of gene symbols to the 
    # CHEA3 API and return TF enrichment results from the named library
    
    getCHEA3 <- function(gL, qN, url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"){
      # This function submits a list of gene names to CHEA3 and
      # returns all enrichment results
      # Adapted from R example at https://maayanlab.cloud/chea3/
      
      toCHEA3 = list(query_name = qN, 
                     gene_set = gL)
      
      response = POST(url = url, body = toCHEA3, encode = "json")
      json = content(response, "text")
      
      #results as list of R dataframes
      results <- try(jsonlite::fromJSON(json))
      return(results)
    }
    
    filterToLibrary <- function(rCHEA3, lN) {
      # The purpose of this function is to subset CHEA3
      # enrichment results to the ReMap ChIP-seq library
      resultsCHEA3Filt <- rCHEA3[[lN]]
      
      return(resultsCHEA3Filt)
    }
    
    x <- getCHEA3(gL = geneList, qN = queryName)
    x <- filterToLibrary(x, 
                         lN = libraryName)
    
    return(x)
  }
  
  libraryResultsList <- mapply(getCHEA3Library,
                               geneList = geneSetList,
                               queryName = names(geneSetList),
                               libraryName = libraryName,
                               SIMPLIFY = FALSE)
  
  names(libraryResultsList) <- names(geneSetList)
  
  return(libraryResultsList)
}

writeListToExcel <- function(x, fName) {
  # The purpose of this script is to write a named list to an excel file.
  # list names will be used to name sheets
  write.xlsx(x = x,
             file = fName,
             rowNames = FALSE)
}

writeListToCSV <- function(x, fName) {
  # The purpose of this script is to collate a list of data.frames
  # into one long data.drame, then write to a CSV file.
  xLong <- Reduce(bind_rows, x)
  
  write_csv(xLong, 
            file = fName)
}

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

###############################################################################
# Submitting RNA features to CHEA3 API

x <- importRNAFeatures(integratedModuleExcelFile)
remap <- getCHEA3LibraryList(x, libaryName)

###############################################################################
# Writing enrichment tables to CSV/excel
if (!dir.exists(outDir)) {dir.create(outDir)}

writeListToExcel(remap, 
                 fName = sprintf("%s/MDD_CHEA3_%s.xlsx", 
                                 outDir, str_extract(libaryName,
                                                     "[:alnum:]+")))
writeListToCSV(remap, 
               fName = sprintf("%s/MDD_CHEA3_%s.csv", 
                               outDir, str_extract(libaryName,
                                                   "[:alnum:]+")))

