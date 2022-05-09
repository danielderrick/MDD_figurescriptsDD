# The purpose of this script is to plot the top-ranked
# transcription factors from CHEA3 analysis of the
# RNAseq features in the MDD integrated modules.
# 
#
###############################################################################

# Libraries
library(tidyverse)
library(cowplot)

# Paths
remapCSVFile <- "../misc/CHEA3_ReMap_TFEnrichments/MDD_CHEA3_ReMap.csv"
outDirPlots  <- "../plots/MDD_manuscript_figures"

libaryName    <- "ReMap--ChIP-seq"
q_thresh      <- 0.2
n_rank        <- 5
###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

###############################################################################
# reading in CHEA3 enrichments from the ReMap library
remapLong <- read_csv(remapCSVFile)

###############################################################################
# Figure 5C
remapLong <- remapLong %>%
  dplyr::rename(Module = `Query Name`)

Figure5C_point1_binary_top5 <-
  remapLong %>% 
  mutate(significant = as.numeric(FDR) < q_thresh) %>% 
  mutate(Rank = as.numeric(Rank)) %>% 
  mutate(topFive = Rank <= 5) %>% 
  group_by(TF) %>% 
  filter(any(Rank <= n_rank)) %>% 
  ungroup() %>% 
  mutate(Module = str_remove(Module, "module_")) %>% 
  mutate(Module = fct_rev(fct_inorder(as.factor(Module)))) %>% 
  arrange(Module, FDR) %>% 
  mutate(TF = (fct_inorder(as.factor(TF)))) %>% 
  dplyr::select(Module, TF, significant, topFive, Rank) %>% 
  mutate(TF = fct_drop(TF)) %>% 
  # complete(Module,TF, fill = list("significant" = 0)) %>%
  data.frame %>% 
  arrange((Module), significant, Rank) %>% 
  filter(significant == 1 |
           topFive) %>% 
  mutate(TF = fct_inorder(TF)) %>% 
  complete(Module,TF, fill = list("topFive" = FALSE,
                                  significant = FALSE)) %>%
  arrange(Module, Rank) %>% 
  mutate(TF = fct_inorder(TF)) %>% 
  mutate(FDR2 = case_when(significant ~ "significant",
                          !significant & topFive ~ "NSTF",
                          !significant & !topFive ~ "nothing")) %>% 
  ggplot(aes(x = TF, y = Module, fill = topFive, color = FDR2)) +
  geom_tile(height = .7, width = .7, size = .6) +
  scale_color_manual("", values = c("significant" = "firebrick4",
                                    "NSTF" = "pink",
                                    "nothing" = "gray95"),
                     labels = c("FDR < 0.2", "", "")) +
  scale_fill_manual("",
                    values = c("TRUE" = "pink", "FALSE" = "gray95"),
                    labels = c("Top Rank", "")) +
  xlab("Transcription Factor") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

if (!dir.exists(outDirPlots)) {dir.create(outDirPlots, recursive = TRUE)}

pdf(sprintf("%s/MDD_F5C.pdf", outDirPlots), height = 4.25, width = 9)
Figure5C_point1_binary_top5
dev.off()

