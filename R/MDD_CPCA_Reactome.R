library("fgsea", quietly = TRUE)

## Perform biological annotation-based comparison 
## INPUTS:
# factorizations = already computed factirizations
# path.database = path to a GMT annotation file
# pval.thr = p-value threshold (default to 0.05)
## OUPUTS: a list containing output values
# selectivity = Selectivity (fraction of significant annotations per all significant factors)
# nonZeroFacs = Number of unique factors that were significant at least once
# total_pathways = Number of clinical annotations with at least one significant factor component

pval.thr<-0.05
# Annotation databases used for biological enrichment
path.database <- "../Data/bio_annotations/ReactomePathways.gmt" #REACTOME
# Load annotation database
pathways <- gmtPathways(path.database)

# Number of factors
num.factors <- 11
metagenes <- moa@loading[,1:num.factors]

# Rename columns
colnames(metagenes)<-paste0("PC",1:num.factors)
# Rename rows to remove "|" characters and keep only the gene name before
rownames(metagenes)<-str_replace(rownames(metagenes), "_.*", "")
# Remove duplicated gene names that could confuse fgsea
duplicated_names <- unique(rownames(metagenes)[duplicated(rownames(metagenes))])
metagenes <- metagenes[!(rownames(metagenes) %in% duplicated_names), ]


# Calculate biological annotation enrichment.
# For each factor,
fgsea_results <- lapply(1:num.factors, function(i){
  #Order the metagenes by the factorization coefficients
  rnk <- setNames(as.matrix(metagenes[,i]), rownames(metagenes)) %>%
    sort(decreasing = TRUE) 
  # Compute fgsea
  fgseaRes <- fgseaMultilevel(pathways, rnk, minSize=15, maxSize=500, nPermSimple = 100) %>%
    tibble() %>%
    filter(padj<pval.thr) %>%
    mutate(PC = i)
  
}) %>%
  bind_rows() %>%
  mutate(PC = factor(PC, levels = as.character(1:num.factors), ordered = TRUE)) %>%
  select(PC, pathway, padj, NES, size, leadingEdge)

df <- fgsea_results %>%
  filter(NES > 2.0) %>%
  mutate(pathway = str_trunc(pathway, width = 30))

p <- ggplot(df, aes(x = pathway, y = PC, size = NES, color = padj)) +
  geom_point() +
  scale_colour_gradient2(midpoint = .025) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
p

pdf("MCF10A_MDD_revisions_reactome.pdf",height = 7, width = 14)
print(p)
res <- dev.off()