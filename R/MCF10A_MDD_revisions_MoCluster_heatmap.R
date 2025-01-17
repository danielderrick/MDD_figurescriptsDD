#testing moCluster methods on MCF10A MDD integrated clusters
#manually load combined18_features_dm
library(mogsa)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)

keep_PBS <- FALSE
PBS_pointer <- "_noPBS2"
if(keep_PBS) PBS_pointer <- ""

#Read in consensus clustering data
combined14_features_data_ann <- read_csv("../Data/MDD_multiomics_14.csv") 

c14_consensus_dm <- combined14_features_data_ann %>%
  select(matches("_[0-9][0-9]")) %>%
  as.matrix()

#load the raw data asay_pk_data for future processing
load("../Data/assay_pk_data.rda")
data_path <- paste0("../Data/integrated_matrix_lfc_rr",PBS_pointer, ".csv")
p_values_path <- paste0("../Data/integrated_adj_p_values", PBS_pointer, ".csv") 
shared_features_path <- paste0("../Data/integrated_shared_features_rr", PBS_pointer, ".csv") 
unique_features_path <- paste0("../Data/integrated_unique_features_rr", PBS_pointer, "_absDir.csv") 

experimentalTimePoints <- c("_24", "_48")
if(keep_PBS) {
  ligand_order <-  c("PBS", "HGF", "OSM", "EGF","BMP2+EGF", "IFNG+EGF", "TGFB+EGF")
} else {
  ligand_order <-  c("HGF", "OSM", "EGF","BMP2+EGF", "IFNG+EGF", "TGFB+EGF")
}

condition_order <- paste0(rep(ligand_order, each = length(experimentalTimePoints)), rep(c("_24", "_48"), times = length(ligand_order)))

combined_clust_num <- 18
cluster_method <- "pam"

#Read in the integrated data, add EGF to some names, set order
if(!file.exists(data_path)) stop("cannot read ",data_path)
all_features_data <- read_csv(data_path) %>%
  pivot_longer(cols = matches("_24|_48"), names_to = "condition") %>%
  mutate(condition = str_replace(condition, "BMP2","BMP2+EGF"),
         condition = str_replace(condition, "IFNG","IFNG+EGF"),
         condition = str_replace(condition, "TGFB","TGFB+EGF"),
         feature_type =  paste(feature, Type, sep = "_")) %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  select(feature, Type, feature_type, all_of(condition_order)) %>%
  drop_na

#Read in the integrated data, add EGF to some names, set order
if(!file.exists(unique_features_path)) stop("cannot read ",unique_features_path)
set.seed(42) #randomize row positions
unique_features <- read_csv(unique_features_path) %>%
  slice_sample(prop = 1)
unique_features_data <- unique_features%>%
  mutate(ligand = str_replace(ligand, "BMP2","BMP2+EGF"),
         ligand = str_replace(ligand, "IFNG","IFNG+EGF"),
         ligand = str_replace(ligand, "TGFB","TGFB+EGF"),
         feature_type =  paste(feature, Type, sep = "_"),
         Set = ligand) %>%
  select(ligand, feature_type, Abs_Direction, Set) %>%
  rename(Direction = Abs_Direction) %>%
  left_join(all_features_data, by = "feature_type") %>%
  rename(Cluster = ligand) %>%
  mutate(Cluster = paste(Cluster, Direction, sep = "_")) %>%
  drop_na

#Read in the integrated data, add EGF to some names, set order
if(!file.exists(shared_features_path)) stop("cannot read ",shared_features_path)
set.seed(42) #randomize row positions
shared_features <- read_csv(shared_features_path) %>%
  slice_sample(prop = 1)
shared_features_data <- shared_features %>%
  mutate(ligand = str_replace(ligand, "BMP2","BMP2+EGF"),
         ligand = str_replace(ligand, "IFNG","IFNG+EGF"),
         ligand = str_replace(ligand, "TGFB","TGFB+EGF"),
         feature_type =  paste(feature, Type, sep = "_"),
         Set = "shared") %>%
  select(feature_type, Set) %>%
  distinct() %>%
  left_join(all_features_data, by = c("feature_type" ))  %>%
  drop_na

set.seed(42) #randomize row positions
combined18_features_data <- bind_rows(shared_features_data, select(unique_features_data, -Direction, -Cluster)) %>%
  slice_sample(prop = 1)
assay_cols <- structure(RColorBrewer::brewer.pal(8, "Paired")[c(3, 2,6,5, 1)], names = unique(combined18_features_data$Type))

#create a numeric matrix with feature row names
combined18_features_dm <- combined18_features_data %>%
  select(all_of(condition_order)) %>%
  as.matrix() 
rownames(combined18_features_dm) <- combined18_features_data$feature_type
combined18_features_dm <- combined18_features_dm[!is.nan(combined18_features_dm[,1]),]


#cluster all features 
if(cluster_method == "kmeans"){
  #seed the clustering with a fixed, random set of points
  set.seed(42)
  centers <- combined18_features_dm[sample(1:nrow(combined18_features_dm), size = combined18_clust_num, replace = FALSE),]
  clusters <- kmeans(x = combined18_features_dm, iter.max = 20, centers = centers)$cluster
} else if(cluster_method == "HA") {
  clust_obj <- hclust(dist(combined18_features_dm))
  clusters <- cutree(clust_obj, k = combined_clust_num)
} else if(cluster_method == "pam") {
  
  if(file.exists(paste0("../Data/MDD_",cluster_method,"_clusters.rda"))){
    load(paste0("../Data/MDD_",cluster_method,"_clusters.rda"))
  } else {
    set.seed(42)
    medoids <- sample(1:nrow(combined18_features_dm), size = combined_clust_num, replace = FALSE)
    clusters  <- pam(x = combined18_features_dm, k = combined_clust_num, medoids = medoids, cluster.only=TRUE, pamonce = 5)
    save(clusters, file = paste0("../integrated_analysis/MDD_",cluster_method,"_clusters.rda"))
  }
  
} else stop("unsupported cluster method chosen")

#Add the cluster assignments to the data
combined18_features_data_ann <-  tibble(feature_type = names(clusters),
                                        Cluster = clusters) %>%
  right_join(combined18_features_data, by = "feature_type") %>%
  mutate(Cluster = factor(Cluster, ordered = TRUE))

combined18_features_ann <- combined18_features_data_ann %>%
  select(Cluster, Type, Set) %>%
  as.data.frame()

#Create a matrix and heatmap of the mean combined18 values
combined18_features_data_mean <- combined18_features_data_ann %>%
  group_by(Cluster) %>%
  summarise(across(.cols = matches("_24|48"), .fns = mean, .groups = "drop"))

combined18_features_data_mean_dm <- combined18_features_data_mean %>%
  select(-Cluster) %>%
  as.matrix()
rownames(combined18_features_data_mean_dm) <- paste0("Module_",combined18_features_data_mean$Cluster)

combined14_features_data_mean <- combined14_features_data_ann %>%
  group_by(Cluster) %>%
  summarise(across(.cols = matches("_24|48"), .fns = mean, .groups = "drop"))

combined14_features_data_mean_dm <- combined14_features_data_mean %>%
  select(-Cluster)  %>%
as.matrix()
rownames(combined14_features_data_mean_dm) <- paste0("Module_",combined14_features_data_mean$Cluster)

#create a heatmap with treatments on the x axis and module means on the y

hm_means_treatments <- Heatmap(combined14_features_data_mean_dm,
                  column_title = paste("MCF10A MDD module means"),
                  name = "rrscale vale\n(mean)",
                  row_names_gp = gpar(fontsize = 10),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  #split = clusters,
                  column_names_gp = gpar(fontsize = 8),
                  col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                  column_title_gp = gpar(fontsize = 8))
hm_means_treatments
#separate the assays for MoCluster analysis
RNAseq_dm <- combined18_features_dm[str_detect(rownames(combined18_features_dm), "RNAseq"),]
RPPA_dm <- combined18_features_dm[str_detect(rownames(combined18_features_dm), "RPPA"),]
GCP_dm <- combined18_features_dm[str_detect(rownames(combined18_features_dm), "GCP"),]
cycIF_dm <- combined18_features_dm[str_detect(rownames(combined18_features_dm), "cycIF"),]
ATACseq_dm <- combined18_features_dm[str_detect(rownames(combined18_features_dm), "motifs"),]

max_JLV_nbr <- 20
cPCA_sparsity <- .9
moa <- mbpca(list(RNAseq_dm, RPPA_dm, GCP_dm, cycIF_dm,ATACseq_dm), ncomp = max_JLV_nbr, method = "globalScore", k = cPCA_sparsity,
             option = "lambda1", center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast", maxiter = 1000)
df <- tibble(eigenvalue = moa@eig, JLV = 1:length( moa@eig))
p_eigen <- ggplot(df, aes(x = JLV, y = eigenvalue)) +
  geom_path()  +
  geom_point() +
  labs(title = "Pseudoeigenvalues to optimize number of JLVs")
p_eigen

max_JLV_nbr <- 8
moa <- mbpca(list(RNAseq_dm, RPPA_dm, GCP_dm, cycIF_dm,ATACseq_dm), ncomp = max_JLV_nbr, method = "globalScore", k = cPCA_sparsity,
             option = "lambda1", center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast", maxiter = 1000)
scr <- moaScore(moa)
colnames(scr) <- str_replace(colnames(scr), "PC", "JLV")

##find correlations between the modules PCs
cor_mdd_cPCA <- cor(t(combined14_features_data_mean_dm), scr, method = "pearson")

cor_hm <- Heatmap(cor_mdd_cPCA,
                  column_title = paste("Correlation of MultiOmic CPCA\n to MDD modules\n sparsity:",cPCA_sparsity),
                  name = "pearson\ncorrelation",
                  row_names_gp = gpar(fontsize = 10),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  #split = clusters,
                  column_names_gp = gpar(fontsize = 8),
                  col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                  column_title_gp = gpar(fontsize = 8))
cor_hm

pdf(paste0("../plots/MCF10A_MDD_MoCluster_heatmaps_sparsity_",cPCA_sparsity,".pdf"), width = 5,height = 5)
print(p_eigen)
print(Heatmap(t(scr), 
              name = "CPCA Score",
              col = brewer.pal(9,"RdBu")[9:1],
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              column_names_gp = gpar(fontsize = 6)))
print(cor_hm)
print(hm_means_treatments)
res <- dev.off()

#Use PC2 loadings to compare 
#start with the features in the top20% based on weights in PC2
#Compare the top PC2 loading  to the features in cluster 3+4 which has a strong IFNG signal
mdd_34_features <- rownames(combined18_features_dm)[clusters %in% c(3,4)]

PC2_loadings <- PC2_loadings <- moa@loading[,2] %>%
  sort(decreasing = TRUE) %>%
  names %>%
  str_remove_all("_data_1")

PC2_top_features <- PC2_loadings[1:length(mdd_34_features)]

common_features <- PC2_top_features[PC2_top_features %in% mdd_34_features] 
