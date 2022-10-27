# Install the Homo sapiens package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
if (!("multiClust" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("multiClust", update = FALSE)
}
if (!("preprocessCore" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("preprocessCore", update = FALSE)
}
if (!("gplots" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("gplots", update = FALSE)
}
if (!("dendextend" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("dendextend", update = FALSE)
}
if (!("graphics" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("graphics", update = FALSE)
}
if (!("grDevices" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("grDevices", update = FALSE)
}
if (!("amap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("amap", update = FALSE)
}
if (!("survival" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("survival", update = FALSE)
}
if (!("ctc" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ctc", update = FALSE)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  # Install this package if it isn't installed yet
  install.packages("BiocManager")
  BiocManager::install("scater")
}
if(!("matrixStats" %in% installed.packages())){
  BioCManager::install.packages("matrixStats", update = FALSE)
}
if(!("dendextend" %in% installed.packages())){
  BioCManager::install.packages("dendextend", update = FALSE)
}
if (!("pheatmap" %in% installed.packages())) {
  # Install pheatmap
  install.packages("pheatmap", update = FALSE)
}
if (!("DESeq2" %in% installed.packages())) {
  # Install DESeq2
  BiocManager::install("DESeq2", update = FALSE)
}

#Load packages and libraries
library(ggplot2)
library(org.Hs.eg.db)
library(magrittr) # We will need this so we can use the pipe: %>%)
library(matrixStats)
library(dendextend)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(GGally)
library(pheatmap)
library(DESeq2)

#-----Read Data-----#
metaData <- readr::read_csv("Fixed Data 2/metadata_file.csv")
expression_df <- readr::read_csv("Fixed Data 2/dataset_file.csv")
  
#--------------------------------------------SUBSETTING THE DATA: 5000 MOST VARIABLE GENES----------------------------------------------------#

#creates a col(71) called row_var in copy of expression_df
#row_var stores the row variance of each gene 

expressionVar_df <- expression_df
expressionVar_df$row_var <- rowVars(as.matrix(expression_df[,c(2:70)]))

sorted_expressionVar <- expressionVar_df[order(expressionVar_df$row_var, decreasing = TRUE),]
subset_df <- sorted_expressionVar[c(1:5001),c(1:70)] #first 50001 rows 1 for the samples, 50000 for the rows

#--------------------------------------------PERFORMING HIERARCHICAL CLUSTERING----------------------------------------------------#

# create sankey plot table 
# make columns and row dendograms

#------------------------------------H-Clust on Subset Data: 5000 Genes (sample based)------------------------------#

# Create a matrix from table of counts
hclust_sample_matrix <- subset_df %>%
  select(-Gene) %>%
  # coerce to a matrix
  as.matrix() %>%
  # transpose the matrix so that rows = samples and columns = variables
  t()

#Assign column names to matrix
colnames(hclust_sample_matrix) <- subset_df$Gene

#Scale expression of genes in matrix to obtain z-scores (transpose -> scale -> transpose)
hclust_sample_matrix <- hclust_sample_matrix %>%
  scale()

#calculate the distance between each gene row in matrix:
sample_dist <- dist(hclust_sample_matrix)

#Perform the hierarchical clustering based on distance
h5000 <- hclust(sample_dist, method = "complete")
dend <- as.dendrogram(h5000)
plot(dend, main = "Cluster Dendrogram", ylab = "Height")

#-----------------------------------------HEATMAP--------------------------------------#
annotation_df <- metaData %>%
  # select only the columns we need for annotation
  dplyr::select(
    Run,
    Group
  )
# matHM <- t(hclust_sample_matrix)
# matHM <- scale(matHM)
# 
# dds <- DESeqDataSetFromMatrix(
#   countData = matHM, # the counts values for all samples
#   colData = annotation_df, # annotation data for the samples
#   design = ~1 # Here we are not specifying a model
#   # Replace with an appropriate design variable for your analysis
# )
# 
# dds_norm <- rlog(dds)
# 
#   #tibble::column_to_rownames("Run")
# # The `pheatmap()` function requires that the row names of our annotation
# # data frame match the column names of our `DESeaDataSet` object
# 
# heatmap_annotated <- pheatmap(
#     matHM,
#     cluster_rows = TRUE,
#     cluster_cols = TRUE
#     show_rownames = FALSE,
#     # show_colnames = TRUE,
#     annotation_col = annotation_df
#     # main = "Annotated Heatmap for Hierarchical Clustering"
#     # #colorRampPalette(c("blue","red","white"))
#   )

#Plot heatmap
#heatmap(hclust_sample_matrix, main = "Annotated Heatmap for Hierarchical Clustering"))


#-------------K CLUSTER GROUPS--------------# 

# k = 2
cluster5000_k2 <- data.frame(cutree(h5000, k = 2))
cutree(dend, k = 2) # on dendrogram
d2 = color_branches(dend, k=2)
plot2 <- plot(d2, main = "Cluster Dendrogram, k = 2", ylab = "Height")

# rect.hclust(h5000, k = 4, which = NULL, x = NULL, h = NULL,
#             border = 2)

# k = 3
cluster5000_k3 <- data.frame(cutree(h5000, k = 3))
cutree(dend, k = 3) # on dendrogram
d3 = color_branches(dend, k=3)
plot3 <- plot(d3, main = "Cluster Dendrogram, k = 3", ylab = "Height")

# k = 4
cluster5000_k4 <- data.frame(cutree(h5000, k = 4))
cutree(dend, k = 4) # on dendrogram
d4 = color_branches(dend, k=4)
plot4 <- plot(d4, main = "Cluster Dendrogram, k = 4", ylab = "Height")

#k = 5
cluster5000_k5 <- data.frame(cutree(h5000, k = 5))
cutree(dend, k = 5) # on dendrogram
d5 = color_branches(dend, k=5)
plot5 <- plot(d5, main = "Cluster Dendrogram, k = 5", ylab = "Height")

# k = 7
cluster5000_k7 <- data.frame(cutree(h5000, k = 7))
cutree(dend, k = 7) # on dendrogram
d7 = color_branches(dend, k=7)
plot7 <- plot(d7, main = "Cluster Dendrogram, k = 7", ylab = "Height")

#----------------------------------10 Genes--------------------------------------#
sorted_10_df <- subset_df[c(1:11),]

# Create a matrix from table of counts
hclust_sample_matrix <- sorted_10_df %>% 
  select(-Gene) %>%
  # coerce to a matrix
  as.matrix() %>%
  # transpose the matrix so that rows = samples and columns = variables
  t()

#Assign column names to matrix
colnames(hclust_sample_matrix) <- sorted_10_df$Gene

#Scale expression of genes in matrix to obtain z-scores (transpose -> scale -> transpose)
hclust_sample_matrix <- hclust_sample_matrix %>%
  scale()

#calculate the distance between each gene row in matrix:
gene_dist_10 <- dist(hclust_sample_matrix)

#Perform the hierarchical clustering based on distance
h10 <- hclust(gene_dist_10, method = "complete")
dend <- as.dendrogram(h10)
plot(dend)

pairs(gene_dist_10)

#-------------K CLUSTER GROUPS--------------# 

# k = 2
cluster10_k2 <- data.frame(cutree(h10, k = 2))
cutree(dend, k = 2) # on dendrogram
d2 = color_branches(dend, k=2)
plot(d2)

# rect.hclust(h10, k = 4, which = NULL, x = NULL, h = NULL,
#             border = 2)

# k = 3
cluster10_k3 <- data.frame(cutree(h10, k = 3))
cutree(dend, k = 3) # on dendrogram
d3 = color_branches(dend, k=3)
plot(d3, main = "Cluster Dendrogram, k = 3", ylab = "Height")

# k = 4
cluster10_k4 <- data.frame(cutree(h10, k = 4))
cutree(dend, k = 4) # on dendrogram
d4 = color_branches(dend, k=4)
plot(d4)

#k = 5
cluster10_k5 <- data.frame(cutree(h10, k = 5))
cutree(dend, k = 5) # on dendrogram
d5 = color_branches(dend, k=5)
plot(d5)

# k = 7
cluster10_k7 <- data.frame(cutree(h10, k = 7))
cutree(dend, k = 7) # on dendrogram
d7 = color_branches(dend, k=7)
plot(d7)

#----------------------------------100 Genes--------------------------------------#
sorted_100_df <- subset_df[c(1:101),]

# Create a matrix from table of counts
hclust_sample_matrix <- sorted_100_df %>% 
  select(-Gene) %>%
  # coerce to a matrix
  as.matrix() %>%
  # transpose the matrix so that rows = samples and columns = variables
  t()

#Assign column names to matrix
colnames(hclust_sample_matrix) <- sorted_100_df$Gene

#Scale expression of genes in matrix to obtain z-scores (transpose -> scale -> transpose)
hclust_sample_matrix <- hclust_sample_matrix %>%
  scale()

#calculate the distance between each gene row in matrix:
gene_dist_100 <- dist(hclust_sample_matrix)

#Perform the hierarchical clustering based on distance
h100 <- hclust(gene_dist_100, method = "complete")
dend <- as.dendrogram(h100)
plot(dend)


#-------------K CLUSTER GROUPS--------------# 

# k = 2
cluster100_k2 <- data.frame(cutree(h100, k = 2))
cutree(dend, k = 2) # on dendrogram
d2 = color_branches(dend, k=2)
plot(d2)

# rect.hclust(h100, k = 4, which = NULL, x = NULL, h = NULL,
#             border = 2)

# k = 3
cluster100_k3 <- data.frame(cutree(h100, k = 3))
cutree(dend, k = 3) # on dendrogram
d3 = color_branches(dend, k=3)
plot(d3, main = "Cluster Dendrogram, k = 3", ylab = "Height")

# k = 4
cluster100_k4 <- data.frame(cutree(h100, k = 4))
cutree(dend, k = 4) # on dendrogram
d4 = color_branches(dend, k=4)
plot(d4)

#k = 5
cluster100_k5 <- data.frame(cutree(h100, k = 5))
cutree(dend, k = 5) # on dendrogram
d5 = color_branches(dend, k=5)
plot(d5)

# k = 7
cluster100_k7 <- data.frame(cutree(h100, k = 7))
cutree(dend, k = 7) # on dendrogram
d7 = color_branches(dend, k=7)
plot(d7)

#----------------------------------1000 Genes--------------------------------------#
sorted_1000_df <- subset_df[c(1:1001),]

# Create a matrix from table of counts
hclust_sample_matrix <- sorted_1000_df %>% 
  select(-Gene) %>%
  # coerce to a matrix
  as.matrix() %>%
  # transpose the matrix so that rows = samples and columns = variables
  t()

#Assign column names to matrix
colnames(hclust_sample_matrix) <- sorted_1000_df$Gene

#Scale expression of genes in matrix to obtain z-scores (transpose -> scale -> transpose)
hclust_sample_matrix <- hclust_sample_matrix %>%
  scale()

#calculate the distance between each gene row in matrix:
gene_dist_1000 <- dist(hclust_sample_matrix)

#Perform the hierarchical clustering based on distance
h1000 <- hclust(gene_dist_1000, method = "complete")
dend <- as.dendrogram(h1000)
plot(dend)

#-------------K CLUSTER GROUPS--------------# 

# k = 2
cluster1000_k2 <- data.frame(cutree(h100, k = 2))
cutree(dend, k = 2) # on dendrogram
d2 = color_branches(dend, k=2)
plot(d2)

# rect.hclust(h100, k = 4, which = NULL, x = NULL, h = NULL,
#             border = 2)

# k = 3
cluster1000_k3 <- data.frame(cutree(h100, k = 3))
cutree(dend, k = 3) # on dendrogram
d3 = color_branches(dend, k=3)
plot(d3, main = "Cluster Dendrogram, k = 3", ylab = "Height")

# k = 4
cluster1000_k4 <- data.frame(cutree(h100, k = 4))
cutree(dend, k = 4) # on dendrogram
d4 = color_branches(dend, k=4)
plot(d4)

#k = 5
cluster1000_k5 <- data.frame(cutree(h100, k = 5))
cutree(dend, k = 5) # on dendrogram
d5 = color_branches(dend, k=5)
plot(d5)

# k = 7
cluster1000_k7 <- data.frame(cutree(h100, k = 7))
cutree(dend, k = 7) # on dendrogram
d7 = color_branches(dend, k=7)
plot(d7)

#----------------------------------10000 Genes--------------------------------------#
sorted_10000_df <- subset_df[c(1:10001),]

# Create a matrix from table of counts
hclust_sample_matrix <- sorted_10000_df %>% 
  select(-Gene) %>%
  # coerce to a matrix
  as.matrix() %>%
  # transpose the matrix so that rows = samples and columns = variables
  t()

#Assign column names to matrix
colnames(hclust_sample_matrix) <- sorted_10000_df$Gene

#Scale expression of genes in matrix to obtain z-scores (transpose -> scale -> transpose)
hclust_sample_matrix <- hclust_sample_matrix %>%
  scale()

#calculate the distance between each gene row in matrix:
gene_dist_10000 <- dist(hclust_sample_matrix)

#Perform the hierarchical clustering based on distance
h10000 <- hclust(gene_dist_10000, method = "complete")
dend <- as.dendrogram(h10000)
plot(dend)

#-------------K CLUSTER GROUPS--------------# 

# k = 2
cluster10000_k2 <- data.frame(cutree(h100, k = 2))
cutree(dend, k = 2) # on dendrogram
d2 = color_branches(dend, k=2)
plot(d2)

# rect.hclust(h100, k = 4, which = NULL, x = NULL, h = NULL,
#             border = 2)

# k = 3
cluster10000_k3 <- data.frame(cutree(h100, k = 3))
cutree(dend, k = 3) # on dendrogram
d3 = color_branches(dend, k=3)
plot(d3, main = "Cluster Dendrogram, k = 3", ylab = "Height")

# k = 4
cluster10000_k4 <- data.frame(cutree(h100, k = 4))
cutree(dend, k = 4) # on dendrogram
d4 = color_branches(dend, k=4)
plot(d4)

#k = 5
cluster10000_k5 <- data.frame(cutree(h100, k = 5))
cutree(dend, k = 5) # on dendrogram
d5 = color_branches(dend, k=5)
plot(d5)

# k = 7
cluster10000_k7 <- data.frame(cutree(h100, k = 7))
cutree(dend, k = 7) # on dendrogram
d7 = color_branches(dend, k=7)
plot(d7)


#----------------------------------Alluvial Plot--------------------------------------#
library(ggalluvial)
clustersNums <- c(1,2,3)
clustersCategories <- as.data.frame(clustersNums)

sankeyTable <- as.data.frame((cbind(cluster5000_k3, cluster10_k3, cluster100_k3, cluster1000_k3, cluster10000_k3, annotation_df$Run, clustersCategories)))
sankeyTable <- setNames(sankeyTable,c("Clusters from 5000 Genes", "Clusters from 10 Genes","Clusters from 100 Genes", "Clusters from 1000 Genes", "Clusters from 10000 Genes", "Samples", "Cluster Assig"))

ggplot(as.data.frame(sankeyTable),
       aes(axis1 = "Clusters from 10 Genes", axis2 = "Clusters from 100 Genes", axis3 = "Clusters from 1000 Genes", axis4 = "Clusters from 10000 Genes")) +
  geom_alluvium(aes(fill = "Samples"), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Clusters"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Clustering of Samples by Number of Genes")

#----------------------------------Stats--------------------------------------#

# hc <- as.hclust(merge(merge(merge(merge(
#   as.dendrogram(gene_hclust), as.dendrogram(h10)), as.dendrogram(h100)),
#   as.dendrogram(h1000)), as.dendrogram(h10000)))

# #-------Statistics-------#
groupCluster <- annotation_df
groupCluster$Group <-replace(groupCluster$Group,  groupCluster$Group == "OC", 1)
groupCluster$Group <-replace(groupCluster$Group,  groupCluster$Group == "CTRL", 2)
groupCluster$Group <-replace(groupCluster$Group,  groupCluster$Group == "Other", 3)
originalGroups <- as.matrix(groupCluster[1:69, 2])
rownames(originalGroups) <- groupCluster$Run

comparison1

chisq1 <- chisq.test(cluster5000_k3, originalGroups)
p = c(chisq1$p.value)

# cluster5k_k3$Group <-replace(cluster5k_k_3$Group,  cluster5k_k_3$Group == 1, 1)
# cluster5k_k3$Group <-replace(cluster5k_k_3$Group,  cluster5k_k_3$Group == 2, 2)
# cluster5k_k3$Group <-replace(cluster5k_k_3$Group,  cluster5k_k_3$Group == 3, 3)


chisq2 <- chisq.test(cluster10_k3, originalGroups)
p2 = c(chisq2$p.value)

chisq3 <- as.data.frame(chisq.test(cluster100_k3, originalGroups))
p3 = c(chisq3$p.value)

chisq4 <- as.data.frame(chisq.test(cluster1000_k3, originalGroups))
p4 = c(chisq4$p.value)

chisq5 <- as.data.frame(chisq.test(cluster10000_k3, originalGroups))
p5 = c(chisq5$p.value)


#-----------------------Table of p_value-----------------------------------#
pValues = c(chisq1$p.value,
      chisq2$p.value,
      chisq3$p.value,
      chisq4$p.value,
      chisq5$p.value)

adjusted <-p.adjust(pValues, method = "fdr", n = length(p))

chisq1a <- chisq1
chisq1a$p.value <- adjusted[2]

chisq2a <- chisq2
chisq2a$p.value <- adjusted[3]

chisq3a <- chisq3
chisq3a$p.value <- adjusted[4]

chisq4a <- chisq4
chisq4a$p.value <- adjusted[5]

chisq5a <- chisq5
chisq5a$p.value <- adjusted[6]

adj_pValues <- data.frame(
  regular = c(chisq1$p.value,
              chisq2$p.value,
              chisq3$p.value,
              chisq4$p.value,
              chisq5$p.value),
  adjusted = c(chisq1a$p.value,
               chisq2a$p.value,
               chisq3a$p.value,
               chisq4a$p.value,
               chisq5a$p.value)
)

pairs(comparison2)
