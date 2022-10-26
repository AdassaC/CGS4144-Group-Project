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

#Load packages and libraries
library(ggplot2)
library(org.Hs.eg.db)
library(magrittr) # We will need this so we can use the pipe: %>%)
library(matrixStats)
library(dendextend)
library(tidyverse)


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

#-----------H-Clust on Subset Data: 5000 Genes (sample based)----------#

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

#Plot using dendogram
plot(h5000, labels = FALSE)
abline(h = 150, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

#Plot heatmap
heatmap(hclust_sample_matrix)

#-----------H-Clust on Subset Data: 10 Genes----------#
sorted_10_df <- sorted_expressionVar[c(1:11),]

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

h10 <- hclust(gene_dist_10, method = "complete")

#Plot using dendogram
plot(h10, labels = FALSE)
#abline(h = 150, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

#-----------100 Genes----------#
sorted_100_df <- sorted_expressionVar[c(1:101),]

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

#Plot using dendogram
plot(h100, labels = FALSE)
# abline(h = 12, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

#-----------1000 Genes----------#
sorted_1000_df <- sorted_expressionVar[c(1:1001),]

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

#Plot using dendogram
plot(h1000, labels = FALSE)
# abline(h = 12, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

#-----------10000 Genes----------#
sorted_10000_df <- sorted_expressionVar[c(1:10001),]

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

plot(h10000, labels = FALSE)
# abline(h = 12, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

# hc <- as.hclust(merge(merge(merge(merge(
#   as.dendrogram(gene_hclust), as.dendrogram(h10)), as.dendrogram(h100)), 
#   as.dendrogram(h1000)), as.dendrogram(h10000)))

