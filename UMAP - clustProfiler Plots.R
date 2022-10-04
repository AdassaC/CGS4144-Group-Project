# Install the Homo sapiens package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}
if (!("umap" %in% installed.packages())) {
  # Install umap package
  BiocManager::install("umap", update = FALSE)
}
if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}
if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}
if (!("ggnewscale" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages('ggnewscale')
}

# Attach the library
library(clusterProfiler)

# Package that contains MSigDB gene sets in tidy format
library(msigdbr)

# Attach the Organism library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

# Attach the tidyverse library
library(tidyverse)

# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

#Attach UMAP library for UMAP plotting
library(umap)

#Attach tidyr library
library(tidyr)

#Attach enrich plot library
library(enrichplot)

#Attach ggnewscale library
library(ggnewscale)

# Set seed for randomization
set.seed(12345)

#-----Define directories-----#

# Define the file path to the plots directory
plots_dir <- "Plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "Results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#-----Read data-----#
metaData <- readr::read_csv("Fixed Data/metadata_file.csv")
expression_df <- readr::read_csv("Fixed Data/dataset_file.csv") %>%
  #IP_expression_df <- readr::read_csv("Fixed Data/IP_data_file.csv") %>%
  
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metaData$Run)

# Check if this is in the same order
all.equal(colnames(expression_df), metaData$Run)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated HUGO SYMBOL IDs
mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)

# Display mapped_list
head(mapped_list)

# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list%>%
  tibble::enframe(name = "Ensembl", value = "Symbol") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Symbol)

#Drop null values
mapped_df[!is.na(mapped_df$Symbol),]

# Display mapped_df
head(mapped_df)

# Use the `summary()` function to show the distribution of HUGO SYMBOL values
# We need to use `as.factor()` here to get the count of unique values
# `maxsum = 10` limits the summary to 10 distinct values
summary(as.factor(mapped_df$Symbol), maxsum = 100)


#------------------------------------------- DENSITY PLOT ---------------------------------------#

log_scaled<-apply(expression_df[,2:70],1,log) 

maxVals<-apply(log_scaled[,2:70],1,max)
plot(maxVals)

minVals<-apply(log_scaled[,2:70],1,min)
plot(minVals) 

#count matrix geneXsample
ranges <- maxVals - minVals
d <- density(ranges) # returns the density data
plot(d) # plots the result


#------------------------------------- PERFORM DESeq2 ANALYSIS ------------------------------------------#
head(metaData$Group)

metaData <- metaData %>%
  # Let's get the cancer status from this variable
  dplyr::mutate(cancer_status = dplyr::case_when(
    stringr::str_detect(Group, "CTRL") ~ "normal",
    stringr::str_detect(Group, "OC") ~ "cancerous",
    stringr::str_detect(Group, "Other") ~ "neither"
    
  ))

dplyr::select(metaData, DESCRIPTION, Group, cancer_status)

str(metaData$cancer_status)

# Make cancer_status a factor and set the levels appropriately
metaData <- metaData %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    cancer_status = factor(cancer_status, levels = c("normal", "cancerous", "neither")),
    DESCRIPTION = factor(DESCRIPTION),
    Group = factor(Group)
  )

levels(metaData$cancer_status)

filtered_expression_df <- (expression_df[,2:70]) %>%
  dplyr::filter(rowSums(.) >= 100)

# round all expression counts
gene_matrix <- round(expression_df[,2:70])

ddset <- DESeqDataSetFromMatrix(
  # Here we supply normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metaData,
  # Supply our experimental variable to `design`
  design = ~cancer_status
)

#Generate DESeq object and DESeq2 analysis results
deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)

head(deseq_results)

# Create dataframe for DESeqResults 
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = pvalue < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

#-------------------------------------------- PCA PLOT ------------------------------------------#
dds_norm <- vst(ddset)
plotPCA(
  dds_norm,
  intgroup = "Group"
)

# We first have to save the results of the `plotPCA()` function for use with `ggplot2`
pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("Group"),
    returnData = TRUE # This argument tells R to return the PCA values
  )

#--------------------------------------------- UMAP PLOT ----------------------------------------#
normalized_counts <- assay(dds_norm) %>%
  t()

# Now perform UMAP on the normalized data
umap_results <- umap::umap(normalized_counts)

umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("Run") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(metaData, by = "Run")

umap_plot_df

# Plot using `ggplot()` function and save to an object
final_annotated_umap_plot <- ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    # plot points with different colors for each `refinebio_treatment` group
    color = cancer_status,
    # plot points with different shapes for each `refinebio_disease` group
  )
) +
  geom_point() # make a scatterplot

# Display the plot that we saved above
final_annotated_umap_plot

# Save plot using `ggsave()` function
# ggsave(
#   file.path(
#     plots_dir,
#     "OC_umap_plot.png" # Replace with a good file name for your plot
#   ),
#   plot = final_annotated_umap_plot
# )

#--------------------------------------------- clustProfiler GSEA ----------------------------------------#
#Display DESeq2 results in dataframe
head(deseq_df)

#Generate geneList
geneList = deseq_df$log2FoldChange
names(geneList) = deseq_df$Gene
geneList = sort(geneList, decreasing = TRUE)


#Display geneList
head(geneList)

#Perform enrichment analysis
CP_ego3 <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "CC",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 #seed = TRUE,
                 verbose      = FALSE)

#Visualize GSEA using directed acyclical graph
CP_network_graph <- goplot(ego3)
CP_network_graph

#Visualize GSEA using enrichment plot --> default
CP_enrichPlotDefault <- gseaplot(
  ego3,
  geneSetID = 1,
  title = "",
  #color.line = "#0d76ff"
)
CP_enrichPlotDefault

CP_dotPlot <- dotplot(ego3, showCategory=30) + ggtitle("dotplot for GSEA")
CP_dotPlot

#data(deseq_df)
xx <- compareCluster(deseq_df, 
                     fun="enrichGO",
                     OrgDb = org.Hs.eg.db, 
                     pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)
CP_clusterPlot <- emapplot(xx)
CP_clusterPlot

# #Visualize GSEA using enrichment plot --> default
# CP_enrichPlot2 <- gseaplot(
#   ego3,
#   geneSetID = 1,
#   by = "runningScore",
#   title = ego3$Description[2], #"HALLMARK_MYC_TARGETS_V2",
#   #color.line = "#0d76ff"
# )
# CP_enrichPlot2
# 
# #Visualize GSEA using enrichment plot --> default
# CP_enrichPlot3 <- gseaplot(
#   ego3,
#   geneSetID = 1,
#   by = "preranked",
#   title = ego3$Description[3], #"HALLMARK_MYC_TARGETS_V2",
#   #color.line = "#0d76ff"
# )
# CP_enrichPlot3

