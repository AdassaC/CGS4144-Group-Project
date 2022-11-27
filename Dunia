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

# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

library(tidyverse)

# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

set.seed(12345)


#-----Read Data-----#
metaData <- readr::read_csv("metadata_file.csv")
expression_df <- readr::read_csv("dataset_file.csv") %>%
  
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

head(mapped_list)

# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Symbol") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Symbol)

head(mapped_df)

# Use the `summary()` function to show the distribution of HUGO SYMBOL values
# We need to use `as.factor()` here to get the count of unique values
# `maxsum = 10` limits the summary to 10 distinct values
summary(as.factor(mapped_df$Symbol), maxsum = 100)


#-----Density Plot-----#

log_scaled<-apply(expression_df[,2:70],1,log) 
#plot(log_scaled)

maxVals<-apply(log_scaled[,2:70],1,max)
plot(maxVals)

minVals<-apply(log_scaled[,2:70],1,min)
plot(minVals) 

#count matrix geneXsample
ranges <- maxVals - minVals
d <- density(ranges) # returns the density data
plot(d) # plots the result


#-----PlotPCA-----#
head(metaData$Group)

metaData <- metaData %>%
  # Let's get the RPL10 mutation status from this variable
  dplyr::mutate(cancer_status = dplyr::case_when(
    stringr::str_detect(Group, "CTRL") ~ "normal",
    stringr::str_detect(Group, "OC") ~ "cancerous",
    #stringr::str_detect(Group, "Other") ~ "neither"
    
  ))

dplyr::select(metaData, Group, cancer_status)

str(metaData$cancer_status)

# Make mutation_status a factor and set the levels appropriately
metaData <- metaData %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    cancer_status = factor(cancer_status, levels = c("normal", "cancerous"))
  )

levels(metaData$cancer_status)
expression_df <- readr::read_csv("dataset_file.csv") %>%
  
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metaData$Run)

# Check if this is in the same order
all.equal(colnames(expression_df), metaData$Run)



filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 100)

# round all expression counts
gene_matrix <- round(expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metaData,
  # Supply our experimental variable to `design`
  design = ~1 #cancer_status
)

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


deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)

head(deseq_results)

deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 3, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

head(deseq_results)

# this is of class DESeqResults -- we want a data frame
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

head(deseq_df)

pcaPlot<-plotCounts(ddset, gene = "ENSG00000000003", intgroup = "cancer_status")

pcaPlot()

# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
volcano_plot

ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "volcano_plot.png")
) # Replace with a plot name relevant to your data

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "SRP123625_volcano_plot.png")
) # Replace with a plot name relevant to your data
readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "SRP123625_diff_expr_results.tsv" # Replace with a relevant output file name
  )
)

