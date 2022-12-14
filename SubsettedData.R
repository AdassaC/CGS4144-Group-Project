if (!requireNamespace("BiocManager", quietly = TRUE)) {
  # Install this package if it isn't installed yet
  install.packages("BiocManager")
  BiocManager::install("scater")
}
if(!("matrixStats" %in% installed.packages())){
  BioCManager::install.packages("matrixStats", update = FALSE)
}

#Load packages and libraries
library(org.Hs.eg.db)
library(magrittr) # We will need this so we can use the pipe: %>%)
library(matrixStats)
library(tidyverse)


#-----Read Data-----#
metaData <- readr::read_csv("Fixed Data 2/metadata_file.csv")
expression_df <- readr::read_csv("Fixed Data 2/dataset_file.csv") %>%
  
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metaData$Run)

# Check if this is in the same order
all.equal(colnames(expression_df), metaData$Run)

#wait to run
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

#--------------------------------------------SUBSETTING THE DATA: 5000 MOST VARIABLE GENES----------------------------------------------------#

#creates a col(71) called row_var in copy of expression_df
#row_var stores the row variance of each gene 

expressionVar_df <- expression_df
expressionVar_df$row_var <- rowVars(as.matrix(expression_df[,c(2:70)]))

sorted_expressionVar <- expressionVar_df[order(expressionVar_df$row_var, decreasing = TRUE),]
subset_df <- sorted_expressionVar[c(1:5001),] #first 50001 rows 1 for the samples, 50000 for the rows
