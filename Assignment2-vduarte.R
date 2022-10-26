results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

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
if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}
if (!("umap" %in% installed.packages())) {
  # Install umap package
  BiocManager::install("umap", update = FALSE)
}
if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicSuperSignature")


browseVignettes("GenomicSuperSignature")


BiocManager::install("bcellViper")

## ----results="hide", message=FALSE, warning=FALSE-----------------------------
library(GenomicSuperSignature)
library(bcellViper)

## ----load_model---------------------------------------------------------------
RAVmodel <- getModel("PLIERpriors", load=TRUE)
RAVmodel

version(RAVmodel)

## ----message=FALSE, warning=FALSE---------------------------------------------
data(bcellViper)
dset

## -----------------------------------------------------------------------------
val_all <- validate(dset, RAVmodel)
head(val_all)

## ----out.height="80%", out.width="80%", message=FALSE, warning=FALSE----------
heatmapTable(val_all, RAVmodel, num.out = 5, swCutoff = 0)

## ----out.height="75%", out.width="75%", plotValidate_function-----------------
plotValidate(val_all, interactive = FALSE)

## -----------------------------------------------------------------------------
validated_ind <- validatedSignatures(val_all, RAVmodel, num.out = 3, 
                                     swCutoff = 0, indexOnly = TRUE)
validated_ind

## ----out.height="60%", out.width="60%"----------------------------------------
set.seed(1) # only if you want to reproduce identical display of the same words
drawWordcloud(RAVmodel, validated_ind[1])
drawWordcloud(RAVmodel, validated_ind[2])
drawWordcloud(RAVmodel, validated_ind[3])

## -----------------------------------------------------------------------------
RAVnum <- validated_ind[2]  # RAV1139
res <- gsea(RAVmodel)[[RAVnum]]   
head(res)

## -----------------------------------------------------------------------------
findSignature(RAVmodel, "Bcell")

## -----------------------------------------------------------------------------
findSignature(RAVmodel, "Bcell", k = 5)

## -----------------------------------------------------------------------------
findKeywordInRAV(RAVmodel, "Bcell", ind = 695)

## -----------------------------------------------------------------------------
subsetEnrichedPathways(RAVmodel, ind = RAVnum, n = 3, both = TRUE)
subsetEnrichedPathways(RAVmodel, ind = 695, n = 3, both = TRUE)
subsetEnrichedPathways(RAVmodel, ind = 1994, n = 3, both = TRUE)

## -----------------------------------------------------------------------------
findStudiesInCluster(RAVmodel, validated_ind[2])

## -----------------------------------------------------------------------------
findStudiesInCluster(RAVmodel, validated_ind[2], studyTitle = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()


BiocManager::install("enrichplot")
# The following initializes usage of Bioc devel
#BiocManager::install(version='devel')

#BiocManager::install("M3C")
#install.packages("Rtsne")

install.packages("gprofiler2")


# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

library(tidyverse)

# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

#library for t-sne plot
library(M3C)

library(Rtsne)

library(gprofiler2)

library(enrichplot)
library(clusterProfiler)
library(msigdbr)




set.seed(12345)


#-----Read Data-----#
metaData <- readr::read_csv("data/metadata_file.csv")
expression_df <- readr::read_csv("data/dataset_file.csv") %>%
  
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

expression_df <- readr::read_csv("data/dataset_file.csv") %>%
  
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

#----pca Plot----#
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

#-----t_NSE plot-----#

set.seed(12345)
my_data <- expression_df

tsne_results <- Rtsne(pca_results,
                      perplexity = 10,
                      eta = 1000,
                      max_iter = 1000,
                      check_duplicates = FALSE)

Y <- as.data.frame(tsne_results$Y)
# load your omic data here as mydata
#library(M3C)
cancer_col <- metaData$cancer_status
gene_col <- expression_df$Gene

ggplot(Y,
       aes(x = V1, y = V2, col = cancer_col)) +
  geom_point()+
  labs(x = "tSNE = 1",
       y = "tSNE = 2",
       color = "cell types")+
  theme_classic()




#-----Volcano Plot-----#

head(metaData$Group)

metaData <- metaData %>%
  # Let's get the RPL10 mutation status from this variable
  dplyr::mutate(cancer_status = dplyr::case_when(
    stringr::str_detect(Group, "CTRL") ~ "normal",
    stringr::str_detect(Group, "OC") ~ "cancerous",
    stringr::str_detect(Group, "Other") ~ "neither"
    
  ))

dplyr::select(metaData, Group, cancer_status)

str(metaData$cancer_status)

# Make mutation_status a factor and set the levels appropriately
metaData <- metaData %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    cancer_status = factor(cancer_status, levels = c("normal", "cancerous","neither"))
  )

levels(metaData$cancer_status)

expression_df <- readr::read_csv("data/dataset_file.csv") %>%
  
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
  design = ~cancer_status
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

readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "data_diff_expr_results.tsv" # Replace with a relevant output file name
  )
)

# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)



#-----gprofiler2----#
gene_names<-mapped_df[,2]

gostres <- gost(query = c(gene_names), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

# The result is a named list where “result” is a data.frame with the enrichment analysis results
# and “meta” containing a named list with all the metadata for the query.

names(gostres)

head(gostres$result, 3)
names(gostres$meta)

gostres2 <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL)

head(gostres2$result, 3)
gostres_link <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                     as_short_link = TRUE)

multi_gostres1 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340", 
                                                 "GO:0005005", "ENSG00000156103", "NLRP1"),
                                    "chromY" = c("Y:1:10000000", "rs17396340", 
                                                 "GO:0005005", "ENSG00000156103", "NLRP1")), 
                       multi_query = FALSE)



head(multi_gostres1$result, 3)

multi_gostres2 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340",
                                                 "GO:0005005", "ENSG00000156103", "NLRP1"),
                                    "chromY" = c("Y:1:10000000", "rs17396340", 
                                                 "GO:0005005", "ENSG00000156103", "NLRP1")), 
                       multi_query = TRUE)


head(multi_gostres2$result, 3)

gostplot(gostres, capped = TRUE, interactive = TRUE)




head(gostres$result)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p


