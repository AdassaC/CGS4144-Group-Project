# Install the Homo sapiens package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

metaData <- readr::read_csv("data/metadata_file.csv")
expression_df <- readr::read_csv("data/dataset_file.csv") %>%
  
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

library(tidyverse)

my_data<-apply(expression_df[2:70],2,log)
plot(my_data)
  
# my_data <- expression_df
# my_data %>% select(2:70)
# maxVals<-apply(my_data,2,max)
# 
# plot(maxVals)

# data %>% ggplot(aes(x=Gene)) + geom_density()

# expression_df %>%
#   pivot_longer(cols=2:70) %>%
#   #filter(ID == "ABCA8") %>%
#   ggplot(aes(x=Gene, y = SRR)) +
#   geom_point()



# plot(head(expression_df$SRR12705190, 57736))
# plot(head(log(expression_df$SRR12705190), 57736), head(log(expression_df$SRR12705235), 57736))
# plot(head(log(expression_df$SRR12705190), 57736), head(log(expression_df$SRR12705235), 57736))


#plot(expression_df$Gene, expression_df)

#Initialize x and y
#x <- 1:69
#y <- 1:57736

#Plot x and y
#plot(x,y,)

#Log scale
#plot(x, y, log="y")
