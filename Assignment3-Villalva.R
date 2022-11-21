
# Install the Homo sapiens package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

## Install matrixStats package

#if(!("matrixStats" %in% installed.packages())){
#  install.packages("matrixStats", update = FALSE)
  
#}
install.packages("matrixStats")

# Installing Packages
install.packages("ClusterR")
install.packages("cluster")
install.packages("ggalluvial")


# Loading package
library(ClusterR)
library(cluster)

library(matrixStats)


# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

library(tidyverse)

library(factoextra)

library(pheatmap) 


library(RColorBrewer)

library(ggplot2)

library(GGally)


#-----Read Data-----#
metaData <- readr::read_csv("data/metadata_file.csv")
expression_df <- readr::read_csv("data/dataset_file.csv") 

#-------Variance-------#

#creates a col(71) called row_var in copy of expression_df
#row_var stores the row variance of each gene 

copy_expression_df <- expression_df
copy_expression_df$row_var = rowVars(as.matrix(expression_df[,c(2:70)]))
copy_expression_df$row_var

copy_expression_df <- copy_expression_df[order(-copy_expression_df$row_var),]
copy_expression_df <- copy_expression_df[c(1:5001),] #first 5001 rows 1 for the samples, 5000 for the rows

#-------K_means Heatmap-------#

# Structure 
str(copy_expression_df)
set.seed(240) 

kmData <- copy_expression_df[,-1]
kmData <- kmData[,1:69] #take away row var
kmData <-t(kmData)

#3 clusters
km = kmeans(kmData, 3)

kmData2 <- cbind(kmData,km$cluster)

o <- order(kmData2[, 70])
o
kmData2 <- kmData2[o, ]

fviz_cluster(km, data =kmData)


heatmap(kmData)#, cluster_rows=F,cluster_cols=F, col=brewer.pal(12,"Set3"),border_color=NA, fontsize = 4)


#-------KMeans 5000 genes ------#
kmData <- copy_expression_df[,-1] #remove gene labels 
kmData <- kmData[,1:69]

#transpose matrix
kmData <-t(kmData)

#3 clusters for kmeans 
km = kmeans(kmData, 3)

#show kmeans heatmap
#add dendograms
fviz_cluster(km, data=kmData, main ="5000 genes with 3 clusters")

#-------KMeans 10 genes-------#

copy_expression_df <- expression_df

#do row variance to find top genes
copy_expression_df$row_var = rowVars(as.matrix(expression_df[,c(2:70)]))
copy_expression_df$row_var

#first 11 rows 1 for the samples, 10 for the rows
copy_expression_df <- copy_expression_df[order(-copy_expression_df$row_var),]
copy_expression_df <- copy_expression_df[c(1:11),] 


kmData10 <- copy_expression_df[,-1]
kmData10 <- kmData10[,1:69]#remove row variance

#transpose matrix
kmData10<-t(kmData10)

km10 = kmeans(kmData10, 3)

fviz_cluster(km10, data=kmData10, main = "10 genes with 3 clusters")
             

#-------KMeans 100 genes ------#

copy_expression_df <- expression_df
copy_expression_df$row_var = rowVars(as.matrix(expression_df[,c(2:70)]))
copy_expression_df$row_var

copy_expression_df <- copy_expression_df[order(-copy_expression_df$row_var),]
copy_expression_df <- copy_expression_df[c(1:101),] #first 11 rows 1 for the samples, 10 for the rows

kmData100 <- copy_expression_df[,-1]
kmData100 <- kmData100[,1:69]#remove row variance

kmData100 <- t(kmData100) #transpose matrix

km100 = kmeans(kmData100, 3)

fviz_cluster(km100, data =kmData100,main = "100 genes with 3 clusters")

#-------KMeans 1000 genes ------#

copy_expression_df <- expression_df
copy_expression_df$row_var = rowVars(as.matrix(expression_df[,c(2:70)]))
copy_expression_df$row_var

copy_expression_df <- copy_expression_df[order(-copy_expression_df$row_var),]
copy_expression_df <- copy_expression_df[c(1:1001),] #first 11 rows 1 for the samples, 10 for the rows


kmData1k <- copy_expression_df[,-1]
kmData1k <- kmData1k[,1:69]#remove row variance
kmData1k <- t(kmData1k) #transpose matrix

km1k = kmeans(kmData1k, 3)

fviz_cluster(km1k, data =kmData1k,main = "1000 genes with 3 clusters")

#-------KMeans 10000 genes ------#

copy_expression_df <- expression_df
copy_expression_df$row_var = rowVars(as.matrix(expression_df[,c(2:70)]))
copy_expression_df$row_var

copy_expression_df <- copy_expression_df[order(-copy_expression_df$row_var),]
copy_expression_df <- copy_expression_df[c(1:10001),] #first 11 rows 1 for the samples, 10 for the rows

kmData10k <- copy_expression_df[,-1]
kmData10k <- kmData10k[,1:69] # remove row variance
kmData10k <- t(kmData10k)

km10k = kmeans(kmData10k, 3)

fviz_cluster(km10k, data =kmData10k, main = "10000 genes with 3 clusters")

#-------Sankey Plot-------#
#s2 <- km10$cluster

cdata <- data.frame(
  genes10 = km10$cluster, 
  genes100 = km100$cluster,
  genes1k = km1k$cluster,
  genes10k = km10k$cluster,
  genes5000 = km$cluster
  
)
ggplot(as.data.frame(cdata),
       aes(y = cdata[,1],
           axis1 = genes10, axis2 = genes100, axis3 = genes1k, axis4 = genes10k, axis5 = genes5000)) +
  geom_alluvium(aes(fill = Class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("genes10", "genes100", "genes1k", "genes10k", "genes5000")) +
  coord_flip() +
  
  ggtitle("Kmeans clustering by different number of genes")


#-------Statistics-------#
copy_expression_df <- expression_df
entire <- copy_expression_df[,-1]
entire <- t(entire)

kmEntire = kmeans(entire, 3)

chisq1 <- chisq.test(table(kmEntire$cluster, km$cluster))
chisq1

chisq2 <- chisq.test(table(kmEntire$cluster,km10$cluster))
chisq2

chisq3 <- chisq.test(table(kmEntire$cluster,km100$cluster))
chisq3

chisq4 <- chisq.test(table(kmEntire$cluster,km1k$cluster))
chisq4

chisq5 <- chisq.test(table(kmEntire$cluster,km10k$cluster))
chisq5

#---Adjust p_value---#
#get all pvalues of the chi tests 
  p = c(chisq$p.value,
        chisq1$p.value,
        chisq2$p.value,
        chisq3$p.value,
        chisq4$p.value,
        chisq5$p.value)
  


adjusted <-p.adjust(p, method = "fdr", n = length(p))

chisqa <- chisq
chisqa$p.value <- adjusted[1]

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

#tests<-data.frame(reg = c(chisq,chisq1) )

pvalues = data.frame(
  regular = c(chisq$p.value,
      chisq1$p.value,
      chisq2$p.value,
      chisq3$p.value,
      chisq4$p.value,
      chisq5$p.value),
  adjusted = c(chisqa$p.value,
              chisq1a$p.value,
              chisq2a$p.value,
              chisq3a$p.value,
              chisq4a$p.value,
              chisq5a$p.value)
)

pairs(pvalues[1:2],pch=21)
pairs(cdata)
