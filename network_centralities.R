#!/usr/bin/env Rscript
# usage: Rsript network_centralities.R input.csv sample_name outdir/.

suppressPackageStartupMessages({
    
    library(igraph)
    library(tidyverse)

})
options(warn=-1)
set.seed(1234)

# Necessary arguments: gene-gene association file, samplename
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  
  # Output a message
  stop("At least three input arguments must be supplied. <input> <sample> [output]", call.=FALSE)

} else if (length(args)==1) {
  
  # Output file name not informed will be assigned as out.txt
  args[2] <- "output"
  args[3] <- "./"
}

inputfile <- args[1]
filename <- args[2]
outputdir <- args[3]



# -----------------------------------------------------------------------
# Program start

# Read input
infile <- read_csv(inputfile,
                  col_names = c("node1", "node2", "cdi"),
                  col_types = "ccd", skip = 1)

# Convert the list of gene-gene associations to igraph unweighted and undirected
ingraph <- graph_from_edgelist(as.matrix(infile[, 1:2]), directed = FALSE)

# Then, add weight to the network
edge.attributes(ingraph)$weight <- infile$cdi



# -----------------------------------------------------------------------
# Calculate all the centralities

# closeness centrality
cm_closeness <- as.data.frame(igraph::closeness(ingraph), make.names = TRUE)
cm_closeness$genes <- row.names(cm_closeness)
colnames(cm_closeness) <- c("cloneness", "genes")
cm_closeness_ord <- cm_closeness[order(cm_closeness$cloneness, decreasing = TRUE),]

# Degree centrality
cm_degree <- as.data.frame(igraph::degree(ingraph))
cm_degree$genes <- row.names(cm_degree)
colnames(cm_degree) <- c("degree", "genes")
cm_degree_ord <- cm_degree[order(cm_degree$degree, decreasing = TRUE),]

# Betweenness centrality
cm_betweenness <- as.data.frame(igraph::betweenness(ingraph))
cm_betweenness$genes <- row.names(cm_betweenness)
colnames(cm_betweenness) <- c("betweenness", "genes")
cm_betweenness_ord <- cm_betweenness[order(cm_betweenness$betweenness, decreasing = TRUE),]

# Pagerank centrality
cm_pagerank <- as.data.frame(igraph::page_rank(ingraph)$vector)
cm_pagerank$genes <- row.names(cm_pagerank)
colnames(cm_pagerank) <- c("pagerank", "genes")
cm_pagerank_ord <- cm_pagerank[order(cm_pagerank$pagerank, decreasing = TRUE),]

# Eigen centrality
cm_eigen <- as.data.frame(igraph::eigen_centrality(ingraph)$vector)
cm_eigen$genes <- row.names(cm_eigen)
colnames(cm_eigen) <- c("eigenvector", "genes")
cm_eigen_ord <- cm_eigen[order(cm_eigen$eigenvector, decreasing = TRUE),]



# -----------------------------------------------------------------------
# Bind all centralities already sorted

# Genes
centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank")
CM_genes <- cbind(cm_degree_ord$genes, cm_betweenness_ord$genes, cm_closeness_ord$genes, cm_eigen_ord$genes, cm_pagerank_ord$genes)
colnames(CM_genes) <- centralities
#head(CM_genes)

# CM values
centralities <- c("degree", "betweenness", "closeness", "eigen", "pagerank")
CM_values <- cbind(cm_degree_ord$degree, cm_betweenness_ord$betweenness, cm_closeness_ord$cloneness, cm_eigen_ord$eigenvector, cm_pagerank_ord$pagerank)
colnames(CM_values) <- centralities
#head(CM_values)



# -----------------------------------------------------------------------
# Export output

outputname1 <- paste0(outputdir, "/", filename, "_CM_genes.csv")
write.csv(CM_genes, outputname1, quote=FALSE, row.names=FALSE, col.names=TRUE)

outputname2 <- paste0(outputdir, "/", filename, "_CM_values.csv")
write.csv(CM_values, outputname2, quote=FALSE, row.names=FALSE, col.names=TRUE)

print("Files exported to ", outputdir)



# -----------------------------------------------------------------------
# TODO: PLOTS 

#1 Correlation of centralitites: CM X Degree
#2 Correlation of all centralities: correlation heatmap
#3 VennDiagram of all centralities
#4 Histogram of CM values
#5 Degree distribution scatter plot to show kind of powerlaw distribution


# Then, a second approach is to merge all datasets for a horizontal comparison
#1 Hub multiplicity: shows how much a gene is unique
#2 GO plots for time-series datasets




# save the environment
#save.image(file='result_CDI_centralities.RData')

# quit
quit(save="no")
