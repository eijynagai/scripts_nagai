#module load python/3.6
#source activate science

library(Seurat)
S3 <- readRDS(file = "S3.rds")

# Save as ps
postscript('plots.ps')

#This use the normalized values
#S3@assays$SCT@scale.data


# Getting the variable genes
var_fea <- VariableFeatures(S3)
length(var_fea)
summary(factor(S3@meta.data$orig.ident))

#This use the log scaled normalized value of the 5000 most variable genes
write.table(S3@assays$SCT@scale.data[var_fea,], file='Gene_NormCount_per_Cell_nonlog.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# filter out only the variable genes of 300 stage
x <- subset(S3, orig.ident == "300")
data <- x@assays$SCT@scale.data[var_fea,]
dim(data)
length(data)
head(data)
class(data)


# Normalize the data for each gene to mean 0 and standard deviation of 1
library(dplyr)
mydata <- apply(data, 1, scale)


# Calculating PCA based on normalized gene expression levels
mydata.pca <- prcomp(mydata, scale = TRUE)
plot(mydata.pca)
head(mydata.pca$x[,1:6])

# Calculating K-means clustering on PCA data






#parei aqui






# Calculating umap and plot
install.packages("devtools")
devtools::install_github("ropenscilabs/umapr")

library(umapr)
library(tidyverse)


embedding <- umap(as.data.frame(data.pca$x))


# plot the result
embedding %>%
  mutate( = iris$Species) %>%
  ggplot(aes(UMAP1, UMAP2, color = Species)) + geom_point()










##test this
# Get principal component vectors using prcomp instead of princomp
pc <- prcomp(data)
# First for principal components
comp <- data.frame(pc$x[,1:4])
# Plot
plot(comp, pch=16, col=rgb(0,0,0,0.5))
# end test this






data <- as.data.frame(data)
head(data)


## Ploting the PCA of 300 stage
#install.packages("ggfortify")
library(ggfortify)
pca.plot <- autoplot(data.pca, data = data)
pca.plot

#k-means on PCA data
for (i in c(10)){
    kmeans_x <- kmeans(x = data.frame(data.pca$x), centers = i, iter.max = 10, nstart = 1)
    print(barplot(table(kmeans_x$cluster)))

}



# Since PCA didn't work well I want to try UMAP

library(umap)

data.umap <- umap(data)

# visualize the data
data.umap

# visualize the content of knn
data.umap$knn


# Using the parameters
umap(
  data,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = NULL,
  learning_rate = 1,
  scale = FALSE,
  init = "spectral",
  init_sdev = NULL,
  spread = 1,
  min_dist = 0.01,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  nn_method = NULL,
  n_trees = 50,
  search_k = 2 * n_neighbors * n_trees,
  approx_pow = FALSE,
  y = NULL,
  target_n_neighbors = n_neighbors,
  target_metric = "euclidean",
  target_weight = 0.5,
  pca = 10,
  pca_center = TRUE,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  ret_model = FALSE,
  ret_nn = FALSE,
  ret_extra = c(),
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  tmpdir = tempdir(),
  verbose = getOption("verbose", TRUE)
)


# running and plotting umap
umap <- umap(data[!duplicated(data), -5])

df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 Species = data[!duplicated(data), 5])

ggplot(df, aes(x, y) +
  geom_point()















kmeans_10 <- kmeans(x = data, centers = 10, iter.max = 10, nstart = 1)$cluster
kmeans_10


# Plot T-sne of the kmeans

# intall packages
#install.packages("caret")
#install.packages("Rtsne")

## Rtsne function may take some minutes to complete...
set.seed(9)
tsne_model_1 = Rtsne(as.matrix(kmeans_10), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)






library(factoextra)
aggregate(data, by=list(cluster=kmeans_10), mean)
dd <- cbind(data, cluster = kmeans_10)
head(dd)

#### until here ok



S2$kmeans_15 <- kmeans(x = data, centers = 15)$cluster
S2$kmeans_20 <- kmeans(x = data, centers = 20)$cluster


## Scaling (modify)
scaling <- data.pca$sdev[1:2] * sqrt(nrow(pilots))

pc1 <- rowSums(t(t(sweep(pilots[,2:7], 2 ,colMeans(pilots[,2:7]))) * s.eigen$vectors[,1] * -1) / scaling[1])
pc2 <- rowSums(t(t(sweep(pilots[,2:7], 2, colMeans(pilots[,2:7]))) * s.eigen$vectors[,2]) / scaling[2])
Collect the PCs in a data.frame and plot using ggplot (loaded when ggfortify was loaded).

df <- data.frame(pc1, pc2, c(rep('Apprentice', 20), rep('Pilot', 20)))
colnames(df) <- c('PC1', 'PC2', 'Group')

ggplot(df, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point()
