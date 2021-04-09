#!/usr/bin/env Rscript

suppressMessages(require(Seurat))
suppressMessages(require(dplyr))
suppressMessages(require(cowplot))
suppressMessages(require(ggplot2))
suppressMessages(require(dplyr))
suppressMessages(require(patchwork))
suppressMessages(require(argparse))

################################################
#running in parallel using FUTURE package
library(future)
availableCores()

# this method used the available cores specified in qsub
#plan(multicore)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 8 * 1024^3)
################################################

# Create parse parameters
parser <- ArgumentParser()
parser$add_argument('--outdir', type='character', default='./',
    help='Directory to be used for output files')
parser$add_argument('--resolution', type='double', default=4.0,
    help='Resolution for granularity')
args <- parser$parse_args()
dirname <- args$outdir
RES <- args$resolution


# In case of running R directly
#dirname <- '/home/nagailae/projects/sc_ciona/seurat3_v1/result/ciona_res0.09'
#RES <- 0.09

# Parameters load file
source(paste0(dirname,'/cfg_params.txt'))

# Output plot file
postscript(paste0(dirname, '/plots.ps'), horizontal=T)
par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)


# Pre-process Seurat Object
#
sample.info <- read.delim("script/aggregate_libraries.csv", header = T, sep=',', stringsAsFactors = F)

sample.list <- list()

for (r in 1:nrow(sample.info)) {
    dat.10x <- Read10X_h5(paste0('/home/nagailae/projects/sc_ciona/seurat3_v1/data/',sample.info$file_name[r]))
    dat.10x <- CreateSeuratObject(counts = dat.10x, project = sample.info$library_id[r],
        min.cells = qtd_min_cells, min.features = qtd_min_features)
    sample.list[[sample.info$library_id[r]]] <- dat.10x
}


# SCTrasnform normalization for each stage
for (i in 1:length(sample.list)) {
    sample.list[[i]] <- SCTransform(sample.list[[i]], verbose = FALSE)
}


# Selecting features for downstream integration including Pearson resuduals
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)

sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = sample.features,
    verbose = FALSE)


# Identify anchors and integrate the datasets
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT",
    anchor.features = sample.features, verbose = FALSE)

sample.integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = "SCT",
    verbose = FALSE)


# now proceed with downstram analysis (visualization and clustering)
sample.integrated <- RunPCA(sample.integrated, verbose = FALSE)
sample.integrated <- RunUMAP(sample.integrated, dims = 1:10)
sample.integrated <- RunTSNE(sample.integrated, dims = 1:10)

plots <- DimPlot(sample.integrated, combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
    byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)


# cluster
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:10)
sample.integrated <- FindClusters(sample.integrated, resolution = RES)


# plots
VlnPlot(sample.integrated, features="nCount_RNA")
VlnPlot(sample.integrated, features="nFeature_RNA")

p1 <- DimPlot(sample.integrated, reduction = "umap", label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='UMAP')
p2 <- DimPlot(sample.integrated, reduction='tsne', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='TSNE')
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
(p1 + p2) & NoLegend()

S3.markers <- FindAllMarkers(sample.integrated, only.pos=TRUE
                            ,min.pct=0
                            ,test.use=findmarkers_test
                            ,logfc.threshold=logFC
                            ,return.thres=0.01
                            ,min.cells.feature=qtd_min_cells
                            ,min.cells.group=qtd_min_cells)

sprintf('test using %s finished', findmarkers_test)

# Make the list of existent clusters
clusters <- unique(S3.markers$cluster)
clusters



#saving the R object
saveRDS(sample.list, file = paste0(dirname, '/sample.list.rds'))

dev.off()
q()
