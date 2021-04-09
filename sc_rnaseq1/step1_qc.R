#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(patchwork)
library(ggplot2)
library(cowplot)
library(argparse)

################################################
#running in parallel using FUTURE package
library(future)
availableCores()

# this method used the available cores specified in qsub
plan(multicore)
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


# Parameters load file
source(paste0(dirname,'/cfg_params.txt'))

# Output plot file
postscript(paste0(dirname, '/plots.ps'), horizontal=T)
#pdf(paste0(dirname, '/plots.pdf'),
#        paper = 'A4r',
#        height = 0.01,
#        width = 0.05,
#        useDingbats=FALSE)
par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)

# Pre-process Seurat Object
#
# STEP1: Read in 10X Cell Ranger matrix and create Seurat object
counts_matrix_filename = cr_dir
CR <- Read10X(data.dir = counts_matrix_filename)

# remove the '.1.1' ... suffixes
colnames(CR) <- sub('.[12].[123]', '', colnames(CR))

# Check for duplicated barcodes
summary(duplicated(colnames(CR)))

# Metrics
dim(CR)
object.size(CR) # size in bytes
object.size(as.matrix(CR)) # size in bytes


## Seurat Object creation
S1 <- CreateSeuratObject(counts = CR[,which(!duplicated(colnames(CR)))]
                            , project = proj_name
                            , assay='RNA'
                            , min.cells = qtd_min_cells
                            , min.features = qtd_min_features
                            , names.delim = '-'
                            , names.field = 2
                            , meta.data = NULL
)

# mitochondrial gene list
mt_gene_list = c('nduo-6', 'atp-6', 'nduo-2', 'ctb-1', 'ctc-3', 'nduo-4', 'ctc-1', 'ctc-2','nduo-3', 'nduo-5')

#check if mitochondrial genes are in the matrix
print('All mitochondrial genes in the data matrix?')
all(mt_gene_list %in% rownames(S1))


S1 <- PercentageFeatureSet(object=S1 , features=mt_gene_list, col.name='Percent_MT')

######## Collect metrics?
metrics <- FALSE

if (metrics) {
    # Data metrics collection
    counts_per_cell <- Matrix::colSums(CR)
    counts_per_gene <- Matrix::rowSums(CR)
    genes_per_cell <- Matrix::colSums(CR>0)
    cells_per_gene <- Matrix::rowSums(CR>0) # only count cells where the gene is expressed

    # Data metrics plots
    p1 <- hist(log10(counts_per_cell+1),main='counts per cell',col=colorplots)
    p2 <- hist(log10(genes_per_cell+1), main='genes per cell', col=colorplots)
    p3 <- plot(counts_per_cell, genes_per_cell, log='xy', col=colorplots)+title('counts vs genes per cell')
    p4 <- hist(log10(counts_per_gene+1), main='counts per gene', col=colorplots)
    p5 <- plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
    p6 <- VlnPlot(S1, features = c('nFeature_RNA', 'nCount_RNA', 'Percent_MT'), ncol = 3, pt.size=0.1, log=TRUE)

    # Apoptotic cells before filtering
    APOPT_S1 <- FeatureScatter(S1, feature1 = 'nCount_RNA', feature2 = 'Percent_MT', pt.size=0.4, cols=colorplots)+NoLegend()
    # Droplets before filtering
    DLET_S1 <- FeatureScatter(S1, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', pt.size=0.4, cols=colorplots)+NoLegend()

    ########## DATA metrics PLOTS
    p1 # Histogram of counts per cell
    p2 # Histogram pf genes per cell
    p3 # Scatter counts per gene per cell
    p4 # Histogram Log10 counts per gene
    p5 # Number of genes per cell ordered
    print(p6) # Violin plot
    print(APOPT_S1)
    print(DLET_S1)

    #FREE memory
    save(CR, file=paste0(dirname, '/CR.RData'))
    save(S1, file=paste0(dirname, '/S1.RData'))
    remove(p1,p2,p3,p4,p5,p6, APOPT_S1, DLET_S1)
    remove(CR)
print('Metrics collected succesfully.')
}
# Summary (attributes(S1))
attr(S1, 'assays')$RNA


### STEP2: Filtering unwanted cells
S2 <- subset(S1, subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & nCount_RNA < 50000 & Percent_MT < max_perc_mt)

# Apoptosis after filtering
p9 <- FeatureScatter(S2, feature1 = 'nCount_RNA', feature2 = 'Percent_MT', pt.size=0.5, cols=colorplots)+NoLegend()

# Droplets after filtering
p10 <- FeatureScatter(S2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', pt.size=0.5, cols=colorplots)+NoLegend()

# Summary (attributes(S2))
attr(S2, 'assays')$RNA

# Free memory
remove(S1)

print('Subset done')
# Run SCTransform? Else run standard pipeline
sctransf <- TRUE

if (sctransf){
    S2 <- SCTransform(S2, vars.to.regress = "Percent_MT"
                        , variable.features.n = nFeatures_var
                        , verbose = TRUE)
    print('SCTraforme done')
} else {

    ### STEP3: Pre-process Seurat Object
    #3.1 Normalization
    S2 <- NormalizeData(S2, normalization.method='LogNormalize', scale.factor = norm_scale, verbose=FALSE)

    #3.2 Selecting variable genes
    S2 <- FindVariableFeatures(S2, selection.method = 'vst', nfeatures = nFeatures_var, verbose=FALSE)

    # These lines will convert NAs to 0
    ##################
    sum(S2@assays$RNA@meta.features$vst.variance == 0)
    table(is.na(S2@assays$RNA@data@x))
    S2@assays$RNA@data@x[is.na(S2@assays$RNA@data@x)] <- 0
    ###################

    #3.3 Scaling
    all.genes <- rownames(S2)
    S2 <- ScaleData(S2, vars.to.regress='Percent_MT', features = all.genes, verbose=FALSE)

}


VlnPlot(S2, features = c('nFeature_RNA', 'nCount_RNA', 'Percent_MT'),ncol = 3, pt.size=0, log=TRUE)


# Getting the variable genes
var_fea <- VariableFeatures(S2)
length(var_fea)
summary(factor(S2@meta.data$orig.ident))

# Identify the N most highly variable genes
topN <- head(VariableFeatures(S2), topx)
topN

# plot variable features with and without labels
p11 <- VariableFeaturePlot(S2)
p12 <- LabelPoints(plot = p11, points = topN, repel = TRUE, xnudge=0, ynudge=0)

# Plot percentage of mitochondrial genes by clusters
VlnPlot(S2, features = c('nFeature_RNA', 'nCount_RNA', 'Percent_MT'),ncol = 3, pt.size=0, log=TRUE)


# PCA with nVF
S2 <- RunPCA(S2, features = VariableFeatures(object = S2)
            , npcs=nPCS
            , verbose=TRUE)
print('PCA done')

names(S2@graphs)

kval <- c(10, 15, 20)

for (k in kval){

    kv <- paste0('kmeans_', k)
    print(kv)
    #k-means on PCA data
    S2$kv <- kmeans(x = S2@reductions$pca@cell.embeddings, centers = k)$cluster
    mytable <- table(S2$kv)
    myframe <- as.data.frame(mytable)
    myframe
    p <- ggplot(myframe, aes(x=Var1, y=Freq), ) +
        geom_point(size=5) +
        ggtitle(kv) +
        xlab('Cluster') +
        ylab('Number of cells')

    p + theme(
        text = element_text(size=18),
        plot.title = element_text(size=20, face='bold'),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)
    ) + ggtitle(kv) +
    xlab('Cluster') +
    ylab('Number of cells')
}



# UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
S2 <- RunUMAP(S2, dims=1:10, verbose=TRUE)
S2 <- RunTSNE(S2, dims=1:10, verbose=TRUE)

for (j in kval){

    jv <- paste0('kmeans_', k)
    DimPlot(S2, reduction = "umap", group.by = jv)+ggtitle(jv)
    DimPlot(S2, reduction = "tsne", group.by = jv)+ggtitle(jv)
}

print('umap and tsne done')



## GO terms of each cluster

















#----------------------- Plots
#print(p9)
#print(p10)
#print(p12) # variablefeatureplots

#p15 <- FeatureScatter(S2, feature1 = 'nCount_RNA', feature2 = 'Percent_MT', pt.size=0.2, cols=colorplots)
#p16 <- FeatureScatter(S2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', pt.size=0.2, cols=colorplots)
#p17 <- DimPlot(S2, reduction='pca', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='PCA')
#p18 <- DimPlot(S2, reduction='umap', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='UMAP')
#p19 <- DimPlot(S2, reduction='tsne', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='TSNE')

#p15 # featurescatter('ncount_rna, percent_mt')
#p16 # featurescatter('ncount_rna, nfeature_rna')
#p17 # dimplot('PCA')
#p18 # dimplot('umap')
#p19 # dimplot('tsne')



#remove(p9)
#remove(p10)
#remove(p12)
#remove(p15)
#remove(p16)
#remove(p17)
#remove(p18)
#remove(p19)
dev.off()

# Store current Seurat object
#save(S2, file=paste0(dirname, '/S2.RData'))
saveRDS(S2, file = paste0(dirname, '/S2.rds'))
#save.image(file = paste0(dirname, '/workspace.RData'))

remove(S2)

print('Script step1_qc.R finished.')

q()
