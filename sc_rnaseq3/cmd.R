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
options(future.globals.maxSize = 12 * 1024^3)
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

# Metrics
dim(CR)
object.size(CR) # size in bytes
object.size(as.matrix(CR)) # size in bytes


## Seurat Object creation
S1 <- CreateSeuratObject(counts = CR, project = proj_name
                            , assay='RNA'
                            , min.cells = qtd_min_cells
                            , min.features = qtd_min_features)


##### testing
orig.levels <- levels(S1)
Idents(S1) <- gsub(pattern = " ", replacement = "-", x = Idents(S1))
orig.levels <- gsub(pattern = " ", replacement = "-", x = orig.levels)
levels(S1) <- orig.levels

# mitochondrial gene list
mt_gene_list = c('nduo-6', 'atp-6', 'nduo-2', 'ctb-1', 'ctc-3', 'nduo-4', 'ctc-1', 'ctc-2','nduo-3', 'nduo-5')

#check if mitochondrial genes are in the matrix
print('All mitochondrial genes in the data matrix?')
all(mt_gene_list %in% rownames(S1))

S1 <- PercentageFeatureSet(object=S1 , features=mt_gene_list, col.name='Percent_MT')

######## Collect metrics?
metrics <- TRUE

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

    # Capture abnormal cells
    APOPT_S1 <- FeatureScatter(S1, feature1 = 'nCount_RNA', feature2 = 'Percent_MT', pt.size=0.4, cols=colorplots)+NoLegend()

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
}


# Summary (attributes(S1))
attr(S1, 'assays')$RNA


### STEP2: Filtering unwanted cells
S2 <- subset(S1, subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & nCount_RNA < 50000 & Percent_MT < max_perc_mt)


# Apoptosis
p9 <- FeatureScatter(S2, feature1 = 'nCount_RNA', feature2 = 'Percent_MT', pt.size=0.5, cols=colorplots)+NoLegend()

# Droplets
p10 <- FeatureScatter(S2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', pt.size=0.5, cols=colorplots)+NoLegend()

# Summary (attributes(S2))
attr(S2, 'assays')$RNA

# Free memory
remove(S1)


# Run SCTransform? Else run standard pipeline
sctransf <- TRUE

if (sctransf){

    S3 <- SCTransform(S2, vars.to.regress = "Percent_MT", verbose = FALSE)

} else {

    ### STEP3: Pre-process Seurat Object
    #3.1 Normalization
    S3 <- NormalizeData(S2, normalization.method='LogNormalize', scale.factor = norm_scale, verbose=FALSE)

    #3.2 Selecting variable genes
    S3 <- FindVariableFeatures(S3, selection.method = 'vst', nfeatures = nFeatures_var, verbose=FALSE)

    # These lines will convert NAs to 0
    ##################
    sum(S3@assays$RNA@meta.features$vst.variance == 0)
    table(is.na(S3@assays$RNA@data@x))
    S3@assays$RNA@data@x[is.na(S3@assays$RNA@data@x)] <- 0
    ###################

    #3.3 Scaling
    all.genes <- rownames(S3)
    S3 <- ScaleData(S3, vars.to.regress='Percent_MT', features = all.genes, verbose=FALSE)

}

# Identify the N most highly variable genes
topN <- head(VariableFeatures(S3), topx)
topN

# plot variable features with and without labels
p11 <- VariableFeaturePlot(S3)
p12 <- LabelPoints(plot = p11, points = topN, repel = TRUE, xnudge=0, ynudge=0)

# Plot percentage of mitochondrial genes by clusters
VlnPlot(S3, features = c('nFeature_RNA', 'nCount_RNA', 'Percent_MT'),ncol = 3, pt.size=0, log=TRUE)


# PCA with nVF
S3 <- RunPCA(S3, features = VariableFeatures(object = S3)
            , ndims.print=1:3
            , nfeatures.print=5
            , npcs=nPCS
            , verbose=FALSE)

# Define automatically the number of PCs to consider or choose in cfg_file
#PCx <- length(attributes(S3[['pca']])$stdev) #automatically
#print(paste0('Number of PCs used:',PCx))


# Determine the dimentionality of the dataset (cannot do on SCTransform)
if (!(sctransf)){
    summary(attributes(S3[['pca']]))
    S3 <- JackStraw(S3, dims=PCx, num.replicate=100, verbose=FALSE)
    S3 <- ScoreJackStraw(S3, dims=1:PCx)
    print(JackStrawPlot(S3, dims=1:PCx, reduction='pca'))
    print(ElbowPlot(S3, ndims=nPCS, reduction='pca') + geom_vline(xintercept=PCx, colour='red'))
}

sprintf('Number of total PCs is %s, and you chose dimension of %s PCs. The second must smaller than the first.', nPCS, PCx)

# UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
S3 <- RunUMAP(S3, dims=1:PCx, verbose=FALSE)
S3 <- RunTSNE(S3, dims=1:PCx)


#----------------------- Plots
#p7
#p8
print(p9)
print(p10)
#p11
print(p12) # variablefeatureplots
#p13 # heatmap top10
#p14 # heatmap top20

#p15 <- FeatureScatter(S3, feature1 = 'nCount_RNA', feature2 = 'Percent_MT', pt.size=0.2, cols=colorplots)
#p16 <- FeatureScatter(S3, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', pt.size=0.2, cols=colorplots)
#p17 <- DimPlot(S3, reduction='pca', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='PCA')
#p18 <- DimPlot(S3, reduction='umap', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='UMAP')
#p19 <- DimPlot(S3, reduction='tsne', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='TSNE')

#p15 # featurescatter('ncount_rna, percent_mt')
#p16 # featurescatter('ncount_rna, nfeature_rna')
#p17 # dimplot('PCA')
#p18 # dimplot('umap')
#p19 # dimplot('tsne')



######## test this
#library(clustree)
#clustree(S3, prefix = 'k', node_colour = 'sc3_stability')



#remove(p7)
#remove(p8)
remove(p9)
remove(p10)
#remove(p11)
remove(p12)
#remove(p13)
#remove(p14)
#remove(p15)
#remove(p16)
#remove(p17)
#remove(p18)
#remove(p19)






#not using it anymore
if(FALSE){
    # for checking marker genes from a list of genes
    source(paste0(dirname,'marker_genes.txt'))

    for (i in 1:11){
        genes_list <- assign(paste0("markers", i))

        # standard
        VlnPlot(S3, features = genes_list, pt.size = 0.2,  combine = FALSE) 

        # log10
        #VlnPlot(S3, features = markers[i] , slot = "counts", log = TRUE,  pt.size = 0.1,  combine = FALSE)

        # umap position plot
        FeaturePlot(S3, features = genes_list, reduction = "umap", sort.cell = TRUE, pt.size = 0.2, label.size = 2, combine = FALSE)
    }
}



#Renaming the clusters
#new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(S3)
#S3 <- RenameIdents(S3, new.cluster.ids)
#DimPlot(S3, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()





dev.off()

# Store current Seurat3 object
save(S3, file=paste0(dirname, '/S3.RData'))
#save(S3.markers, file=paste0(dirname, '/S3.markers.RData'))
saveRDS(S3, file = paste0(dirname, '/S3.rds'))
#save.image(file = paste0(dirname, '/workspace.RData'))


remove(S3)
#remove(S3.markers)

print('Script cmd.R finished.')

q()
