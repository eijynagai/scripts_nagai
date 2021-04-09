#!/usr/bin/env Rscript

library(Seurat)
library(argparse)
library(ggplot2)
library(dplyr)
library(patchwork)
################################################
#running in parallel using FUTURE package
library(future)
#availableCores()

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
postscript(paste0(dirname, '/plots2.ps'), horizontal=T)
#pdf(paste0(dirname, '/plots2.pdf'),
#        paper = 'A4r',
#        height = 0.01,
#        width = 0.05,
#        useDingbats=FALSE)
par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)


# input S3 data
S2 <- readRDS(paste0(dirname, '/S2.rds'))


## Clustering based on PCA distance
S3 <- FindNeighbors(S2, dims=1:PCx, verbose=FALSE, nn.eps = 0.5)
S3 <- FindClusters(S3, resolution = RES, verbose=FALSE, n.start = 10)

# [QC] Expression levels in each cluster
VlnPlot(S3,features="nCount_RNA")
VlnPlot(S3,features="nFeature_RNA")

# [QC] Percentage of largest gene in each cluster
apply(
  S3@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x)
) -> S3$Percent.Largest.Gene
VlnPlot(S3,features="Percent.Largest.Gene")


# PLOT TSNE and UMAP together
p1 <- DimPlot(S3, reduction='umap', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='UMAP')
p2 <- DimPlot(S3, reduction='tsne', label=TRUE, pt.size=0.1, label.size=6)+ggtitle(label='TSNE')
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
(p1 + p2) & NoLegend()


# Find markers for every cluster compared to all remaining cells.
# logFC >0; highly expressed at the 1st group
# p_val_adj; bonferroni correction
# pct1.; percentage of cells where the gene is detected in the 1st group
# pct2.; percentage of cells where the gene is detercted in the 2nd group
S3.markers <- FindAllMarkers(S3, only.pos=TRUE
                             ,min.pct=0
                             ,test.use=findmarkers_test
                             ,logfc.threshold=logFC
                             ,return.thres=0.01
                             ,min.cells.feature=qtd_min_cells
                             ,min.cells.group=qtd_min_cells)

sprintf('test using %s finished', findmarkers_test)

if (findmarkers_test == 'wilcox') {
    idx <-S3.markers[, 'p_val_adj'] < qvThres
    S3.markers <- S3.markers[idx,]
    #S3.markers
}


# Make the list of existent clusters
clusters <- unique(S3.markers$cluster)
clusters


### Number of genes in each cluster
# tabulate cells by cluster ID
cat('## Number of cells per cluster\n', file=paste0(dirname, '/clusters.txt'), append = FALSE)
write.table(table(Idents(S3)), file=paste0(dirname, '/clusters.txt'), quote=F, col.names=T, sep='\t', row.names=T, append=T)

cat('\n\n## Distribution of cells in each cluster\n', file=paste0(dirname, '/clusters.txt'), append=TRUE)
write.table(prop.table(table(Idents(S3))), file=paste0(dirname, '/clusters.txt'), quote=F, col.names=T, sep='\t', row.names=T, append=T)


# marker genes path included in cfg file
source(markers_path)

print('go ok and markers loaded')   ###test

#empty column
column_values <- list()

for (i in 1:length(clusters)) {
  S3.markers.cluster <- S3.markers %>% group_by(cluster) %>% filter(cluster == clusters[[i]])
  S3.markers.cluster.genes <- S3.markers.cluster$gene

  #empty row
  row_values = list()
  for (j in 1:11){

    genes_list <- eval(parse(text=paste("markers", j, sep='')))
    match_genes <- S3.markers.cluster.genes[S3.markers.cluster.genes %in% genes_list]

    #bind to the table
    row_values[j] <- length(unique(match_genes))

    # if there is any match, plot
    if (length(unique(match_genes))) {
      # umap position plot
      p3 <- FeaturePlot(S3, features = match_genes, reduction = "umap", order = TRUE, pt.size = 0.1, label.size = 2, combine = FALSE)
      p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
      wrap_plots(p3)
    }
  }

  column_values[[i]] <- row_values

}

# gather all results and make the dataframe
df <- do.call(rbind, column_values)
colnames(df) <- c(paste0("marker", 1:11))
row.names(df) <- c(clusters)

print(df)



# -------- Printing the table to the plot

library(grid)
library(gridExtra)


grid.ftable <- function(d, padding = unit(4, "mm"), ...) {

  nc <- ncol(d)
  nr <- nrow(d)

  ## character table with added row and column names
  extended_matrix <- cbind(c("", rownames(d)),
                           rbind(colnames(d),
                                 as.matrix(d)))

  ## string width and height
  w <- apply(extended_matrix, 2, strwidth, "inch")
  h <- apply(extended_matrix, 2, strheight, "inch")

  widths <- apply(w, 2, max)
  heights <- apply(h, 1, max)

  padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)

  x <- cumsum(widths + padding) - 0.5 * padding
  y <- cumsum(heights + padding) - padding

  rg <- rectGrob(x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 width = unit(widths + padding, "in"),
                 height = unit(heights + padding, "in"))

  tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 just = "center")

  g <- gTree(children = gList(rg, tg), ...,
             x = x, y = y, widths = widths, heights = heights)

  grid.draw(g)
  invisible(g)
}


grid.newpage()
print(grid.ftable(head(df, length(clusters)), gp = gpar(fill = rep(c("grey90", "grey95"), each = (length(clusters)+1)))))

#saving the R object
save(S3.markers, file=paste0(dirname, '/S3.markers.RData'))
saveRDS(S3, file = paste0(dirname, '/S3.rds'))

print('Script step2_clustering.R finished.')

dev.off()
q()
