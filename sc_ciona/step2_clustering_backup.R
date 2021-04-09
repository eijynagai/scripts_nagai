#!/usr/bin/env Rscript

library(Seurat)
library(argparse)
library(ggplot2)
library(dplyr)
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
S3 <- readRDS(paste0(dirname, '/S3.rds'))


## Clustering based on PCA distance
S3 <- FindNeighbors(S3, dims=1:PCx, verbose=FALSE, nn.eps = 0.5)
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


#------------ GO terms enriched per cluster
library(topGO)
topGOterms = function( fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Mouse",
                       ontology.use = "BP",
                       stats.use = "fisher",
                       algorithm.use = "weight01",
                       topnodes.print=20,
                       num.char=100){

  if (is.null(fg.genes) | is.null(bg.genes)){
    stop("Error : Both gene lists are empty")
  }

  require(topGO)
  if (organism == "Mouse"){
    mapping.use = "org.Mm.eg.db"
    library(org.Mm.eg.db)
  } else if (organism == "Human"){
    mapping.use = "org.Hs.eg.db"
    library(org.Hs.eg.db)
  } else if (organism == "Celegans"){
    mapping.use = "org.Ce.eg.db"
    library(org.Ce.eg.db)
  } else {
    stop("Error : Organisms other than mouse not supported currently")
  }

  n = length(bg.genes)
  geneList = integer(n)
  names(geneList) = bg.genes
  geneList[intersect(names(geneList), fg.genes)]=1
  print(paste0("Total ", length(geneList), " genes. ", sum(geneList), " genes in the foreground"))
  geneList = factor(geneList)

  if (ontology.use %in% c("BP", "CC", "MF")){
    print(paste0("Using Ontology : ", ontology.use))
  } else {
    stop("Error: Ontology not available. Should be one of BP, CC or MF")
  }
  # Make GO object
  GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = ontology.use,
                allGenes = geneList,
                annot = annFUN.org,
                mapping = mapping.use,
                ID = "SYMBOL",
                nodeSize = 10)
  print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
  res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)
  to.return = list()
  to.return$GOdata = GOdata
  to.return$res.table <- GenTable(GOdata, pval = res.result, topNodes = topnodes.print, numChar = num.char)
  return(to.return)
}


# delete previous outputfile
#Check its existence
if (file.exists(paste0(dirname, '/GO_cluster_ALL.txt')))
  #Delete file if it exists
  file.remove(paste0(dirname, '/GO_cluster_ALL.txt'))



# marker genes path included in cfg file
source(markers_path)

print('go ok and markers loaded')   ###test

#empty column
column_values <- list()

for (i in 1:length(clusters)) {
  S3.markers.cluster <- S3.markers %>% group_by(cluster) %>% filter(cluster == clusters[[i]])
  S3.markers.cluster.genes <- S3.markers.cluster$gene
  GOterms.S3.markers.cluster.genes = topGOterms(fg.genes = S3.markers.cluster.genes, bg.genes = rownames(S3), organism = 'Celegans')
  #print(GOterms.S3.markers.cluster.genes)


  print('go calculated ok') ###test

  # Plot Go terms p-values
  bp_plot <- GOterms.S3.markers.cluster.genes$res.table
  bp_plot$pval <- as.numeric(as.character(bp_plot$pval))

  #in case of NA values in pval
  bp_plot$pval[is.na(bp_plot$pval)] <- 1

  pp <- ggplot(bp_plot, aes(x = Term, y = -log10(as.numeric(pval)))) +
                geom_col() +
                ylab("Enrichment") +
                xlab("Biological process") +
                ggtitle("GO term enrichment") +
                scale_y_continuous(breaks = round(seq(0, max(-log10(as.numeric(bp_plot$pval))), by = 2), 1)) +
                scale_x_discrete(limits=bp_plot$Term) +
                theme_bw(base_size=24) +
                theme(
                    legend.position='none',
                    legend.background=element_rect(),
                    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
                    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
                    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
                    axis.title=element_text(size=24, face="bold"),
                    legend.key=element_blank(),     #removes the border
                    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
                    legend.text=element_text(size=18),  #Text size
                    title=element_text(size=18)) +
                guides(colour=guide_legend(override.aes=list(size=2.5))) +
                coord_flip()
  print(pp)  # plot GO terms




  # export individual cluster file
  capture.output(GOterms.S3.markers.cluster.genes, file=paste0(dirname, '/GO_cluster', clusters[[i]],'.txt'), append=FALSE)

  # making a consolidate result only with GO terms
  cat('\n\n## Cluster ', clusters[[i]], '\n', file=paste0(dirname, '/GO_cluster_ALL.txt'), append=TRUE)
  capture.output(GOterms.S3.markers.cluster.genes$res.table, file=paste0(dirname, '/GO_cluster_ALL.txt'), append = TRUE)

  if (findmarkers_test == 'roc'){
     # export POWER classification generated using ROC
     write.table(S3.markers.cluster, file=paste0(dirname, '/ROC_cluster',clusters[[i]],'.txt'), quote=F, col.names=T, sep='\t', row.names=F, append=F)
     # ploting
     x <- arrange(S3.markers.cluster, desc(power)) %>% top_n(9)
  } else {
     x <- S3.markers.cluster %>% top_n(9)
  }
  best.markers <- x$gene
  print(best.markers)

  print(VlnPlot(S3, features = best.markers, log=FALSE))
  #print(FeaturePlot(S3, features = best.markers, reduction = 'tsne'))
  print(FeaturePlot(S3,features=best.markers, reduction = 'umap'))
  #print(RidgePlot(S3, features = best.markers))

  print('go tables and marker plots ok') ###test

  #empty row
  row_values = list()
  for (j in 1:11){

    genes_list <- eval(parse(text=paste("markers", j, sep='')))
    match_genes <- S3.markers.cluster.genes[S3.markers.cluster.genes %in% genes_list]

    #markers name
    print(paste("markers", j, sep=''))

    #bind to the table
    row_values[j] <- length(unique(match_genes))
    # bind the genes instead of number
    #row_values[j] <- unique(match_genes)

    # if there is any match, plot
    if (length(unique(match_genes))) {
      VlnPlot(S3, features = match_genes, pt.size = 0.2,  combine = TRUE)

      # log10
      #VlnPlot(S3, features = match_genes, slot = "counts", log = TRUE,  pt.size = 0.1,  combine = TRUE)

      # umap position plot
      FeaturePlot(S3, features = match_genes, reduction = "umap", sort.cell = TRUE, pt.size = 0.2, label.size = 2, combine = TRUE)
    }
  }

  column_values[[i]] <- row_values

} # from the main

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







print('Script cluster.R finished.')

dev.off()
q()
