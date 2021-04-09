#!/usr/bin/env perl

require("/home/park/PERL/park_basic.pl");

### narrow downling 
if(@ARGV != 3){
    print "%>perl " . $0 . " [Seurat3 DIR] [OUT_DIR] [Cluster IDs: 0,1,2,...]\n";
    exit;
}

$SR_INPUT_DIR = shift(@ARGV);      #Input Seurat3
$SR_OUTPUT_DIR = shift(@ARGV);  #Output Seurat3

@CLSs = split(/\,/, shift(@ARGV) );   # Given Cluster IDs of Input S3 object
if( scalar(@CLSs) < 1){
    print "Given Cluster iDs were zero \n";
    exit;
}
$given_clusters =  join(", ", @CLSs);


if( !defined($SR_INPUT_DIR = &Chk_Dir($SR_INPUT_DIR)) ){ exit; } 

$this_output_dir = $SR_OUTPUT_DIR; ## . $this;
$this_output_dir = &myMKDIR($this_output_dir);

$param_file = $SR_INPUT_DIR  . "cfg_params.txt";
$Rdata_S3   = $SR_INPUT_DIR . "Seurat3.RData";

if( &Chk_File($param_file) == 0 ){ exit; }
if( &Chk_File($Rdata_S3) == 0 ){ exit; }

#marker genes from Seurat3
$this_cls_dir = &myMKDIR($this_output_dir . "ClsMarkers");
#marker genes from Seurat3 and Monocle3 modules (co-regulated genes) in UMAP clusters
$this_Umap_module_dir = &myMKDIR($this_output_dir . "Modules_UMAP");
#marker genes for Monocle3 trajectory cells
$this_Traj_module_dir = &myMKDIR($this_output_dir . "Modules_Traj");

#run once again, these must be clean up
if(-e $this_cls_dir){
    system("rm -fr $this_cls_dir");
}
if(-e $this_Umap_module_dir){
    system("rm -fr $this_Umap_module_dir");
}
if(-e $this_Traj_module_dir){
    system("rm -fr $this_Traj_module_dir");
}

$this_cls_dir = &myMKDIR($this_output_dir . "ClsMarkers");
$this_Umap_module_dir = &myMKDIR($this_output_dir . "Modules_UMAP");
$this_Traj_module_dir = &myMKDIR($this_output_dir . "Modules_Traj");



$stat_file = $this_output_dir . "numbers.txt";
$ps_file   = $this_output_dir . "plots.ps";
$Rdata_file       = $this_output_dir . "Seurat3_subCLS.RData";
$S3marker_file = $this_output_dir . "Seurat3_subCLS.markers.RData";
$CDS_file = $this_output_dir . "CDS_subCLS.RData";
$M3marker_file = $this_output_dir . "Monocle3_subCLS.markers.RData";

system("cp $param_file $this_output_dir");

unlink($ps_file);
unlink($Rdata_file);
unlink($S3marker_file);
unlink($CDS_file);
unlink($M3marker_file);

#----------------------------------
open(OUT, ">$stat_file");
print OUT "#Selected Clusters\n";
for($i=0; $i <scalar(@CLSs); $i++){
    print OUT "# Cluster" . $CLSs[$i] . "\n";
}
close(OUT);
#----------------------------------


&Seurat3_SubCLS_R();

exit;

sub Seurat3_SubCLS_R{
    my $r_cmd = $this_output_dir . "R.cmd";
    my $r_out   = $this_output_dir . "R.out";

## Prepare R command script 
open(OUT, ">$r_cmd");
print OUT <<EOF;

library(Seurat)
library(monocle3)
library(dplyr)
library(cowplot)
library(ggplot2)

# subroutine for handling Monocle3 modules
myClsModule <- function(title, output_dir, agg_mat, gene_module, seurat.markers ){

    #Plot heatmap
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", 
		       main=title, fontsize=12)
    outfile <- paste0(output_dir, "Matrix_module_cluster.txt")
    write.table( round(agg_mat,6) , outfile, quote=F, col.names=T, append=F, sep="\\t", row.names=T)
     
    # NOTE: one module includes one cluster
    ### Write a link of Cluster to Raw Module
    for( Cls in colnames(agg_mat) ){
	cID <- as.numeric( unlist( strsplit(Cls, " ") )[2] )
    
            # get Max module score
	    Mod <- names( sort( agg_mat[, Cls], decreasing=T )[1] )
	    mID <- as.numeric( unlist( strsplit(Mod, " ") )[2] )

	    idx <- gene_module\$module == mID
	    if( sum(idx) < 1){ next }
	    this_genes <- sort( names( gene_module\$module[ idx ] ) )

            # ClusterX_markerModuleY.txt; one cluster and one module genes (=marker genes)
	    tbl <- seurat.markers[this_genes, c("cluster", "p_val", "p_val_adj", "avg_logFC", "gene")]
	    notNAidx <- !is.na( tbl[,1] )
	    if(sum(notNAidx) < 1){ next; }
	    tbl <- tbl[notNAidx,]
	    idx <- tbl[, "cluster"] == cID
	    tbl <- tbl[idx,]

	    # Cluster0_marker <- Cluster--Module link based on maximal scored moduel
	    outfile <- gsub( " ", "", paste0(output_dir, Cls, "_marker", Mod, ".txt") )
	    idx <- sort( tbl[, "avg_logFC"], decreasing=TRUE, index.return=TRUE)\$ix 
	    write.table( unique(tbl[idx,]), outfile, quote=F, col.names=T, append=F, sep="\\t", row.names=F)

	    remove(notNAidx)
	    remove(idx)
	    remove(this_genes)
    }

    ### Write modules on Trajectory
    for(mID in levels(gene_module\$module) ){
	idx <- gene_module\$module == mID
	if( sum(idx) < 1){ next }
	this_genes <- sort( names( gene_module\$module[ idx ] ) )

	    # raw module -> module ovl with Cluster ->

	    outfile <- paste0(output_dir, "raw_Module_", mID, ".txt")
	    tbl <- cds_pr_test[ this_genes, ]
	    notNAidx <- !is.na( tbl[,1] )
	    if(sum(notNAidx) < 1){ next; }
	    tbl <- tbl[notNAidx,]
	    tbl <- tbl[ , c("p_value","morans_I","q_value","gene_short_name") ]
	    colnames(tbl) <- c("p_value","morans_I","q_value","gene")
	    write.table( unique(tbl), outfile, quote=F, col.names=T, append=F, sep="\\t", row.names=F)
	    remove(tbl)
	    remove(notNAidx)

	    # ClusterX_markerModuleY.txt; one cluster and one module genes (=marker genes)
	    # ovlClsMarker_ModuleY.txt: All module Y genes also present at Cluster marker (cluster IDs are no matter)
	    # So question; why a module is DEGs in several clusters?
	    outfile <- paste0(output_dir, "ovlClsMarker_Module_", mID, ".txt")
	    tbl <- seurat.markers[this_genes, c("cluster", "p_val", "p_val_adj", "avg_logFC", "gene")]
	    notNAidx <- !is.na( tbl[,1] )
	    if(sum(notNAidx) < 1){ next; }
	tbl <- tbl[notNAidx,]
	    remove(notNAidx)

	    # Cluster and genes in this module
	    # header	
	    cnt <-  tbl %>% group_by(cluster) %>% summarize(count=n())
	    cnt <- rbind(as.numeric( t(cnt[,"cluster"])), as.numeric( t(cnt[,"count"]) ))
	    colnames(cnt) <- c( paste0("#Cluster", cnt[1,]) )

    idx <- sort( cnt[2,], decreasing=T, index.return=TRUE)\$ix
    cnt <- cnt[ ,idx]
    cnt <- rbind( colnames(cnt), cnt)
    #rownames(cnt)<- c("#ID", "Cluster", "Count")
    cat("#Cluster_ID\\tCluster\\tCount\\n", file=outfile, append=FALSE)
    write.table( t(cnt), outfile, quote=F, col.names=F, append=T, sep="\\t", row.names=F)

    idx <- sort( tbl[, "avg_logFC"], decreasing=TRUE, index.return=TRUE)\$ix 
    write.table( unique(tbl[idx,]), outfile, quote=F, col.names=T, append=T, sep="\\t", row.names=F)

    remove(tbl)
    remove(idx)
    remove(this_genes)
    remove(cnt)
    }

#----------------------------
}



source("$param_file")
load("$Rdata_S3")

S3_subCLS <- subset(S3, idents=c($given_clusters) )
## Pre-process Seurat object
if( SCTransf ){
    S3_subCLS <- SCTransform(S3_subCLS 
		,variable.features.rv.th = 1.3  #rv: residual variance cutoff
		,variable.features.n = NULL
		,vars.to.regress="Percent_MT", 
		return.only.var.genes = TRUE,
		verbose=FALSE)
#warnings()
}else{

    # Step 1: Normalizing
    S3_subCLS <- NormalizeData(S3_subCLS, normalization.method = "LogNormalize", scale.factor = Norm_Scale, verbose=FALSE)
    # Step 2: Selecting variable genes
    S3_subCLS <- FindVariableFeatures(S3_subCLS, selection.method = "vst", nfeatures = Variable_Feature, verbose=FALSE)

    # Step 3: Scaling
    ## ------- scaling with all genes
    if(scaleAll){
	all.genes <- rownames(S3_subCLS)
        S3_subCLS <- ScaleData( S3_subCLS, vars.to.regress = "Percent_MT", features = all.genes, verbose=FALSE)
    }else{
        ##------- scaling with selected nVF only
	S3_subCLS <- ScaleData( S3_subCLS, vars.to.regress = "Percent_MT", verbose=FALSE )
    }

}


postscript("$ps_file", horizontal=T)
op <- par(cex=1.2)

TOP_Gene <- head(VariableFeatures(S3_subCLS), TOPx)
TOP_Gene <- c(TOP_Gene, "Foxc1", "Pitx2", "Dkk2", "Tfap2b", "Lmx1b", "Emx2")

p1 <- VariableFeaturePlot(S3_subCLS)
p2 <- LabelPoints(plot=p1, points=TOP_Gene, repel=TRUE, xnudge=0,  ynudge=0)

S3_subCLS <- RunPCA( S3_subCLS, ndims.print = 1:3, nfeatures.print=5, npcs=nPCS, verbose=FALSE)
PCx <- length( attributes( S3_subCLS[["pca"]])\$stdev)
(PCx)

p4 <- ElbowPlot( S3_subCLS, ndims = PCx, reduction = "pca")

# UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
S3_subCLS <- RunUMAP( S3_subCLS, dims = 1:PCx, verbose=FALSE)
S3_subCLS <- RunTSNE( S3_subCLS, dims=1:PCx)

### Clustering based on PCA distance
S3_subCLS <- FindNeighbors( S3_subCLS, dims = 1:PCx, verbose=FALSE)
S3_subCLS <- FindClusters( S3_subCLS, resolution = Resolution, verbose=FALSE)

### write original cluster info
S3_subCLS\$orig.seurat_clusters <- S3_subCLS\$seurat_clusters
S3_subCLS\$orig.seurat_clusters <- S3\$seurat_clusters[ Cells(S3_subCLS) ]
remove(S3)

### find markers at a cluster compared to others. Check the rownames
### The row.names are not genename. gene_name+auto_increased_numbers
## FindMarkers (1 vs. 1)
## FindConservedMarkers (2 vs. others)
## 


S3_subCLS.markers <- FindAllMarkers( S3_subCLS, only.pos = TRUE 
			      ,min.pct = 0.25
			      ,test.use = "wilcox"
			      ,logfc.threshold = logFC
			      ,return.thres = 0.01
			      ,min.cells.feature=minCell
			      ,min.cells.group=minCell)
idx <- S3_subCLS.markers[, "p_val_adj"] < qvThres
S3_subCLS.markers <- S3_subCLS.markers[idx,]



# Store current Seurat3 object
save(S3_subCLS, file="$Rdata_file")
save(S3_subCLS.markers, file="$S3marker_file")

#----------------- Plot

# at here, after clustering
p5 <- DimPlot( S3_subCLS, reduction = "pca", label=TRUE, pt.size=2, label.size=8) + NoLegend() + theme(aspect.ratio=1)
p6 <- DimPlot( S3_subCLS, reduction = "umap", label=TRUE, pt.size=2,label.size=8) + NoLegend() + theme(aspect.ratio=1)
p9 <- DimPlot( S3_subCLS, reduction = "tsne", label=TRUE, pt.size=2,label.size=8) + NoLegend() + theme(aspect.ratio=1)

hmTOPx <- max( summary(S3_subCLS.markers\$cluster) )
topGenes <- S3_subCLS.markers %>% group_by(cluster) %>% top_n(n = hmTOPx, wt = avg_logFC)
p10 <-  DoHeatmap(S3_subCLS, features=topGenes\$gene) + NoLegend()

#------------ Plot
p1 + p2 # Average Expression and Standarized Variance

p4 #PC loading
p5 #PCA
p6 #UMAP

#######
remove(p1)
remove(p2)
remove(p4)
remove(p5)
remove(p6)
#######



#----------------------------------------------------
# Now, we run Monocle3
#----------------------------------------------------

gene_annotation <- as.data.frame(rownames(S3_subCLS\@reductions[["pca"]]\@feature.loadings),
		row.names = rownames(S3_subCLS\@reductions[["pca"]]\@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(S3_subCLS\@assays[["RNA"]]\@counts\@Dimnames[[2]],
	      row.names = S3_subCLS\@assays[["RNA"]]\@counts\@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
cat( "Number of Genes\\t", nrow(gene_annotation), "\\n")
cat("Number of Cells\\t", nrow(cell_metadata), "\\n")


New_matrix <- S3_subCLS\@assays[["RNA"]]\@counts #this is raw counts
###New_matrix <- S3_subCLS\@assays[["RNA"]]\@scale.data

New_matrix <- New_matrix[rownames(S3_subCLS\@reductions[["pca"]]\@feature.loadings), ]
expression_matrix <- New_matrix
cat(" Expression Matrix\\t", nrow(expression_matrix), "x", ncol(expression_matrix), "\\n")

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#summary( attributes(cds))

recreate.partition <- c(rep(1, length(cds\@colData\@rownames)))
names(recreate.partition) <- cds\@colData\@rownames
recreate.partition <- as.factor(recreate.partition)

cds\@clusters\@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info
list_cluster <- S3_subCLS\@meta.data[["seurat_clusters"]]
names(list_cluster) <- S3_subCLS\@assays[["RNA"]]\@data\@Dimnames[[2]]
cds\@clusters\@listData[["UMAP"]][["clusters"]] <- list_cluster

cds\@clusters\@listData[["UMAP"]][["louvain_res"]] <- "NA"
### Assign UMAP coordinate
cds\@int_colData\@listData\$reducedDims\@listData[["UMAP"]] <-S3_subCLS\@reductions[["umap"]]\@cell.embeddings

### Assign feature loading for downstream module analysis
cds\@preprocess_aux\$gene_loadings <- S3_subCLS\@reductions[["pca"]]\@feature.loadings
cds <- learn_graph(cds, use_partition = T)
save(cds, file="$CDS_file")

# Trajectory plot
p7 <- plot_cells(cds, reduction_method="UMAP", group_label_size=8, graph_label_size=5, cell_size=2, show_trajectory_graph=TRUE)+theme(aspect.ratio=1)
p7 #UMAP + Trajectory
remove(p7)

#----------------------------------------------
# Markers from Monocle3 "top_markers" using Jenson-Shannon
# A gene will be included in multiple cluster as a marker
#----------------------------------------------
geneTest <- 100
if( length(summary(S3_subCLS.markers\$cluster)) > 20 ){
  geneTest <- 50
}
(geneTest)

M3.markers <- top_markers( cds,
			   group_cells_by = "cluster",
			   genes_to_test_per_group = geneTest,
			   reduction_method = "UMAP",
			   marker_sig_test = TRUE,
			   reference_cells = NULL,
			   cores = 2,
			   verbose = FALSE )

M3.markers <- M3.markers %>% filter(fraction_expressing >= 0.10) %>% filter(marker_test_q_value <qvThres)
save(M3.markers, file="$M3marker_file")

# Rank top 10
top_M3 <- M3.markers %>% group_by(cell_group) %>% top_n(10, pseudo_R2)
top_M3_ids <- unique(top_M3 %>% pull(gene_id))

# Plot
plotM <- plot_genes_by_group(cds, top_M3_ids, group_cells_by="cluster", max.size=5, ordering_type="cluster_row_col")
plotM + labs( title="Monocle Markers (Top10)" ) + theme(axis.text.\y=element_text(size=7), title=element_text(size=15) )
remove(plotM)

for(cls in unique( M3.markers[,"cell_group"] ) ){
    idx <- M3.markers[,"cell_group"] == cls
    if( sum(idx) < 1){ next }

    outfile <- paste0("$this_cls_dir", "M3_markers_Cls_", cls, ".txt")
    tbl <- M3.markers[idx, c("cell_group", "marker_score", "mean_expression", "fraction_expressing", "specificity", "pseudo_R2", "marker_test_p_value", "marker_test_q_value", "gene_id") ]

    idx <- sort(tbl[, "pseudo_R2"], decreasing=TRUE, index.return=TRUE)\$ix 
    write.table( unique(tbl[idx,]), outfile, quote=F, col.names=T, append=F, sep="\\t", row.names=F)
}
#-------------------------------------------




#-----------------------------------------------------------------
# DEGs based on the UMAP clusters (Knn)
#-----------------------------------------------------------------
cds_pr_test <- graph_test( cds, neighbor_graph="knn", cores=1, verbose=FALSE)
colnames( cds_pr_test )

cds_pr_deg_ids <- row.names( subset( cds_pr_test, q_value < 0.05))
# co-regulated genes using Louvain community analysis
gene_module_umap <- find_gene_modules( cds[cds_pr_deg_ids,], resolution=c(10^seq(-6,-1)), verbose=FALSE )

# Building module matrix
Cell <- names( clusters(cds) )
Cell_Clusters <- stringr::str_c("Cluster ", clusters(cds))
names(Cell_Clusters) <- Cell
cell_group_df <- tibble::tibble( cell=Cell, cell_group=Cell_Clusters)

# used for plotting a heatmap: Module (row) x Cluster (col)
agg_mat <- aggregate_gene_expression(cds, gene_module_umap, cell_group_df)
myClsModule("Based on Monocle Knn", "$this_Umap_module_dir", agg_mat, gene_module_umap, S3_subCLS.markers )

remove(cell_group_df)
remove(Cell)
remove(Cell_Clusters)
remove(agg_mat)
remove(cds_pr_deg_ids)
remove(cds_pr_test)


#-----------------------------------------------------------------
# DEGs based on the Trajectory graph (Principal_graph)
# to test whether cells at similar positions on the trajectory have correlated expression
#-----------------------------------------------------------------
cds_pr_test <- graph_test( cds, neighbor_graph="principal_graph", cores=1, verbose=FALSE)
colnames( cds_pr_test )

cds_pr_deg_ids <- row.names( subset( cds_pr_test, q_value < 0.05))
# co-regulated genes using Louvain community analysis
gene_module_traj <- find_gene_modules( cds[cds_pr_deg_ids,], resolution=c(10^seq(-6,-1)), verbose=FALSE )
# Refer gene_module_df\$module -> "gene" module id

# Building module matrix
Cell <- names( clusters(cds) )
Cell_Clusters <- stringr::str_c("Cluster ", clusters(cds))
names(Cell_Clusters) <- Cell
cell_group_df <- tibble::tibble( cell=Cell, cell_group=Cell_Clusters)

# used for plotting a heatmap: Module (row) x Cluster (col)
agg_mat <- aggregate_gene_expression(cds, gene_module_traj, cell_group_df)
myClsModule("Based on Monocle Trajectory", "$this_Traj_module_dir", agg_mat, gene_module_traj, S3_subCLS.markers )

#----------------------------------------
# End of Finding Modules (co-regulated genes) based on Monocle Trajectory
#----------------------------------------
remove(cell_group_df)
remove(Cell)
remove(Cell_Clusters)
remove(agg_mat)
remove(cds_pr_deg_ids)
remove(cds_pr_test)

remove(cds)
#----------------------- END of Monocle3 --------------------

p9 # TSNE
p10 # Heatmap for Clustering (Seurat3)

remove(p9)
remove(p10)


#===========================
## Write Seurat Marker Information
#===========================
no_cells_eachCluster <- S3_subCLS\@meta.data %>% group_by(seurat_clusters) %>% summarize(count=n())
colnames(no_cells_eachCluster)[1] <- "clusters"

no_markers_eachCluster <- summary(S3_subCLS.markers\$cluster)
total_markers <- length(S3_subCLS.markers\$cluster)

#----------- Write to numbers. txt
#cat("##Number of cells\\n", file="$stat_file", append=FALSE)
cat("##Number of cells\\n", file="$stat_file", append=TRUE)
write.table( no_cells_eachCluster, file="$stat_file", quote=F, col.names=T, append=T, sep="\\t", row.names=F)
cat("\\n", file="$stat_file", append=TRUE)

cat("##Number of cluster markers (Seurat3)\\n", file="$stat_file", append=TRUE)
write.table( no_markers_eachCluster, file="$stat_file", quote=F, col.names=F, append=T, sep="\\t", row.names=T)
cat("\\n", file="$stat_file", append=TRUE)

cat("## Knn module informaion (Monocle3)\\n", file="$stat_file", append=TRUE)
cat("## Module\\tGenes\\n", file="$stat_file", append=TRUE)
write.table( summary(gene_module_umap\$module), file="$stat_file", quote=F, col.names=F, append=T, sep="\\t", row.names=T)
cat("\\n", file="$stat_file", append=TRUE)

cat("## Trajectory module informaion (Monocle3)\\n", file="$stat_file", append=TRUE)
cat("## Module\\tGenes\\n", file="$stat_file", append=TRUE)
write.table( summary(gene_module_traj\$module), file="$stat_file", quote=F, col.names=F, append=T, sep="\\t", row.names=T)
cat("\\n", file="$stat_file", append=TRUE)


### Look at Pitx2 Foxc1 Foxc2 Tfap2b Ngfr
genes <- c("Pitx2", "Foxc1", "Foxc2", "Tfap2b", "Lmx1b", "Emx2", "Ngfr", "Rax", "Vsx2")
p11 <- VlnPlot(S3_subCLS, features=genes,
	    ncol=3,  log=T, pt.size=0.1,
            same.\y.lims=TRUE,  combine=T) + NoLegend()
plot( p11 + labs( caption="Markers") )
remove(p11)

### Write marker genes and plot as violin (top 12) Why top12?? (3x4 plots / page)
for(cls in 0: (length(no_markers_eachCluster)-1) ){
    idx <- S3_subCLS.markers[, "cluster"] == cls
    if( sum(idx) < 1){ next }

    # make output filename
    outfile <- paste0("$this_cls_dir", "S3_markers_Cls_", cls, ".txt")
    # make table and sort
    tbl <- S3_subCLS.markers[idx, c("cluster", "p_val", "p_val_adj", "avg_logFC", "gene")]
    idx <- sort( tbl[, "avg_logFC"], decreasing=TRUE, index.return=TRUE )\$ix 
    # write it
    write.table( unique(tbl[idx,]), outfile, quote=F, col.names=T, append=F, sep="\\t", row.names=F )


    ### plot violin with TOP12 (why 12? 3x4 plots / page)
    if( length(idx) > 11 ){
       genes <- tbl[ idx[1:12], "gene"]
    }else{
       genes <- tbl[ idx, "gene"]
    }
    
    p12 <- VlnPlot(S3_subCLS, features=genes,
		   ncol=3,
		   log=T, pt.size=0.1,
		   same.\y.lims=TRUE, combine=T) + NoLegend()

    plot( p12 + labs(caption=paste0("Cluster ", cls, " (DEGs)")) )

    remove(p12)
    remove(tbl)
    remove(idx)
}


remove(S3_subCLS)
remove(S3_subCLS.markers)

dev.off()
q()

# Raw Counts
##VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

EOF
close(OUT);
# End of R command lines

# Run it
system("R CMD BATCH --no-save $r_cmd $r_out");
undef(@data);
@data = `tail -3 $r_out`;

#system("cat $r_out");
if( !($data[0] =~ /^\>\sproc\.time.*/) ){
    system("cat $r_out");
    exit;
}

#unlink($r_cmd);
#unlink($r_out);

}
