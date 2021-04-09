#!/usr/bin/env perl

require("/home/park/PERL/park_basic.pl");

# Input is CellRanger "cellranger_outputs/filtered_feature_bc_matrix/"
# Create Output directory
if(@ARGV != 5){
    print "%>perl " . $0 . " [Input: filtered_feature_bc_matrix] [Output Directory] [Resolution] [Hbb?0|1] [Tubb?0|1]\n";
    exit;
}

$CR_INPUT_DIR = shift(@ARGV);
$SR_OUTPUT_DIR = shift(@ARGV); 
$RESOLUTION = shift(@ARGV);
$HBB   = shift(@ARGV);
$TUBB = shift(@ARGV);

if( !defined($CR_INPUT_DIR = &Chk_Dir($CR_INPUT_DIR)) ){ exit; } 
$SR_OUTPUT_DIR = &myMKDIR($SR_OUTPUT_DIR);

#---------- Configuration Sets
undef(%PARAM);  #paramter sets
undef(%DESC);     #description for each paramter

# UMI : total RNAs sequenced, nFeature: number of genes
$DESC{ "CR_DIR" }    = "Input dir (Cell Ranger)";
$PARAM{ "CR_DIR" } = "\"" . $CR_INPUT_DIR . "\"";

$DESC{ "Prj" }    = "Project Name";
$PARAM{ "Prj" } = "\"mNCC\"";

$DESC{ "UMI_Cutoff" } = "Gene expression threshold";
$PARAM{ "UMI_Cutoff" } = 1; 

$DESC{ "minCell" } = "Minimal number of cells that a gene is expressed";
$PARAM{ "minCell" } = 5; 

$DESC{ "minFeatures" } = "Minimal number of genes";
$PARAM{ "minFeatures" } = 200; 

$DESC{ "MT_Pattern" }   = "Capturing mitochondrial genes";
## mouse case
$PARAM{ "MT_Pattern" } = "\"^mt-\"";

$DESC{ "MT_Max" }   = "Apoptosis: Maximal percentage of mitochondrial gene expression";
$PARAM{ "MT_Max" } = 25;

$DESC{ "MT_Min" }   = "Apoptosis: Minimal Percentage threshold of mitochondrial gene expression";
$PARAM{ "MT_Min" } = 0.1;

$DESC{ "Gene_Cutoff" }  = "Cell doublets, multiplets: nFeature_RNA = Number of genes ";
$PARAM{ "Gene_Cutoff" } = 8000;

$DESC{ "Norm_Scale" } = "Scaling factor for normalization of gene expression";
$PARAM{"Norm_Scale"} = 10000;

$DESC{ "TOPx"} = "Top x labels of diff. genes to be shown";
$PARAM{ "TOPx" } = 10;

$DESC{ "Resolution" }  = "Granularity of cluster resolution";
$PARAM{ "Resolution" } = $RESOLUTION;

$DESC{ "nPCS" }  = "Number of PCs at RunPCA";
$PARAM{ "nPCS" } = 100;

$DESC{"logFC"} = "Expression fold-change in log-scale";
$PARAM{"logFC"} = 0.25;  #<<--- log(1.28), 1.28-fold change

$DESC{"qvThres"} = "Threshold of p_adj_val for DEGs";
$PARAM{"qvThres"} = 0.001;

$DESC{"Variable_Feature"} = "Number of diff. genes to be used for clustering";
$PARAM{"Variable_Feature"} = 3000;

$DESC{"SCTransf"} = "Use SCTransform";
$PARAM{"SCTransf"} = 0;

$DESC{"scaleAll"} = "Scaling uses all genes";
$PARAM{"scaleAll"} = 0;

$DESC{"removeHbb"} = "Cells not expressing Hbb-y";
$PARAM{"removeHbb"} = 0;
if( $HBB ){
    $PARAM{"removeHbb"} = 1;
}

$DESC{"removeTubb"} = "Cells not expressing Tubb";
$PARAM{"removeTubb"} = 0;
if( $TUBB ){
    $PARAM{"removeTubb"} = 1;
}

$this_output_dir = $SR_OUTPUT_DIR; ## . $this;
$this_output_dir = &myMKDIR($this_output_dir);

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




$param_file = $this_output_dir . "cfg_params.txt";
$stat_file = $this_output_dir . "numbers.txt";

$ps_file  = $this_output_dir . "plots.ps";

$Rdata_file     = $this_output_dir . "Seurat3.RData";
$S3marker_file = $this_output_dir . "Seurat3.markers.RData";
$CDS_file = $this_output_dir . "CDS.RData";
$M3marker_file = $this_output_dir . "Monocle3.markers.RData";

# Write parameter file
open(OUT, ">$param_file");
foreach $key (keys %PARAM){
    print OUT "## " . $DESC{ $key }  . "\n";
    print OUT $key  . " <- " . $PARAM{ $key } . "\n" . "\n";
}
close(OUT);

unlink($ps_file);
unlink($Rdata_file);
unlink($S3marker_filej);
unlink($CDS_file);
unlink($M3marker_file);


#unlink($res_file);
&Seurat3_R();

exit;

sub Seurat3_R{
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

#Smaple1:Mouse E12.5 eye Neural Crest Cells (E12.5 p75+/ITGA4- fraction)
#Sample2:Human iPSC-derived Neural Crest Cells (Day26 p75+/ITGA4- fraction)
# p75 == Ngfr
# Itga4 <- Blood cells
#眼神経堤細胞のマーカーとしてはPITX2, FOXC1, TFAP2B
# Lmx1B, Emx2
postscript("$ps_file", horizontal=T)
op <- par(cex=1.2)

## Pre-process Seurat object 
## STEP 1: Read 10X Cell Ranger Data and create Seurat Object
CR <- Read10X(data.dir = CR_DIR)
S1 <- CreateSeuratObject(counts=CR, project=Prj, assay="RNA"
			 ,min.cells=minCell
			 ,min.features=minFeatures)
S1 <- PercentageFeatureSet(S1, pattern=MT_Pattern, col.name="Percent_MT")
# For capture abnormal cells
APOPT_S1 <- FeatureScatter(S1, feature1 = "nCount_RNA", feature2 = "Percent_MT", pt.size=0.5, cols="gray")+NoLegend()
DLET_S1 <- FeatureScatter(S1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.5, cols="gray")+NoLegend()
#summary( attributes( S1 ) )
attr( S1, "assays")\$RNA

### Filtering unwanted cells
#S2 <- subset( S1, subset= Itga4 < UMI_Cutoff & Ngfr > UMI_Cutoff)

S2 <- S1
if( removeHbb ){
  S2 <- subset( S2, subset=`Hbb-y` <UMI_Cutoff)
}
if( removeTubb ){
  S2 <- subset( S2, subset=Tubb3 <UMI_Cutoff)
}

S2 <- PercentageFeatureSet(S2, pattern=MT_Pattern, col.name="Percent_MT")
APOPT_S2 <- FeatureScatter(S2, feature1 = "nCount_RNA", feature2 = "Percent_MT", pt.size=0.5, cols="gray")+NoLegend()
DLET_S2 <- FeatureScatter(S2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.5, cols="gray")+NoLegend()
#summary( attributes( S2 ) )
attr( S2, "assays")\$RNA


S3 <- subset(S2, subset=nFeature_RNA > minFeatures & nFeature_RNA < Gene_Cutoff & nCount_RNA < 50000 & Percent_MT < MT_Max & Percent_MT > MT_Min)
S3 <- PercentageFeatureSet(S3, pattern=MT_Pattern, col.name="Percent_MT")
#S3[["Percent_MT"]] <- PercentageFeatureSet(S3, pattern=MT_Pattern)
APOPT_S3 <- FeatureScatter(S3, feature1 = "nCount_RNA", feature2 = "Percent_MT", pt.size=0.5, cols="black")+NoLegend()
DLET_S3 <- FeatureScatter(S3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.5, cols="black")+NoLegend()
#summary( attributes( S3 ) )
attr( S3, "assays")\$RNA

remove(CR)
remove(S1)
remove(S2)

## Pre-process Seurat object
if( SCTransf ){
    S3 <- SCTransform(S3 
		,variable.features.rv.th = 1.3  #rv: residual variance cutoff
		,variable.features.n = NULL
		,vars.to.regress="Percent_MT", 
		return.only.var.genes = TRUE,
		verbose=FALSE)
#warnings()
}else{

    # Step 1: Normalizing
    S3 <- NormalizeData(S3, normalization.method = "LogNormalize", scale.factor = Norm_Scale, verbose=FALSE)
    # Step 2: Selecting variable genes
    S3 <- FindVariableFeatures(S3, selection.method = "vst", nfeatures = Variable_Feature, verbose=FALSE)

    # Step 3: Scaling
    ## ------- scaling with all genes
    if(scaleAll){
	all.genes <- rownames(S3)
        S3 <- ScaleData( S3, vars.to.regress = "Percent_MT", features = all.genes, verbose=FALSE)
    }else{
        ##------- scaling with selected nVF only
	S3 <- ScaleData( S3, vars.to.regress = "Percent_MT", verbose=FALSE )
    }
    #summary( attributes( S3 ) )
    #summary( attributes( S3[["RNA"]]))
}

TOP_Gene <- head(VariableFeatures(S3), TOPx)
TOP_Gene <- c(TOP_Gene, "Foxc1", "Pitx2", "Dkk2", "Tfap2b", "Lmx1b", "Emx2")

p1 <- VariableFeaturePlot(S3)
p2 <- LabelPoints(plot=p1, points=TOP_Gene, repel=TRUE, xnudge=0,  ynudge=0)

#   PCA with nVF
S3  <- RunPCA(S3
#	      ,features = VariableFeatures(object = S3)
	      ,ndims.print = 1:3, nfeatures.print=5, npcs=nPCS, verbose=FALSE)
PCx <- length( attributes(S3[["pca"]])\$stdev)
(PCx)
#summary(attributes(S3[["pca"]]))

#S3 <- JackStraw( S3, dims=PCx, num.replicate = 100)
#S3 <- ScoreJackStraw( S3, dims = 1:PCx)
#p3 <- JackStrawPlot( S3, dims = 1:PCx,reduction="pca")
#p3 # I don't need to plot it
p4 <- ElbowPlot( S3, ndims = PCx, reduction = "pca")

# UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
S3 <- RunUMAP( S3, dims = 1:PCx, verbose=FALSE)
S3 <- RunTSNE( S3, dims=1:PCx)

### Clustering based on PCA distance
S3 <- FindNeighbors( S3, dims = 1:PCx, verbose=FALSE)
S3 <- FindClusters( S3, resolution = Resolution, verbose=FALSE)

# find markers for every cluster compared to all remaining cells,
# logFC >0; highly expressed at the 1st group
# p_val_adj; bonferroni correction
# pct1.; percentage of cells where the gene is detected in the 1st group
# pct2.; percentage of cells where the gene is detected in the 2nd group
S3.markers <- FindAllMarkers( S3, only.pos = TRUE 
			      ,min.pct = 0.25
			      ,test.use = "wilcox"
			      ,logfc.threshold = logFC
			      ,return.thres = 0.01
			      ,min.cells.feature=minCell
			      ,min.cells.group=minCell)
idx <- S3.markers[, "p_val_adj"] < qvThres
S3.markers <- S3.markers[idx,]

## output looks like
#                 p_val                  avg_logFC  pct.1  pct.2  p_val_adj                cluster  gene
#Tfap2b     1.013074e-158 0.9849337 0.956 0.428 1.801853e-154       0          Tfap2b
#Nbl1         2.484890e-106 0.6987402 0.878 0.461 4.419625e-102       0          Nbl1
#Fst            5.046267e-104 0.7543293 0.863 0.425 8.975290e-100       0          Fst

# Store current Seurat3 object
save(S3, file="$Rdata_file") #load("S3.RData")
save(S3.markers, file="$S3marker_file")

#----------------- Plot
# Plots after clustering
APOPT_S4 <- FeatureScatter(S3, feature1 = "nCount_RNA", feature2 = "Percent_MT", pt.size=0.5)+NoLegend()
DLET_S4    <- FeatureScatter(S3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.5)+NoLegend()

# when you dont need legend
#DimPlot( S3, reduction = "umap", label=TRUE, pt.size=0.7,label.size=8) + NoLegend()
#p1<-DimPlot( S3, reduction = "umap", label=TRUE, pt.size=0.7,label.size=8) + ggtitle(label="UMAP")
#p1 <- AugmentPlot(plot = p1)
#p2 <- AugmentPlot(plot = p2)
#(p1 + p2) & NoLegend()

# at here, after clustering
p5 <- DimPlot( S3, reduction = "pca", label=TRUE, pt.size=1, label.size=7) + NoLegend() + theme(aspect.ratio=1)
p6 <- DimPlot( S3, reduction = "umap", label=TRUE, pt.size=1,label.size=7) + NoLegend() + theme(aspect.ratio=1)
p9 <- DimPlot( S3, reduction = "tsne", label=TRUE, pt.size=1,label.size=7) + NoLegend() + theme(aspect.ratio=1)

############## Example for finding marker genes
##VlnPlot( S1, features = c("MS4A1", "CD79A"))
# Raw Counts
#VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
##Cls1.markers <- FindMarkers( S1, ident.1 = 1, min.pct = 0.25)

hmTOPx <- max( summary(S3.markers\$cluster) )
topGenes <- S3.markers %>% group_by(cluster) %>% top_n(n = hmTOPx, wt = avg_logFC)
p10 <-  DoHeatmap(S3, features=topGenes\$gene) + NoLegend()

#------------ Plot
p1 + p2 # Average Expression and Standarized Variance

lb <- c( 'S1 (Apoptosis)', 'S1 (Doublets)' , 'S3 (Apoptosis)', 'S3 (Doublets)', 'S3 (Apoptosis)', 'S3 (Doublets)' )
plot_grid(
     APOPT_S1, DLET_S1 , APOPT_S3, DLET_S3, APOPT_S4, DLET_S4

    ,labels=lb
    ,label_size= 12
    ,nrow=3, ncol=2
)
#------------------------- 

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

gene_annotation <- as.data.frame(rownames(S3\@reductions[["pca"]]\@feature.loadings),
		row.names = rownames(S3\@reductions[["pca"]]\@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(S3\@assays[["RNA"]]\@counts\@Dimnames[[2]],
	      row.names = S3\@assays[["RNA"]]\@counts\@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
cat( "Number of Genes\\t", nrow(gene_annotation), "\\n")
cat("Number of Cells\\t", nrow(cell_metadata), "\\n")


New_matrix <- S3\@assays[["RNA"]]\@counts #this is raw counts
###New_matrix <- S3\@assays[["RNA"]]\@scale.data

New_matrix <- New_matrix[rownames(S3\@reductions[["pca"]]\@feature.loadings), ]
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
list_cluster <- S3\@meta.data[["seurat_clusters"]]
names(list_cluster) <- S3\@assays[["RNA"]]\@data\@Dimnames[[2]]
cds\@clusters\@listData[["UMAP"]][["clusters"]] <- list_cluster

cds\@clusters\@listData[["UMAP"]][["louvain_res"]] <- "NA"
### Assign UMAP coordinate
cds\@int_colData\@listData\$reducedDims\@listData[["UMAP"]] <-S3\@reductions[["umap"]]\@cell.embeddings

### Assign feature loading for downstream module analysis
cds\@preprocess_aux\$gene_loadings <- S3\@reductions[["pca"]]\@feature.loadings
cds <- learn_graph(cds, use_partition = T)
save(cds, file="$CDS_file")

# Trajectory plot
p7 <- plot_cells(cds, reduction_method="UMAP", group_label_size=7, graph_label_size=4, cell_size=1, show_trajectory_graph=TRUE) + theme(aspect.ratio=1)
p7 #UMAP + Trajectory
remove(p7)

# NOTE: partitions in Monocle obejct is only one, but clusters exist
## End of building cds ====================


# Differentially expressed genes across a single-cell trajectory
# Spatial correlation analysis, the Moran’s I test.
# Moran’s I is a measure of multi-directional and multi-dimensional spatial autocorrelation
# The statistic tells you whether cells at nearby positions on a trajectory will have similar (or dissimilar) expression levels for the gene being tested
# Pearson correlation and Moran’s I ranges from -1 to 1,
# +1 means that nearby cells will have perfectly similar expression;
# 0 represents no correlation, and -1 means that neighboring cells will be anti-correlated.

#----------------------------------------------
# Markers from Monocle3 "top_markers" using Jenson-Shannon
# A gene will be included in multiple cluster as a marker
#----------------------------------------------
geneTest <- 100
if( length(summary(S3.markers\$cluster)) > 20 ){
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
myClsModule("Based on Monocle Knn", "$this_Umap_module_dir", agg_mat, gene_module_umap, S3.markers )

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
myClsModule("Based on Monocle Trajectory", "$this_Traj_module_dir", agg_mat, gene_module_traj, S3.markers )

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
no_cells_eachCluster <- S3\@meta.data %>% group_by(seurat_clusters) %>% summarize(count=n())
colnames(no_cells_eachCluster)[1] <- "clusters"

no_markers_eachCluster <- summary(S3.markers\$cluster)
total_markers <- length(S3.markers\$cluster)

#----------- Write to numbers. txt
cat("##Number of cells\\n", file="$stat_file", append=FALSE)
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
p11 <- VlnPlot(S3, features=genes,
	    ncol=3,  log=T, pt.size=0.1,
            same.\y.lims=TRUE,  combine=T) + NoLegend()
plot( p11 + labs( caption="Markers") )
remove(p11)

### Write marker genes and plot as violin (top 12)
for(cls in 0: (length(no_markers_eachCluster)-1) ){
    idx <- S3.markers[, "cluster"] == cls
    if( sum(idx) < 1){ next }

    outfile <- paste0("$this_cls_dir", "S3_markers_Cls_", cls, ".txt")

    tbl <- S3.markers[idx, c("cluster", "p_val", "p_val_adj", "avg_logFC", "gene")]
    idx <- sort(tbl[, "avg_logFC"], decreasing=TRUE, index.return=TRUE)\$ix 
    write.table(unique(tbl[idx,]), outfile, quote=F, col.names=T, append=F, sep="\\t", row.names=F)
 
    if( length(idx) > 11 ){
       genes <- tbl[ idx[1:12], "gene"]
    }else{
       genes <- tbl[ idx, "gene"]
    }

    p12 <- VlnPlot(S3, features=genes,
	    ncol=3,
	    log=T, pt.size=0.1,
            same.\y.lims=TRUE, combine=T) + NoLegend()

    plot( p12 + labs(caption=paste0("Cluster ", cls, " (DEGs)")) )

    remove(p12)
    remove(tbl)
    remove(idx)
}


remove(S3)
remove(S3.markers)

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
