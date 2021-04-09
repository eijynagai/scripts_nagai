#!/usr/bin/env perl
 
require("/home/park/PERL/park_basic.pl");

if(@ARGV != 2){
    print "%>perl " . $0 . " [Input_Dir] [Output_Dir]\n";
    exit;
}

$indir = shift(@ARGV);
$outdir = shift(@ARGV);

if( !defined($indir = &Chk_Dir($indir)) ){ exit; }

## Check the number 
$number_file = $indir . "numbers.txt";
if( &Chk_File($number_file) == 0){ exit; }

#mkdir and cleanup and mkdir
if(-e $outdir){
    $tmp = $outdir . "/Cluster*";
    system("rm -f $tmp");
    $tmp = $outdir . "/Freq_*";
    system("rm -f $tmp");
}
$outdir = &myMKDIR($outdir);

$ClsMarker_dir = $indir . "ClsMarkers/";
$Module_UMAP_dir = $indir . "Modules_UMAP/";
$Module_Traj_dir = $indir . "Modules_Traj/";

if( !defined($ClsMarker_dir = &Chk_Dir($ClsMarker_dir)) ){ exit; }
if( !defined($Module_UMAP_dir = &Chk_Dir($Module_UMAP_dir)) ){ exit; }
if( !defined($Module_Traj_dir = &Chk_Dir($Module_Traj_dir)) ){ exit; }

#----------------------------------------
#Get Cluster IDs in this directory
undef(@CLSs); #keep this order
&Find_ClusterIDs();
#print join("\t", @CLSs) . "\n";

undef(%MARKERS); # $MARKERS{ gene }[cluster_id] = value

$target_dir = $ClsMarker_dir;
&Get_ClusterMarker1();

$target_dir = $Module_UMAP_dir;
&Get_ClusterMarker2();

$target_dir = $Module_Traj_dir;
&Get_ClusterMarker2();

### Write and Seperate
$Freq_Matrix = $outdir . "Freq_clsMarkers_Matrix.txt";
&Print_Markers($Freq_Matrix);
print $Freq_Matrix . "\n";
if( &Chk_File($Freq_Matrix) == 0){ exit; }

$flag = 0;
$subCLS_mode = 0;
$S3_file = $indir . "Seurat3.RData";
if( &Chk_File($S3_file) == 0){ $flag++; }

if($flag != 0){
    $flag = 0;
    $subCLS_mode = 1;
    $S3_file = $indir . "Seurat3_subCLS.RData";

    if( &Chk_File($S3_file) == 0){
	$flag++;
    }else{
	print $S3_file . " found\n";
    }
}
if($flag){ exit; }


$flag = 0;
$CDS_file = $indir . "CDS.RData";
if( &Chk_File($CDS_file) == 0){ $flag++; }
if($flag != 0){ #not found, but is it subCLS?
    $flag = 0;
    $CDS_file = $indir . "CDS_subCLS.RData";
    if( &Chk_File($CDS_file) == 0){
	$flag++;
    }else{
	print $CDS_file . " found\n";
    }
}
if($flag){ exit; }


$thres   = 2; #>2 only
#$thres   = 3;
$ps_file = $outdir . "Plot_cells.ps";
&R_PlotCell( $CDS_file, $S3_file, $Freq_Matrix,  $thres, $outdir, $ps_file);

exit;


sub R_PlotCell{
    my ($CDS, $S3, $input_file, $thres, $this_outdir, $ps_file) = @_;

    my $r_cmd = $this_outdir . "R.cmd";
    my $r_out = $this_outdir . "R.out";
    
## Prepare R command script 
    open(OUT, ">$r_cmd");
print OUT <<EOF;
    
library(Seurat)
library(monocle3)
library(dplyr)
library(cowplot)
library(ggplot2)

library(cluster)
library(reshape2)
library(scales)
library(plyr)


load("$CDS")
load("$S3")

psize=1; lsize= 7;
if( $subCLS_mode ){
    psize=2; lsize= 8;
    S3 <- S3_subCLS
}

tbl   <- as.matrix( read.table("$input_file", header=T, sep="\\t", row.names=1) )

tbl [ tbl <= $thres ] <- 0 
idx <- apply(tbl, 1, sum) > 0
tbl <- tbl[ idx, ]

idx <- apply(tbl, 2, sum) > 0
tbl <- tbl[, idx]

#sort based on cluster 0
# because of the ordered based on clustering, can skip
idx <- sort(tbl[,1], decreasing=TRUE,index.return=TRUE)\$ix
tbl <- tbl[idx,]

## For clustering
my.hclust<- function(d) hclust(d, method="average")
my.dist  <- function(x) dist(x, method="euclidean")

A  <- tbl
A  <- t(scale(t(A), center=TRUE, scale=TRUE))

#------------------------#
HC  <- my.hclust( my.dist(A) )
#ORDER by rowname (from bottom to top in heatmap)
CLSROW <- HC\$label[ HC\$order]

CM <- cor(A)
HC2 <- my.hclust( as.dist( 1 - CM) ) 
#ORDER by colname (from ???)
CLSCOL <- HC2\$label[ HC2\$order]
#------------------------#


M <- as.data.frame( tbl )
M <- cbind(rownames(M), M)

colnames(M)[1] <- "Gene"
M.\m <- melt(M)
colnames(M.\m)[2] <- "Cluster"

# reordering M
#M.m\$Cluster <- factor(M.m\$Cluster, levels=as.character(CLSCOL) )
M.m\$Cluster <- factor(M.\m\$Cluster, levels=rev( levels(M.\m\$Cluster) ) )
M.m\$Gene <- factor(M.m\$Gene, levels=as.character(CLSROW) )

base_size  <- 6
label_size <- 22

#----------- Plot heatmap with Detection frequency
# this plot shows only over thres.freq. and may be overlaped between clusters (one gene can be included at several clusters)
p <- ggplot(M.\m, aes(\y=Cluster, x=Gene, fill=value)) + geom_tile() + scale_fill_gradient( low="white",high="red", name="Freq") 
p <- p + theme_gray(base_size=label_size)
p <- p + labs(y="", x="Consensus DEGs") + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0) )

p <- p + theme(axis.ticks=element_blank(),
    axis.text.x=element_text(size=base_size*1.0, angle=90, hjust=1, vjust=0.5, colour="gray"),
    axis.text.y=element_text(size=base_size*1.7, colour="black", hjust=1), 
    legend.title=element_text(size=12), legend.text=element_text(size=12)  )

p1 <- DimPlot( S3, reduction = "umap", label=TRUE, pt.size=psize,label.size=lsize) + NoLegend() + theme(aspect.ratio=1)

postscript("$ps_file", horizontal=T)
op <- par(cex=1.2)

p1
remove(p1)

p
remove(p)



#-------------- This is for ClusterX_consMakers.txt
#-------------- only a gene for only one clusters

tbl [ tbl <= $thres ] <- 0 
idx <- apply(tbl, 1, sum) > 0
tbl <- tbl[ idx, ]
idx <- apply(tbl, 2, sum) > 0
tbl <- tbl[, idx]

## tbl has been already cutted by thres
## convert to binary (i.e. anyway detected genes )
tbl_binary <- tbl
tbl_binary[ tbl_binary != 0] <- 1

MD <- c()
CLSs <- colnames(tbl) 
all_found_genes <- c()

for(i in 1:length(CLSs)){

    ID  <- CLSs[i] #cluster ID from colnames
    this_cls_id <- as.numeric( gsub("Cluster","", ID) ) # cluster ID as numeric
    this_cls_all_genes <- c()
    this_cls_specific_genes <- c()
    partner_cluster_IDs <- c()
    this_cls_other_marker_genes <- c()
    genes_multiple_cls <-  c()

    # pick up all genes found for this clusters (value == 1)
    # Is there any possibility that this is empty?
    # No.. a column includes at least 1.
    # If the selected rows from a matrix is one, the new object becomes a vector
    idx <- tbl_binary[, ID] == 1
    if(sum(idx) > 0){
           this_cls_all_genes <- tbl_binary[ idx, ] #matrix or vector
    }
if(is.vector( this_cls_all_genes) ){ #means idx == 1
    this_cls_all_genes <- t(as.matrix( this_cls_all_genes ) )
    rownames(this_cls_all_genes) <- rownames(tbl_binary)[idx]
}

    # pick up this cluster specific genes
    idx <- apply( this_cls_all_genes, 1, sum) == 1
    if(sum(idx) > 0){
       # gene_name vector!!
       this_cls_specific_genes <-  rownames( this_cls_all_genes)[idx]
    }

    # multiple cluster genes
    #remove this cluster column
     tmp <- this_cls_all_genes[ !idx, -i]
    if( is.vector(tmp) ){
       tmp <- t( as.matrix(tmp) )
       rownames(tmp) <- rownames( this_cls_all_genes)
    }

if(nrow(tmp) > 0 && ncol(tmp) > 0){
    idx <- apply(tmp, 2, sum) != 0 #remove col. having all zero 
    genes_multiple_cls <-  tmp[,idx] # all partner of this_cls_id
    partner_cluster_IDs <- as.numeric( gsub("Cluster","", colnames(genes_multiple_cls)) )
}
    remove(tmp)

    # markers for this and other cls
    if( length(partner_cluster_IDs) > 0 ){

       for( j in 1:length(partner_cluster_IDs) ){
	   # do not use "j" directly
	   partner_ID <- colnames(genes_multiple_cls)[j] ###"Clusterxxx"
	   partner_cls_id <- partner_cluster_IDs[j]             ### "xxx" only

           idx <- genes_multiple_cls[, partner_ID] != 0 # pickup genes observed from cluster i and j
           if(sum(idx) < 1){ next; }

	   # to be confirmed
           target_genes  <- rownames( genes_multiple_cls)[idx] #target genes of pair this_cls_id -- j
	   Test_TwoCls   <- FindMarkers(S3,
				      this_cls_id, partner_cls_id
				      ,grouping.var="seurat_clusters"
				      ,only.pos=FALSE
				      ,test.use="wilcox"
				      ,min.pct = 0.25
				      ,logfc.threshold = logFC
				      ,return.thres = 0.01
				      ,min.cells.feature=minCell
				      ,min.cells.group=minCell
				      ,verbose=FALSE, assay="RNA")
	   # If positive logFC means high in the first, negative FC means high in the second
	   # If not found (<- ""only.pos==false"" important!!!!), these target_genes are markers for BOTH
	   found_genes <- Test_TwoCls[ target_genes, "avg_logFC"]
           names( found_genes ) <- target_genes

	   # for both  clusters
	   tmp <- names(found_genes)[ is.na(found_genes) ]
	   if( length( tmp ) > 0){
	       this_cls_other_marker_genes <- c( this_cls_other_marker_genes, tmp)
           }
           remove(tmp)

	   # find positive FC. for this cluster
	   found_genes[ found_genes < 0 ] <- NA
	   found_genes <- names(found_genes)[ !is.na(found_genes) ]
	   if( length( found_genes) < 1){ next; }
	   # Add these genes to 
	   this_cls_specific_genes <- c( this_cls_specific_genes, found_genes)
	   remove(found_genes)
      }

   }

   if( length(this_cls_specific_genes) > 0){
       tmp <- cbind( this_cls_specific_genes, ID) # gene_name, cluster_id
       MD  <- rbind( MD, tmp)
       remove(tmp)
   }
   if( length(this_cls_other_marker_genes) > 0){
        tmp <- cbind( this_cls_other_marker_genes, ID)
        MD  <- rbind( MD, tmp) # gene_name, cluster_id
        remove(tmp)
   }

   tmp <- c()
   if( length(this_cls_specific_genes) > 0){
       tmp <- cbind( this_cls_specific_genes, 2) # gene_name, cluster_id
   }
   if( length(this_cls_other_marker_genes) > 0){
        tmp2 <- cbind( this_cls_other_marker_genes, 1) # gene 2 or gene 1 (2: specific for thie cluster, 1: by testing with partners)
        tmp <- rbind(tmp, tmp2)
        remove(tmp2)
   }

#if( length(tmp) < 1 ){ next }
if( nrow(tmp) < 1){ next }

    colnames(tmp) <- c("Gene", ID);
#idx <- sort(tmp[,2], decreasing=TRUE, index.return=TRUE)\$ix

    outfile   <- paste0("$this_outdir", ID, "_consMarkers.txt")
    sink(file=outfile, append=FALSE)
    cat ("## These genes are detected by methods at least > X times"); cat("\n")
    cat ("## These genes are specific for this cluster if the col. value is 2"); cat("\n")
    cat ("## These genes are for multiple clusters if the col. value is 1"); cat("\n")
    sink()
#    write.table( tmp[idx,], outfile, quote=F, col.names=T, append=T, sep="\\t", row.names=F)
    write.table( tmp, outfile, quote=F, col.names=T, append=T, sep="\\t", row.names=F)

    all_found_genes <- c(all_found_genes, tmp[, "Gene"])

    remove(tmp)
    remove( this_cls_specific_genes )
    remove( this_cls_other_marker_genes )
    remove(idx)
    remove( genes_multiple_cls )
    remove( partner_cluster_IDs )
}

# Heatmap
p3 <- DoHeatmap(S3, features= unique(all_found_genes), hjust=0, angle=0, size=2)
p3 <- p3 + theme(text=element_text(size=3), legend.title=element_text(size=10), legend.text=element_text(size=7) )
##  + NoLegend()
p3
remove(p3)

#MD[ MD == 0] <- NA
plot_cells( cds, genes=MD, show_trajectory_graph=FALSE,  reduction_method = "UMAP",
	    alpha=1, 
	    label_cell_groups=TRUE, cell_size=0.5,
	    group_label_size=4, 
#	    min_expr=3 ) + scale_colour_gradient2( low="orange", mid="cyan", high="brown", name="Expression")
	    min_expr=3 ) + scale_colour_gradient( na.value="white", low="cyan", high="brown", name="Expression") + theme( aspect.ratio=1)

remove(MD)
dev.off()
q()

EOF
close(OUT);
# End of R command lines

# Run it
  system("R CMD BATCH --no-save $r_cmd $r_out");
  undef(@data);
  @data = `tail -3 $r_out`;
    
#    system("cat $r_out");
  if( !($data[0] =~ /^\>\sproc\.time.*/) ){
      system("cat $r_out");
      exit;
  }

  unlink($r_cmd);
  unlink($r_out);
}




sub Get_ClusterMarker1{
    my(@files, $file, $gene, $cluster, @data, $i);

    undef(@files);
    @files = glob($target_dir . "*_markers_Cls_*.txt"); #M3 and S3, markers_Cls_xx.txt
    if(scalar(@files < 1)){
	print "Zero M3 or S3 files in " . $target_dir . "\n";
	exit;
    }
    
    foreach $file (@files){
	$file =~ /.*(M3|S3)_.*markers_Cls_(\d+)\.txt/;
	$method  = $1;
	$cluster = $2;
	if( &Chk_File($file) == 0){ exit; }
	
	# last field is gene id
	$line = `head -1 $file`;
	chop($line);
	undef(@data);
	@data = split(/\t/, $line);
	$last_field = scalar(@data);
	
	undef(@data);
	@data = `cut -f $last_field $file | grep -v gene`;
	chop(@data);
	
	foreach $gene (@data){
	    if( !defined($MARKERS{$gene})){
		for($i=0; $i<scalar(@CLSs); $i++){
		    $MARKERS{ $gene }[$CLSs[$i]] = 0;
		}
	    }
	    $MARKERS{ $gene }[$cluster] += 1;

	    #if($gene eq "Pitx2"){
	#	print $gene . "\t" . $cluster . "\t" . $file . "\n";
	 #   }
	}
    }

}



sub Get_ClusterMarker2{
    my(@files, $file, $gene, $cluster, @data, $i);

    undef(@files);
    @files = glob($target_dir . "*Cluster*.txt");
    if(scalar(@files < 1)){
	print "Zero M3 or S3 files in " . $target_dir . "\n";
	exit;
    }

    foreach $file (@files){
	$file =~ /.*Cluster(\d+)_markerModule.*\.txt/;
	$cluster = $1;
	if( &Chk_File($file) == 0){ exit; }

	# last field is gene id
	$line = `head -1 $file`;
	chop($line);
	undef(@data);
	@data = split(/\t/, $line);
	$last_field = scalar(@data);

	undef(@data);
	@data = `cut -f $last_field $file | grep -v gene`;
	chop(@data);
	
	foreach $gene (@data){
	    if( !defined($MARKERS{$gene})){
		for($i=0; $i<scalar(@CLSs); $i++){
		    $MARKERS{ $gene }[$CLSs[$i]] = 0;
		}
	    }
	    $MARKERS{ $gene }[$cluster] += 1;

	    #if($gene eq "Pitx2"){
	#	print $gene . "\t" . $cluster . "\t" . $file . "\n";
	 #   }

	}

    }

}




sub Print_Markers{
    my $outfile = $_[0];
    my ($i, $gene );

    open(PM, ">$outfile");
    print PM "Gene";
    for($i=0; $i<scalar(@CLSs); $i++){
	print PM "\t" . "Cluster" . $CLSs[$i];
    }
    print PM "\n";
    
    foreach $gene (sort {$a cmp $b} keys %MARKERS){
	$outline = $gene;
	for($i=0; $i<scalar(@CLSs); $i++){
	    $outline .= "\t" . $MARKERS{ $gene }[ $CLSs[$i] ];
	}
	$outline .= "\n";
	print PM $outline;

	#if($gene eq "Pitx2"){
	 #   print $outline . "\n";
	#}
    }
    close(PM);
}


#&Find_ClusterIDs();

sub Find_ClusterIDs{
    undef(@CLSs);

    undef(@files);
    @files = glob($ClsMarker_dir . "*_markers_Cls_*.txt"); #M3 and S3, markers_Cls_xx.txt
    if(scalar(@files < 1)){
	print "Zero M3 or S3 files in " . $ClsMarker_dir . "\n";
	exit;
    }

    foreach $file (@files){
	$file =~ /.*(M3|S3)_markers_Cls_(\d+).txt/;
	$method  = $1;
	$cluster = $2;
	push(@CLSs, $cluster);
    }
    @CLSs = &Get_Unique_Elements_Sort_Num( @CLSs);
    
    undef(@files);
    @files = glob($Module_UMAP_dir . "Cluster*_markerModule*.txt");
    if(scalar(@files < 1)){
	print "Zero Cluster marker files in " . $Module_UAMP_dir . "\n";
	exit;
    }
    foreach $file (@files){
	$file =~ /.*Cluster(\d+)_markerModule.*\.txt/;
	$cluster = $1;
	push(@CLSs, $cluster);
    }
    @CLSs = &Get_Unique_Elements_Sort_Num( @CLSs);
    
    
    undef(@files);
    @files = glob($Module_Traj_dir . "Cluster*_markerModule*.txt");
    if(scalar(@files < 1)){
	print "Zero Cluster marker files in " . $Module_Traj_dir . "\n";
	exit;
    }
    foreach $file (@files){
	$file =~ /.*Cluster(\d+)_markerModule.*\.txt/;
	$cluster = $1;
	push(@CLSs, $cluster);
    }
    @CLSs = &Get_Unique_Elements_Sort_Num( @CLSs);
}