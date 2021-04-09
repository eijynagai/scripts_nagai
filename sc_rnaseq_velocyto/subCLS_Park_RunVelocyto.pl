#!/usr/bin/env perl

require("/home/park/PERL/park_basic.pl");

# Create Output directory
if(@ARGV != 4){
    print "%>perl " . $0 . " [Seurat3 DIR] [DropEst Rds] [Output Directory] [Fresh_Run? 1|0]\n";
    exit;
}

$SR_INPUT_DIR     = shift(@ARGV);  #Seurat3
$DR_INPUT_FILE    = shift(@ARGV);  #DropEst (matrices.rds)
$VR_OUTPUT_DIR = shift(@ARGV); # Velocyto 
$New_RUN            = shift(@ARGV); # 1= creat new new output directory, 0=load Veloctyo.Rdata"

#----------- Check
if( !defined($SR_INPUT_DIR = &Chk_Dir($SR_INPUT_DIR)) ){ exit; }
if( &Chk_File($DR_INPUT_FILE) == 0 ){ exit; }

if( $New_RUN){
    if(-e $VR_OUTPUT_DIR){
	system("rm -fr $VR_OUTPUT_DIR");
    }
}
$this_output_dir = &myMKDIR($VR_OUTPUT_DIR);

$param_file = $SR_INPUT_DIR   . "cfg_params.txt";
$Rdata_S3   = $SR_INPUT_DIR   . "Seurat3_subCLS.RData";
$ConsMarker = $SR_INPUT_DIR . "ConsensusMarkers/Freq_clsMarkers_Matrix.txt";

if( &Chk_File($param_file) == 0 ){ exit; }
if( &Chk_File($Rdata_S3) == 0 ){ exit; }

#------ Create
$pdf_file  = $this_output_dir . "velo_plot.pdf";
$resi_pdf_file  = $this_output_dir . "resid_plot.pdf";
$Rdata_Vel = $this_output_dir . "Velocyto_subCLS.Rdata"; #S4, if $NEW_RUN = 1
$Rdata_Vel_marker = $this_output_dir . "Velocyto_subCLS_marker.Rdata";

#check it
if(!$New_RUN){
    if( &Chk_File($Rdata_Vel) == 0 ){
	print $Rdata_Vel . " not found. Run it\n";
	$New_RUN = 1; # run velocyto
    }
}

#unlink($ps_file);
unlink($pdf_file);
unlink($resi_pdf_file);
&Velocyto_R();

exit;

sub Velocyto_R{
    my $r_cmd = $this_output_dir . "R.cmd";
    my $r_out   = $this_output_dir . "R.out";

## Prepare R command script 
open(OUT, ">$r_cmd");
print OUT <<EOF;

library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(ggplot2)

source("$param_file")

# estimate using top/bottom 2% of cells
fit.quantile <- 0.02
kCells <- 5
cell.alpha <- 0.8
cell.cex <- 1
arrow.lwd <- 1.5
cons_marker <- 2 #over this cons_marker (max, 4)

if($New_RUN == 1){

    load("$Rdata_S3") ## load Seurat3 object

   DrpEst <- readRDS(file = "$DR_INPUT_FILE")

   S4 <- S3_subCLS ## Copy. This is not required indeed.
   S5 <- S3_subCLS ## reduced to marker gene only
   remove(S3_subCLS)

#summary( attributes(S4) )
#attributes(S4)\$assays

### Intersect. All DrpEst genes
   cells <- intersect(  intersect( colnames( attributes(S4)\$assays\$RNA ), colnames( DrpEst\$exon ) ), colnames( DrpEst\$intron ) )
   genes <- intersect( intersect( rownames(S4\@assays\$RNA), rownames( DrpEst\$exon ) ),  rownames( DrpEst\$intron ) )

   S4[["spliced"]]     <- CreateAssayObject( counts = DrpEst\$exon[genes, cells] )
   S4[["unspliced"]] <- CreateAssayObject( counts = DrpEst\$intron[genes, cells])
   S4[["spanning"]]  <- CreateAssayObject( counts = DrpEst\$spanning[, cells])

#summary( attributes(S4) )
#attributes(S4)\$assays
   #remove(genes)
   #remove(cells)

#S4 <- SCTransform(object = bm, assay = "spliced")
   S4 <- RunVelocity(object = S4, deltaT = 1, kCells = kCells, fit.quantile = fit.quantile)
#summary( attributes(S4) )
#summary( attributes(S4)\$tools )
#summary( attributes(S4)\$tools\$RunVelocity)
   save(S4, file="$Rdata_Vel")

   ## Loading consensus markers already detected by Seurat3 and Monocle3
   #$ConsMarker = $SR_INPUT_DIR . "ConsensusMarkers/Freq_clsMarkers_Matrix.txt";
   ## gene cluster0 cluster1 ....
    marker <- read.table("$ConsMarker",  header=T, sep="\\t", row.names=1)
    idx <- apply(marker, 1, sum) > cons_marker
    marker_genes <- rownames(marker)[idx]
    
    #pickup intersect genes
    genes<- intersect( marker_genes, genes)

    S5[["spliced"]]     <- CreateAssayObject( counts = DrpEst\$exon[genes, cells] )
    S5[["unspliced"]] <- CreateAssayObject( counts = DrpEst\$intron[genes, cells])
    S5[["spanning"]]  <- CreateAssayObject( counts = DrpEst\$spanning[, cells])
    S5 <- RunVelocity(object = S5, deltaT = 1, kCells = kCells, fit.quantile = fit.quantile)
    save(S5, file="$Rdata_Vel_marker")

    remove(genes)
    remove(cells)
    remove(marker_genes)
    remove(marker)

}else{

  #load previous S4 that velocity has been performed
  load("$Rdata_Vel")
  load("$Rdata_Vel_marker")
}

plot_velocity <- function (S, arrow.scale) {

  ## Pre-processing for plotting Velocity
  ### Clusters
   ident.colors <- (scales::hue_pal())(n = length(x = levels(x = S)))
   names(x = ident.colors) <- levels(x = S)

   ### Cell colors
   cell.colors <- ident.colors[Idents(object = S)]
   names(x = cell.colors) <- colnames(x = S)

   emb = Embeddings(object = S, reduction = "umap")
   emat <- S\$spliced
   nmat <- S\$unspliced
##smat <- S4\$spanning #S4$ambiguous

   vel   = Tool(object = S, slot = "RunVelocity")
   #summary ( vel )

## Plot requires MEMORY !!!
# correlation-based transition probability matrix within the kNN graph
# n: number of neighbors

#this_scale=c("sqrt", "log", "rank", "linear")
    this_scale=c( "sqrt" )
    for (i in 1:length(this_scale)){
        show.velocity.on.embedding.cor( emb, vel, n=100, scale = this_scale[i],
				    cell.colors = ac(x = cell.colors,  alpha =cell.alpha),
				    cex = cell.cex, arrow.scale = arrow.scale,  arrow.lwd = arrow.lwd,
				    show.grid.flow = TRUE,
				    min.grid.cell.mass = 0.5,  grid.n = 40,
				    xlab="UMAP_1", ylab="UMAP_2",
				    do.par = FALSE,  expression.scaling = TRUE,
				    cell.border.alpha = 0.2) 

#   show.velocity.on.embedding.cor( emb, vel, n=100, scale = "sqrt", cell.colors = ac(x = cell.colors,  alpha =cell.alpha), cex = cell.cex, arrow.scale = arrow.scale,  arrow.lwd = arrow.lwd, show.grid.flow = TRUE, min.grid.cell.mass = 0.5,  grid.n = 40, xlab="UMAP_1", ylab="UMAP_2", do.par = FALSE,  expression.scaling = TRUE, cell.border.alpha = 0.2)
   }

}


#postscript("$ps_file", horizontal=T)
pdf(file="$pdf_file")
op <- par(cex=1.2)

DimPlot( S4, reduction = "umap", label=TRUE, pt.size=2,label.size=8) + NoLegend() + theme(aspect.ratio=1)

plot_velocity( S4, 4)
remove(S4)

plot_velocity( S5, 2)
remove(S5)

par(op)
dev.off()

### If you need to see each gene
#PITX2, FOXC1, TFAP2B
# Lmx1B, Emx2
##rvel <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,deltaT2 = 1,kCells = kCells)
#pdf(file="$resi_pdf_file",onefile=TRUE, width=11,height=3)
#op <- par(mfrow=c(1,4), cex=0.8)
#plot.gene <- c("Pitx2", "Foxc1", "Tfap2b", "Lmx1b")

#for (i in 1:length(plot.gene) ){
#   gene <- plot.gene[i]
#   if( sum(rownames(emat) == gene)){
#       gene.relative.velocity.estimates(emat,nmat,show.gene=gene,cell.emb=emb,cell.colors=cell.colors,
#				 kCells=kCells,  do.par=FALSE, fit.quantile=fit.quantile)
#  }
#}
#par(op)
#dev.off()

q()

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
