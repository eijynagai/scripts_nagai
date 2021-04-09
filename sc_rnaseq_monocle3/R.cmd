#module use /usr/local/package/modulefiles
#module load R/3.6

library(monocle3)
library(ggplot2)
library(dplyr)

CR <- "/home/park/Hayashi/CellRanger/DrHayashi_No174_Sample1_cellranger/"
Genome <- "mm10"
UMI_Cutoff <- 0
PCx <- 100

## Cell Data Set
cds <- load_cellranger_data( pipestance_path = CR, genome = Genome,
     barcode_filtered = TRUE,
     umi_cutoff = UMI_Cutoff)

cds <- preprocess_cds(cds, method="PCA", num_dim = PCx, norm_method = "log")
#cds <- align_cds( cds, alignment_group = "batch")

cds <- reduce_dimension(cds,  preprocess_method="PCA", reduction_method="UMAP")
cds <- reduce_dimension(cds,  preprocess_method="PCA", reduction_method="tSNE")

cds <- cluster_cells(cds
       ,reduction_method="UMAP"
       ,cluster_method="louvain"   #louvain,  K and resoution are ignored
       ,num_iter = 10
       ,k=20  #bigger K results lower resolution
)    

cds<-learn_graph(cds)

postscript("plots.ps", horizontal=T)
op<-par(cex=1.2)

plot_cells(cds, reduction_method="UMAP", group_label_size=4, cell_size=1,show_trajectory_graph=FALSE)
plot_cells(cds, reduction_method="UMAP" , group_label_size=4,  cell_size=1,show_trajectory_graph=TRUE)

plot_cells(cds, reduction_method="tSNE" , group_label_size=4,  cell_size=1,show_trajectory_graph=FALSE)

par(op)
dev.off()
q()
