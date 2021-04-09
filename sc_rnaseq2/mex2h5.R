library(Seurat)
library(Matrix)
library(DropletUtils)

# Convert MEX 10x files (matrix.mtx, genes.tsv, barcodes.tsv) to H5
wd <- '/home/nagailae/projects/sc_celegans/seurat3/'
outputdir <- paste0(wd, 'data/h5_raw/')
outputname <- 'celegans.h5'
inputdir <- paste0(wd, 'data/raw/')
genomeversion <- "unknown" #"WS260"

# Read in count matrix
CR <- Read10X(data.dir = inputdir)

# Transpose the matrix to be cell on row and gene on columns
RC <- t(CR)

# Save to H5 file
write10xCounts(path=paste0(outputdir,outputname),
    x=RC,
    barcodes = colnames(RC),
    gene.id = rownames(RC),
    gene.symbol = rownames(RC),
    gene.type = "Gene Expression",
    overwrite = TRUE,
    type = "HDF5",
    version="2",
    genome = genomeversion)

remove(RC)
remove(CR)
q()
