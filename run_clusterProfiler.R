# Script to run GO analysis in tab delimited files using clusterProfiler


# How to run
# Rscript run_clusterProfiler_eijy.v1.R --i DEGs --o GOterms


# TODO




# Libraries
suppressPackageStartupMessages(library("clusterProfiler"))
library(org.Hs.eg.db)
suppressPackageStartupMessages(library("argparse"))
#library(enrichplot)
if(FALSE){
# Create parse parameters
parser <- ArgumentParser()
parser$add_argument('--indir', type='character', default='./',
                    help='Indicate the input directory.')
parser$add_argument('--outdir', type='character', default='./',
                    help='Indicate the output directory.')
args <- parser$parse_args()
input_path <- args$indir
output_path <- args$outdir


# Check if input and output exist.
if (!file.exists(input_path)) {
  cat("Input directory does not exists. Please check again.\n")
} else if (!file.exists(output_path)) {
  cat("Output directory does not exist, creating...")
  dir.create(file.path(output_path))
  cat(" done!\n")
}
}## end if


# All DEGs list
allDEGs <- read.csv(file = "AllDEGs.txt", header = FALSE)$V1
head(allDEGs)
length(allDEGs)
class(allDEGs) #must be character

# Imput gene list
setwd("Documents/_work/nakato_DEGs_GO_2021/")
data <-read.csv(file = "DEGs/AllDEGs.1e-4.Zscore.cluster0.tsv", sep = "\t")
head(data)
deg.list <- data[,1]   #using symbol, set 2 for ensemblID
#head(deg.list)
#length(deg.list)
#class(deg.list) #must be character

ego <- enrichGO(gene          = deg.list,
                universe      = allDEGs,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
#head(ego)

# Bar plot
barplot(ego, showCategory=20)

# Dot plot
p1 <- dotplot(ego, showCategory=30) + ggtitle("dotplot")
p1

# Gene-concept Netword
p2 <- cnetplot(ego, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
p2

# Heatmap-like functional classification
p3 <- heatplot(ego)
p3

# UpSet Plot
p4 <- upsetplot(ego)
p4


