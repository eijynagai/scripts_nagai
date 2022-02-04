# Script to run GO analysis in tab delimited files using clusterProfiler


# How to run
# Rscript run_clusterProfiler_eijy.R --i DEGs --o GOterms

# TODO
# make better plot including p-value, number of significant
# reverse the GO term values
#



# Libraries
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("enrichplot"))
suppressPackageStartupMessages(library("ggplot2"))


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


# for testing
if(FALSE){
  setwd("Documents/_work/nakato_DEGs_GO_2021/")
  input_path <- "DEGs"
  output_path <- "GOterms"
}



# All DEGs list
allDEGs <- read.csv(file = paste(input_path,"AllDEGs.txt", sep = "/"), 
                    header = FALSE)$V1
#head(allDEGs)
#length(allDEGs)
#class(allDEGs) #must be character




# Input gene list
files <- list.files(path=input_path, pattern = "\\.tsv$")
#file <- files[7]
#file

for (file in files){
  # Initiate loop to all files in the directory

  # Create the pdf file
  fname <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file))

  # deleting previous outputfile
  #Check its existence
  if (file.exists(paste0(output_path, fname,"_clusterProfiler.pdf")))
    #Delete file if it exists
    file.remove(paste0(output_path, fname,"_clusterProfiler.pdf"))

  # Modify the size of page if necessary
  #pdf(paste0(output_path, "/", fname,"_clusterProfiler.pdf"), width=10, height=7) # include the path with name

  cat("Preparing the file: ", fname)
  cat("\n")
  cat("It may take a while...\n")

  # Read in each file separated
  data <-read.csv(file = paste(input_path, file, sep="/"), sep = "\t", header = TRUE) 
  deg.list <- data[,1]  #using symbol, set 2 for ensemblID

  deg.list

  ego <- enrichGO(gene        = deg.list,
                universe      = allDEGs,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH", # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
  dim(ego)[1] == 0
  #tail(ego)
  #length(ego)

  # Test if no enrichment term
  if (dim(ego)[1] == 0){
    cat("O-oh! No enrichment for this cluster. Skipping...\n\n")
    next
  }

  # Modify the size of page if necessary
  #pdf(paste0(output_path, "/", fname,"_clusterProfiler.pdf"), width=15, height=7) # include the path with name

  # Bar plot
  p1 <- barplot(ego,
                title = "GO term enrichment by clusterProfiler",
                showCategory=30)
  #print(p1)
  #p1
  ggsave(filename = paste0(output_path, "/", fname,"_clusterProfiler.pdf"), plot = p1, device = NULL)
  
  
  # Dot plot
  #p2 <- dotplot(ego, showCategory=30) + ggtitle("dotplot")
  #print(p2)

  # Gene-concept Netword
  #p3 <- cnetplot(ego, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
  #print(p3)

  # Heatmap-like functional classification
  #p4 <- heatplot(ego)
  #print(p4)

  # UpSet Plot
  #p5 <- upsetplot(ego)
  #print(p5)

  #dev.off()
  cat("Done.\n\n")

}
