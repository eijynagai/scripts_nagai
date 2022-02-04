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




# All DEGs list
allDEGs <- read.csv(file = paste(input_path,"AllDEGs.txt", sep = "/"), 
                    header = FALSE)$V1
#head(allDEGs)
#length(allDEGs)
#class(allDEGs) #must be character




# input files
files <- list.files(path = input_path, pattern = "\\.tsv$")
files


my_list <- list()

for (i in 1:length(files)){
    
    #pattern of the files
    fname <- sub(pattern = "updownDEGsimilarity.kmeans.updown.top1000.k20.", replacement = '', basename(files[[i]]))
    fname <- sub(pattern = ".tsv", replacement = '', basename(fname))
    clustern <- sub(pattern = "cluster", replacement = '', basename(fname))
    #print(fname)
    
    #read in each file
    d <- read.csv(file = paste(input_path, files[[i]], sep = '/'), header = TRUE, sep = '\t')
    #print(head(d))
    
    # convert and take only entrezID because cluster compare can only accept entrezID
    e <- bitr(d$X1, fromType = "SYMBOL",
            toType = c("ENSEMBL", "ENTREZID"),
            OrgDb = org.Hs.eg.db)
    #print(head(e))
    
    # get only entrezID
    my_list[[paste0("C", clustern)]] <- e$ENTREZID
    #my_list <- list(my_list, list(e$ENTREZID))
}
#my_list


lapply(my_list, head)

# GO term enrichment
cg <- compareCluster(geneCluster = my_list, fun = "enrichGO", 
                     OrgDb='org.Hs.eg.db',
                     pvalueCutoff=0.01
                    )
head(as.data.frame(cg))



#options(repr.plot.width = 16, repr.plot.height = 20)
fig1 <- dotplot(cg, showCategory=10, color="pvalue")
#fig1
ggsave(filename = paste0(output_path, "/clusterProfiler_compare_BP.pdf"),plot = fig1, device = NULL,  width=12, height=14)



# KEGG pathway enrichment
ck <- compareCluster(geneCluster = my_list, fun = "enrichKEGG",
                     pvalueCutoff=0.01)
head(as.data.frame(ck))


#options(repr.plot.width = 12, repr.plot.height = 15)
#options(repr.plot.width = 13, repr.plot.height = 15)
fig2 <- dotplot(ck, showCategory=10, color="pvalue")
#fig2
ggsave(filename = paste0(output_path, "/clusterProfiler_compare_pathways.pdf"),plot = fig2, device = NULL, width=11, height=14)
