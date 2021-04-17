# Script to run GO analysis in tab delimited files using topGO


# How to run
# Rscript run_topGO_eijy.v2.R --i DEGs --o GOterms

# TODO
# make better plot including p-value, number of significant
# reverse the GO term values
# 



# Libraries
suppressPackageStartupMessages(library("topGO"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("argparse"))



# Create parse parameters
parser <- ArgumentParser()
parser$add_argument('--indir', type='character', default='./',
                    help='Indicate the input directory.')
parser$add_argument('--outdir', type='character', default='./',
                    help='Indicate the output directory.')
args <- parser$parse_args()
input_path <- args$indir
output_path <- args$outdir

if(FALSE){
setwd("Documents/_work/nakato_DEGs_GO_2021/")
input_path <- "DEGs"
output_path <- "GOterms"
}

# Check if input and output exist.
if (!file.exists(input_path)) {
  cat("Input directory does not exists. Please check again.\n")
} else if (!file.exists(output_path)) {
  cat("Output directory does not exist, creating...")
  dir.create(file.path(output_path))
  cat(" done!\n")
}




#topGO function
# Please modify here the number nodes, statistics, and ontology
topGOterms = function( fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Human",
                       ontology.use = "BP",  #"BP", "CC", "MF"
                       stats.use = "fisher",  #fisher, ks, t, globaltest, sum
                       algorithm.use = "classic", #classic, elim, weight, weight01, lea, parentChild
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
    stop("Error : Organisms other than human, mouse, and worm are not supported currently")
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
                annot = annFUN.org, #annFUN.GO2genes, GO2genes=allGO2genes
                mapping = mapping.use,
                ID = "SYMBOL",
                nodeSize = 10)
  print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
  res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use, scoreOrder="increasing")
  to.return = list()
  to.return$GOdata = GOdata
  to.return$res.table <- GenTable(GOdata, pval = res.result, orderBy="pval", ranksOf = "pval", topNodes = topnodes.print, numChar = num.char)
  return(to.return)
}



#########################
# Background gene list

# remove iffalse if want to use whole gene annotation from biomaRt
if (FALSE){
# Downloading background gene list...it may take a while.
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#create vector of chromosomes
my_chr <- c(1:22,'X','Y')

#fetch data
all.genes = getBM(attributes = "hgnc_symbol",
                  filters = 'chromosome_name',
                  values = my_chr,
                  mart = ensembl)

# how many entries
dim(all.genes)

# convert from dataframe to list
bg.list <- as.list(as.data.frame(t(all.genes)))

# check if there are duplicated entries
cat("Is there any duplicated entry in background genes?\n")
table(duplicated(bg.list))
} #end of iffalse


## Creating the universe genes list based on input files
#system(paste0("Rscript run_create_background_gene_list.R --i ", input_path))


# All DEGs list
bg.list <- read.csv(file = paste(input_path,"AllDEGs.txt", sep = "/"), 
                    header = FALSE)$V1


#########################
# Input gene list
files <- list.files(path=input_path, pattern = "\\.tsv$")
#files
for (file in files){
  
  #file <- files[[1]]
  
  # Initiate loop to all files in the directory
  
  # Create the pdf file
  fname <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file))
  
  
  # deleting previous outputfile
  #Check its existence
  if (file.exists(paste0(output_path, fname,"_topGO.pdf")))
  #Delete file if it exists
    file.remove(paste0(output_path, fname,"_topGO.pdf"))
  
  # Modify the size of page if necessary
  pdf(paste0(output_path, "/", fname,"_topGO.pdf"), width=7, height=4) # include the path with name
  
  cat("Preparing the file: ", fname)
  cat("\n\n")
  
  # Read in each file separated
  data <-read.csv(file = paste(input_path, file, sep="/"), sep = "\t", header = TRUE) 
  deg.list <- data[,1]  #using symbol, set 2 for ensemblID
  
  #deg.list
  
  #########################
  # GO analysis
  
  # Create topGO object
  GOterms = topGOterms(fg.genes = deg.list, bg.genes = bg.list, organism = 'Human')
  #GOterms
  
  # Plot Go terms p-values
  bp_plot <- GOterms$res.table
  bp_plot$pval <- as.numeric(as.character(bp_plot$pval))
  #bp_plot
  #in case of NA values in pval
  bp_plot$pval[is.na(bp_plot$pval)] <- as.numeric(1e-30)
  
  # create the plot
  
  pp <- ggplot(bp_plot, 
          aes(x = Term, y = Significant, fill = -log10(as.numeric(pval)))) +
          geom_col() +
          ylab("Significant genes") +
          xlab("Biological process") +
          ggtitle("GO term enrichment by topGO") + 
          scale_x_discrete(limits=rev(bp_plot$Term)) +
          scale_fill_continuous(low = 'blue', high = 'red') +
          theme_bw(base_size=11) +  #size of font
      theme(
          legend.position='bottom') +
          labs(fill="-log10(P-val)") +
      coord_flip() 
  
  print(pp)  
  #pp
  
  dev.off()
  cat("Done.\n\n\n")

}



