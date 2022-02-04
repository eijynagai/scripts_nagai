# Function to create a file that contains all unique genes sorted from 
# a set of several files in a directory.

# Input: a directory containing all cluster genes
# Output: a file containing all a unique set of all genes, sorted.

suppressPackageStartupMessages(library("argparse"))

# Create parse parameters
parser <- ArgumentParser()
parser$add_argument('--indir', type='character', default='./',
                    help='Indicate the input directory.')
args <- parser$parse_args()
input_path <- args$indir


# Check if input and output exist.
if (!file.exists(input_path)) {
  stop("Input directory does not exists. Please check again.\n")
}

# Delete previous output file
#Check its existence
if (file.exists(paste0(input_path, "/AllDEGs.txt"))){
  #Delete file if it exists
  cat("There is a 'AllDEGs.txt' file, we will replace it with a updated one!\n")
  invisible(file.remove(paste0(input_path, "/AllDEGs.txt")))
}

# Create a empty list
finallist <- list()

# Loop all files in a directory
# Input gene list
files <- list.files(path=input_path, pattern = "\\.tsv$")
for (file in files){
  file
  data <-read.csv(file = paste(input_path, file, sep="/"), sep = "\t") 
  deg.list <- data[,1]  #using symbol, set 2 for ensemblID
  deg.list
# Store the gene column
  finallist <- append(finallist, list(deg.list))

}


files <- list.files(path=input_path, pattern = "\\.tsv$", full.names = TRUE)

final.df <- lapply(files, function(x){
  df <- read.csv(x, sep = "\t")
  df[,1]#using symbol, set 2 for ensemblID
})

#finallist<- as.list(as.data.frame(t(final.df))) #as list
finallist <- unlist(final.df) #as character

#class(finallist)
#length(finallist)

# Sort-unique values
#length(unique(sort(finallist)))
finallist <- unique(sort(finallist))
#finallist

# Save to a file
write.table(finallist, file = paste(input_path, "AllDEGs.txt", sep = "/"), 
            quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)

# Finished!
cat("\nDone!\n")
q()
