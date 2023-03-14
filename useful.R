# List of useful packages and commands in R

# Load libraries and install if not installed
require2 <- function(x) { 
  if (!base::require(x, character.only = TRUE)) {
  install.packages(x, dep = TRUE) ; 
  base::require(x, character.only = TRUE)
  }
  base::library(x, character.only = TRUE)
}
require2("Seurat")
require2("tidyverse")
require2("patchwork")
require2("BiocManager")
BiocManager::install("dittoSeq", update = FALSE)



# str_c() concatenate strings
# dir() returns the file/directories contained in a path, similar to ls command
# str_c(dir()) returns the files and you can attach the path to them
# output is a list of paths to the respective files
str_c("path/to/dir", dir("path/do/dir"), sep = '/')


# str_replace to replace parts of string
# you can select a pattern and replace to another one
# example: remove the extension from files
str_replace(pattern = ".csv", replacement = "")


# map() will apply in <some_data>, <some_command>
# it is particular useful for list of lists or similar
map(so, plot())


# commandArgs() allows arguments to your Rscript
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
