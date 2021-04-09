setwd("/Users/nagai/Documents/vandenbon/")

#available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))
#installed.genomes(splitNameParts=FALSE)
library(BiocManager)
#install("BSgenome.Celegans.UCSC.ce10")   #WS220
#getBSgenome("BSgenome.Celegans.UCSC.ce11", masked=FALSE)
library(BSgenome.Celegans.UCSC.ce10)
library(GenomicRanges)

## Number of sequences in this genome:
length(Celegans) 

## Display a summary of the sequences:
Celegans

## Index of single sequences:
seqnames(Celegans)

## Lengths (i.e. number of nucleotides) of the single sequences:
seqlengths(Celegans)



# input table common gene name,gene name,chr,start position,strand,distance to the start codon,other gene names (or overlapped genes)
annot <- read.table(file = "TSS/TableS2_representativeTSS.input.csv", sep = ",", quote="", fill=FALSE )
annot2 <- read.table(file = "TSS/TableS3_representativeTSS.input.csv", sep = ",", quote="", fill=FALSE)

# check if excel data issue happen
table(annot$V1=='01-Apr') #bad
table(annot$V1=='apr-1') #good

annot[1:10,]
annot2[1:10,]
dim(annot)
dim(annot2)

# Merge the two files
merged <- rbind(annot,annot2)
dim(merged)

# how many overlpas between embryo and adults?
table(annot2$V1 %in% annot$V1)
table(annot$V1 %in% annot2$V1)
table(duplicated(merged$V1))

# example of overlpas
example <- head(annot[annot$V1 %in% annot2$V1,]$V1, 10)
annot[annot$V1%in%example,][,1:6]
annot2[annot2$V1%in%example,][,1:6]

# keeping the first entry
library(dplyr)
ordered <- merged %>%
            arrange(V1) %>% 
            group_by(V1) %>% 
            distinct()
table(duplicated(ordered$V1))

#keep only first entry of the duplicated different position values
dedup <- ordered[!duplicated(ordered$V1),]
dim(dedup)
table(duplicated(dedup$V1))


# make it as dataframe
class(dedup)
dedup <- as.data.frame(dedup)

# which is the highest value of chromosome2 ?
head(dedup)
max(dedup[dedup$V3=='chrII',]$V4) #15273954
min(dedup[dedup$V3=='chrII',]$V4) #15273954

dim(dedup)
dedup <- dedup[!((dedup$V3=='chrII') & (dedup$V4==878)),]
dim(dedup)

# converting to GRanges 
gr <- GRanges(seqnames = Rle(dedup[,3]), ranges =   #chromosome
                IRanges(start = dedup[,4]  #TSS start
                        , width = 1 #
                        , names = dedup[,2])  # use 2 or 1
                        , strand = Rle(dedup[,5]))
gr
length(gr)
size = promoters(gr, upstream = 1000, downstream = 1000)  # +-1kb


# solving the out of range chromosome problem
size
seqs <- getSeq(Celegans, size)
seqs


# finding overlaps regions (stoped here)
findOverlaps(gr,)






######################
##  TFBSTools
######################
BiocManager::install("TFBSTools")
#BiocManager::install("JASPAR2020") cannot install in my R version
BiocManager::install("JASPAR2018")

suppressMessages(library(JASPAR2018))

# Search a motif in a sequence
searchSeq() 

# plot a sequence logo 
seqLogo()

# Container for storing a set of putative TFBS on a nucletide sequence
SiteSet()

# Motifset is a container for storing the negerated motif from Motif identification tool such as MEME
MotifSet()






