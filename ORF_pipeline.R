# clean global enviroment
rm(list=ls())

# packages
library(traits)
library(seqinr)
library(ORFik)


# Import from local files
sequences <- read.fasta(file="A:/Praktikum_Chris/data/Covid_ref/sequence.fasta",as.string=TRUE)

# extract specific fasta file from list
sequences <- sequences[[1]]

# extract the nucleotide sequence
seq <- sequences[1]

orfs<- findORFs(seq, minimumLength = 98)

orf_matrix <- as.matrix(orfs@unlistData)

orf_matrix <- cbind(orf_matrix,orf_matrix[,1]+orf_matrix[,2])

colnames(orf_matrix)<- c("start","length","end")
