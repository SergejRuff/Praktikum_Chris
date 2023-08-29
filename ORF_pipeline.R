# clean global enviroment
rm(list=ls())

# packages
library(traits)
library(seqinr)
library(ORFik)
library(dplyr)

# Import from local files
sequences <- read.fasta(file="A:/Praktikum_Chris/data/Covid_ref/sequence.fasta",as.string=TRUE)

# extract specific fasta file from list
sequences <- sequences[[1]]

# extract the nucleotide sequence
seq <- sequences[1]

# find the ORFs in + direction with length of 300 or more. Using ORFik-package
orfs<- findORFs(seq, minimumLength = 98)

print(orfs)

# Convert the Data to a data frame
orf_df <- as.data.frame(orfs@unlistData)

# extract last position of 5´-UTR
firstposition <- min(orf_df$start)-1

# extract first position of 3´-UTR
lastposition <- max(orf_df$end)+1

cat(paste("5´-UTR goes from position 1-", firstposition,
            "\n3´-UTR starts at position:", lastposition,
            "and ends with the last nucleotide"))


# extract 5´- and 3´-UTR
five_UTR <- substring(seq, first = 1, last=firstposition)
three_UTR <- substring(seq, first =lastposition )

#RNAstructure requires nucleotides to be upper case
five_UTR <- toupper(five_UTR)

three_UTR <- toupper(three_UTR)


