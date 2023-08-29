# clean global enviroment
rm(list=ls())

# packages
library(traits)
library(seqinr)

setwd("A:/Praktikum_Chris/R/preperation")
current_date <- Sys.Date()



# import of files
sequences <- ncbi_byid(ids="NC_045512.2",verbose=TRUE)
seq <- sequences$sequence # extract sequence only

# Import from local files
#sequences <- read.fasta(file="A:/Praktikum_Chris/R/preperation/sequence.fasta",as.string=TRUE)
#sequences <- sequences[[1]]
#seq <- sequences[1]

#NameofTaxon <- getAnnot(sequences)
# returns as a list - can import multiple fasta



# data preprocessing

firstposition <-265
lastposition <-29534


five_UTR <- substring(seq, first = 1, last=firstposition)
three_UTR <- substring(seq, first =lastposition )

#RNAstructure requires nucleotides to be upper case
five_UTR <- toupper(five_UTR)

three_UTR <- toupper(three_UTR)

# analysis

# Export



basedirectory <- getwd()
folder_name <- paste0("fasta_files_",Sys.Date())
folder_path <- file.path(basedirectory, folder_name)

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}

setwd(folder_path)

write.fasta(five_UTR,names=paste("5_UTR of",sequences$taxon),file.out = paste("5_URT",current_date,".fasta"))

write.fasta(three_UTR,names=paste("3_UTR of",sequences$taxon),file.out = paste("3_URT",current_date,".fasta"))
