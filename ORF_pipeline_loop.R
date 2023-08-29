# clean global enviroment
rm(list=ls())

###########
# packages
###########

library(traits)
library(seqinr)
library(ORFik)
library(gt)
library(webshot2)

##########
# Import
###########

fasta_files <- "A:/Praktikum_Chris/data/nido_roniviruses_n6_2908/nido_roniviruses_n6.fasta"
current_date <- Sys.Date()
# Import from local files

sequences <- read.fasta(file=fasta_files,as.string=TRUE)


for (i in 1:length(sequences)){

  ##############
  # preprocessing
  ###############

  # extract specific fasta file from list
  sequence_ <- sequences[[i]]

  # extract the nucleotide sequence
  seq <- sequence_[1]

  seq<- toupper(seq)

  # get name of taxa
  NameofTaxon <- getName(sequence_)

  # find the ORFs in + direction with length of 300 or more. Using ORFik-package
  orfs<- findORFs(seq, minimumLength = 98,startCodon = "ATG")      # !!!!!!!! might need to change length later.
                                                                   # !!!! Startcodon changes results. Default uses alternative codons as well.


  # Convert the Data to a data frame
  orf_df <- as.data.frame(orfs@unlistData)

  # extract last position of 5´-UTR
  firstposition <- min(orf_df$start)-1

  # extract first position of 3´-UTR
  lastposition <- max(orf_df$end)+1


  # extract 5´- and 3´-UTR
  five_UTR <- substring(seq, first = 1, last=firstposition)
  three_UTR <- substring(seq, first =lastposition )

  #RNAstructure requires nucleotides to be upper case
  five_UTR <- toupper(five_UTR)

  three_UTR <- toupper(three_UTR)

  ###################
  # print-statements
  ###################

  cat(paste("for virus",NameofTaxon))
  print(orfs)

  # Create a data frame with the information
  data <- data.frame(
    "UTRs" = c("5´-UTR Position", "3´-UTR Position"),
    "position" = c(paste("1-",firstposition), paste(lastposition,"-",nchar(seq))),
    "length"= c(firstposition,nchar(seq)-lastposition+1)
  )

  # Create a gt table
  table <- gt(data)%>%
    tab_header(title ="UTR Positions and Lengths ",
               subtitle=md(paste("**virus**: ",NameofTaxon)))%>%
    opt_align_table_header(align="left")

  print(table)

  #########
  # export
  #########

  basedirectory <- "A:/Praktikum_Chris/output"
  folder_name <- paste0("fasta_files_",Sys.Date())
  folder_path <- file.path(basedirectory, folder_name)

  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }

  folder_path <- file.path(folder_path, NameofTaxon)

  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }


  setwd(folder_path)

  write.fasta(five_UTR,names=paste("5_UTR of",NameofTaxon),file.out = paste0("5_URT_",NameofTaxon,current_date,".fasta"))

  write.fasta(three_UTR,names=paste("3_UTR of",NameofTaxon),file.out = paste0("3_URT_",NameofTaxon,current_date,".fasta"))

  gtsave(data=table,filename = paste("table_",NameofTaxon,".pdf"))

  setwd("A:/Praktikum_Chris/R/code")




}
