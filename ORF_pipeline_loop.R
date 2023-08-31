# clean global enviroment
rm(list=ls())

###########
# packages
###########


library(seqinr) # package for read.fasta
library(ORFik) # package for ORF prediction
library(gt)
library(LncFinder) # runRNA_fold
library(ncRNAtools) # package for
library(RRNA) # good import ct-function
library(RNAsmc)

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
  orfs<- findORFs(seq, minimumLength = 98,startCodon = "ATG")      # !!!!!!!! might need to change length later. 248 for 750
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

  #cat(paste("for virus",NameofTaxon))
  #print(orfs) # print info about available ORFs.

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

  print(table) # print info about chosen ORFs

  #########
  # export
  #########

  # exprort for each virus in individual folders and than collect each fasta file in one folder.
  # Meaning fasta files are exported twice.

  basedirectory <- "A:/Praktikum_Chris/output"
  folder_name <- paste0("fasta_files_",Sys.Date())
  folder_path_ <- file.path(basedirectory, folder_name)

  if (!dir.exists(folder_path_)) {
    dir.create(folder_path_, recursive = TRUE)
  }

  folder_path <- file.path(folder_path_, NameofTaxon)

  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }


  setwd(folder_path) # save all fasta files and position-info in each folder per virus.

  write.fasta(five_UTR,names=paste("5_UTR of",NameofTaxon),file.out = paste0("5_URT_",NameofTaxon,current_date,".fasta"))

  write.fasta(three_UTR,names=paste("3_UTR of",NameofTaxon),file.out = paste0("3_URT_",NameofTaxon,current_date,".fasta"))

  gtsave(data=table,filename = paste("table_",NameofTaxon,".pdf"))



  folder_path2 <- file.path(folder_path_, "Alle_FASTAFiles") # check for existing folder and creat if not there.

  if (!dir.exists(folder_path2)) {
    dir.create(folder_path2, recursive = TRUE)
  }

  setwd(folder_path2)

  # save all fasta files in one folder called Alle_Fastafiles

  write.fasta(five_UTR,names=paste("5_UTR of",NameofTaxon),file.out = paste0("5_URT_",NameofTaxon,current_date,".fasta"))

  write.fasta(three_UTR,names=paste("3_UTR of",NameofTaxon),file.out = paste0("3_URT_",NameofTaxon,current_date,".fasta"))

  setwd("A:/Praktikum_Chris/R/code")




}

##############################################################################
########### Run RNA Fold #####################################################
##############################################################################

fasta_files <- list.files(path=folder_path2,pattern="*.fasta",full.names = TRUE)

results_list <- list()  # Create an empty list to store the results

# perform RNA-fold calculations for all fasta-files in a directory.
for (i in 1:length(fasta_files)){


  sequences2 <- read.fasta(file=fasta_files[i],as.string=TRUE) # import 5 and 3UTR files in directory

  name_virus <- getAnnot(sequences2) # get virus_name

  print(paste("performing RNAfold for virus:",name_virus))

  results <- run_RNAfold(sequences2,RNAfold.path = "A:/Praktikum_Chris/RNAfold/RNAfold.exe",parallel.cores = -1 ) # perform RNAfold using LncFinder-package

  results_list[[i]] <- results # save results into a list. dataframes containing sequence, dot bracket notation and energy.

  names(results_list)[i]<- name_virus
}


# Export as CT-file

# creat a new folder called CT_Fold_Ergebnisse and save all CT_Files in that folder.
for (i in 1:length(results_list)){

  # creat a CT-results folder if not already created
  CT_folder <- "CT_Fold_Ergebnisse"
  cT_path <- file.path(folder_path_, CT_folder)

  if (!dir.exists(cT_path)) {
    dir.create(cT_path, recursive = TRUE)
  }

  # extract name of virus from results_list
  virus_nam <-names(results_list)[i]

  # remove > otherwise file can´t be created
  virus_nam <- sub("^>", "", virus_nam )
  CT_name <- paste0(cT_path,"/",virus_nam,".ct")

  # save files as Connectivity table format using ncRNAtools´ writeCT-function
  writeCT(filename=CT_name,sequence=results_list[[i]][1,],secondaryStructure = results_list[[i]][2,],sequenceName = virus_nam)

}


#####################################
### plot secondary structure ########
#####################################


# plotting with VARNA and a bash-script.
# here we are importing the


# create dendogram for each UTR

ctUTR_files <- list.files(path=cT_path,pattern="*.ct",full.names = TRUE)

fiveUTR_list <- list()

threeUTR_list <- list()


# import only 3UTR and save them in list
for (file in ctUTR_files){
  if (grepl("^3_UTR", basename(file))){
    three_ct <- loadCt(file)
    name_v1 <- sub(".ct", "", basename(file) )
    threeUTR_list[[name_v1]]<- three_ct
  }

}

# import only 5UTR and save them in list
for (file in ctUTR_files){
  if (grepl("^5_UTR", basename(file))){
    five_ct <- loadCt(file)
    name_v2 <- sub(".ct", "", basename(file) )
    fiveUTR_list[[name_v2]]<- five_ct
  }

}


RNAstrCluster(threeUTR_list)



RNAstrCluster(fiveUTR_list)

for (i in 1:length(fiveUTR_list)){
  ct_circle <- RNAcirPlot (fiveUTR_list[[i]],cex=0.6)
  print(ct_circle)
}


#for (i in 1:length(threeUTR_list)){
#  ct_circle3 <- RNAcirPlot (threeUTR_list[[i]],cex=0.6)
#  print(ct_circle3)
#}


