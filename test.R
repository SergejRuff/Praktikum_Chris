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
library(gridExtra)
library(cowplot)
library(magick)
library(patchwork)
library(ggplot2)
library(Biostrings)
library(ape)
library(dendextend)
library(phytools)
library(phangorn)
library(TreeDist)
library(stringr)


fasta_files <- list.files(path="A:/Praktikum_Chris/Alle_FASTAFiles",pattern="*.fasta",full.names = TRUE)

cT_path <- "A:/Praktikum_Chris/Alle_FASTAFiles/Neuer Ordner"


results_list <- list()  # Create an empty list to store the results

# perform RNA-fold calculations for all fasta-files in a directory.
for (i in 1:length(fasta_files)){


  sequences2 <- read.fasta(file=fasta_files[i],as.string=TRUE) # import 5 and 3UTR files in directory

  name_virus <- getAnnot(sequences2) # get virus_name

  print(paste("performing RNAfold for virus:",name_virus))

  results <- run_RNAfold(sequences2,RNAfold.path = "A:/Praktikum_Chris/viennarna/RNAfold.exe",parallel.cores = -1 ) # perform RNAfold using LncFinder-package

  results_list[[i]] <- results # save results into a list. dataframes containing sequence, dot bracket notation and energy.

  names(results_list)[i]<- name_virus
}


# Export as CT-file

# creat a new folder called CT_Fold_Ergebnisse and save all CT_Files in that folder.
for (i in 1:length(results_list)){





  # extract name of virus from results_list
  virus_nam <-names(results_list)[i]

  # remove > otherwise file can´t be created
  virus_nam <- sub("^>", "", virus_nam )
  CT_name <- paste0(cT_path,"/",virus_nam,".dot")

  # save files as Connectivity table format using ncRNAtools´ writeCT-function

  writeDotBracket (filename=CT_name,sequence=results_list[[i]][1,],secondaryStructure = results_list[[i]][2,],sequenceName = virus_nam)

}

collect_dot_files <- function(folder_path) {
  # List all files in the folder
  all_files <- list.files(path = folder_path, pattern = ".*", full.names = TRUE)

  # Filter only dot-files
  dot_files <- all_files[file.info(all_files)$isdir == FALSE]

  # Initialize an empty character vector to store content
  file_contents <- character(0)

  # Loop through dot-files and collect their content
  for (file in dot_files) {
    content <- readLines(file, warn = FALSE)
    file_contents <- c(file_contents, content)
  }

  # Create a combined text file
  combined_file <- file.path(folder_path, "combined.txt")
  writeLines(file_contents, combined_file)

  cat("Dot-files content has been collected and saved to 'combined.txt' in the same folder.\n")
}

# Usage example:
collect_dot_files(cT_path )
