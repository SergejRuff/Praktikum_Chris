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
library(gridExtra)
library(cowplot)
library(magick)
library(patchwork)
library(ggplot2)
library(Biostrings)


##########
# Import
###########

fasta_files <- "A:/Praktikum_Chris/data/nido_roniviruses_n6_2908/nido_roniviruses_n6.fasta"
current_date <- Sys.Date()

basedirectory <- "A:/Praktikum_Chris/output" # where should the output be placed. ! Important.
# All other folders will be generated in basedirectory
# Import from local files

sequences <- read.fasta(file=fasta_files,as.string=TRUE)

all_data <- list()


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

  # Create a data frame with the information for this run
  data <- data.frame(
    "Virus" = NameofTaxon,
    "UTRs" = c("5´-UTR Position", "3´-UTR Position"),
    "position" = c(paste("1-",firstposition), paste(lastposition,"-",nchar(seq))),
    "length"= c(firstposition, nchar(seq) - lastposition + 1)
  )

  # Add the data frame to the list
  all_data[[i]] <- data

  #print(table) # print info about chosen ORFs

  #########
  # export
  #########

  # exprort for each virus in individual folders and than collect each fasta file in one folder.
  # Meaning fasta files are exported twice.

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

# Combine all data frames into one
combined_data <- do.call(rbind, all_data)

# Create a gt table using the combined data
combined_table <- gt(combined_data) %>%
  tab_header(title = "UTR Positions and Lengths") %>%
  opt_align_table_header(align = "left")

# Print the combined table
print(combined_table)

gtsave(data=combined_table,filename ="table_Gesamt.pdf",path=folder_path_)

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

  results <- run_RNAfold(sequences2,RNAfold.path = "A:/Praktikum_Chris/viennarna/RNAfold.exe",parallel.cores = -1 ) # perform RNAfold using LncFinder-package

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

# Excecute bash script from R to plot fold-CTs

# run command in Ubuntu to convert windows syntax to linux.
#otherwise code doesnt run.
#sudo apt-get install dos2unix
#dos2unix visualise_ct_files_naview.sh
## bash files dont accept windows_path ways. Need to convert to linux_path.

# Replace only the uppercase drive letter 'A' with 'a'
wsl_path <- gsub("^A:", "a", cT_path )

# Replace the Windows drive letter with the WSL path
wsl_path <- paste0("/mnt/",wsl_path)

# Replace only the uppercase drive letter 'A' with 'a'
plot_path <- gsub("^A:", "a", folder_path_ )

# Replace the Windows drive letter with the WSL path
plot_path <- paste0("/mnt/",plot_path)

new_PLOTS_DIR <- paste0(plot_path,"/plots/")



# Path to your Bash script (including the script file name)
bash_script_path <- "/mnt/a/Praktikum_Chris/visualise_ct_files_naview.sh"

# Build the command with arguments
command <- sprintf("bash \"%s\" \"%s\" \"%s\"", shQuote(bash_script_path), shQuote(wsl_path), shQuote(new_PLOTS_DIR))

# Run the Bash script with the specified paths as arguments using Git Bash
system(command)


####################################################
# combine all pngs generated by VARNA into one pdf #
####################################################


# Path to the directory containing plots
plots_directory <- paste0(folder_path_,"/","plots")

# List of plot file names
plot_files <- list.files(plots_directory, pattern = "\\.png$", full.names = TRUE)

# Separate plot files into those containing "3_UTR" and "5_UTR"
utr3_plot_files <- plot_files[grep("3_UTR", plot_files)]
utr5_plot_files <- plot_files[grep("5_UTR", plot_files)]

# Function to import and create a plot for a list of files
import_and_plot <- function(files) {
  plots <- list()
  for (file in files) {
    plots[[basename(file)]] <- ggdraw() + draw_image(image_read(file))
  }
  return(plots)
}

# Import and plot "3_UTR" images
imported_plots_utr3 <- import_and_plot(utr3_plot_files)

# add title to plots
plots_list_utr3 <- lapply(names(imported_plots_utr3), function(title) {
  plot = imported_plots_utr3[[title]]
  title2 <- gsub("\\.png$", "", title)
  plot + ggtitle(title2) + theme(plot.title = element_text(size = 10, hjust = 0.5))
})

# Import and plot "5_UTR" images
imported_plots_utr5 <- import_and_plot(utr5_plot_files)

# add title to plots
plots_list_utr5 <- lapply(names(imported_plots_utr5), function(title) {
  plot = imported_plots_utr5[[title]]
  title2 <- gsub("\\.png$", "", title)
  plot + ggtitle(title2) + theme(plot.title = element_text(size = 10, hjust = 0.5))
})

# Combine plots using patchwork for "3_UTR" images
combined_plot_utr3 <- wrap_plots(plots = plots_list_utr3, ncol = 1)

# Combine plots using patchwork for "5_UTR" images
combined_plot_utr5 <- wrap_plots(plots = plots_list_utr5, ncol = 1)

# Display the combined plots
print(combined_plot_utr3)
print(combined_plot_utr5)

setwd(plots_directory)

# Export plots as PDF files
ggsave("combined_plot_utr3.pdf", combined_plot_utr3, width = 8, height = 10,dpi=320)
ggsave("combined_plot_utr5.pdf", combined_plot_utr5, width = 8, height = 10,dpi=320)

# set path back to code directory
setwd("A:/Praktikum_Chris/R/code")

##########################################################
###### Perform global and local Allignment ###############
##########################################################

# Initialize empty lists to store the separated objects
list_gt3_UTR <- list()
list_lt5_UTR <- list()



for(i in 1:length(results_list)){
  if (grepl("^>3_UTR", names(results_list)[i])) {
    list_gt3_UTR[[names(results_list)[i]]] <- results_list[[names(results_list)[i]]]
  } else if (grepl(">5_UTR", names(results_list)[i])) {
    list_lt5_UTR[[names(results_list)[i]]] <- results_list[[names(results_list)[i]]]
  }

}



allignment_gl <- function(UTR_list,start=1){

  n <- length(UTR_list)
  alignment_scores_local_UTR <- matrix(NA, n, n)
  alignment_scores_global_UTR <- matrix(NA, n, n)

  alignment_results_local_UTR <- list()
  alignment_results_global_UTR <- list()

  lev_distance_matrix_UTR <- matrix(0, n, n)
  lev_distance_percentage_matrix_UTR <- matrix(0, n, n)

  # Create a vector of sequence names
  sequence_names <- names(UTR_list)

  # Assign row and column names
  rownames(lev_distance_matrix_UTR) <- sequence_names
  colnames(lev_distance_matrix_UTR) <- sequence_names

  rownames(lev_distance_percentage_matrix_UTR) <- sequence_names
  colnames(lev_distance_percentage_matrix_UTR) <- sequence_names

  # Loop through UTR_list
  for (i in 1:n) {
    for (j in 1:n) {
      # Example dot-bracket notations from UTR_list
      sequence_1 <- UTR_list[[i]][2, ]
      sequence_2 <- UTR_list[[j]][2, ]

      # Calculate the Levenshtein distance
      lev_distance <- stringdist::stringdist(sequence_1, sequence_2)

      # Length of the longer sequence
      max_length <- max(nchar(sequence_1), nchar(sequence_2))

      # Calculate the Levenshtein distance as a percentage
      lev_distance_percentage <- (lev_distance / max_length) * 100

      # Store results in matrices
      lev_distance_matrix_UTR[i, j] <- lev_distance
      lev_distance_percentage_matrix_UTR[i, j] <- lev_distance_percentage
    }
  }


  # Loop through results_list for 5_UTR
  for (i in start:length(UTR_list)) {


    for (j in start:length(UTR_list)) {

      sequence_1 <- UTR_list[[i]][2,]
      sequence_2 <- UTR_list[[j]][2,]

      #print(sequence_1)
      #print(sequence_2)

      alignment <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2, type = "local")
      alignment_g <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2, type = "global")

      #print(paste("score for global",score(alignment_g)))
      #print(paste("score for global",score(alignment)))
      alignment_results_local_UTR[[paste("alignment_", names(UTR_list[i]), "_", names(UTR_list[j]), sep = "")]] <- alignment
      alignment_results_global_UTR[[paste("alignment_", names(UTR_list[i]), "_", names(UTR_list[j]), sep = "")]] <- alignment_g

      alignment_scores_local_UTR[i, j] <- score(alignment)
      alignment_scores_global_UTR[i, j] <- score(alignment_g)


    }
  }

  # Add row and column names
  rownames(alignment_scores_local_UTR) <- names(UTR_list)
  colnames(alignment_scores_local_UTR) <- names(UTR_list)

  rownames(alignment_scores_global_UTR) <- names(UTR_list)
  colnames(alignment_scores_global_UTR) <- names(UTR_list)

  return(list(alignment_scores_local_UTR=alignment_scores_local_UTR,alignment_scores_global_UTR=alignment_scores_global_UTR,
              alignment_results_local_UTR=alignment_results_local_UTR,alignment_results_global_UTR=alignment_results_global_UTR,
              lev_distance_matrix_UTR=lev_distance_matrix_UTR,lev_distance_percentage_matrix_UTR=lev_distance_percentage_matrix_UTR))

}


fiveUTR_allignment <- allignment_gl(UTR_list =list_lt5_UTR,start=1)
threeUTR_allignment <- allignment_gl(UTR_list =list_gt3_UTR,start=2) # change start for future !

# alvinovirus is too short. it was ignored for allignment. We remove its rows and coloumns to remove NAs
threeUTR_allignment[["alignment_scores_global_UTR"]] <- threeUTR_allignment[["alignment_scores_global_UTR"]][-1, -1]                          # !!!!!!! remove for future virus-data.
threeUTR_allignment[["alignment_scores_local_UTR"]]<- threeUTR_allignment[["alignment_scores_local_UTR"]][-1, -1]                            # !!!!!!! remove for future virus-data.

###################### ######################
# Export Allignment and Levenstein_matrices #
#############################################

# creat a score-results folder if not already created
score_folder <- "allignmentscores"
scores_path <- file.path(folder_path_, score_folder)

if (!dir.exists(scores_path)) {
  dir.create(scores_path, recursive = TRUE)
}


# Write the matrix to a CSV file
write.csv(threeUTR_allignment[["alignment_scores_local_UTR"]], file = paste0(folder_path_,"/allignmentscores/","alignment_scores_local_3UTR.csv"))

# Write the matrix to a CSV file
write.csv(threeUTR_allignment[["alignment_scores_global_UTR"]], file = paste0(folder_path_,"/allignmentscores/","alignment_scores_global_3UTR.csv"))

# Write the matrix to a CSV file
write.csv(fiveUTR_allignment[["alignment_scores_local_UTR"]], file = paste0(folder_path_,"/allignmentscores/","alignment_scores_local_5UTR.csv"))
# Write the matrix to a CSV file
write.csv(fiveUTR_allignment[["alignment_scores_global_UTR"]], file = paste0(folder_path_,"/allignmentscores/","alignment_scores_global_5UTR.csv"))

### Lvenstein

# Write the matrix to a CSV file
write.csv(fiveUTR_allignment[["lev_distance_matrix_UTR"]], file = paste0(folder_path_,"/allignmentscores/","lev_distance_matrix_5UTR.csv"))

# Write the matrix to a CSV file
write.csv(fiveUTR_allignment[["lev_distance_percentage_matrix_UTR"]], file = paste0(folder_path_,"/allignmentscores/","lev_distance_percentage_matrix_5UTR.csv"))

# Write the matrix to a CSV file
write.csv(threeUTR_allignment[["lev_distance_matrix_UTR"]], file = paste0(folder_path_,"/allignmentscores/","lev_distance_matrix_3UTR.csv"))
# Write the matrix to a CSV file
write.csv(threeUTR_allignment[["lev_distance_percentage_matrix_UTR"]], file = paste0(folder_path_,"/allignmentscores/","lev_distance_percentage_matrix_3UTR.csv"))




