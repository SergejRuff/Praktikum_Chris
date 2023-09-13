# clean global enviroment
rm(list=ls())

set.seed(700) # change before running code.

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


##########
# Import
###########

fasta_files <- "A:/Praktikum_Chris/data/nido_roniviruses_n6_2908/nido_invertebrate_refseq.fasta"
current_date <- Sys.Date()
current_time <- format(Sys.time(), format = "%H-%M-%S")
current_date <- paste0(current_date,"_",current_time)

range = 100

basedirectory <- "A:/Praktikum_Chris/output" # where should the output be placed. ! Important.
# All other folders will be generated in basedirectory
# Import from local files

# Path to your Bash script (including the script file name)
bash_script_path <- "/mnt/a/Praktikum_Chris/visualise_ct_files_naview.sh"

sequences <- read.fasta(file=fasta_files,as.string=TRUE)

all_data <- list()

robinsonfould_run_results <- matrix(NA, nrow=range, ncol=3)
colnames(robinsonfould_run_results) <- c("run","Rob_5UTR","Rob_3UTR")

##################
# Functions ######
##################



allignment_gl <- function(UTR_list,start=1){

  n <- length(UTR_list)

  alignment_scores_global_UTR <- matrix(NA, n, n)

  alignment_results_global_UTR <- list()

  lev_distance_matrix_UTR_global <- matrix(0, n, n)
  lev_distance_percentage_matrix_UTR_global <- matrix(0, n, n)



  # Create a vector of sequence names
  sequence_names <- names(UTR_list)

  # Assign row and column names
  rownames(lev_distance_matrix_UTR_global) <- sequence_names
  colnames(lev_distance_matrix_UTR_global) <- sequence_names

  rownames(lev_distance_percentage_matrix_UTR_global) <- sequence_names
  colnames(lev_distance_percentage_matrix_UTR_global) <- sequence_names





  # Loop through results_list for 5_UTR
  for (i in start:length(UTR_list)) {


    for (j in start:length(UTR_list)) {

      sequence_1 <- UTR_list[[i]][2,]
      sequence_2 <- UTR_list[[j]][2,]

      #print(sequence_1)
      #print(sequence_2)


      alignment_g <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2, type = "global")

      #print(paste("score for global",score(alignment_g)))
      #print(paste("score for global",score(alignment)))
      alignment_results_global_UTR[[paste("alignment_", names(UTR_list[i]), "_", names(UTR_list[j]), sep = "")]] <- alignment_g


      alignment_scores_global_UTR[i, j] <- score(alignment_g)

      # Example dot-bracket notations from UTR_list
      subject_global <- toString(subject(alignment_g ))
      pattern_global <- toString(pattern(alignment_g ))

      #print(paste("subject:",subject_global))
      #print(paste("pattern:",pattern_global))

      # Calculate the Levenshtein distance for global
      lev_distance_global <- stringdist::stringdist(subject_global, pattern_global)



      # Length of the longer sequence
      max_length_g <- max(nchar(subject_global), nchar(pattern_global))


      # Calculate the Levenshtein distance as a percentage
      lev_distance_percentage_g <- (lev_distance_global / max_length_g) * 100

      # Store results in matrices global
      lev_distance_matrix_UTR_global[i, j] <-  lev_distance_global
      lev_distance_percentage_matrix_UTR_global[i, j] <- lev_distance_percentage_g




    }

  }

  # Add row and column names

  rownames(alignment_scores_global_UTR) <- names(UTR_list)
  colnames(alignment_scores_global_UTR) <- names(UTR_list)

  # replace NaNs by 0 if error occurs, where 0 are replaced by NaN
  lev_distance_percentage_matrix_UTR_global[is.na(lev_distance_percentage_matrix_UTR_global)]<-0

  return(list(alignment_scores_global_UTR=alignment_scores_global_UTR,
              alignment_results_global_UTR=alignment_results_global_UTR,
              lev_distance_matrix_UTR_global=lev_distance_matrix_UTR_global,
              lev_distance_percentage_matrix_UTR_global=lev_distance_percentage_matrix_UTR_global))

}


utr_nj <- function(matrix,filename,test){


  distance_matrix <- as.dist(matrix)

  njtree <- nj(distance_matrix)

  #set midpoint.root
  njtree <- midpoint.root(njtree)


  #plot results.
  png(filename, width = 800, height = 600)

  # Create a larger plot area
  par(mfrow = c(1, 1))
  plot(njtree,main = paste("njplot for",test), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")



  dev.off()

  return(njtree)


}


tanglegram_plot <- function(x,y,all_x,all_y,filename){

  # Check if x and y are valid phylogenetic trees
  if (!is(x, "phylo") || !is(y, "phylo")) {
    stop("Input x and y must be valid phylogenetic trees.")
  }

  # Check if all_x and all_y are character vectors
  if (!is.character(all_x) || !is.character(all_y)) {
    stop("all_x and all_y must be character vectors.")
  }

  pattern <- "NC_[[:alnum:]_]{6}"

  tiplabel1<- x$tip.label
  tiplabel2<- y$tip.label




  association_matrix <- matrix(0,nrow=length(tiplabel1),ncol=2)

  # Iterate through tiplabel1 and find matches in tiplabel2
  for (i in 1:length(tiplabel1)) {
    # Check if tiplabel1[i] matches the pattern
    if (grepl(pattern, tiplabel1[i])) {
      # Extract the NC_number from tiplabel1
      nc_number1 <- sub(".*NC_([[:alnum:]_]{6}).*", "\\1", tiplabel1[i])

      # Find matching elements in tiplabel2 based on the NC_number
      matching_indices <- which(grepl(paste0("NC_", nc_number1), tiplabel2))

      # Check if there are matching elements in tiplabel2
      if (length(matching_indices) > 0) {
        matching_tip2 <- tiplabel2[matching_indices]
        association_matrix[i, 1] <- tiplabel1[i]
        association_matrix[i, 2] <- matching_tip2
      }
    }
  }


  # Print the association_matrix
  #print(association_matrix)


  tangleplot<- cophylo(x,y,assoc=association_matrix)

  #plot results.
  png(filename, width = 800, height = 600)

  # Create a larger plot area
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot(tangleplot)
  title(main = paste("tangleplot for",all_x,"and",all_y))


  dev.off()
}

# function to calculate the Robinson-Foulds (RF) distance between two trees
calculate_rf_distance <- function(tree1, tree2) {
  if (!is(tree1, "phylo") || !is(tree2, "phylo")) {
    stop("Input tree1 and tree2 must be valid phylogenetic trees.")
  }
  return(RobinsonFoulds(tree1, tree2))
}

#######################
# creat directory #####
#######################

folder_name <- paste0("random_files_",current_date)
folder_path1 <- file.path(basedirectory, folder_name)

if (!dir.exists(folder_path1)) {
  dir.create(folder_path1, recursive = TRUE)
}

###############################

for (k in 1:range){

  print(paste("run number:",k))

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

    three_UTR_ <- three_UTR
    five_UTR_ <- five_UTR


    # Convert the string to a character vector
    five_UTR <- unlist(strsplit(five_UTR, ""))

    # Shuffle the character vector
    five_UTR <- sample(five_UTR)

    # Convert the shuffled character vector back to a string
    five_UTR <- paste(five_UTR, collapse = "")

    # Convert the string to a character vector
    three_UTR <- unlist(strsplit(three_UTR, ""))

    # Shuffle the character vector
    three_UTR <- sample(three_UTR)

    # Convert the shuffled character vector back to a string
    three_UTR <- paste(three_UTR, collapse = "")


    print(paste("random 5_UTR and og 5UTR are identical?",identical(five_UTR,five_UTR_)))
    print(paste("random 3_UTR and og 3UTR are identical?",identical(three_UTR,three_UTR_)))




    ###################
    # print-statements
    ###################

    #cat(paste("for virus",NameofTaxon))
    #print(orfs) # print info about available ORFs.



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

    folder_name_r <- paste0("samplerun_",k)
    folder_path_ <- file.path(folder_path1, folder_name_r)

    if (!dir.exists(folder_path_)) {
      dir.create(folder_path_, recursive = TRUE)
    }





    folder_path2 <- file.path(folder_path_, "Alle_FASTAFiles") # check for existing folder and creat if not there.

    if (!dir.exists(folder_path2)) {
      dir.create(folder_path2, recursive = TRUE)
    }

    setwd(folder_path2)

    # save all fasta files in one folder called Alle_Fastafiles

    write.fasta(five_UTR,names=paste0("5_UTR_of_",NameofTaxon),file.out = paste0("5_URT_",NameofTaxon,current_date,".fasta"))

    write.fasta(three_UTR,names=paste0("3_UTR_of_",NameofTaxon),file.out = paste0("3_URT_",NameofTaxon,current_date,".fasta"))

    setwd("A:/Praktikum_Chris/R/code")




  }

  # Combine all data frames into one
  combined_data <- do.call(rbind, all_data)

  # Create a gt table using the combined data
  combined_table <- gt(combined_data) %>%
    tab_header(title = "UTR Positions and Lengths") %>%
    opt_align_table_header(align = "left")

  # Print the combined table
  #print(combined_table)

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



  fiveUTR_allignment <- allignment_gl(UTR_list =list_lt5_UTR,start=1)
  threeUTR_allignment <- allignment_gl(UTR_list =list_gt3_UTR,start=1) # change start for future !

  # alvinovirus is too short. it was ignored for allignment. We remove its rows and coloumns to remove NAs
  #threeUTR_allignment[["alignment_scores_global_UTR"]] <- threeUTR_allignment[["alignment_scores_global_UTR"]][-1, -1]                          # !!!!!!! remove for future virus-data.
  #threeUTR_allignment[["alignment_scores_local_UTR"]]<- threeUTR_allignment[["alignment_scores_local_UTR"]][-1, -1]

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
  write.csv(threeUTR_allignment[["lev_distance_matrix_UTR_global"]], file = paste0(folder_path_,"/allignmentscores/","threeUTR_global_lev_distance_matrix.csv"))

  # Write the matrix to a CSV file
  write.csv(threeUTR_allignment[["lev_distance_percentage_matrix_UTR_global"]], file = paste0(folder_path_,"/allignmentscores/","threeUTR_global_percentage_lev_distance_matrix.csv"))


  # Write the matrix to a CSV file
  write.csv(fiveUTR_allignment[["lev_distance_matrix_UTR_global"]], file = paste0(folder_path_,"/allignmentscores/","fiveUTR_global_lev_distance_matrix.csv"))

  # Write the matrix to a CSV file
  write.csv(fiveUTR_allignment[["lev_distance_percentage_matrix_UTR_global"]], file = paste0(folder_path_,"/allignmentscores/","fiveUTR_global_percentage_lev_distance_matrix.csv"))


  gc()  # Trigger garbage collection to release memory

  ###########################
  # Perfom Neighbor Joining #
  ###########################

  protein_sequences <- read.tree(
    file="A:/Praktikum_Chris/data/nido_roniviruses_n6_2908/nido_invertebrate_refseq_trimmed.phy_phyml_SH_tree_annotated_rooted.nwk")

  setwd(scores_path)



  # Call the function with allignment matrix
  tree_3_global <- utr_nj(threeUTR_allignment[["lev_distance_percentage_matrix_UTR_global"]],filename="threeUTR_njtree_global.png",test="Global Alignment of 3´-UTR")

  tree_5_global <- utr_nj(fiveUTR_allignment[["lev_distance_percentage_matrix_UTR_global"]],filename="fiveUTR_njtree_global.png",test="Global Alignment of 5´-UTR")


  # call functions to creat tangleplots.
  tanglegram_plot(tree_5_global,protein_sequences ,all_x = "5´-UTR Global",all_y = "Proteinsequences",filename="tangle_5UTR-Global-Protein.png")
  tanglegram_plot(tree_3_global,protein_sequences,all_x = "3´-UTR Global",all_y = "Proteinsequences",filename="tangle_3UTR-Global-Protein.png")




  # Update tip labels
  tree_5_global$tip.label <- sub(".*NC_([[:alnum:]_]{6}).*", "NC_\\1", tree_5_global$tip.label)

  # Update tip labels
  tree_3_global$tip.label <- sub(".*NC_([[:alnum:]_]{6}).*", "NC_\\1", tree_3_global$tip.label)

  # Update tip labels
  protein_sequences$tip.label <- sub(".*NC_([[:alnum:]_]{6}).*", "NC_\\1", protein_sequences$tip.label)




  # Calculate RF distances for your trees
  rf_distance_5_global_protein <- calculate_rf_distance(tree_5_global, protein_sequences)
  rf_distance_3_global_protein <- calculate_rf_distance(tree_3_global, protein_sequences)

  # Print RF distances
  #cat("RF Distance between 5´-UTR Global and Protein sequence:", rf_distance_5_global_protein , "\n")
  #cat("RF Distance between 3´-UTR Global and Protein sequence:", rf_distance_3_global_protein , "\n")


  distance_df <- data.frame(
    rf_distance_5_global_protein,
    rf_distance_3_global_protein
  )




  # Create a gt table
  table_Robinson <- gt(distance_df)

  # Change column labels
  table_Robinson <- table_Robinson %>%
    cols_label(
      rf_distance_5_global_protein = "5´-UTR and Protein",
      rf_distance_3_global_protein = "3´-UTR and Protein"
    )

  # Add table header
  table_Robinson <- table_Robinson %>%
    tab_header(
      title = "Distance between UTRs and Protein",
      subtitle = "Comparison of 5´UTR and 3´UTR with Protein"
    )%>%
    tab_source_note(
      source_note = "value represents the number of bipartitions (splits) that differ between two phylogenetic trees.A smaller RF distance implies a higher degree of similarity."

    )

  #print(table_Robinson)

  # Save the modified table with the title and column labels to a PDF
  gtsave(table_Robinson, filename = "table_Robinson.pdf")

  # Set path back to the code directory
  setwd("A:/Praktikum_Chris/R/code")

  robinsonfould_run_results[k,1]<-k
  robinsonfould_run_results[k,2]<-rf_distance_5_global_protein
  robinsonfould_run_results[k,3]<-rf_distance_3_global_protein


}

setwd(folder_path1)
write.csv(x=robinsonfould_run_results,"robinson_run_results.csv")
setwd("A:/Praktikum_Chris/R/code")
