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
library(ape)
library(dendextend)
library(phytools)
library(phangorn)


##########
# Import
###########

fasta_files <- "A:/Praktikum_Chris/data/nido_roniviruses_n6_2908/nido_invertebrate_refseq.fasta"
current_date <- Sys.Date()
current_time <- format(Sys.time(), format = "%H-%M-%S")
current_date <- paste0(current_date,"_",current_time)

basedirectory <- "A:/Praktikum_Chris/output" # where should the output be placed. ! Important.
# All other folders will be generated in basedirectory
# Import from local files

# Path to your Bash script (including the script file name)
bash_script_path <- "/mnt/a/Praktikum_Chris/visualise_ct_files_naview.sh"

sequences <- read.fasta(file=fasta_files,as.string=TRUE)

# Define the sequences to remove
sequences_to_remove <- c("Empoasca_onukii_nidovirus_2", "Parasteatoda_nidovirus_1")

# Remove the specified sequences
sequences <- sequences[!names(sequences) %in% sequences_to_remove]

all_data <- list()

has_poly_a_tail <- function(sequence) {
  poly_a_pattern <- "A{5,}$"
  grepl(poly_a_pattern, sequence)
}



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

  has_poly_a <- has_poly_a_tail(three_UTR)


  ###################
  # print-statements
  ###################

  #cat(paste("for virus",NameofTaxon))
  #print(orfs) # print info about available ORFs.

  # Create a data frame with the information
  data <- data.frame(
    "UTRs" = c("5´-UTR Position", "3´-UTR Position"),
    "position" = c(paste("1-",firstposition), paste(lastposition,"-",nchar(seq))),
    "length"= c(firstposition,nchar(seq)-lastposition+1),
    "Poly-A Tail" = has_poly_a  # Add the 'Poly-A Tail' column
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
    "length"= c(firstposition, nchar(seq) - lastposition + 1),
    "Poly-A Tail" = has_poly_a  # Add the 'Poly-A Tail' column
  )

  # Add the data frame to the list
  all_data[[i]] <- data

  #print(table) # print info about chosen ORFs

  #########
  # export
  #########

  # exprort for each virus in individual folders and than collect each fasta file in one folder.
  # Meaning fasta files are exported twice.

  folder_name <- paste0("fasta_files_",current_date)
  folder_path_ <- file.path(basedirectory, folder_name)

  if (!dir.exists(folder_path_)) {
    dir.create(folder_path_, recursive = TRUE)
  }

  folder_path <- file.path(folder_path_, NameofTaxon)



  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }

  setwd(folder_path) # save all fasta files and position-info in each folder per virus.

  write.fasta(five_UTR,names=paste0("5_UTR_of_",NameofTaxon),file.out = paste0("5_URT_",NameofTaxon,current_date,".fasta"))

  write.fasta(three_UTR,names=paste0("3_UTR_of_",NameofTaxon),file.out = paste0("3_URT_",NameofTaxon,current_date,".fasta"))

  gtsave(data=table,filename = paste("table_",NameofTaxon,".pdf"))




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



# Build the command with arguments
command <- sprintf("bash \"%s\" \"%s\" \"%s\"", shQuote(bash_script_path), shQuote(wsl_path), shQuote(new_PLOTS_DIR))

# Run the Bash script with the specified paths as arguments using Git Bash
system(command)


####################################################
# combine all pngs generated by VARNA into one pdf #
####################################################


# Path to the directory containing plots
plots_directory <- paste0(folder_path_,"/","plots")

if (!dir.exists(plots_directory)) {
  dir.create(plots_directory, recursive = TRUE)
}

setwd(plots_directory)


# Function to combine and save plots with unique filenames
combine_and_save_plots <- function(plot_files, prefix) {
  num_plots <- length(plot_files)
  num_per_combined <- 3  # Number of plots to combine in each image
  num_combined <- ceiling(num_plots / num_per_combined)

  timestamp <- Sys.Date()


  for (i in seq_len(num_combined)) {
    start_idx <- (i - 1) * num_per_combined + 1
    end_idx <- min(i * num_per_combined, num_plots)

    current_plot_files <- plot_files[start_idx:end_idx]
    imported_plots <- import_and_plot(current_plot_files)

    # Add titles to plots
    plots_list <- lapply(names(imported_plots), function(title) {
      plot <- imported_plots[[title]]
      title2 <- gsub("\\.png$", "", title)
      plot + ggtitle(title2) + theme(plot.title = element_text(size = 10, hjust = 0.5))
    })

    # Combine plots using patchwork
    combined_plot <- wrap_plots(plots = plots_list, ncol = 1)

    # Generate a unique output filename with a timestamp
    output_filename <- paste0(prefix, "_", timestamp, "_", i, ".pdf")

    # Export the combined plot as a PNG file
    ggsave(output_filename, combined_plot, width = 8, height = 10, dpi = 320)
    gc()  # Trigger garbage collection to release memory
  }
}

# Path to the directory containing plots
plots_directory <- paste0(folder_path_, "/", "plots")

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

# Combine and save "3_UTR" plots with unique filenames
combine_and_save_plots(utr3_plot_files, "combined_plot_utr3")
gc()  # Trigger garbage collection to release memory

# Combine and save "5_UTR" plots with unique filenames
combine_and_save_plots(utr5_plot_files, "combined_plot_utr5")
gc()  # Trigger garbage collection to release memory

rm("imported_plots_utr5","imported_plots_utr3","plots_list_utr5","plots_list_utr3")
gc()  # Trigger garbage collection to release memory

# Set path back to the code directory
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
  print(association_matrix)


  tangleplot<- cophylo(x,y,assoc=association_matrix)

  #plot results.
  png(filename, width = 1920, height = 1080)

  # Create a larger plot area
  par(par(mfrow=c(1,2), mar = c(5, 4, 4, 2) + 0.1))
  plot(tangleplot, use.edge.length = TRUE, cex = 0.8)
  title(main = paste("tangleplot for", all_x, "and", all_y), line = -1)




  dev.off()
}

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



# function to calculate the Robinson-Foulds (RF) distance between two trees
calculate_rf_distance <- function(tree1, tree2) {
  if (!is(tree1, "phylo") || !is(tree2, "phylo")) {
    stop("Input tree1 and tree2 must be valid phylogenetic trees.")
  }
  return(RobinsonFoulds(tree1, tree2))
}

# Calculate RF distances for your trees
rf_distance_5_global_protein <- calculate_rf_distance(tree_5_global, protein_sequences)
rf_distance_3_global_protein <- calculate_rf_distance(tree_3_global, protein_sequences)

# Print RF distances
cat("RF Distance between 5´-UTR Global and Protein sequence:", rf_distance_5_global_protein , "\n")
cat("RF Distance between 3´-UTR Global and Protein sequence:", rf_distance_3_global_protein , "\n")


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

print(table_Robinson)

# Save the modified table with the title and column labels to a PDF
gtsave(table_Robinson, filename = "table_Robinson.pdf")

# Set path back to the code directory
setwd("A:/Praktikum_Chris/R/code")


#########################
# test for significance #
#########################

# import 100 sample runs.

test <- read.csv("A:/Praktikum_Chris/output/random_files_2023-09-13_08-44-00/robinson_run_results.csv")

five_random <- test[[3]]
three_random  <- test[[4]]

################################
# check for normal destribution#
################################

test1 <- shapiro.test(test[[3]])
test2 <- shapiro.test(test[[4]])



message("shapiro-wilk-test for 5_UTR-random runs\n")
print(test1)
cat("\n")
message("shapiro-wilk-test for 3_UTR-random runs\n")
print(test2)
cat("\n")

# p> 0.05 -> normal destribution. p<0.05 -> not normaly dest.

####################
# not normaly dist #
####################

# Perform the Wilcoxon rank-sum test for 5 utr
result_5 <- wilcox.test(five_random, y =  rf_distance_5_global_protein)

# Perform the Wilcoxon rank-sum test for 3 utr
result_3 <- wilcox.test(three_random, y =  rf_distance_3_global_protein)


setwd(scores_path)
# Open a text file for writing
sink("wilcoxon_test_results.txt")

# Capture and write the result_5 to the file
cat("Wilcoxon rank sum test with continuity correction\n")
cat("data: five_random and rf_distance_5_global_protein\n")
print(result_5)
cat("\n\n")  # Adding a newline for separation

# Capture and write the result_3 to the file
cat("Wilcoxon rank sum test with continuity correction\n")
cat("data: three_random and rf_distance_3_global_protein\n")
print(result_3)

# Close the file connection#
sink()

# Print a message indicating where the file was saved
cat("Results exported to wilcoxon_test_results.txt\n")

setwd("A:/Praktikum_Chris/R/code")


