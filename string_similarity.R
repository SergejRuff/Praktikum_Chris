# Load the required packages
library(Biostrings)
library(stringdist)
library(DECIPHER)
library(msa)
library(seqinr)

alignment_results_local <- list()
alignment_results_global <- list()
for (i in 1:length(results_list)){
  for(j in 1:length(results_list)){
    if(i!=j){
      # Example dot-bracket notations from your results_list
      sequence_1 <- results_list[[i]][2,]
      sequence_2 <- results_list[[j]][2,]

      # Perform pairwise alignment
      alignment <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2,type="local")

      lokal_seq <- c(alignedPattern(alignment), alignedSubject(alignment))
      #BrowseSeqs(lokal_seq)
      alignment_results_local[[paste("alignment_", i, "_", j, sep = "")]] <- alignment

      alignment_g <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2,type="global")
      global_seq <- c(alignedPattern(alignment_g ), alignedSubject(alignment_g ))
      #BrowseSeqs(global_seq)
      alignment_results_global[[paste("alignment_", i, "_", j, sep = "")]] <- alignment_g
    }

  }

}


#####################################################################


# Create an empty data frame to store results
results_df_lev <- data.frame(name1 = character(0),
                         name2 = character(0),
                         lev_distance = numeric(0),
                         lev_distance_percentage = numeric(0))
for (i in 1:length(results_list)){
  for(j in 1:length(results_list)){
    if(i!=j){

      # Example dot-bracket notations from your results_list
      sequence_1 <- results_list[[i]][2,]
      sequence_2 <- results_list[[j]][2,]

      name1 <- names(results_list)[i]
      name2 <- names(results_list)[j]

      # Calculate the Levenshtein distance
      lev_distance <- stringdist::stringdist(sequence_1, sequence_2)

      cat("Levenshtein Distance:", lev_distance, "\n")
      # shows how many edits need to be taken to transform one sequence into another

      # Length of the longer sequence
      max_length <- max(nchar(sequence_1), nchar(sequence_2))

      # Calculate the Levenshtein distance as a percentage
      lev_distance_percentage <- (lev_distance / max_length) * 100

      # Add data to the data frame
      results_df_lev <- rbind(results_df_lev, data.frame(name1, name2, lev_distance, lev_distance_percentage))
    }

  }

}

# Print the results data frame
View(results_df_lev)

##############################################################################
# Extract the second row from each data frame
#dot_sequences <- lapply(results_list, function(df) df[1, ])


#aa_sequences <- AAStringSet(unlist(dot_sequences))
#testmsa <- msa(aa_sequences)
#print(testmsa, show="complete")



###################################################
#####################################################################

# Example dot-bracket notations from your results_list
#sequence_1 <- results_list[[5]][2,]
#sequence_2 <- results_list[[6]][2,]

# Perform pairwise alignment
#alignment <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2,type="local")

#lokal_seq <- c(alignedPattern(alignment), alignedSubject(alignment))
#BrowseSeqs(lokal_seq)

#alignment_g <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2,type="global")
#global_seq <- c(alignedPattern(alignment_g ), alignedSubject(alignment_g ))
#BrowseSeqs(global_seq)

###############
###################################
# Print the results data frame
#View(results_df_lev)

# Calculate the Levenshtein distance
#lev_distance <- stringdist::stringdist(sequence_1, sequence_2)

#cat("Levenshtein Distance:", lev_distance, "\n")
# shows how many edits need to be taken to transform one sequence into another

# Length of the longer sequence
#max_length <- max(nchar(sequence_1), nchar(sequence_2))

# Calculate the Levenshtein distance as a percentage
#lev_distance_percentage <- (lev_distance / max_length) * 100

#cat("Levenshtein Distance Percentage:", lev_distance_percentage, "%\n")
# what percentage of dots and points needs to be changed in the longer sequence to match teh short sequence.





# Initialize matrices
n <- length(list_gt3_UTR)
lev_distance_matrix <- matrix(0, n, n)
lev_distance_percentage_matrix <- matrix(0, n, n)

# Create a vector of sequence names
sequence_names <- names(list_gt3_UTR)

# Assign row and column names
rownames(lev_distance_matrix) <- sequence_names
colnames(lev_distance_matrix) <- sequence_names

rownames(lev_distance_percentage_matrix) <- sequence_names
colnames(lev_distance_percentage_matrix) <- sequence_names

# Loop through list_gt3_UTR
for (i in 1:n) {
  for (j in 1:n) {
    # Example dot-bracket notations from list_gt3_UTR
    sequence_1 <- list_gt3_UTR[[i]][2, ]
    sequence_2 <- list_gt3_UTR[[j]][2, ]

    # Calculate the Levenshtein distance
    lev_distance <- stringdist::stringdist(sequence_1, sequence_2)

    # Length of the longer sequence
    max_length <- max(nchar(sequence_1), nchar(sequence_2))

    # Calculate the Levenshtein distance as a percentage
    lev_distance_percentage <- (lev_distance / max_length) * 100

    # Store results in matrices
    lev_distance_matrix[i, j] <- lev_distance
    lev_distance_percentage_matrix[i, j] <- lev_distance_percentage
  }
}
