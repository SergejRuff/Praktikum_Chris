# Load the required packages
library(Biostrings)
library(stringdist)

# Example dot-bracket notations from your results_list
sequence_1 <- results_list[[3]][2,]
sequence_2 <- results_list[[4]][2,]

# Perform pairwise alignment
alignment <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2,type="local")


alignment_g <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2,type="global")

###############

# Calculate the Levenshtein distance
lev_distance <- stringdist::stringdist(sequence_1, sequence_2)

cat("Levenshtein Distance:", lev_distance, "\n")
# shows how many edits need to be taken to transform one sequence into another

# Length of the longer sequence
max_length <- max(nchar(sequence_1), nchar(sequence_2))

# Calculate the Levenshtein distance as a percentage
lev_distance_percentage <- (lev_distance / max_length) * 100

cat("Levenshtein Distance Percentage:", lev_distance_percentage, "%\n")
