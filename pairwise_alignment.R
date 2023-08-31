library(Biostrings)

alignment_results_local <- list()
alignment_results_global <- list()

# Loop through results_list
for (i in 1:length(results_list)) {
  name_i <- results_list[[i]]  # Get the name of sequence i
  prefix_i <- substring(names(name_i), 1, 6)  # Get the first 6 characters (5_UTR or 3_UTR)

  for (j in 1:length(results_list)) {
    if (i != j) {
      name_j <- results_list[[j]]  # Get the name of sequence j
      prefix_j <- substring(names(name_j), 1, 6)  # Get the first 6 characters (5_UTR or 3_UTR)

      # Perform pairwise alignment only if prefixes are the same
      if (prefix_i == prefix_j) {
        sequence_1 <- results_list[[i]][2,]
        sequence_2 <- results_list[[j]][2,]

        #print(paste("using", prefix_j,"and",prefix_i))

        alignment <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2, type = "local")

        lokal_seq <- c(alignedPattern(alignment), alignedSubject(alignment))
        alignment_results_local[[paste("alignment_", i, "_", j, sep = "")]] <- alignment

        alignment_g <- pairwiseAlignment(pattern = sequence_1, subject = sequence_2, type = "global")
        global_seq <- c(alignedPattern(alignment_g), alignedSubject(alignment_g))
        alignment_results_global[[paste("alignment_", i, "_", j, sep = "")]] <- alignment_g
      }
    }
  }
}
