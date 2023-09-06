######################
# compare to 6 files #
######################

threeutr_ofsix <- read.csv("A:/Praktikum_Chris/output/fasta_files_2023-09-06_10-03-48/allignmentscores/alignment_scores_global_3UTR.csv")

threeutr_ofsix <- threeutr_ofsix[, -1]

fiveutr_ofsix <-read.csv("A:/Praktikum_Chris/output/fasta_files_2023-09-06_10-03-48/allignmentscores/alignment_scores_global_5UTR.csv")

fiveutr_ofsix <- fiveutr_ofsix[, -1]


distances <- as.matrix(dist(threeutr_ofsix,threeUTR_allignment[["alignment_scores_global_UTR"]], method = "euclidean"))

distances2 <- as.matrix(dist(fiveutr_ofsix,fiveUTR_allignment[["alignment_scores_global_UTR"]], method = "euclidean"))



ranks <- apply(distances, 1, rank)
ranks2 <- apply(distances2, 1, rank)

# Find the top 1 most similar reference sequence for each of your 6 sequences
get_sim <- function(matrix, ranks) {
  top_reference_indices <- apply(ranks, 1, function(x) which.min(x))
  most_similar_references <- rownames(matrix)[top_reference_indices]
  return(most_similar_references)
}

sim_3 <- get_sim(threeUTR_allignment[["alignment_scores_global_UTR"]],ranks)
sim_5 <-get_sim(fiveUTR_allignment[["alignment_scores_global_UTR"]],ranks2)


# Print the results for threeUTR
for (i in 1:length(sim_3)) {
  cat("For Sequence", colnames(threeutr_ofsix)[i], "using threeUTR:", "\n")
  cat("Most Similar Reference Sequences:", sim_3[[i]], "\n")
}

# Print the results for fiveUTR
for (i in 1:length(sim_5)) {
  cat("For Sequence", colnames(fiveutr_ofsix)[i], "using fiveUTR:", "\n")
  cat("Most Similar Reference Sequences:", sim_5[[i]], "\n")
}
