# Load the required packages
library(dplyr)
library(gt)

# Load data
threeutr_ofsix <- read.csv("A:/Praktikum_Chris/output/fasta_files_2023-09-06_10-03-48/allignmentscores/alignment_scores_global_3UTR.csv")
threeutr_ofsix <- threeutr_ofsix[, -1]

fiveutr_ofsix <- read.csv("A:/Praktikum_Chris/output/fasta_files_2023-09-06_10-03-48/allignmentscores/alignment_scores_global_5UTR.csv")
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

sim_3 <- get_sim(threeUTR_allignment[["alignment_scores_global_UTR"]], ranks)
sim_5 <- get_sim(fiveUTR_allignment[["alignment_scores_global_UTR"]], ranks2)

# Create data frames
threeUTR_table_data <- data.frame(
  Sequence = colnames(threeutr_ofsix),
  Most_Similar_Reference_Sequences = sim_3
)

fiveUTR_table_data <- data.frame(
  Sequence =  colnames(fiveutr_ofsix),
  Most_Similar_Reference_Sequences = sim_5
)

# Create gt tables
threeUTR_table <- gt(threeUTR_table_data) %>%
  tab_header(title = "Most Similar Reference Sequences for three UTR")

fiveUTR_table <- gt(fiveUTR_table_data) %>%
  tab_header(title = "Most Similar Reference Sequences for five UTR")

# Print the tables
print(threeUTR_table)
print(fiveUTR_table)
