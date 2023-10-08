protein_sequences <- read.tree(
  file="A:/Praktikum_Chris/data/nido_roniviruses_n6_2908/nido_invertebrate_new.nwk")

setwd(scores_path)

utr_nj <- function(matrix,filename,test){


  distance_matrix <- as.dist(matrix)

  njtree <- nj(distance_matrix)

 # Check for negative branch lengths and adjust them
  neg_branch_lengths <- which(njtree$edge.length < 0)
  if (length(neg_branch_lengths) > 0) {
    cat("Negative branch lengths found and adjusted.\n")
    for (i in neg_branch_lengths) {
      diff_length <- abs(njtree$edge.length[i])
      njtree$edge.length[i] <- 0
      parent_node <- njtree$edge[i, 1]
      njtree$edge.length[parent_node] <- njtree$edge.length[parent_node] + diff_length
    }
  }

  #set midpoint.root
  njtree <- midpoint.root(njtree)


  #plot results.
  png(filename, width = 1920, height = 1080)

  # Create a larger plot area
  par(mfrow = c(1, 1))
  plot(njtree,main = paste("njplot for",test), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,use.edge.length = TRUE)



  dev.off()

  return(njtree)


}


tanglegram_plot <- function(x, y, all_x, all_y, filename) {

  # Check if x and y are valid phylogenetic trees
  if (!is(x, "phylo") || !is(y, "phylo")) {
    stop("Input x and y must be valid phylogenetic trees.")
  }

  # Check if all_x and all_y are character vectors
  if (!is.character(all_x) || !is.character(all_y)) {
    stop("all_x and all_y must be character vectors.")
  }

  pattern <- ">(5|3)_UTR_of_"

  tiplabel1 <- x$tip.label

  tiplabel1 <- tiplabel1[!grepl("segment[2-9]$", tiplabel1)]

  print("Length of tiplabel1:")
  print(length(tiplabel1))

  tiplabel2 <- y$tip.label

  print("Length of tiplabel2:")
  print(length(tiplabel2))

  association_matrix <- matrix(0, nrow = length(tiplabel1), ncol = 2)

  # Iterate through tiplabel1 and find matches in tiplabel2
  for (i in 1:length(tiplabel1)) {
    # Check if tiplabel1[i] matches the pattern
    if (grepl(pattern, tiplabel1[i])) {
      # Extract the NC_number from tiplabel1
      test_name <- sub(">(5|3)_UTR_of_", "", tiplabel1[i])

      # Find matching elements in tiplabel2 based on the NC_number
      matching_indices <- which(str_detect(test_name, tiplabel2))

      # Check if there are matching elements in tiplabel2
      if (length(matching_indices) > 0) {
        matching_tip2 <- tiplabel2[matching_indices]
        association_matrix[i, 1] <- tiplabel1[i]
        association_matrix[i, 2] <- paste(matching_tip2, collapse = ", ")  # Combine matching labels
      } else {
        print(paste("No match found for", tiplabel1[i]))
      }
    }
  }

  # Print the association_matrix
  print("Association matrix:")
  print(association_matrix)

  tangleplot <- cophylo(x, y, assoc = association_matrix)

  # Plot results with adjusted parameters
  png(filename, width = 1920, height = 1080)

  # Create a larger plot area
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot(tangleplot, use.edge.length = TRUE, cex = 0.8)  # Adjust cex for label size
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
tree_5_global$tip.label <- gsub(">5_UTR_of_", "", tree_5_global$tip.label)

# Update tip labels
tree_3_global$tip.label <- gsub(">3_UTR_of_", "", tree_3_global$tip.label)

# Update tip labels
protein_sequences$tip.label



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
