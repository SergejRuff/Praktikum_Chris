# Define a sequence of numbers
sequence <- c(1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15)

# Calculate the length of the sequence
L <- length(sequence)

# Nested loops to iterate over pairs (i, j) and print the result
for (i in 1:(L - 1)) {
  for (j in (i + 1):L) {
    # Print the values of i and j
    cat("i =", i, "j =", j, "\n")

    # Access the elements at positions i and j in the sequence
    element_i <- sequence[i]
    element_j <- sequence[j]


  }
}
