
sequence <- c(1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15)


n <- length(sequence)

for (l in 2:n) {
  for (i in 1:(n - l + 1)) {
    j <- i + l - 1

    cat("i =", i, "j =", j, "\n")
  }
}
