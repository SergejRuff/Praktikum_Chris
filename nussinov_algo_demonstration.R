
sequence <- c(1, 2, 3, 4, 5)
n <- length(sequence)
for (l in 2:n) {
  for (i in 1:(n - l + 1)) {
    j <- i + l - 1

    cat("i =", i, "j =", j, "\n")
  }
}





###############################
# Zuker #
#########

