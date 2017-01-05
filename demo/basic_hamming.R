
# Example:

library(MABC)

mabc <- MABC(n_dimensional_input_space = 10)

f <- function(x) {
	return(sum(x))
}

results <- run(mabc, f, n_employed_bees = 8, n_max_iterations = 100)