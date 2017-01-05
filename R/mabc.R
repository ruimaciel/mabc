# Maciel's implementation of artificial bee colony algorithm for optimization problems with a binary search space
#
# The algorithm is implemented as a S4 class, it's developed for minimization problems, and doesn't apply a cross-over operator

library("parallel")

MABC <- setClass(
	Class="MABC",
	
	slots = list(
		exploitation_limit = "numeric",
		n_iteration = "numeric",
		n_dimensional_input_space = "numeric"
	),
	
	prototype=list(
		exploitation_limit=10,
		n_iteration = 0,
		n_dimensional_input_space = 0	 # must be set
	)
)


setGeneric(	name="initialize",
			def <- function(theObject, n_bees, n_dimensional) {
				standardGeneric("initialize")
			}
)

setMethod(f="initialize",
		  signature = "MABC",
		  definition = function(theObject, n_bees, n_dimensional)
		  {
		  		bees <- lapply(seq(1,n_bees), function(i) generate_random_bee(theObject, n_dimensional) )
		  		
		  		theObject@n_iteration <- 1;
		  		
		  		return(bees)
		  }
)


setGeneric(name="generate_random_bee",
		   def <- function(theObject, n_dimensional) {
		   	standardGeneric("generate_random_bee")
		   }
)


setMethod(f="generate_random_bee",
		  signature="MABC",
		  definition = function(theObject, n_dimensional)
		  {
		  	
		  	bee <- list(
		  		iteration = theObject@n_iteration,
		  		x = sample(c(TRUE,FALSE), n_dimensional, replace = TRUE),
		  		fitness_value = Inf
		  	)
		  	
		  	return(bee)
		  }
)


setGeneric(name="mutate",
		   def <- function(theObject, bee) {
		   	standardGeneric("mutate")
		   }
)


setMethod(f="mutate",
		  signature = "MABC",
		  definition <- function(theObject, bee) {
		  		# applies a bit-flip on a random bit
		  		i = sample(1:length(bee$x),1)
		  		
		  		bee$x[[i]] <- !bee$x[[i]]
		  		
		  		return(bee)
		  }
)


setGeneric(name="employed_bees_stage",
		   def <- function(theObject, bees, objective_function, cl) {
		   		standardGeneric("employed_bees_stage")
		   }
)


setMethod(f="employed_bees_stage",
		  signature = "MABC",
		  definition <- function(theObject, bees, objective_function, cl) {
		  		# apply mutator operator to all bees
		  		
		  		mutated_bees <- lapply(bees, function(bee) mutate(theObject,bee))
		  		mutated_bees <- evaluate(theObject, objective_function, mutated_bees, cl)
		  		
		  		# check which mutated bee improves
		  		for(i in 1:length(mutated_bees)) {
		  			if(bees[[i]]$fitness_value > mutated_bees[[i]]$fitness_value) {
		  				bees[[i]] <- mutated_bees[[i]]
		  				bees[[i]]$iteration <- theObject@n_iteration;
		  			}
		  		}
		  	
		  		return(bees)
		  }
)


setGeneric(name="onlooker_bees_stage",
		   def <- function(theObject, bees, objective_function, n_employed_bees, cl) {
		   	standardGeneric("onlooker_bees_stage")
		   }
)


setMethod(f="onlooker_bees_stage",
		  signature = "MABC",
		  definition <- function(theObject, bees, objective_function, n_employed_bees, cl) {
		  		# apply selection operator
		  		
		  		fitness_score_calculation <- function(b) {
		  			f <- b$fitness_value
		  			if( f > 0) {
		  				return( 1/(1+f))
		  			}else {
		  				return(1+abs(f))
		  			}
		  		}
		  		
		  		fitness_score <- sapply(bees, fitness_score_calculation)
		  		
		  		for(i in 2:length(fitness_score)) {
		  			fitness_score[i] = fitness_score[i] + fitness_score[i-1];
		  		}
		  		fitness_score = fitness_score/max(fitness_score)
		  		
		  		selected_idx = findInterval( runif(n_employed_bees), fitness_score)+1

		  		selected_bees = bees[selected_idx]
		  		
		  		# apply mutator operator
		  		mutated_bees <- lapply(selected_bees, function(bee) mutate(theObject, bee))
		  		mutated_bees <- evaluate(theObject, objective_function, mutated_bees, cl)
		  	
		  		# update current bees
		  		for(i in selected_idx) {
		  			original_idx = selected_idx[i]
		  			
		  			#updates the old record if there is an update
		  			if( bees[[original_idx]]$fitness_value > mutated_bees[[i]]$fitness_value) {
		  				bees[[original_idx]] <- mutated_bees[[i]]
		  				bees[[original_idx]]$iteration <- theObject@n_iteration
		  			}
		  		}
		  		
		  		return(bees)
		  }
)


setGeneric(name="scout_bees_stage",
		   def <- function(theObject, bees, objective_function, cl) {
		   	standardGeneric("scout_bees_stage")
		   }
)


setMethod(f="scout_bees_stage",
		  signature = "MABC",
		  definition <- function(theObject, bees, objective_function, cl) {
		  		exhausted <- sapply(bees, function(b) theObject@n_iteration-b$iteration > theObject@exploitation_limit)
		  		
		  		if( any(exhausted)) {
		  			exhausted_bees <- bees[exhausted]
		  			replacement_bees <- lapply(exhausted_bees, function(b) generate_random_bee(theObject,length(b$x)))
		  			replacement_bees <- evaluate(theObject, objective_function, replacement_bees, cl)
		  			bees[exhausted] <- replacement_bees
		  		}
		  		
		  		return(bees)
		  }
)


setGeneric(	name="run",
	def <- function(theObject, objective_function, n_employed_bees, n_max_iterations, cl = NULL) {
		standardGeneric("run")
	}
)

setMethod(f="run",
		  signature = "MABC",
		  definition = function(theObject, objective_function, n_employed_bees, n_max_iterations, cl = NULL)
		  {
		  		theObject@n_iteration <- 0
		  		n_dimensional <- theObject@n_dimensional_input_space
		  		
		  		n_initial_bees <- n_employed_bees;
		  		
		  		employed_bees_list <- initialize(theObject, n_initial_bees, n_dimensional)
		  		employed_bees_list <- evaluate(theObject, objective_function, employed_bees_list, cl)
		  		
		  		best_bees = list()
		  		
		  		
		  		while(theObject@n_iteration < n_max_iterations) {
		  			
		  			#TODO apply artificial bee colony algorithm
		  			bees <- employed_bees_stage(theObject, employed_bees_list, objective_function, cl)
		  			
		  			bees <- onlooker_bees_stage(theObject, bees, objective_function, n_employed_bees, cl)
		  			
		  			# store best point from each iteration
		  			best_idx = which.min( sapply(bees, function(b) b$fitness_value))
		  			best_bees[[length(best_bees)+1]] = bees[[best_idx]]
		  			
		  			bees <- scout_bees_stage(theObject, bees, objective_function, cl)
		  			
		  			employed_bees_list <- bees
		  			
		  			# update iteration variables
		  			theObject@n_iteration <- theObject@n_iteration + 1
		  		}
		  		
		  		# produce output data structure
		  		best_idx = which.min( lapply(best_bees, function(b) b$fitness_value))
		  		best_bee = best_bees[[best_idx]]
		  		
		  		output = list(
		  			best_per_iteration = best_bees,
		  			absolute_best = best_bee
		  		)
		  		
		  		return(output)
		  }
)


setGeneric(	name="evaluate",
			def <- function(theObject, f, bee_list, cl = NULL) {
				standardGeneric("evaluate")
			}
)

setMethod(f="evaluate",
		  signature = "MABC",
		  definition = function(theObject, f, bee_list, cl = NULL)
		  {
			  	if (is.null(cl))
			  	{
			  		fitness_value <- lapply(bee_list, function(b) f(b$x))
			  	}
		  		else
			  	{
			  		fitness_value <- parLapply(cl, bee_list, function(b) f(b$x))
			  	}
		  		
		  		for(i in 1:length(bee_list)) {
		  			bee_list[[i]]$fitness_value <- fitness_value[[i]]
		  		}
		  	
		  		return(bee_list)
		  }
)


