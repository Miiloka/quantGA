checkFun <- function(residuals, tau) {
    if (tau <= 0 || tau >= 1) {
        stop("tau must be between 0 and 1")
    }

  pos <- residuals >= 0
  sum(residuals * (pos * tau + (!pos) * (tau - 1)))
}

objFun <- function(beta, X, y, tau) {
    if (length(beta) != ncol(X)) {
        return(.Machine$double.xmax)
    }

    tryCatch(
        {
            fitted_values <- as.numeric(X %*% beta)
            residuals <- y - fitted_values

            if (any(!is.finite(fitted_values))) {
                return(.Machine$double.xmax)
            }

            checkFun(residuals, tau)
        },
        error = function(e) {
            .Machine$double.xmax
        }
    )
}

popFitness <- function(pop, X, y, tau) {
  resid_mat <- - tcrossprod(pop, X)
  resid_mat <- sweep(resid_mat, 2, y, "+")

  pos  <- resid_mat >= 0
  w    <- ifelse(pos,  tau, tau - 1)
  matrixStats::rowSums2(resid_mat * w) # plus rapide que la fonction R
}








initPop <- function(pop_size, n_params, bounds, X = NULL, y = NULL) {
    lower <- bounds$lower
    upper <- bounds$upper

    # Gestion des bornes scalaires
    if (length(lower) == 1) lower <- rep(lower, n_params)
    if (length(upper) == 1) upper <- rep(upper, n_params)

    population <- matrix(0, nrow = pop_size, ncol = n_params)

    # Example initialization using an "hybrid strategy" for demonstration
    n_uniform <- ceiling(pop_size * 0.6)
    n_gaussian <- ceiling(pop_size * 0.3)
    n_lm_based <- pop_size - n_uniform - n_gaussian

    # ----- Initialisation U(lower, upper)
    for (i in seq_len(n_uniform)) {
        population[i, ] <- runif(n_params, lower, upper)
    }

    # ----- Initialisation N(0, (upper - lower) / 4)
    start_idx <- n_uniform + 1
    end_idx <- n_uniform + n_gaussian
    if (end_idx >= start_idx) {
        for (i in start_idx:end_idx) {
            population[i, ] <- rnorm(
                n_params, 0,
                pmin(abs(upper - lower) / 4, 1)
            )
            population[i, ] <- pmax(pmin(population[i, ], upper), lower)
        }
    }

    # ----- Initialization by linear regression (Beta = (X'X)^(-1) * X'y)
    if (!is.null(X) && !is.null(y) && n_lm_based > 0) {
        tryCatch(# fallback if X'X is not invertible
            {
                lm_coef <- solve(crossprod(X), crossprod(X, y))
                if (all(is.finite(lm_coef))) {
                    start_idx <- n_uniform + n_gaussian + 1
                    for (i in start_idx:pop_size) {
                        # add noise for more diversity
                        noise <- rnorm(n_params, 0, abs(lm_coef) * 0.1 + 0.1)
                        population[i, ] <- lm_coef + noise
                        population[i, ] <- pmax(pmin(population[i, ], upper), lower)
                    }
                }
            },
            error = function(e) {
                for (i in (n_uniform + n_gaussian + 1):pop_size) {
                    population[i, ] <- runif(n_params, lower, upper)
                }
            }
        )
    }

    population
}

tournoi <- function(fitness, tournament_size = 3, n_selections = 1) {
    n_individuals <- length(fitness)

    if (n_individuals == 0) {
        return(integer(0))
    }

    valid_indices <- which(is.finite(fitness))
    if (length(valid_indices) == 0) {
        return(sample.int(n_individuals, n_selections, replace = TRUE))
    }

    selected <- integer(n_selections)

    for (i in seq_len(n_selections)) {
        tournament_indices <- sample(valid_indices,
            min(tournament_size, length(valid_indices)),
            replace = TRUE
        )

        selected[i] <- tournament_indices[which.min(fitness[tournament_indices])]
    }

    selected
}

crossoverUniform <- function(parent1, parent2, crossover_prob) {
    if (runif(1) >= crossover_prob) {
        return(rbind(parent1, parent2))
    }

    # Masque de croisement
    crossover_mask <- runif(length(parent1)) < 0.5
    child1 <- ifelse(crossover_mask, parent1, parent2)
    child2 <- ifelse(crossover_mask, parent2, parent1)
    rbind(child1, child2)
}

mutGaussAdapt <- function(individual, mutation_prob, mutation_strength,
                          bounds, generation = 1, max_generations = 100) {
    # Adaptation de la force de mutation au cours du temps
    adaptation_factor <- 1 - (generation / max_generations) * 0.5
    adapted_strength <- mutation_strength * adaptation_factor

    # Application de la mutation
    mutation_mask <- runif(length(individual)) < mutation_prob

    if (any(mutation_mask)) {
        mutations <- rnorm(sum(mutation_mask), 0, adapted_strength)
        individual[mutation_mask] <- individual[mutation_mask] + mutations
        individual <- pmax(pmin(individual, bounds$upper), bounds$lower)
    }

    individual
}

GA.QR <- function(
    X, y,
    tau = 0.5,
    pop_size = 50,
    max_generations = 200,
    crossover_prob = 0.8,
    mutation_prob = 0.1,
    mutation_strength = 0.5,
    bounds = list(lower = -100, upper = 100),
    tournament_size = 3,
    elitism_rate = 0.1,
    convergence_tolerance = 1e-6,
    convergence_generations = 50,
    verbose = TRUE,
    seed = NULL) {

    timeStart <- Sys.time()

    # Validation, intialisation
    if (!is.null(seed)) set.seed(seed)
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data.frame")
    }
    X <- as.matrix(X)
    y <- as.numeric(y)
    if (nrow(X) != length(y)) {
        stop("Number of rows of X must equal length of y")
    }
    if (tau <= 0 || tau >= 1) stop("tau must be between 0 and 1")
    if (pop_size < 4) stop("pop_size must be at least 4")
    if (max_generations < 1) stop("max_generations must be positive")
    n_params <- ncol(X)
    n_elite <- max(1, floor(pop_size * elitism_rate))
    if (verbose) {
        cat(sprintf(
            "GA start: %d parameters, %d individuals, %d generations max\n",
            n_params, pop_size, max_generations
        ))
    }

    population <- initPop(pop_size, n_params, bounds, X, y)
    fitness <- popFitness(population, X, y, tau)
    if (!any(is.finite(fitness))) {
        stop("No valid fitness in the initial population. Check parameters.")
    }




    best_fitness_history <- numeric(max_generations)
    mean_fitness_history <- numeric(max_generations)
    best_individual <- population[which.min(fitness), ]
    best_fitness <- min(fitness[is.finite(fitness)])
    generations_without_improvement <- 0

    for (generation in seq_len(max_generations)) {
        # Sort population by fitness (best to worst)
        fitness_order <- order(fitness)
        population <- population[fitness_order, , drop = FALSE]
        fitness <- fitness[fitness_order]

        # Elitism: keep the best
        elite_population <- population[seq_len(n_elite), , drop = FALSE]
        elite_fitness <- fitness[seq_len(n_elite)]

        # nouvelle population
        offspring_population <- matrix(0, nrow = pop_size - n_elite, ncol = n_params)

        for (i in seq(1, pop_size - n_elite, by = 2)) {
            # Parent selection
            parent_indices <- tournoi(fitness, tournament_size, 2)
            parent1 <- population[parent_indices[1], ]
            parent2 <- population[parent_indices[2], ]

            # Crossover
            offspring <- crossoverUniform(parent1, parent2, crossover_prob)

            # Mutation
            for (j in seq_len(nrow(offspring))) {
                offspring[j, ] <- mutGaussAdapt(
                    offspring[j, ], mutation_prob, mutation_strength, bounds,
                    generation, max_generations
                )
            }

            # Add to population
            if (i < pop_size - n_elite) {
                offspring_population[i, ] <- offspring[1, ]
            }
            if (i + 1 <= pop_size - n_elite) {
                offspring_population[i + 1, ] <- offspring[2, ]
            }
        }

        # New population
        population <- rbind(elite_population, offspring_population)

        # new fitness
        offspring_fitness <- popFitness(offspring_population, X, y, tau)
        fitness <- c(elite_fitness, offspring_fitness)

        # Update best individual
        current_best_fitness <- min(fitness[is.finite(fitness)])
        if (current_best_fitness < best_fitness) {
            best_fitness <- current_best_fitness
            best_individual <- population[which.min(fitness), ]
            generations_without_improvement <- 0
        } else {
            generations_without_improvement <- generations_without_improvement + 1
        }

        best_fitness_history[generation] <- best_fitness
        mean_fitness_history[generation] <- mean(fitness[is.finite(fitness)])

        # verbose de debug
        if (verbose && (generation %% 10 == 0 || generation == max_generations)) {
            cat(sprintf(
                "Gen %3d | Best: %.6f | Mean: %.6f | No improve: %d\n",
                generation, best_fitness, mean_fitness_history[generation],
                generations_without_improvement
            ))
        }

        if (generations_without_improvement >= convergence_generations) {
            if (verbose) {
                cat(sprintf("Stopped by convergence after %d generations\n", generation))
            }
            best_fitness_history <- best_fitness_history[seq_len(generation)]
            mean_fitness_history <- mean_fitness_history[seq_len(generation)]
            break
        }
    } # end of generation loop

    final_generation <- min(generation, max_generations)
    fitted_values <- as.numeric(X %*% best_individual)
    residuals <- y - fitted_values

    result <- list(
        coefficients = best_individual,
        fitted_values = fitted_values,
        residuals = residuals,
        loss = best_fitness,
        tau = tau,
        convergence = list(
            converged = generations_without_improvement >= convergence_generations,
            iterations = final_generation,
            final_loss = best_fitness
        ),
        history = list(
            best_fitness = best_fitness_history,
            mean_fitness = mean_fitness_history
        ),
        algorithm_settings = list(
            pop_size = pop_size,
            max_generations = max_generations,
            crossover_prob = crossover_prob,
            mutation_prob = mutation_prob,
            mutation_strength = mutation_strength,
            tournament_size = tournament_size,
            elitism_rate = elitism_rate,
            bounds = bounds
        ),
        call = match.call()
    )

    class(result) <- "quantileGA"
    timeStop <- Sys.time()
    cat(sprintf("\n\n Execution time: %.2f seconds\n", as.numeric(difftime(timeStop, timeStart, units = "secs"))))
    result
}

print.quantileGA <- function(x, ...) {
    cat("Quantile Regression by Genetic Algorithm\n")
    cat("===========================================\n\n")
    cat(sprintf("Target quantile: %.3f\n", x$tau))
    cat(sprintf("Final loss: %.6f\n", x$loss))
    cat(sprintf(
        "%s after %d generations\n",
        ifelse(x$convergence$converged, "Convergence", "No convergence"),
        x$convergence$iterations
    ))
    cat("\nEstimated coefficients:\n")
    print(x$coefficients)
}







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exemple
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(42)
n <- 500
x1 <- rnorm(n)
x2 <- runif(n, -1, 1)
X <- cbind(1, x1, x2) # include intercept
beta_true <- c(2, -1, 0.5)
errors <- rt(n, df = 3) # heavy-tailed
y <- as.numeric(X %*% beta_true + errors)

result <- GA.QR(
    X, y,
    tau = 0.5,
    pop_size = 100,
    tournament_size = 3, 
    elitism_rate = 0.1,
    max_generations = 500,
    mutation_prob = 0.5,
    mutation_strength = 0.5,
    convergence_generations = 75,
    verbose = FALSE,
    seed = 42
)

qrSOTA <- quantreg::rq(y ~ x1 + x2, tau = 0.5, data = data.frame(y, x1, x2))$coefficients # quasiment pareil !!

# print(result)
# result$convergence

data.frame(
    `True values` = beta_true,
    Estimates = result$coefficients,
    `Absolute error` = beta_true - result$coefficients,
    `Relative error` = (beta_true - result$coefficients) / beta_true * 100,
    `Estimate by quantreg::rq()` = qrSOTA,
    check.names = FALSE # keep spaces in column names
)



# profilage
# library(profvis)
# aaaa <- profvis({
#     result <- GA.QR(
#         X, y,
#         tau = 0.5,
#         pop_size = 200,
#         max_generations = 500,
#         mutation_strength = 0.25,
#         convergence_generations = 60,
#         verbose = FALSE,
#         seed = 42
#     )
# })

# aaaa
