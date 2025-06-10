# Good paper: https://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR09/TR09-02.pdf
# library(Rcpp)
# library(RcppArmadillo)
# sourceCpp("src/helpers.cpp")
# source("R/plot.quantileGA.R")
# source("R/print.quantileGA.R")
# source("R/summary.quantileGA.R")


.initPop <- function(pop_size, n_params, bounds, X = NULL, y = NULL) {
    lower <- bounds$lower
    upper <- bounds$upper

    # if bounds is scalar
    if (length(lower) == 1) lower <- rep(lower, n_params)
    if (length(upper) == 1) upper <- rep(upper, n_params)

    population <- matrix(0, nrow = pop_size, ncol = n_params)

    # Example initialization for demonstration
    n_uniform <- ceiling(pop_size * 0.6)
    n_gaussian <- ceiling(pop_size * 0.3)
    n_lm_based <- pop_size - n_uniform - n_gaussian

    # ----- Using U(lower, upper)
    for (i in seq_len(n_uniform)) {
        population[i, ] <- runif(n_params, lower, upper)
    }

    # ----- Using N(0, (upper - lower) / 4)
    start_idx <- n_uniform + 1
    end_idx <- n_uniform + n_gaussian
    if (end_idx >= start_idx) {
        for (i in start_idx:end_idx) {
            population[i, ] <- rnorm(
                n_params, 0,
                pmin(abs(upper - lower) / 4, 2) # avoid too large sd
            )
            population[i, ] <- pmax(pmin(population[i, ], upper), lower)
        }
    }

    # ----- Initialization by linear regression (beta^(0) = (X'X)^(-1) * X'y)
    if (!is.null(X) && !is.null(y) && n_lm_based > 0) {
        tryCatch( # fallback if X'X is not invertible
            {
                lm_coef <- solve(crossprod(X), crossprod(X, y))
                if (all(is.finite(lm_coef))) {
                    start_idx <- n_uniform + n_gaussian + 1
                    for (i in start_idx:pop_size) {
                        # add noise to diversify initialization
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

.tournoi <- function(fitness, tournament_size = 3, n_selections = 1) {
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

.crossoverUniform <- function(parent1, parent2, crossover_prob) {
    if (runif(1) >= crossover_prob) {
        return(rbind(parent1, parent2))
    }

    # Crossover mask
    crossover_mask <- runif(length(parent1)) < 0.5
    child1 <- ifelse(crossover_mask, parent1, parent2)
    child2 <- ifelse(crossover_mask, parent2, parent1)
    rbind(child1, child2)
}

.mutGaussAdapt <- function(individual, mutation_prob, mutation_strength,
                          bounds, generation = 1, max_generations = 100) {
    # Adapt mutation strength over time
    adaptation_factor <- 1 - (generation / max_generations) * 0.5
    adapted_strength <- mutation_strength * adaptation_factor

    # mutation
    mutation_mask <- runif(length(individual)) < mutation_prob

    if (any(mutation_mask)) {
        mutations <- rnorm(sum(mutation_mask), 0, adapted_strength)
        individual[mutation_mask] <- individual[mutation_mask] + mutations
        individual <- pmax(pmin(individual, bounds$upper), bounds$lower)
    }

    individual
}


.qrVar <- function(X, resid, tau,
                  bw = "Silverman",
                  se = c("Approx", "Powell")) {
    if (!is.matrix(X) && !is.data.frame(X)) {stop("X must be a matrix or data.frame")}
    if (!is.numeric(resid)) {stop("resid must be numeric")}
    if (tau <= 0 || tau >= 1) {stop("tau must be between 0 and 1")}
    if (length(resid) != nrow(X)) {stop("Length of resid must equal the number of rows of X")}

    n <- length(resid)
    p <- ncol(X)

    if (is.numeric(bw)) {
        h <- bw
    } else {
        # https://stats.stackexchange.com/questions/6670/which-is-the-formula-from-silverman-to-calculate-the-bandwidth-in-a-kernel-densi
        bw <- match.arg(bw)
        sigma_r <- min(sd(resid), IQR(resid) / 1.34)
        h <- 1.06 * sigma_r * n^(-1 / 5)
    }

    se <- match.arg(se)

    if (se == "Approx") {
        f0 <- sum(dnorm(resid / h)) / (n * h) # Gaussian KDE estimate of f(0)
        XtX_inv <- solve(crossprod(X)) # (X'X)^-1
        V <- (tau * (1 - tau) / (f0^2)) * XtX_inv
    } else if (se == "Powell") {
        # Powell (1991)
        S_hat <- crossprod(X) / n # (1/n) X'X
        inside <- abs(resid) <= h
        J_hat <- crossprod(X[inside, , drop = FALSE]) / (2 * n * h) # (1/(2n h)) sum{|r_i| <= h} x_i x_i'
        J_hat_inv <- solve(J_hat)
        V <- (tau * (1 - tau) / n) * (J_hat_inv %*% S_hat %*% J_hat_inv)
    }

    colnames(V) <- rownames(V) <- colnames(X)
    list(vcov = V, se = sqrt(diag(V)))
}




#' Fits a linear quantile regression model using a genetic algorithm.
#'
#' @param formula \code{formula}\cr Model formula specifying the response and
#'   predictors. If supplied, \code{X} and \code{y} are ignored. \code{NULL} by
#'   default.
#' @param data \code{data.frame}\cr Data in which to evaluate the formula.
#'   Only used when \code{formula} is provided.
#' @param X \code{matrix}\cr Design matrix for the regression. Defaults to
#'   \code{NULL} and will be created from \code{formula} when provided.
#' @param y \code{numeric}\cr Response vector. Defaults to \code{NULL}.
#' @param tau \code{numeric}\cr Quantile to be modelled (between 0 and 1).
#'   Defaults to \code{0.5}.
#' @param intercept \code{logical}\cr Should an intercept be included when a
#'   formula is supplied? Defaults to \code{TRUE}.
#' @param pop_size \code{integer}\cr Size of the population in the genetic
#'   algorithm. Defaults to \code{50}.
#' @param max_generations \code{integer}\cr Maximum number of generations to
#'   evolve. Defaults to \code{200}.
#' @param crossover_prob \code{numeric}\cr Probability of crossover between two
#'   parents. Defaults to \code{0.8}.
#' @param mutation_prob \code{numeric}\cr Probability that a coefficient is
#'   mutated. Defaults to \code{0.1}.
#' @param mutation_strength \code{numeric}\cr Standard deviation of the
#'   Gaussian mutation. Defaults to \code{0.5}.
#' @param bounds \code{list}\cr List with elements \code{lower} and
#'   \code{upper} giving numeric bounds for the coefficients. Defaults to
#'   \code{list(lower = -100, upper = 100)}.
#' @param tournament_size \code{integer}\cr Number of individuals taking part in
#'   tournament selection. Defaults to \code{3}.
#' @param elitism_rate \code{numeric}\cr Proportion of the best individuals
#'   carried over to the next generation. Defaults to \code{0.1}.
#' @param convergence_tolerance \code{numeric}\cr Not currently used. Intended
#'   relative tolerance for convergence. Defaults to \code{1e-6}.
#' @param convergence_generations \code{integer}\cr Number of successive
#'   generations without improvement before stopping. Defaults to \code{50}.
#' @param seed \code{numeric}\cr Optional seed for reproducibility. Defaults to
#'   \code{NULL}.
#'
#' @returns A list of class \code{quantileGA} containing the fitted coefficients,
#'   standard errors, fitted values, residuals and diagnostics.
#' @export
#' 
#' @examples
#' set.seed(45)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- runif(n, -1, 1)
#' beta_true <- c(2, -1, 0.5)
#' y <- as.numeric(cbind(1, x1, x2) %*% beta_true + rt(n, df = 3))
#' df <- data.frame(y = y, x1 = x1, x2 = x2)
#' 
#' result <- GA.QR(
#'   y ~ x1 + x2,
#'   data = df,
#'   tau = 0.5,
#'   pop_size = 50,
#'   max_generations = 100,
#'   tournament_size = 10,
#'   elitism_rate = 0.5,
#'   mutation_prob = 0.5,
#'   mutation_strength = 0.9,
#'   convergence_generations = 20,
#'   seed = 42
#' )
#' summary(result)
GA.QR <- function(
    formula = NULL,
    data = NULL,
    X = NULL,
    y = NULL,
    tau = 0.5,
    intercept = TRUE,
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
    seed = NULL) {
    # timeStart <- Sys.time()

    if (!is.null(seed)) set.seed(seed)
    if (!is.null(formula)) {
        mf <- model.frame(formula, data)
        y <- model.response(mf)
        X <- model.matrix(formula, mf)

        if (intercept && !("(Intercept)" %in% colnames(X))) {
            X <- cbind(1, X)
            colnames(X)[1] <- "(Intercept)"
        } else if (!intercept && "(Intercept)" %in% colnames(X)) {
            X <- X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
        }
    }

    y <- as.numeric(y)
    if (nrow(X) != length(y)) {stop("Number of rows of X must equal length of y")}
    if (tau <= 0 || tau >= 1) stop("tau must be between 0 and 1")
    if (pop_size < 4) stop("pop_size must be at least 4")
    if (max_generations < 1) stop("max_generations must be positive")
    n_params <- ncol(X)
    n_elite <- max(1, floor(pop_size * elitism_rate))

    population <- .initPop(pop_size, n_params, bounds, X, y)
    fitness <- popFitness(population, X, y, tau)
    if (!any(is.finite(fitness))) {
        stop("No valid fitness in the initial population. Check parameters.")
    }

    # tracking
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

        # Generate the new population
        offspring_population <- matrix(0, nrow = pop_size - n_elite, ncol = n_params)

        for (i in seq(1, pop_size - n_elite, by = 2)) {
            # Parent selection
            parent_indices <- .tournoi(fitness, tournament_size, 2)
            parent1 <- population[parent_indices[1], ]
            parent2 <- population[parent_indices[2], ]

            # Crossover
            offspring <- .crossoverUniform(parent1, parent2, crossover_prob)

            # Mutation
            for (j in seq_len(nrow(offspring))) {
                offspring[j, ] <- .mutGaussAdapt(
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

        if (generations_without_improvement >= convergence_generations) {
            best_fitness_history <- best_fitness_history[seq_len(generation)]
            mean_fitness_history <- mean_fitness_history[seq_len(generation)]
            break
        }
    } # end of generation loop
    final_generation <- min(generation, max_generations)
    fitted_values <- as.numeric(X %*% best_individual)
    residuals <- y - fitted_values
    varobj <- .qrVar(X, residuals, tau, se = "Powell")

    # Names of coefficients
    if (!is.null(colnames(X))) {
        names(best_individual) <- colnames(X)
    } else if (intercept && ncol(X) > 1) {
        names(best_individual) <- c("(Intercept)", paste0("X", seq_len(ncol(X) - 1)))
    } else if (!intercept) {
        names(best_individual) <- paste0("X", seq_len(ncol(X)))
    }

    result <- list(
        coefficients = best_individual,
        vcov = varobj$vcov,
        se = varobj$se,
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
            bounds = bounds,
            intercept = intercept
        ),
        call = match.call()
    )

    class(result) <- "quantileGA"
    # timeStop <- Sys.time()
    # cat(sprintf("\n\n Execution time: %.2f seconds\n", as.numeric(difftime(timeStop, timeStart, units = "secs"))))
    result
}
