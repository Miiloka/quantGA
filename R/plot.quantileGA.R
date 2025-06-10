#' Plot diagnostics for a \code{quantileGA} model.
#'
#' @param x \code{quantileGA}\cr Object returned by \code{GA.QR}.
#' @param type \code{character}\cr Type of plot. Either \code{"convergence"}
#'   to show the evolution of the fitness or \code{"residuals"} to plot
#'   residuals. Defaults to \code{"convergence"}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#'
#' @returns Invisibly returns \code{NULL}. A plot is produced as a side effect.
#' @export
plot.quantileGA <- function(x, type = c("convergence", "residuals"), ...) {
    type <- match.arg(type)
    if (type == "convergence") {
        generations <- seq_along(x$history$best_fitness)
        plot(generations, x$history$best_fitness,
             type = "l", col = "blue", lwd = 2,
             xlab = "Generation", ylab = "Fitness",
             main = "Genetic algorithm convergence", ...)
        lines(generations, x$history$mean_fitness, col = "red", lwd = 1, lty = 2)
        legend("topright", c("Best fitness", "Mean fitness"),
               col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1))
    } else if (type == "residuals") {
        plot(x$fitted_values, x$residuals,
             xlab = "Fitted values", ylab = "Residuals",
             main = "Residual plot", ...)
        abline(h = 0, col = "red", lty = 2)
    }
}


