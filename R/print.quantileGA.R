#' Print method for \code{quantileGA} objects.
#'
#' @param x \code{quantileGA}\cr Object produced by \code{GA.QR}.
#' @param ... Further arguments passed to or from other methods. Currently
#'   ignored.
#'
#' @returns Invisibly returns a matrix of coefficients and related statistics.
#' @export
print.quantileGA <- function(x, ...) {
    cat("Quantile regression via genetic algorithm\n")
    cat("========================================\n\n")
    cat(sprintf("Quantile: %.3f\n", x$tau))
    cat(sprintf("Final loss: %.6f\n", x$loss))
    cat(sprintf(
        "%s after %d generations\n",
        ifelse(x$convergence$converged, "Convergence", "No convergence"),
        x$convergence$iterations
    ))
    z <- x$coefficients / x$se
    p <- 2 * pnorm(-abs(z))
    out <- cbind(
        Estimate = x$coefficients,
        `Std. Error` = x$se,
        `t value` = z,
        `Pr(>|t|)` = p
    )
    printCoefmat(out, P.values = TRUE, has.Pvalue = TRUE)
    invisible(out)
}
