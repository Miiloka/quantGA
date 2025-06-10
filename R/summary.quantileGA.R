#' Summary method for \code{quantileGA} objects.
#'
#' @param object \code{quantileGA}\cr Object returned by \code{GA.QR}.
#' @param ... Additional arguments passed to or from other methods. Currently
#'   ignored.
#'
#' @returns Invisibly returns \code{object} after printing a short summary.
#' @export
summary.quantileGA <- function(object, ...) {
    print(object)
    cat("\n Additional information:\n")
    cat(sprintf("   - Residual variance: %.6f\n", var(object$residuals)))
    invisible(object)
}
