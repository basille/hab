## acf.test
##
##' Test of autocorrelation and partial autocorrelation in SSF models,
##' based on the estimation of autocorrelation functions (ACF).
##'
##' @title Test of autocorrelation in SSFs
##' @param residuals A vector of residuals on which to compute the
##' ACF. The residuals must be sorted chronologically.
##' @param id A vector of corresponding animal IDs.
##' @param type A character string giving the type of ACF to be
##' computed.  Allowed values are \code{correlation} (default),
##' \code{covariance} or \code{partial}.
##' @param ci A numeric giving the confidence interval to be used to
##' test the ACF.
##' @return A list, with the following parameters:
##' \itemize{
##' \item \code{acfk}: a list with the individual autocorrelation
##' functions
##' \item \code{threshold}: a vector with individual thresholds
##' \item \code{lag}: a vector with the resulting individual lags}
##' @seealso See \code{\link[stats]{acf}} for further details on the
##' ACF.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
acf.test <- function(residuals, id, type = c("correlation", "covariance",
    "partial"), ci = 0.95) {
    ## Check the 'type'
    type <- match.arg(type)
    ## Compute the ACF for each individual
    acfk <- lapply(levels(id), function(x) acf(residuals[id ==
        x], type = type, plot = FALSE))
    ## Compute the CIs to find the threshold (see stats:::plot.acf)
    threshold <- unlist(lapply(acfk, function(x) qnorm((1 + ci)/2)/sqrt(x$n.used)))
    ## Compute the lag after which there is no autocorrelation, given
    ## the CI level
    lag <- unlist(lapply(1:length(acfk), function(i) which(acfk[[i]]$acf <
        threshold[i])[1]))
    ## Return a list with individuals ACF (list), thresholds (vector)
    ## and lags (vector)
    return(list(acfk = acfk, threshold = threshold, lag = lag))
}

