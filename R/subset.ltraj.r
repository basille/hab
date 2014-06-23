##' Return subsets of a ltraj which meet conditions (over the descriptive
##' parameters of the ltraj or is infolocs).
##'
##' @title Subsetting a ltraj
##' @param x A ltraj object.
##' @param subset Logical expression indicating elements or rows to keep:
##' missing values are taken as false.
##' @param rec Logical, whether to recompute the ltraj parameters or not
##' (default = FALSE).
##' @S3method subset ltraj
##' @return A ltraj object.
##' @seealso See \code{\link[adehabitatLT]{which.ltraj}} to identify the
##' relocations fullfilling a condition.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Load puechcirc data and add some infolocs:
##' data(puechcirc)
##' info <- list(data.frame(A = rnorm(80), B = runif(80), C = rpois(80,
##'     1)), data.frame(A = rnorm(69, 10), B = runif(69, 10, 11),
##'     C = rpois(69, 10)), data.frame(A = rnorm(66, -10), B = runif(66,
##'     -10, -9), C = rpois(66, 100)))
##' ## Watch the row names:
##' info <- mapply(function(x, y) {
##'     row.names(x) <- row.names(y)
##'     return(x)
##' }, info, puechcirc, SIMPLIFY = FALSE)
##' infolocs(puechcirc) <- info
##'
##' ## Different subsets:
##' (xsub1 <- subset(x, dist > 200, rec = FALSE))
##' (xsub2 <- subset(x, C == 3, rec = TRUE))
##' (xsub3 <- subset(x, C == 3, rec = FALSE))
subset.ltraj <- function (x, subset, rec = FALSE, ...)
{
    ## Check that x is a ltraj
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")
    ## If no 'subset', returns the ltraj
    if (missing(subset))
        return(x)
    ## Converts to a data.frame
    x <- ld(x)
    ## Code from subset.data.frame
    r <- if (missing(subset))
        rep_len(TRUE, nrow(x))
    else {
        e <- substitute(subset)
        r <- eval(e, x, parent.frame())
        if (!is.logical(r))
            stop("'subset' must be logical")
        r & !is.na(r)
    }
    ## Subset the data frame, drop levels that are now absent
    x <- droplevels(x[r, ])
    ## Return the subsetted ltraj, recomputed or not
    return(dl(x, strict = rec))
}
