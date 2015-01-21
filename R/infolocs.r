## infolocs
##
##' Modified version of \code{\link[adehabitatLT]{infolocs}} that
##' returns \code{NULL} if \code{infolocs} exists, but \code{which} is
##' not a colum of it (note that \code{adehabitatLT::infolocs} already
##' returns \code{NULL} if there is no infolocs).
##'
##' @title Infolocs of an Object of Class ltraj
##' @seealso See \code{\link[adehabitatLT]{infolocs}} for further
##' details on the function and all available arguments.
##' @param by Character, replaces \code{perani}. Either \code{"burst"}
##' (identical to \code{perani = FALSE}), \code{"id"} (identical to
##' \code{perani = TRUE}), or \code{"none"} to return only one element
##' for all bursts together.
##' @param simplify Logical. If a single variable is requested, should
##' the function return 1-column data frames (\code{FALSE}, default)
##' or vector (\code{TRUE}).
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Load puechcirc data, and add some random infolocs:
##' data(puechcirc)
##' infolocs(puechcirc)
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
##' puechcirc
##'
##' ## Parameters "by" and "simplify"
##' infolocs(puechcirc, "A", simplify = TRUE)
##' infolocs(puechcirc, "A", by = "id", simplify = TRUE)
##' infolocs(puechcirc, "A", by = "none", simplify = TRUE)
##'
##' ## Try to retrieve the column `toto`:
##' infolocs(puechcirc, "toto")
infolocs <- function(ltraj, which, by = c("burst", "id", "none"),
    simplify = FALSE)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    ## Match the 'by' argument
    by <- match.arg(by)
    if (!is.null(attr(ltraj[[1]], "infolocs"))) {
        if (missing(which))
            which <- names(attr(ltraj[[1]], "infolocs"))
        ## If infolocs exists but which is not in it, returns NULL
        else if (all(!(which %in% names(attr(ltraj[[1]], "infolocs")))))
            return(NULL)
        ## To simplify 1 column data frames
        drop <- ifelse(simplify & length(which) == 1, TRUE, FALSE)
        re <- lapply(ltraj, function(y) {
            res <- attr(y, "infolocs")
            return(res[, names(res) %in% which, drop = drop])
        })
        ## Extract by "none" (all together)
        if (by == "none") {
            if (drop)
                re <- unlist(re)
            else re <- do.call(rbind, re)
        }
        ## Extract by "id"
        if (by == "id") {
            names(re) <- id(ltraj)
            ## Binds unique animals together (bug in case of factors,
            ## returns the numeric values of the levels)
            un <- unique(names(re))
            if (drop)
                re <- lapply(un, function(n) do.call(c, re[names(re) ==
                  n]))
            else
                re <- lapply(un, function(n) do.call(rbind,
                  re[names(re) == n]))
            names(re) <- un
        }
        ## Extract by "burst" or "none"
        if (by == "burst")
            names(re) <- burst(ltraj)
        return(re)
    }
    ## Return NULL if there is no such variable in infolocs
    else {
        return(NULL)
    }
}
