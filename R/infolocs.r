## infolocs
##
##' Modified version of \code{\link[adehabitatLT]{infolocs}} that
##' returns \code{NULL} if \code{infolocs} exists, but \code{which} is
##' not a colum of it.
##'
##' @title Infolocs of an Object of Class ltraj
##' @seealso See \code{\link[adehabitatLT]{infolocs}} for further
##' details on the function and all available arguments.
##' @param perani Logical. Should the function return one element per burst
##' (\code{FALSE}, default) or per animal (\code{TRUE}).
##' @param simplify Logical. If a single variable is requested, should the
##' function return 1-column data frames (\code{FALSE}, default) or vector
##' (\code{TRUE}).
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Load puechcirc data, and add some random infolocs:
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
##' ## Try to retrieve the column `toto`:
##' infolocs(puechcirc, "toto")
infolocs <- function(ltraj, which, perani = FALSE, simplify = FALSE)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
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
        ## Give id/burst names to the result
        if (perani) {
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
        else names(re) <- burst(ltraj)
        return(re)
    }
    else {
        return(NULL)
    }
}
