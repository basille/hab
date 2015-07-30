## plotltr
##
##' Four new arguments to allow a better control by the user (note
##' that arguments \code{pch} and \code{cex} have been removed and are
##' now used in \code{ppar}).
##'
##' Also corrects a bug if the burst contains \code{infolocs}, in
##' which case other attributes (notably \code{burst} and \code{id})
##' were lost.
##'
##' @title Changes in Traject Parameters Over Time
##' @param perani logical.  If \code{FALSE} (default), one plot is
##' drawn for each value of \code{burst}. If \code{TRUE}, one plot is
##' drawn for each value of \code{id}, and the several bursts are
##' superposed on the same plot for a given animal.
##' @param mfrow A vector of the form \code{c(nr, nc)}, which allows
##' the user to define the numbers of rows (\code{nr}) and columns
##' (\code{nc}) in the device (the default uses
##' \code{n2mfrow(length(id))} if \code{length(id) <= 12}, and
##' \code{mfrow = c(3, 4)} otherwhise).
##' @param ppar A list of arguments that allows the user to modify
##' point display, using any argument available to
##' \code{points}. Default is \code{list(pch = 21, col = "black", bg =
##' "white")}.
##' @param lpar A list of arguments that allows the user to modify
##' line display, using any argument available to
##' \code{lines}. Default is \code{list()}, i.e. an empty list.
##' @seealso See \code{\link[adehabitatLT]{plotltr}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##'
##' adehabitatLT::plotltr(puechcirc, "cos(rel.angle)")
##' plotltr(puechcirc, "cos(rel.angle)")
##' \dontrun{
##' plotltr(puechcirc, "cos(rel.angle)", ppar = list(pch = 2, cex = 2),
##'     lpar = list(lty = 2, lwd = 2), mfrow = c(2, 1))}
##'
##' adehabitatLT::plotltr(puechcirc, "dist")
##' plotltr(puechcirc, "dist")
##' plotltr(puechcirc, "dist", ppar = list(pch = 3, col = "blue"),
##'     lpar = list(lty = 3, col = "red"), perani = TRUE)
##'
##' adehabitatLT::plotltr(puechcirc, "dx")
##' plotltr(puechcirc, "dx")
##' plotltr(puechcirc, "dx", ppar = list(col = rep(1:8, each = 6)),
##'     addlines = FALSE)
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
plotltr <- function(x, which = "dist", perani = FALSE, addlines = TRUE,
    mfrow, ppar = list(pch = 16, cex = 0.7), lpar = list(), ...)
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class ltraj")
    ## Allows user interaction with  'perani' and 'mfrow'
    ## opar <- par(mfrow = n2mfrow(length(x)))
    if (perani) {
        x <- bindltraj(x)
    }
    if(missing(mfrow)) {
        if(length(x) > 12)
            mfrow = c(3, 4)
        else
            mfrow = n2mfrow(length(x))
    }
    if(length(x) > prod(mfrow))
        opar <- par(mfrow = mfrow, ask = TRUE)
    else
        opar <- par(mfrow = mfrow)
    ## End of modification
    on.exit(par(opar))
    toto <- lapply(x, function(i) {
        ## Keep 'id' and 'burst' attributes
        ## if (!is.null(attr(i, "infolocs")))
        ##     i <- cbind(i, attr(i, "infolocs"))
        if (!is.null(attr(i, "infolocs"))) {
            id <- attr(i, "id")
            burst <- attr(i, "burst")
            i <- cbind(i, attr(i, "infolocs"))
            attr(i, "id") <- id
            attr(i, "burst") <- burst
        }
        ## End of modification
        ex <- parse(text = which)
        coin <- eval(ex, envir = i)
        ## Allows point modification
        ## plot(i$date, coin, main = attr(i, "burst"), xlab = "Time",
        ##     ylab = which, pch = pch, cex = cex, ...)
        do.call(plot, c(list(i$date, coin, main = attr(i, "burst"),
            xlab = "Time", ylab = which, ...), ppar))
        ## End of modification
        if (addlines)
            ## Allows line modification
            ## lines(i$date, coin)
            do.call(lines, c(list(i$date, coin), lpar))
            ## End of modification
    })
    invisible()
}
