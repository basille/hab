## plotNAltraj
##
##' Five new arguments to allow a better control by the user.
##'
##' @title Highlighting the Patterns in Missing Values in Trajects
##' @param perani logical.  If \code{FALSE} (default), one plot is
##' drawn for each value of \code{burst}. If \code{TRUE}, one plot is
##' drawn for each value of \code{id}, and the several bursts are
##' superposed on the same plot for a given animal.
##' @param addlines Logical.  Indicates whether lines should be added
##' to the plot.
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
##' @seealso See \code{\link[adehabitatLT]{plotNAltraj}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##' adehabitatLT:::plotNAltraj(puechcirc)
##' plotNAltraj(puechcirc, perani = TRUE, addlines = FALSE, mfrow = c(1,
##'     2), ppar = list(pch = 15, cex = 0.5))
plotNAltraj <- function(x, perani = FALSE, addlines = TRUE, mfrow,
    ppar = list(pch = 16, cex = 0.3), lpar = list(), ...)
{
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
    lapply(x, function(i) {
        ## Allows point modification
        ## plot(i$date, is.na(i$x), ylim = c(0, 1), pch = 16, cex = 0.3,
        ##     xlab = "Time", ylab = "Missing values", main = attr(i,
        ##         "burst"), ...)
        do.call(plot, c(list(i$date, is.na(i$x), ylim = c(0,
            1), xlab = "Time", ylab = "Missing values", main = attr(i,
            "burst"), ...), ppar))
        ## End of modification
        if (addlines)
            ## Allows line modification
            ## lines(i$date, is.na(i$x))
            do.call(lines, c(list(i$date, is.na(i$x)), lpar))
            ## End of modification
    })
    invisible()
}
