## rec
##
##' Modified version of \code{\link[adehabitatLT]{rec}} that keeps the
##' original \code{row.names}. Also throws an error if the \code{row.names}
##' are different between the ltraj and its infolocs.
##'
##' @title Recalculates the descriptive parameters of a ltraj
##' @seealso See \code{\link[adehabitatLT]{rec}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##' (bla <- rec(puechcirc))
##' head(puechcirc[[2]])
##' head(bla[[2]])
rec <- function (x, slsp = c("remove", "missing"))
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")
    lif <- infolocs(x)
    if (!is.null(lif)) {
        ## Simply stop if row.names are different between ltraj and its
        ## infolocs
        for (i in 1:length(x)) {
        ##     if (!all(row.names(x[[i]]) %in% row.names(lif[[i]]))) {
        ##         x[[i]] <- x[[i]][row.names(x[[i]]) %in%
        ##           row.names(lif[[i]]), drop = FALSE]
        ##         attr(x[[i]], "infolocs") <- lif[[i]]
        ##     }
        ##     if (!all(row.names(lif[[i]]) %in% row.names(x[[i]]))) {
        ##         lif[[i]] <- lif[[i]][row.names(lif[[i]]) %in%
        ##           row.names(x[[i]]), drop = FALSE]
        ##         attr(x[[i]], "infolocs") <- lif[[i]]
        ##     }
	if (!all(row.names(lif[[i]]) == row.names(x[[i]])))
             stop("The infolocs component should have the same row.names
as the ltraj object")
        }
    }
    slsp <- match.arg(slsp)
    if (attr(x, "typeII")) {
        y <- adehabitatLT:::.traj2df(adehabitatLT:::.ltraj2traj(x))
        if (!is.null(lif)) {
            infol <- do.call("rbind", lif)
            al <- as.ltraj(xy = y[, c("x", "y")], date = y$date,
                id = y$id, burst = y$burst, slsp = slsp, typeII = TRUE,
                infolocs = infol)
        }
        else {
            al <- as.ltraj(xy = y[, c("x", "y")], date = y$date,
                id = y$id, burst = y$burst, slsp = slsp, typeII = TRUE,
                infolocs = lif)
        }
        ## Get the row names from the initial ltraj
        for (i in 1:length(al))
            row.names(al[[i]]) <- row.names(x[[i]])
        return(al)
    }
    else {
        attr(x, "typeII") <- TRUE
        y <- adehabitatLT:::.traj2df(adehabitatLT:::.ltraj2traj(x))
        if (!is.null(infolocs(x))) {
            infol <- do.call("rbind", infolocs(x))
            al <- as.ltraj(xy = y[, c("x", "y")], id = y$id,
                burst = y$burst, slsp = slsp, typeII = FALSE,
                infolocs = infol)
        }
        else {
            al <- as.ltraj(xy = y[, c("x", "y")], id = y$id,
                burst = y$burst, slsp = slsp, typeII = FALSE)
        }
        ## Get the row names from the initial ltraj
        for (i in 1:length(al))
            row.names(al[[i]]) <- row.names(x[[i]])
        return(al)
    }
}
