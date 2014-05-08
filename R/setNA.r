## setNA
##
##' The function does not allow a numeric for \code{date.ref}
##' anymore. In addition, a warning is issued if the length of
##' \code{date.ref} is not 1 or the length of the \code{ltraj}
##' object. See \code{\link[adehabitatLT]{setNA}} for further details
##' on the function and all available arguments.
##'
##' @title Place Missing Values in Objects of Class \code{ltraj}
##' @seealso See \code{\link[adehabitatLT]{setNA}} for further details
##' on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
setNA <- function(ltraj, date.ref, dt, tol = dt/10, units = c("sec",
    "min", "hour", "day"), ...)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!attr(ltraj, "typeII"))
        stop("ltraj should be of type II (time recorded)")
    ## Does not allow for numeric date.ref + check inheritance from
    ## POSIXt (instead of POSIXct) first + check if length of date.ref
    ## is 1 or the length of the ltraj object (issue a warning if not)
    ## if (is.numeric(date.ref)) {
    ##     class(date.ref) <- c("POSIXct", "POSIXt")
    ##     attr(date.ref, "tzone") <- attr(ltraj[[1]]$date, "tzone")
    ## }
    if (!inherits(date.ref, "POSIXt"))
        stop("date.ref should be of class \"POSIXt\"")
    if (inherits(date.ref, "POSIXlt"))
        date.ref <- as.POSIXct(date.ref)
    if (!(length(date.ref) %in% c(1, length(ltraj))))
        warning("the date.ref used is not a single value or of the length of the ltraj object and will be recycled accordingly")
    ## if (!inherits(date.ref, "POSIXct"))
    ##     stop("date.ref should be of class \"POSIXct\"")
    ## End of modification
    units <- match.arg(units)
    ## When the main function is edited outside of the package,
    ## 'adehabitatLT:::.' is needed to call hidden functions
    ## dt <- .convtime(dt, units)
    ## tol <- .convtime(tol, units)
    dt <- adehabitatLT:::.convtime(dt, units)
    tol <- adehabitatLT:::.convtime(tol, units)
    ## End of modification
    if (length(date.ref) == 1)
        date.ref <- rep(date.ref, length(ltraj))
    ## The function is named to allow for finer debugging
    ## res <- lapply(1:length(ltraj), function(oo) {
    setNAburst <- function(oo) {
        ## End of modification
        x <- ltraj[[oo]]
        infol <- attr(x, "infolocs")
        date.refp <- date.ref[oo]
        dc <- x$date
        da <- as.numeric(x$da) - as.numeric(date.refp)
        glou <- round(da/dt, 0) * dt + as.numeric(date.refp)
        if (any(abs(as.numeric(dc) - as.numeric(glou)) > tol))
            stop("ltraj contains irregular data (time lag > or < tol)")
        laou <- as.integer(round(da/dt, 0))
        mlaou <- min(laou)
        laou <- laou - mlaou + 1
        xx <- rep(NA, max(laou))
        yy <- rep(NA, max(laou))
        da <- (((1:max(laou)) - 1) + mlaou) * dt + as.numeric(date.refp)
        if (!is.null(infol)) {
            infol <- do.call("data.frame", lapply(infol, function(y) {
                ll <- rep(NA, max(laou))
                ll[laou] <- y
                return(ll)
            }))
        }
        xx[laou] <- x$x
        yy[laou] <- x$y
        da[laou] <- x$date
        class(da) <- c("POSIXct", "POSIXt")
        attr(da, "tzone") <- attr(ltraj[[1]]$date, "tzone")
        return(as.ltraj(data.frame(xx, yy), da, id = attr(x,
            "id"), burst = attr(x, "burst"), typeII = TRUE, infolocs = infol,
            ...))
    }
    ## And the function is then applied:
    res <- lapply(1:length(ltraj), setNAburst)
    ## End of modification
    return(do.call("c.ltraj", res))
}
