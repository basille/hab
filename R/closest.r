##' Find closest relocation from ltraj objects.
##'
##' Distances are computed between the relocations of two ltraj
##' objects within a temporal window. If only one ltraj is provided,
##' distances are computed for relocations of other individuals or
##' bursts (see parameter \code{by}).
##'
##' \code{dt} allows to define the time window using \code{Inf}
##' (e.g. \code{dt = c(-Inf, 0)} corresponds to all previous
##' locations), as well as positive values (e.g. \code{dt = c(-3600,
##' 3600)} corresponds to all locations from one hour before to one
##' hour after).
##'
##' @title Closest relocation
##' @param from A ltraj object.
##' @param to A ltraj object, supposed to be different from
##' \code{from}.
##' @param dt A numeric of length 2, giving the boundaries of the
##' temporal window. See Details.
##' @param units A character string indicating the time units for
##' \code{dt}.
##' @param prefix A character string attached to the names of the
##' variables returned.
##' @param by Character. Only if \code{to == NULL}, either \code{"id"}
##' to exclude the relocations from the same individual in the
##' computation of distances, or \code{"burst"} to exclude relocations
##' from the same burst.
##' @param range A character string indicating the range type, using a
##' representation with square and round brackets: a square bracket
##' means "inclusive" and a round bracket means "exclusive". For
##' example, \code{"[)"} includes the beginning, but not the end of
##' the interval.
##' @param protect A character string indicating other variables to
##' keep from \code{to} in \code{infolocs}.
##' @note The function assumes that both ltraj are projected, and in
##' the same projection.
##' @return The same ltraj as \code{from}, with additional variables
##' giving x and y coordinates, the date, id and burst of the closest
##' relocation, together with the distance to this relocation.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @examples
##' data(puechcirc)
##'
##' toto <- closest(puechcirc)
##' tail(infolocs(toto)[[3]])
##'
##' tata <- closest(puechcirc, dt = c(-3, 2), units = "hour", by = "burst")
##' tail(infolocs(tata)[[3]])
##'
##' titi <- closest(puechcirc, dt = c(-5, -2), units = "hour", protect = "abs.angle")
##' tail(infolocs(titi)[[3]])
##' @export
closest <- function(from, to = NULL, dt = c(-3600, 0), units = c("sec",
    "min", "hour", "day"), prefix = "to.", by = c("id", "burst"),
    range = c("[]", "[)", "(]", "()"), protect = NULL)
{
    ## Check objects and arguments
    if (!inherits(from, "ltraj"))
        stop("from should be of class ltraj")
    same <- FALSE
    if (is.null(to)) {
        same <- TRUE
        to <- from
    }
    if (!inherits(to, "ltraj"))
        stop("to should be of class ltraj")
    units <- match.arg(units)
    if (!(is.numeric(dt) & length(dt) == 2))
        stop("dt should be an integer of length 2")
    if (!(dt[2] >= dt[1]))
        stop("dt should defined an interval equal or greater than zero (i.e. dt[2] >= dt[1])")
    dt <- adehabitatLT:::.convtime(dt, units)
    by <- match.arg(by)
    range <- match.arg(range)

    ## Convert ltraj to data frames
    from <- ld(from, strict = FALSE)
    to <- ld(to, strict = FALSE)

    ## Function for one location
    distLoc <- function(loc, to, dt) {
        ## If `loc` is a missing loc, return a 1-row df filled with
        ## NAs
        if (is.na(loc[["x"]])) {
            tmp <- subset(to, FALSE, select = c("x", "y", "date",
                "id", "burst", protect))
            tmp[1, ] <- NA
            tmp$distloc <- NA
            return(tmp)
        }
        ## Logical for locs of `to` in the temporal window
        if (range == "[]")
            cond <- to[["date"]] >= loc[["date"]] + dt[1] &
                    to[["date"]] <= loc[["date"]] + dt[2]
        if (range == "[)")
            cond <- to[["date"]] >= loc[["date"]] + dt[1] &
                    to[["date"]] <  loc[["date"]] + dt[2]
        if (range == "(]")
            cond <- to[["date"]] >  loc[["date"]] + dt[1] &
                    to[["date"]] <= loc[["date"]] + dt[2]
        if (range == "()")
            cond <- to[["date"]] >  loc[["date"]] + dt[1] &
                    to[["date"]] <  loc[["date"]] + dt[2]
        ## Logical for locs of `to` with different id/burst
        if (same)
            condid <- to[[by]] != loc[[by]]
        else condid <- TRUE
        ## Subset every loc from `to` in the temporal window
        tmp <- subset(to, cond & condid, select = c("x", "y",
            "date", "id", "burst", protect))
        ## If at least one loc is selected, compute the distance
        ## between loci and all locs
        if (nrow(tmp) > 0) {
            tmp$distloc <- sqrt((loc[["x"]] - tmp[["x"]])^2 +
                                (loc[["y"]] - tmp[["y"]])^2)
            ## If all distances are NAs, return a 1-row df filled with
            ## NAs
            if (all(is.na(tmp$distloc))) {
                tmp <- tmp[1, ]
                tmp[1, ] <- NA
                return(tmp)
            }
            ## Otherwise return the row with minimum distance
            else return(tmp[which.min(tmp$distloc), ])
        }
        ## Otherwise return a 1-row df filled with NAs
        else {
            tmp[1, ] <- NA
            tmp$distloc <- NA
            return(tmp)
        }
    }

    ## Compute the closest distance for each row of `from`
    dists <- do.call(rbind, lapply(1:nrow(from), function(i)
        distLoc(loc = from[i, ], to = to, dt = dt)))
    ## Add prefix to colum names
    names(dists) <- paste0(prefix, names(dists))
    ## Match row names
    row.names(dists) <- row.names(from)
    ## Bind and convert back to ltraj
    suppressWarnings(from <- dl(cbind(from, dists), strict = FALSE))
    return(from)
}
