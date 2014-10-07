## dl
##
##' @rdname ld
##' @export
dl <- function(x, strict = TRUE) {
    if (!inherits(x, "data.frame"))
        stop("x should be of class data.frame")
    ## 'strict = TRUE' corresponds to the regular 'dl' function
    if (strict) {
    ## End of modification
        nareq <- c("x", "y", "date")
        sapply(nareq, function(j) {
            if (!any(names(x) == j))
                stop(paste("No variable named", j))
        })
        if (!any(names(x) == "id")) {
            ## Does not need 'bur <-' here, 'bur' is assigned later
            ## idd <- bur <- rep(paste(sample(letters, 5), collapse = ""),
            ##     nrow(x))
            idd <- rep(paste(sample(letters, 5), collapse = ""),
                nrow(x))
            ## End of modification
        }
        else {
            idd <- x$id
        }
        if (!any(names(x) == "burst")) {
            bur <- idd
        }
        else {
            bur <- x$burst
        }
        ## No need to store 'xy' and 'dat'
        ## xy <- x[, c("x", "y")]
        ## dat <- x$date
        ## End of modification
        ## Clearer code using 'ifelse'
        ## type2 <- TRUE
        ## if (!inherits(dat, "POSIXct"))
        ##     type2 <- FALSE
        type2 <- ifelse(inherits(x$date, "POSIXct"), TRUE, FALSE)
        ## End of modification
        pasinf <- c("x", "y", "date", "dx", "dy", "dist", "dt",
            "R2n", "abs.angle", "rel.angle", "id", "burst")
        ## Condition on 'names', not on 'x'
        ## x <- x[, !(names(x) %in% pasinf), drop = FALSE]
        ## if (ncol(x) > 0) {
        ##     inf <- x
        ## }
        ## else {
        ##     inf <- NULL
        ## }
        if (all(names(x) %in% pasinf))
            inf <- NULL
        else inf <- x[, !(names(x) %in% pasinf), drop = FALSE]
        ## End of modification
        ## 'xy' and 'date' are only created now:
        ## as.ltraj(xy, dat, idd, bur, type2, infolocs = inf)
        as.ltraj(x[, c("x", "y")], x$date, idd, bur, type2, infolocs = inf)
        ## End of modification
    }
    ## 'strict = FALSE' does not use 'as.ltraj' but get all the values
    ## from the data frame
    else {
        warning("Parameters of the trajectory are not recomputed.")
        trajnam <- c("x", "y", "date", "dx", "dy", "dist", "dt",
            "R2n", "abs.angle", "rel.angle")
        ## Extract unique IDs
        idd <- tapply(as.character(x$id), x$burst, unique)
        ## Split the data frame by burst
        traj <- split(x[, names(x) %in% trajnam], x$burst)
        ## 'ltraj' names, class and attributes
        names(traj) <- NULL
        class(traj) <- c("ltraj", "list")
        attr(traj, "typeII") <- TRUE
        attr(traj, "regular") <- is.regular(traj)
        ## In case of 'infolocs' data
        if (any(!(names(x) %in% c(trajnam, "id", "burst")))) {
            ## Split the infolocs by burst
            inf <- split(x[, !(names(x) %in% c(trajnam, "id",
                "burst")), drop = FALSE], x$burst)
            ## Loop in the ltraj to add 'id', 'burst' and 'infolocs'
            for (i in (1:length(traj))) {
                attr(traj[[i]], "id") <- as.character(idd[i])
                attr(traj[[i]], "burst") <- names(idd[i])
                attr(traj[[i]], "infolocs") <- inf[[i]]
            }
        }
        ## If no infolocs, loop in the ltraj to add 'id' and 'burst'
        else for (i in (1:length(traj))) {
            attr(traj[[i]], "id") <- as.character(idd[i])
            attr(traj[[i]], "burst") <- names(idd[i])
        }
        return(traj)
    }
}


## ld
##
##' Faster versions of \code{\link[adehabitatLT]{ld}} and
##' \code{\link[adehabitatLT]{dl}}.
##'
##' In \code{\link[adehabitatLT]{ld}}, \code{strict = FALSE} can be up
##' to 10 times faster, but assumes that the \code{ltraj} is well
##' structured (i.e. not modified by the user). In
##' \code{\link[adehabitatLT]{dl}}, \code{strict = FALSE} can be up to
##' 20 times faster, but assumes that the trajectory parameters in the
##' data frame (x/y increments, angles, etc.) are still valid (e.g. no
##' locations have been removed).
##' @title Quick Conversion of Objects of Class ltraj from and to
##' Dataframes
##' @param strict Logical. Whether to use the regular
##' \code{\link[adehabitatLT]{ld}} or \code{\link[adehabitatLT]{dl}}
##' functions, which enforce more verifications.
##' @seealso See \code{\link[adehabitatLT]{ld}} for further details on
##' the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##' puechcirc ## class ltraj
##'
##' ## ld
##' df1 <- adehabitatLT:::ld(puechcirc)
##' df2 <- ld(puechcirc, strict = FALSE)
##' all.equal(df1, df2)
##' ## Note a difference in row names:
##' attr(df1, "row.names")
##' attr(df2, "row.names")
##'
##' ## dl
##' all.equal(dl(df2), adehabitatLT:::dl(df2))
##' dl(df2, strict = FALSE)
##' ## Comparison regarding 'strict'
##' all.equal(dl(df2), dl(df2, strict = FALSE))
##' ## Differences in row.names (numeric in regular 'dl', characters using
##' ## 'strict = FALSE') + NAs in R2n (for a reason, 'puechcirc[[2]]'
##' ## starts by a sequence of missing values, but has several 'R2n'
##' ## values. As a result, 'strict = FALSE' keeps the 'R2n' values)
ld <- function(ltraj, strict = TRUE)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    ## 'strict = TRUE' corresponds to the regular 'ld' function
    if (strict) {
        ## Note the use of 'adehabitatLT::id' to avoid conflicts with
        ## other packages (e.g. with 'plyr::id')
        iidd <- rep(adehabitatLT::id(ltraj), sapply(ltraj, nrow))
        bur <- rep(burst(ltraj), sapply(ltraj, nrow))
        inf <- infolocs(ltraj)
        tr <- do.call("rbind", ltraj)
        if (!is.null(inf))
            return(data.frame(tr, id = iidd, burst = bur, do.call("rbind",
                inf)))
        return(data.frame(tr, id = iidd, burst = bur))
    }
    ## 'strict = FALSE' builds the data frame without calls to 'rbind'
    else {
        inf <- infolocs(ltraj)
        df <- data.frame(
            x = unlist(lapply(ltraj, function(x) x$x)),
            y = unlist(lapply(ltraj, function(x) x$y)),
            date = unlist(lapply(ltraj, function(x) x$date)),
            dx = unlist(lapply(ltraj, function(x) x$dx)),
            dy = unlist(lapply(ltraj, function(x) x$dy)),
            dist = unlist(lapply(ltraj, function(x) x$dist)),
            dt = unlist(lapply(ltraj, function(x) x$dt)),
            R2n = unlist(lapply(ltraj, function(x) x$R2n)),
            abs.angle = unlist(lapply(ltraj, function(x) x$abs.angle)),
            rel.angle = unlist(lapply(ltraj, function(x) x$rel.angle)),
            ## Note the use of 'adehabitatLT::' to ensure the use of
            ## the 'id' function from this package, and avoid
            ## conflicts (e.g. with 'plyr::id')
            id = rep(adehabitatLT::id(ltraj), sapply(ltraj, nrow)),
            burst = rep(burst(ltraj), sapply(ltraj, nrow)))
        class(df$date) <-  c("POSIXct", "POSIXt")
        attr(df$date, "tzone") <- attr(ltraj[[1]]$date, "tzone")
        if (!is.null(inf)) {
            nc <- ncol(inf[[1]])
            infdf <- as.data.frame(matrix(nrow = nrow(df), ncol = nc))
            names(infdf) <- names(inf[[1]])
            for(i in 1:nc)
                infdf[[i]] <- unlist(lapply(inf, function(x) x[[i]]))
            df <- cbind(df, infdf)
        }
        return(df)
    }
}
