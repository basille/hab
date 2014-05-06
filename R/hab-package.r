##' A set of functions related to habitat selection and movement
##' analyses. Also includes a few patches for functions from the
##' \code{adehabitatXY} package family. For a list of documented
##' functions, use \code{library(help = "hab")}
##'
##' @title Habitat selection and movement analyses
##' @docType package
##' @name hab
##' @author Mathieu Basille \email{basille@@ase-research.org}
NULL


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


## as.ltraj
##
##' Faster version of \code{\link[adehabitatLT]{as.ltraj}}, which can
##' be up to 5-10 times faster.
##'
##' @title Working with Trajectories in 2D Space: the Class ltraj
##' @seealso See \code{\link[adehabitatLT]{ld}} for further details on
##' the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechabonsp)
##' locs <- puechabonsp$relocs
##' xy <- coordinates(locs)
##' df <- as.data.frame(locs)
##' id <- df[,1]
##' da <- as.POSIXct(strptime(as.character(df$Date), "%y%m%d"))
##' ltr1 <- adehabitatLT:::as.ltraj(xy, da, id = id)
##' ltr2 <- as.ltraj(xy, da, id = id)
##' all.equal(ltr1, ltr2)
as.ltraj <- function (xy, date = NULL, id, burst = id, typeII = TRUE, slsp = c("remove",
    "missing"), infolocs = data.frame(pkey = paste(id, date,
    sep = ".")))
{
    if (typeII) {
        if (!inherits(date, "POSIXct"))
            stop("For objects of type II,\n date should be of class \"POSIXct\"")
    }
    else {
        date <- 1:nrow(xy)
    }
    if (length(date) != nrow(xy))
        stop("date should be of the same length as xy")
    slsp <- match.arg(slsp)
    if (!is.null(infolocs)) {
        if (nrow(infolocs) != nrow(xy))
            stop("infolocs should have the same number of rows as xy")
    }
    if (length(id) == 1)
        id <- factor(rep(as.character(id), nrow(xy)))
    if (length(id) != nrow(xy))
        stop("id should be of the same length as xy, or of length 1")
    if (min(table(id)) == 0)
        stop("some id's are not present in the data")
    if (length(burst) == 1)
        burst <- factor(rep(as.character(burst), nrow(xy)))
    if (length(burst) != nrow(xy))
        stop("burst should be of the same length as xy, or of length 1")
    if (min(table(burst)) == 0)
        stop("some bursts are not present in the data")
    id1 <- factor(id)
    burst1 <- factor(burst)
    if (!all(apply(table(id1, burst1) > 0, 2, sum) == 1))
        stop("one burst level should belong to only one id level")
    x <- xy[, 1]
    y <- xy[, 2]
    res <- split(data.frame(x = x, y = y, date = date), burst)
    ## Test of duplicated dates now:
    ## Condition tested on the fly and not stored + use of
    ## 'duplicated'
    ## rr <- any(unlist(lapply(res, function(x) (length(unique(x$date)) !=
    ##     length(x$date)))))
    ## if (any(unlist(lapply(res, function(x) (length(unique(x$date)) !=
    ##     length(x$date))))))
    if (any(unlist(lapply(res, function(x) any(duplicated(x$date))))))
        stop("non unique dates for a given burst")
    ## x/y not needed in resbb
    ## resbb <- split(data.frame(x = x, y = y, date = date), id1)
    resbb <- split(date, id1)
    ## Condition tested on the fly and not stored + use of
    ## 'duplicated'
    ## rr <- any(unlist(lapply(resbb, function(x) (length(unique(x$date)) !=
    ##     length(x$date)))))
    ## if (any(unlist(lapply(resbb, function(x) (length(unique(x)) != length(x))))))
    if (any(unlist(lapply(resbb, function(x) any(duplicated(x))))))
        stop("non unique dates for a given id")
    ## EOM
    ## Not needed now, see below
    ## liid <- split(id, burst)
    ## EOM
    ## Check the condition only once
    ## if (!is.null(infolocs))
    if (!is.null(infolocs)) {
        linfol <- split(infolocs, burst)
    ## if (!is.null(infolocs))
    ## EOM
        linfol <- lapply(1:length(linfol), function(j) linfol[[j]][order(res[[j]]$date),
            , drop = FALSE])
    }
    res <- lapply(res, function(y) y[order(y$date), , drop = FALSE])
    ## Duplication of code
    ## x <- xy[, 1]
    ## y <- xy[, 2]
    ## EOM
    foo <- function(x) {
        x1 <- x[-1, ]
        x2 <- x[-nrow(x), ]
        dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2), NA)
        R2n <- (x$x - x$x[1])^2 + (x$y - x$y[1])^2
        dt <- c(unclass(x1$date) - unclass(x2$date), NA)
        dx <- c(x1$x - x2$x, NA)
        dy <- c(x1$y - x2$y, NA)
        abs.angle <- ifelse(dist < 1e-07, NA, atan2(dy, dx))
        so <- cbind.data.frame(dx = dx, dy = dy, dist = dist,
            dt = dt, R2n = R2n, abs.angle = abs.angle)
        return(so)
    }
    speed <- lapply(res, foo)
    res <- lapply(1:length(res), function(i) cbind(res[[i]],
        speed[[i]]))
    ang.rel <- function(df, slspi = slsp) {
        ang1 <- df$abs.angle[-nrow(df)]
        ang2 <- df$abs.angle[-1]
        if (slspi == "remove") {
            dist <- c(sqrt((df[-nrow(df), "x"] - df[-1, "x"])^2 +
                (df[-nrow(df), "y"] - df[-1, "y"])^2), NA)
            wh.na <- which(dist < 1e-07)
            if (length(wh.na) > 0) {
                no.na <- (1:length(ang1))[!(1:length(ang1)) %in%
                  wh.na]
                for (i in wh.na) {
                  indx <- no.na[no.na < i]
                  ang1[i] <- ifelse(length(indx) == 0, NA, ang1[max(indx)])
                }
            }
        }
        res <- ang2 - ang1
        res <- ifelse(res <= (-pi), 2 * pi + res, res)
        res <- ifelse(res > pi, res - 2 * pi, res)
        return(c(NA, res))
    }
    rel.angle <- lapply(res, ang.rel)
    res <- lapply(1:length(res), function(i) data.frame(res[[i]],
        rel.angle = rel.angle[[i]]))
    ## Assign id/burst attributes to each burst
    ## res <- lapply(1:length(res), function(i) {
    ##     x <- res[[i]]
    ##     attr(x, "id") <- as.character(liid[[i]][1])
    ##     attr(x, "burst") <- levels(factor(burst))[i]
    ##     return(x)
    ## })
    ## Create a list of c(idX, burstX)
    lidburst <- tapply(as.character(id1), burst1, function(x) x[1])
    lidburst <- mapply(c, as.list(lidburst), as.list(names(lidburst)),
        SIMPLIFY = FALSE)
    ## Function that assigns id/burst attributes
    setidburst <- function(x, y) {
        attr(x, "id") <- y[1]
        attr(x, "burst") <- y[2]
        x
    }
    ## Using 'mapply' instead (big speed enhancement)
    res <- mapply(setidburst, res, lidburst, SIMPLIFY = FALSE)
    ## EOM
    if (!is.null(infolocs)) {
        res <- lapply(1:length(res), function(i) {
            x <- res[[i]]
            y <- linfol[[i]]
            row.names(y) <- row.names(x)
            attr(x, "infolocs") <- y
            return(x)
        })
    }
    class(res) <- c("ltraj", "list")
    attr(res, "typeII") <- typeII
    attr(res, "regular") <- is.regular(res)
    return(res)
}


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


## infolocs
##
##' Modified version of \code{\link[adehabitatLT]{infolocs}} that
##' returns \code{NULL} if \code{infolocs} exists, but \code{which} is
##' not a colum of it.
##'
##' @title Infolocs of an Object of Class ltraj
##' @seealso See \code{\link[adehabitatLT]{infolocs}} for further
##' details on the function and all available arguments.
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
infolocs <- function(ltraj, which) {
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    if (!is.null(attr(ltraj[[1]], "infolocs"))) {
        if (missing(which))
            which <- names(attr(ltraj[[1]], "infolocs"))
        ## If infolocs exists but which is not in it, returns NULL
        else if (!(which %in% names(attr(ltraj[[1]], "infolocs"))))
            return(NULL)
        re <- lapply(ltraj, function(y) {
            res <- attr(y, "infolocs")
            return(res[, names(res) %in% which, drop = FALSE])
        })
        return(re)
    }
    else {
        return(NULL)
    }
}


## kerneloverlap
##
##' Modified version of \code{\link[adehabitatHR]{kerneloverlap}},
##' which now allows to select several methods at once, and works on
##' \code{SpatialPointsDataFrame} or \code{estUDm} in the same
##' function.
##'
##' @title Spatial Interaction between Animals Monitored Using Radio-Tracking
##'
##' @param x an object of class \code{SpatialPointsDataFrame}
##' containing only one column (which is a factor indicating the
##' identity associated to the relocations) or an object of class
##' \code{estUDm} containing several home-ranges for which the overlap
##' is to be calculated
##' @param method the desired method(s)for the estimation of overlap
##' (several methods can be chosen, defaults to all at once)
##' @return An array of class \code{kerneloverlap} with as many
##' matrices as chosen methods (if only one method is chosen, a matrix
##' of class \code{kerneloverlap}); each matrix gives the value of
##' indices of overlap for all pairs of animals.
##' @seealso See \code{\link[adehabitatHR]{kerneloverlap}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Prepare the data:
##' data(puechabonsp)
##' loc <- puechabonsp$relocs
##' elev <- puechabonsp$map
##'
##' ## Kernel overlap using various approaches:
##' kerneloverlap(puechabonsp$relocs[,1], grid = 200, meth = "HR", conditional = TRUE)
##' kerneloverlap(puechabonsp$relocs[,1], grid = elev, meth = c("HR", "VI"), conditional = TRUE)
##' kerneloverlap(puechabonsp$relocs[,1], grid = 200, conditional = TRUE)
##'
##' ## Summarizing the results:
##' summary(kerneloverlap(puechabonsp$relocs[,1], grid = 200))
##' summary(kerneloverlap(puechabonsp$relocs[,1], grid = 200, percent = 50))
kerneloverlap <- function(x, method = c("HR", "PHR", "VI", "BA",
    "UDOI", "HD"), percent = 95, conditional = FALSE, ...)
{
    ## Verifications
    ## method <- match.arg(method)
    method <- match.arg(method, several.ok = TRUE)
    ## x must inherits SPDF or estUDm
    if (!inherits(x, c("SpatialPoints", "estUDm")))
        stop("x should inherit the class SpatialPoints, or be of class estUDm")
    ## if (ncol(coordinates(xy))>2)
    ##     stop("xy should be defined in two dimensions")
    if (percent > 100)
        stop("percent should not be >100")
    ## In case of a SPDF, get the kernelUD estimate
    if (inherits(x, c("SpatialPoints")))
        x <- kernelUD(x, same4all = TRUE, ...)
    ## Other checks...
    if (length(x) == 1)
        stop("several animals are needed for this function")
    if (slot(x[[1]], "vol"))
        stop("x should not be a volume under UD")
    ## Check the grid (must be the same for all)
    for (i in 2:length(x)) if (!identical(attr(x[[i]], "grid"),
        attr(x[[1]], "grid")))
        stop("kernelUD must use the same grid for all animals")
    ## Back to the function
    vol <- getvolumeUD(x)
    ## Checks that the UD is distributed over several pixels
    tmp <- lapply(1:length(vol), function(i) {
        y <- vol[[i]]
        su <- sum(y[[1]] < 95)
        if (su < 5)
            warning(paste("For animal ", names(vol)[i], ", the most of the UD is distributed ",
                "in less than 5 pixels.\n The results will probably be inconsistent.\n",
                "Try to increase the parameter grid.\n"))
    })
    x <- lapply(x, function(y) {
        coo <- coordinates(y)
        y[order(coo[, 1], coo[, 2]), ]
    })
    vol <- lapply(vol, function(y) {
        coo <- coordinates(y)
        y[order(coo[, 1], coo[, 2]), ]
    })
    gp <- gridparameters(vol[[1]])
    ## Array of results
    ## res <- matrix(0, ncol=length(x), nrow=length(x))
    res <- array(0, dim = c(length(x), length(x), length(method)),
        dimnames = list(names(x), names(x), method))
    ## Loop for each animal
    for (i in 1:length(x)) {
        for (j in 1:i) {
            vi <- x[[i]][[1]]
            vj <- x[[j]][[1]]
            ai <- vol[[i]][[1]]
            aj <- vol[[j]][[1]]
            ai[ai <= percent] <- 1
            ai[ai > percent] <- 0
            aj[aj <= percent] <- 1
            aj[aj > percent] <- 0
            if (conditional) {
                vi <- vi * ai
                vj <- vj * aj
            }
            ## if (method=="HR") {
            if ("HR" %in% method) {
                ak <- ai * aj
                ## res[i, j] <- sum(vk) / sum(vi)
                ## res[j, i] <- sum(vk) / sum(vj)
                res[i, j, "HR"] <- sum(ak)/sum(ai)
                res[j, i, "HR"] <- sum(ak)/sum(aj)
            }
            ## if (method=="PHR") {
            if ("PHR" %in% method) {
                ## res[j,i] <- sum(vi*aj)*(gp[1,2]^2)
                ## res[i,j] <- sum(vj*ai)*(gp[1,2]^2)
                res[j, i, "PHR"] <- sum(vi * aj) * (gp[1, 2]^2)
                res[i, j, "PHR"] <- sum(vj * ai) * (gp[1, 2]^2)
            }
            ## if (method=="VI") {
            if ("VI" %in% method) {
                ## res[i,j] <- res[j,i] <- sum(pmin(vi, vj))*(gp[1,2]^2)
                res[i, j, "VI"] <- res[j, i, "VI"] <- sum(pmin(vi,
                  vj)) * (gp[1, 2]^2)
            }
            ## if (method=="BA") {
            if ("BA" %in% method) {
                ## res[j,i] <- res[i,j] <- sum(sqrt(vi)*sqrt(vj))*(gp[1,2]^2)
                res[j, i, "BA"] <- res[i, j, "BA"] <- sum(sqrt(vi) *
                  sqrt(vj)) * (gp[1, 2]^2)
            }
            ## if (method=="UDOI") {
            if ("UDOI" %in% method) {
                ak <- sum(ai * aj) * (gp[1, 2]^2)
                ## res[j,i] <- res[i,j] <- ak * sum(vi*vj)*(gp[1,2]^2)
                res[j, i, "UDOI"] <- res[i, j, "UDOI"] <- ak *
                  sum(vi * vj) * (gp[1, 2]^2)
            }
            ## if (method=="HD") {
            if ("HD" %in% method) {
                ## res[j,i] <- res[i,j] <- sqrt(sum((sqrt(vi) - sqrt(vj))^2*(gp[1,2]^2)))
                res[j, i, "HD"] <- res[i, j, "HD"] <- sqrt(sum((sqrt(vi) -
                  sqrt(vj))^2 * (gp[1, 2]^2)))
            }
        }
    }
    ## rownames(res) <- names(x)
    ## colnames(res) <- names(x)
    ## Simplification to matrix if only one method
    if (length(method) == 1)
        res <- res[, , 1]
    attr(res, "method") <- method
    class(res) <- "kerneloverlap"
    return(res)
}
##' @rdname kerneloverlap
##' @param object An array or matrix of class \code{kerneloverlap}
##' @S3method summary kerneloverlap
##' @export
summary.kerneloverlap <- function(object, ...)
{
    ## Check the class of object
    if (!inherits(object, "kerneloverlap"))
        stop("object should be of class kerneloverlap")
    ## Store the method(s) used
    method <- attr(object, "method")
    ## If object is a matrix (1 method only), make it an array
    if (length(method) == 1)
        object <- array(object, dim = c(4, 4, 1), dimnames = list(rownames(object),
            colnames(object), method))
    ## Prepare a list for the output
    ll <- list()
    ## If HR/PHR, results are not symmetrical: need to take both sides
    ## into account
    if ("HR" %in% method) {
        mat <- object[, , "HR"]
        ll$HR12 <- t(mat)[lower.tri(t(mat))]
        ll$HR21 <- mat[lower.tri(mat)]
    }
    if ("PHR" %in% method) {
        mat <- object[, , "PHR"]
        ll$PHR12 <- t(mat)[lower.tri(t(mat))]
        ll$PHR21 <- mat[lower.tri(mat)]
    }
    ## If VI/BA/UDIO/HD, results are symmetrical
    if ("VI" %in% method) {
        mat <- object[, , "VI"]
        ll$VI <- mat[lower.tri(mat)]
    }
    if ("BA" %in% method) {
        mat <- object[, , "BA"]
        ll$BA <- mat[lower.tri(mat)]
    }
    if ("UDOI" %in% method) {
        mat <- object[, , "UDOI"]
        ll$UDOI <- mat[lower.tri(mat)]
    }
    if ("HD" %in% method) {
        mat <- object[, , "HD"]
        ll$HD <- mat[lower.tri(mat)]
    }
    ## Make a data.frame from the list, with the combinations as row
    ## names
    data.frame(ll, row.names = apply(combn(dimnames(object)[[1]],
        2), 2, paste, collapse = "-"))
}


## kernelUD
##
##' Modified version of \code{\link[adehabitatHR]{kernelUD}}, which
##' silences the use of \code{same4all} when working on a
##' \code{SpatialPixels} as a \code{grid}.
##'
##' @title Estimation of Kernel Home-Range
##' @seealso See \code{\link[adehabitatHR]{kernelUD}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Prepare the data:
##' data(puechabonsp)
##' loc <- puechabonsp$relocs
##' elev <- puechabonsp$map
##'
##' ## Compute the kernel estimation
##' ker1 <- kernelUD(puechabonsp$relocs[,1], grid = elev, same4all = TRUE)
##'
##' ## Summarizing the overlaps
##' summary(kerneloverlap(ker1, conditional = TRUE))
kernelUD <- function (xy, h = "href", grid = 60, same4all = FALSE, hlim = c(0.1,
    1.5), kern = c("bivnorm", "epa"), extent = 0.5, boundary = NULL)
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class SpatialPoints")
    if (ncol(coordinates(xy)) > 2)
        stop("xy should be defined in two dimensions")
    pfs1 <- proj4string(xy)
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy) != 1) {
            warning("xy should contain only one column (the id of the animals)\nid ignored")
            id <- rep("a", nrow(as.data.frame(xy)))
            m <- 2
        }
        else {
            id <- xy[[1]]
            m <- 1
        }
    }
    else {
        id <- rep("a", nrow(as.data.frame(xy)))
        m <- 2
    }
    if (!is.null(boundary)) {
        if (!inherits(boundary, "SpatialLines"))
            stop("the boundary should be an object of class SpatialLines")
    }
    if (min(table(id)) < 5)
        stop("At least 5 relocations are required to fit an home range")
    id <- factor(id)
    xy <- as.data.frame(coordinates(xy))
    ## if (same4all) {
    ##     if (inherits(grid, "SpatialPoints"))
    ##         stop("when same4all is TRUE, grid should be a number")
    ##     grid <- .makegridUD(xy, grid, extent)
    ## }
    if (!inherits(grid, "SpatialPixels"))
        if (same4all)
            ## When the main function is edited outside of the package,
            ## 'adehabitatHR:::.' is needed to call hidden functions
            grid <- adehabitatHR:::.makegridUD(xy, grid, extent)
    lixy <- split(xy, id)
    res <- lapply(1:length(lixy), function(i) {
        if (is.list(grid)) {
            grida <- grid[names(lixy)[i]][[1]]
        }
        else {
            grida <- grid
        }
        x <- lixy[names(lixy)[i]][[1]]
        if (!is.null(boundary)) {
            ## When the main function is edited outside of the package,
            ## 'adehabitatHR:::.' is needed to call hidden functions
            bdrk <- adehabitatHR:::.boundaryk(SpatialPoints(x, proj4string = CRS(as.character(pfs1))),
                boundary, h)
            if (!is.null(bdrk)) {
                sigg <- attr(bdrk, "sign")
                bdrk <- as.data.frame(coordinates(bdrk))
                names(bdrk) <- names(x)
                x <- rbind(x, bdrk)
            }
        }
        ## When the main function is edited outside of the package,
        ## 'adehabitatHR:::.' is needed to call hidden functions
        ud <- adehabitatHR:::.kernelUDs(SpatialPoints(x, proj4string = CRS(as.character(pfs1))),
            h = h, grid = grida, hlim = hlim, kern = kern, extent = extent)
        if (!is.null(boundary)) {
            if (!is.null(bdrk)) {
                ## When the main function is edited outside of the package,
                ## 'adehabitatHR:::.' is needed to call hidden functions
                ud <- adehabitatHR:::.fbboun(ud, boundary, sigg, h)
                slot(ud, "data")[, 1] <- slot(ud, "data")[, 1]/(sum(slot(ud,
                  "data")[, 1]) * gridparameters(ud)[1, 2]^2)
            }
        }
        if (!is.na(pfs1))
            proj4string(ud) <- CRS(pfs1)
        return(ud)
    })
    names(res) <- names(lixy)
    class(res) <- "estUDm"
    if (m == 2) {
        res <- res[[1]]
    }
    return(res)
}


## kfold
##
##' Cross-validation for regression models.
##'
##' Note: needs complete names in the coxph/clogit call, and not
##' 'cos(var)' or 'I(var*10)', except for 'strata()' and
##' 'cluster()'.
##'
##' Also needs complete case in the model (e.g. it fails if
##' there is at least one NA in the observed steps for any variable of
##' the model). Returns an error if some stratas have no case.
##' @title kfold
##' @param mod A fitted model for which there exists a \code{kfold}
##' method (currently \code{coxph} and \code{clogit} models).
##' @param k The number of equal size subsamples of the partition.
##' @param nrepet The number of repetitions.
##' @param jitter Logical, whether to add some random noise to the
##' predictions (useful when the model is fitted on categorical
##' variables, which can produces error in the ranking process).
##' @param reproducible Logical, whether to use a fixed seed for each
##' repetition.
##' @param details Logical, whether to return details of each
##' repetition (useful for debugging).
##' @return A data frame with the correlations (\code{cor}) and the
##' type of value (\code{type}).
##' @author Mathieu Basille \email{basille@@ase-research.org}, with
##' the help of Terry Therneau and Guillaume Bastille-Rousseau
##' @export
##' @S3method kfold coxph
kfold <- function(mod, k = 5, nrepet = 100, jitter = FALSE,
    reproducible = TRUE, details = FALSE)
    UseMethod("kfold")
## kfold.glm <- function(mod, k = 5, nrepet = 100, nbins = 10, jitter = FALSE,
##     random = TRUE, reproducible = TRUE)
## {
## kfold.mer <- function(mod, k = 5, nrepet = 100, nbins = 10, jitter = FALSE,
##     random = TRUE, reproducible = TRUE)
## {
##     if (!inherits(mod, c("glm", "mer")))
##         stop("Model of class '(g)lm' or '(g)lmer' expected")
##     if (inherits(mod, c("mer")))
##         require(lme4)
##     dt <- model.frame(mod)
##     ## dt <- mod$data
##     kfold <- rd <- numeric(length = nrepet)
##     resp <- as.character(attr(terms(mod), "variables"))[attr(terms(mod),
##         "response") + 1]
##     ## resp <- all.vars(terms(mod))[attr(terms(mod), "response")]
##     for (i in 1:nrepet) {
##         dt$sets <- "train"
##         if(reproducible)
##             set.seed(i)
##         dt$sets[sample(which(dt[, resp] == 1), sum(dt[, resp] ==
##             1)/k)] <- "test"
##         reg <- update(mod, data = subset(dt, sets == "train"))
##         if (inherits(mod, "glm"))
##             ## predall <- exp(predict(reg, newdata = dt)) ## does the same
##             predall <- exp(as.numeric(model.matrix(terms(reg),
##                 dt) %*% coef(reg)))
##         else if (inherits(mod, "mer"))
##             predall <- exp(as.numeric(model.matrix(terms(reg),
##                 dt) %*% lme4::fixef(reg)))
##         if (jitter) {
##             if(reproducible)
##                 set.seed(i)
##             predall <- jitter(predall)
##         }
##         quant <- quantile(predall[dt[, resp] == 0], probs = seq(from = 0,
##             to = 1, length.out = nbins + 1))
##         quant[1] <- 0
##         quant[length(quant)] <- Inf
##         int <- factor(findInterval(predall[dt$sets == "test"],
##             quant), levels = 1:nbins)
##         kfold[i] <- cor(1:nbins, table(int), method = "spearman")
##         if (random) {
##             if (reproducible)
##                 set.seed(i)
##             dt$sets[sample(which(dt[, resp] == 0), sum(dt[, resp] ==
##                 1)/k)] <- "rd"
##             int <- factor(findInterval(predall[dt$sets == "rd"],
##                 quant), levels = 1:nbins)
##             rd[i] <- cor(1:nbins, table(int), method = "spearman")
##         }
##     }
##     if (random)
##         return(data.frame(kfold = c(kfold, rd), type = rep(c("obs",
##             "rand"), each = nrepet)))
##     else return(kfold)
## }
##
## debug(kfoldRSF)
## (bli <- kfoldRSF(reg1, k = 4, nrepet = 4, jitter = TRUE))
## (bla <- kfoldRSF(reg1, k = 4, nrepet = 4, random = TRUE))
## by(bla$kfold, bla$type, summary)
## t.test(bla$kfold ~ bla$type)
## boxplot(bla$kfold ~ bla$type)
## abline(h = 0)
##
## Error with dummy variables (use 'jitter = TRUE'):
## Error in findInterval(predall[dt$sets == "test"], quant) :
##   'vec' must be sorted non-decreasingly
kfold.coxph <- function(mod, k = 5, nrepet = 100, jitter = FALSE,
    reproducible = TRUE, details = FALSE)
{
    ## Check the class of the model (should be "coxph" or "clogit")
    if (!inherits(mod, "coxph"))
        stop("Model of class 'coxph' expected.")
    ## Load survival
    require(survival)
    ## Try to retrieve the data
    dt <- try(model.frame(mod), silent = TRUE)
    ## If it failed, stop and give a solution
    if (class(dt) == "try-error")
        stop("'model.frame' was unable to retrieve the data.",
          "Use 'model = TRUE' in the 'coxph' or 'clogit' call.")
    ## The first column is named 'srv' instead of 'Surv(faketime,
    ## case)'
    names(dt)[1] <- "srv"
    ## Which column is the strata?
    nstr <- attr(terms(mod), "specials")$strata
    ## Ugly regexp to extract and apply the strata variable name
    names(dt)[nstr] <- namestr <- sub("strata\\((.*)\\)", "\\1",
        names(dt)[nstr])
    ## If there is a cluster term...
    if (!is.null(attr(terms(mod), "specials")$cluster)) {
        ## Which column is the cluster?
        nclu <- attr(terms(mod), "specials")$cluster
        ## Ugly regexp to extract and apply the cluster variable name
        names(dt)[nclu] <- sub("cluster\\((.*)\\)", "\\1",
            names(dt)[nclu])
    }
    ## Is it really a problem?
    ## ncase <- table(tapply(dt$srv[, 2], dt[, nstr], function(x) sum(x == 1)))
    ## if (any(names(ncase) == "0"))
    ##     stop(paste("Some stratas had no case.",
    ##       "It is likely that NAs were present in the variables for some cases."))
    ## Prepare the 'kfold', 'rd' and 'warn' objects
    kfold <- rd <- warn <- numeric(length = nrepet)
    ## 'dbg' object for debugging when 'details = TRUE'
    if (details)
        dbg <- list()
    ## The core of the kfold, each repetition
    for (i in 1:nrepet) {
        ## Create a 'set' column, which defaults to "train"
        dt$sets <- "train"
        ## Allows for reproducibility
        if (reproducible)
            set.seed(i)
        ## Sample the "test" data set
        dt$sets[dt[, namestr] %in% sample(unique(dt[, namestr]),
            length(unique(dt[, namestr]))/k)] <- "test"
        ## Update the regression using the training data
        reg <- update(mod, srv ~ ., data = subset(dt, sets ==
            "train"), model = TRUE)
        ## Extract the "test" data set
        dtest <- droplevels(subset(dt, sets == "test"))
        ## And compute the predictions associated to this data set
        ## using the training regression
        dtest$predall <- exp(predict(reg, type = "lp", newdata = dtest,
            reference = "sample"))
        ## In case of equality among predictions (e.g. categorical
        ## variable), add some noise to the predictions
        if (jitter) {
            ## Allows for reproducibility
            if (reproducible)
                set.seed(i)
            dtest$predall <- jitter(dtest$predall)
        }
        ## The function to compute the rank within a strata
        samplepred <- function(df) {
            ## Number of controls
            nrand <- sum(df$srv[, 2] == 0)
            ## Rank of the case (among case + controls)
            obs <- rank(df$predall)[df$srv[, 2] == 1]
            ## Rank of a random control (among controls only!)
            if (reproducible)
                set.seed(i)
            rand <- sample(rank(df$predall[df$srv[, 2] == 0]),
                1)
            return(data.frame(obs = obs, rand = rand, nrand = nrand))
        }
        ## Compute the ranks for each strata and bind them together
        ranks <- do.call(rbind, by(dtest, dtest[, namestr], samplepred))
        ## Is there the same number of controls per strata?
        nrand <- unique(ranks$nrand)
        ## If no, use the greatest number of controls (and keep track
        ## of it)
        if (length(nrand) != 1) {
            nrand <- max(nrand)
            warn[i] <- 1
        }
        ## Compute the Spearman correlation on the ranks for the cases
        kfold[i] <- cor(1:(nrand+1), table(factor(ranks$obs,
            levels = 1:(nrand+1))), method = "spearman")
        ## Same for the random controls
        rd[i] <- cor(1:(nrand), table(factor(ranks$rand, levels = 1:(nrand))),
            method = "spearman")
        ## Store the ranks for debugging if 'details'
        if (details)
            dbg[[i]] <- ranks
    }
    ## Create a data frame with the correlations and the type of value
    res <- data.frame(cor = c(kfold, rd), type = rep(c("obs",
        "rand"), each = nrepet))
    ## Return the debugging information as an attribute if 'details'
    if (details)
        attr(res, "details") <- dbg
    ## If the number of controls is not the same for each strata, send
    ## a warning
    if (sum(warn) > 0)
        warning(paste("The number of controls was not equal among stratas for",
          sum(warn), "repetitions. Correlations might be biased.",
          "Use 'details = TRUE' to get more details."))
    ## Return the result data frame
    return(res)
}
## kfold.coxme <- function(mod) {
##     ## else if (inherits(mod, "coxme"))
##     ##     predall <- exp(as.numeric(model.matrix(terms(reg),
##     ##         dt) %*% coxme::fixef(reg)))
## }
##
## debug(kfoldSSF)
## (bla <- kfoldSSF(mod, k = 4, nrepet = 4))
## by(bla$kfold, bla$type, summary)
## by(bla$kfold, bla$type, sd)
## t.test(bla$kfold ~ bla$type)
## boxplot(bla$kfold ~ bla$type, names = c("Observed", "Random"), main = "K-fold cross-validation", ylab = "Spearman correlation")
## abline(h = 0)
##
## survival:::predict.coxph
## Pour le subset : https://stat.ethz.ch/pipermail/r-help/2011-December/298434.html


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
        ## Note the use of 'adehabitatLT::' to ensure the use of
        ## the 'id' function from this package, and avoid
        ## conflicts (e.g. with 'plyr::id')
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


## lincircor.test
##
##' Compute Linear-circular correlation, and test it.
##'
##' @title Linear-circular correlation
##' @param xl A numeric, providing the linear variable.
##' @param theta A numeric, providing the circular variable, in
##' radians (\code{rad = deg * pi / 180}).
##' @return \code{lincircor} returns a numeric providing the
##' linear-circular correlation. \code{lincircor.test} returns a list
##' with class \code{"htest"} containing the following components:
##' \itemize{
##' \item \code{statistic}: the value of the test statistic
##' \item \code{parameter}: the degrees of freedom of the test
##' statistic, which follows a F distribution if \code{x}, and
##' \code{theta} are independent and \code{x} is normally distributed
##' \item \code{p.value}: the p-value of the test
##' \item \code{estimate}: the estimated measure of linear-circular
##' correlation
##' \item \code{null.value}: the value of the correlation measure
##' under the null hypothesis, always \code{0}
##' \item \code{alternative}: a character string describing the
##' alternative hypothesis
##' \item \code{method}: a character string indicating how the
##' correlation was measured
##' \item \code{data.name}: a character string giving the names of the
##' data}
##' @references Mardia, K. V. & Jupp, P. E. (2000) Directional
##' statistics. Wiley, 429 pp.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## From Mardia & Jupp (2000), p.246
##' xl <- c(28, 85.2, 80.5, 4.7, 45.9, 12.7, 72.5, 56.6, 31.5, 112,
##'     20, 72.5, 16, 45.9, 32.6, 56.6, 52.6, 91.8, 55.2)
##' theta <- c(327, 91, 88, 305, 344, 270, 67, 21, 281, 8, 204, 86,
##'     333, 18, 57, 6, 11, 27, 84)
##' ## 'xl' is not in radians (but the test is computed anyway)
##' lincircor(xl, theta)
##' ## Conversion to radians
##' lincircor(xl, theta*pi/180)
##' ## Test
##' lincircor.test(xl, theta*pi/180)
lincircor <- function(xl, theta) {
    ## Check the lengths of 'xl' and 'theta'
    if (length(xl) != length(theta))
        stop("x and theta should be vector of same length.")
    ## Quick check of 'theta' in radians (but calculations are done
    ## anyway)
    if (max(theta, na.rm = TRUE) > 2 * pi | min(theta, na.rm = TRUE) <
        -2 * pi)
        warning("theta should be in radians. Calculations done anyway.")
    ## Remove NAs
    nna <- is.na(xl) | is.na(theta)
    if (sum(nna) > 0)
        warning(paste(sum(nna), "NA(s) removed from the data."))
    xl <- xl[!nna]
    theta <- theta[!nna]
    ## Correlation arguments
    rxc <- cor(xl, cos(theta))
    rxs <- cor(xl, sin(theta))
    rcs <- cor(cos(theta), sin(theta))
    ## Return the linear-circular correlation
    return((rxc^2 + rxs^2 - 2 * rxc * rxs * rcs)/(1 - rcs^2))
}
##' @export
##' @rdname lincircor
lincircor.test <- function(xl, theta) {
    ## Computes the linear-circular correlation
    lincir <- lincircor(xl, theta)
    names(lincir) <- "cor"
    ## Df of the test
    df2 <- sum(!is.na(xl) & !is.na(theta)) - 3
    parm <- c(2, df2)
    names(parm) <- c("df1", "df2")
    ## Compute the test statistic
    stat <- df2 * lincir/(1 - lincir)
    names(stat) <- "stat"
    ## P-value of the test
    pv <- 1 - pf(stat, parm[1], parm[2])
    names(pv) <- "p.value"
    ## Null value of the test
    nv <- 0
    names(nv) <- "cor"
    ## Data names
    Call <- match.call()
    nm <- paste(as.character(Call[2]), as.character(Call[3]),
        sep = " and ")
    ## Make a 'htest' test (package stats) for pretty printing
    test <- list(statistic = stat, parameter = parm, p.value = pv,
        estimate = lincir, null.value = nv, method = "Linear-circular correlation",
        data.name = nm, alternative = "greater")
    class(test) <- "htest"
    return(test)
}


## ltraj2spdf & ltraj2sldf
##
##' \code{ltraj2spdf}: Add a \code{proj4string} parameter, and keeps
##' \code{id} and \code{burst}. \code{ltraj2sldf}: Add a
##' \code{proj4string} parameter, and \code{by} is now a character,
##' which accepts \code{"id"}, \code{"burst"} (default) and
##' \code{"step"}. Warning: the conversion to SLDF using \code{by =
##' "step"} is significantly longer, of a factor of ~50, than the
##' other two (it takes about 10 sec for each 10,000 steps on a
##' relatively recent computer).
##'
##' @title Conversion of the class \code{ltraj} to the package \code{sp}
##' @param strict Logical (\code{TRUE} by default). See \code{ld}. For
##' \code{ltraj2sldf}, only applies to \code{by = "step"}.
##' @param proj4string A valid proj4 string (see
##' \code{\link[sp]{proj4string}} for more details).
##' @seealso See \code{\link[adehabitatLT]{ltraj2spdf}} for further
##' details on the function and all available arguments.
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
##' ## Conversion to SPDF:
##' summary(adehabitatLT:::ltraj2spdf(puechcirc))
##' summary(ltraj2spdf(puechcirc, strict = FALSE))
##' summary(ltraj2spdf(puechcirc, strict = FALSE, proj4string = CRS("+init=epsg:27572")))
##'
##' ## Conversion to SLDF:
##' summary(adehabitatLT:::ltraj2sldf(puechcirc, byid = TRUE))
##' summary(ltraj2sldf(puechcirc, by = "id"))
##' summary(ltraj2sldf(puechcirc, proj4string = CRS("+init=epsg:27572")))
##' summary(ltraj2sldf(puechcirc, strict = FALSE, by = "step",
##'     proj4string = CRS("+init=epsg:27572")))
ltraj2spdf <- function (ltr, strict = TRUE, proj4string = CRS(as.character(NA)))
{
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class \"ltraj\"")
    ## Use ld to retrieve id/burst/infolocs at once
    ## tr <- do.call("rbind", ltr)
    ## if (!is.null(infolocs(ltr))) {
    ##     tr2 <- do.call("rbind", infolocs(ltr))
    ##     tr <- cbind(tr, tr2)
    ## }
    ## class(tr) <- "data.frame"
    tr <- ld(ltr, strict = strict)
    ## End of modification
    ## Use subset on x | y to remove NAs
    ## xy <- tr[!is.na(tr$x), c("x", "y")]
    ## tr <- tr[!is.na(tr$x), ]
    ## tr$y <- tr$x <- NULL
    tr <- subset(tr, !(is.na(x) | is.na(y)))
    ## End of modification
    ## Allows for a proj4string in the SPDF, and
    ## returns tr
    ## res <- SpatialPointsDataFrame(xy, tr)
    ## return(res)
    tr <- SpatialPointsDataFrame(tr[, c("x", "y")], tr[, !(colnames(tr) %in%
        c("x", "y"))], proj4string = proj4string)
    return(tr)
    ## End of modification
}
##' @param by Sets the level of aggregation: if \code{id}, one object
##' of class \code{Lines} corresponds to one animal; if \code{burst}
##' (default), it corresponds to one burst; if \code{step}, it
##' corresponds to one step.
##' @rdname ltraj2spdf
##' @export
ltraj2sldf <- function (ltr, by = c("burst", "id", "step"), strict = TRUE,
    proj4string = CRS(as.character(NA)))
{
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class \"ltraj\"")
    ## Match the 'by' argument
    by <- match.arg(by)
    ## End of modification
    ## Embed id/burst in a 'if' statement
    if (by == "id" | by == "burst") {
    ## End of modification
        ## Use subset on x | y to remove NAs
        ## lixy <- lapply(ltr, function(x) Line(as.matrix(x[!is.na(x$x),
        ##     c("x", "y")])))
        lixy <- lapply(ltr, function(x) Line(as.matrix(subset(x,
            !(is.na(x) | is.na(y)), select = c("x", "y")))))
        ## End of modification
        id <- unlist(lapply(ltr, function(x) attr(x, "id")))
        bu <- unlist(lapply(ltr, function(x) attr(x, "burst")))
        ## 'by' is now a character, not a logical ('id')
        ## if (byid) {
        if (by == "id") {
        ## End of modification
            lev <- levels(factor(id))
            ## Replace 're1' by 'res' for consistency and memory
            ## usage, and keep the conversion in SpatialLines for
            ## later
            ## re1 <- lapply(lev, function(x) Lines(lixy[id == x], ID = x))
            ## res <- SpatialLines(re1)
            res <- lapply(lev, function(x) Lines(lixy[id == x],
                ID = x))
            ## End of modification
            df <- data.frame(id = lev)
            row.names(df) <- lev
        }
        ## 'by' is now a character, not a logical ('burst')
        ## else {
        if (by == "burst") {
        ## End of modification
            res <- lapply(1:length(lixy), function(i) Lines(list(lixy[[i]]),
                ID = bu[i]))
            ## Keep the conversion in SpatialLines for later
            ## res <- SpatialLines(res)
            ## End of modification
            df <- data.frame(id = id, burst = bu)
            row.names(df) <- bu
        }
    }
    ## Add the 'step' solution
    if (by == "step") {
        ## We first convert the ltraj to a dataframe
        df <- ld(ltr, strict = strict)
        ## We remove NAs in dist, which actually removes NAs in x/y +
        ## dx/dy (required for steps)
        df <- subset(df, !is.na(dist))
        ## We extract only x/y/dx/dy and create an 'id' numeric field
        ## (it has to be numeric for the apply call, otherwise
        ## everything is converted to character)
        coords <- data.frame(df[, c("x", "y", "dx", "dy")], id = as.numeric(row.names(df)))
        ## Each row of the data frame (i.e. each step) forms a
        ## Lines(Line); note the use of 'format(..., scientific =
        ## FALSE')' to ensure that the textual representation of the
        ## 'id' corresponds exactly to the full numeric representation
        res <- apply(coords, 1, function(dfi) Lines(Line(matrix(c(dfi["x"],
            dfi["y"], dfi["x"] + dfi["dx"], dfi["y"] + dfi["dy"]),
            ncol = 2, byrow = TRUE)), ID = format(dfi["id"],
            scientific = FALSE, trim = TRUE)))
    }
    ## End of modification
    ## Convert to SpatialLines, and allows to set the 'proj4string'
    ## res <- SpatialLinesDataFrame(res, data = df)
    res <- SpatialLinesDataFrame(SpatialLines(res, proj4string = proj4string),
        data = df)
    ## End of modification
    return(res)
}


## makeCluster
##
##' Create independant clusters from a sequence of (observed + random)
##' steps.
##'
##' There must be one and only one case per strata. In addition, the
##' sequences must be orderd in ascending order (there is no check on
##' this), and the case and strata should correspond to that order.
##'
##' It is only necessary to provide \code{nclust} to get at least that
##' many clusters (if possible). The argument \code{nloc} can be set
##' instead, if one wants exactly a number of successive steps in a
##' cluster.
##' @title Create independent clusters
##' @param seq A numeric, indicating the continuous sequences of
##' steps. See Details.
##' @param strata A numeric, indicating the strata of each observed
##' and random steps.
##' @param case A numeric, with 1 for observed steps and 0 for random
##' steps.
##' @param ndrop The number of steps to drop between each
##' cluster. They will be NAs in the resulting vector.
##' @param nclust The ideal number of clusters to obtain. See Details.
##' @param nstep The number of steps consisting in a single
##' cluster. See Details.
##' @return A vector of the same length as \code{seq}, \code{strata}
##' and \code{case}, with the cluster number for each step, and three
##' attributes giving the number of clusters (\code{nclust}), the
##' number of successive steps kept per cluster (\code{nstep}) and the
##' number of successive steps dropped between each cluster
##' (\code{ndrop}).
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## case (1 observed + 9 random)
##' case <- rep(rep(1:0, c(1, 9)), 100)
##' ## 100 stratas of 10 steps (1 observed + 9 random)
##' strata <- rep(1:100, each = 10)
##' ## 5 sequences of 20*10 steps
##' seq <- rep(1:5, each = 200)
##' head(data.frame(case, strata, seq), 22)
##' ## Make at least 15 clusters by dropping 2 steps between each of them
##' makeCluster(seq, strata, case, ndrop = 2, nclust = 15)
##' ## Same result with 'nstep = 7' direcly
##' makeCluster(seq, strata, case, ndrop = 2, nstep = 7)
makeCluster <- function(seq, strata, case, ndrop, nclust = NULL,
    nstep = NULL)
{
    ## Check that there is exactly 1 case for each strata, else stop
    tabcase <- table(case, strata)
    if (any(tabcase["1", ] != 1))
        stop("Some stratas with no or more than one case.")
    ## One and only one of 'nclust' or 'nstep' must be provided
    if ((is.null(nclust) & is.null(nstep)) | (!is.null(nclust) &
        !is.null(nstep)))
        stop("One, and only one, of 'nclust' or 'nstep' must be provided.")
    ## Sub-function to create the clusters
    cluster <- function(seq, strata, case, ndrop, nstep) {
        ## Create the pattern [(0*nstep),(NA*ndrop)]
        pat <- rep(c(0, NA), times = c(nstep, ndrop))
        ## Initiate a counter
        compt <- 0
        ## Initiate the 'cluster' vector
        clust <- numeric()
        ## For each sequence, do:
        for (i in unique(seq)) {
            ## How many steps (observed+random) per strata
            ## corresponding to this sequence?
            tabi <- table(strata[seq == i])
            ## How many clusters can fit in this sequence? ('ceiling'
            ## because it suffices that a cluster starts)
            maxi <- ceiling(sum(seq[case == 1] == i)/(nstep +
                ndrop))
            ## If more than one cluster, repeat the pattern, and set
            ## the corresponding cluster number to each kept step; cut
            ## it to the number of cases
            if (maxi > 1)
                clusti <- (rep(pat, times = maxi) + rep(compt +
                  1:maxi, each = (nstep + ndrop)))[1:sum(seq[case ==
                  1] == i)]
            ## If only one cluster, add 'maxi' (1) to the counter, and
            ## repeat it as many times as there are cases
            else clusti <- rep(compt + maxi, sum(seq[case ==
                1] == i))
            ## Repeat the cluster/case given the number of steps and
            ## add it to the already existing vector of clusters
            clust <- c(clust, rep(clusti, times = tabi))
            ## Increment the counter using the max number of clusters
            compt <- compt + maxi
        }
        ## Return the vector of clusters
        return(clust)
    }
    ## If 'clust' is not provided, use 'nstep' directly:
    if (is.null(nclust))
        clust <- cluster(seq = seq, strata = strata, case = case,
            ndrop = ndrop, nstep = nstep)
    ## Else iterate 'nstep' in 'cluster' until the number of cluster
    ## is lower than
    else {
        ## Automatic approximate 'nstep'
        nstep <- floor(sum(case)/(nclust)) - ndrop
        clust <- cluster(seq = seq, strata = strata, case = case,
            ndrop = ndrop, nstep = nstep)
        ## If the number of clusters is equal or greater than the
        ## required one, iterate through increasing 'nstep'
        if (max(clust, na.rm = TRUE) >= nclust) {
            for (i in (nstep + 1):sum(case)) {
                clusttmp <- cluster(seq = seq, strata = strata,
                  case = case, ndrop = ndrop, nstep = i)
                ## If the number of clusters is no more equal or
                ## greater than the required one, breaks
                if (max(clusttmp, na.rm = TRUE) < nclust) {
                  break
                }
                ## Else store 'clusttmp' and 'nstep' (and go to the
                ## next iteration)
                else {
                  clust <- clusttmp
                  nstep <- i
                }
            }
        }
        ## Else if the number of clusters is smaller than the required
        ## one, iterate through decreasing 'nstep'
        else {
            for (i in (nstep - 1):0) {
                clusttmp <- cluster(seq = seq, strata = strata,
                    case = case, ndrop = ndrop, nstep = i)
                ## If the number of clusters is no more smaller than
                ## the required one, breaks
                if (max(clusttmp, na.rm = TRUE) >= nclust) {
                  break
                }
                ## Else store 'clusttmp' and 'nstep' (and go to the
                ## next iteration)
                else {
                  clust <- clusttmp
                  nstep <- i
                }
            }
        }
        if (max(clust, na.rm = TRUE) < nclust)
            stop("Impossible to return that many clusters. Try with a smaller 'nclust'.")
    }
    ## Set 'nclust', 'nstep' and 'ndrop' attributes
    attr(clust, "nclust") <- max(clust, na.rm = TRUE)
    attr(clust, "nstep") <- nstep
    attr(clust, "ndrop") <- ndrop
    return(clust)
}


## modWeights
##
##' Compute model weights using Information Theory criteria.
##'
##' @title Model weights using Information Theory
##' @param ... All fitted model objects.
##' @param criterion The function to be used as the Information Theory
##' criterion. Only tested with \code{QIC} (default).
##' @param names Logical, whether to assign names to the result.
##' @return A numeric vector with the corresponding weights.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
modWeights <- function(..., criterion = QIC, names = TRUE)
{
    ## Compute the criterion of choice on each model
    CritMod <- unlist(lapply(list(...), criterion))
    ## Compute the difference to the minimum criterion
    dCritMod <- CritMod - min(CritMod)
    ## Compute the weights
    weights <- exp(-0.5 * dCritMod)/sum(exp(-0.5 * dCritMod))
    ## If 'names', assign the names of the models
    if (names)
        names(weights) <- as.character(match.call()[-1L])
    ## Return the weights
    return(weights)
}


## plot.ltraj
##
##' New arguments to allow a better control by the user, plus a
##' computation of xlim/ylim that ensures that the bounding box of the
##' maps are the same for each burst (see parameter
##' \code{center}). Also enables the plot of a
##' \code{SpatialPoints(DataFrame)} in the background.
##'
##' @title Graphical Display of an Object of Class "ltraj"
##' @param spotdf An object of class \code{SpatialPoints}.
##' @param center Logical, whether to center each plot around the
##' current burst (default = \code{FALSE}).
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
##' @param spixdfpar A list of arguments that allows the user to
##' modify SpatialPixelsDataFrame display, using any argument
##' available to \code{image}. Default is \code{list(col =
##' gray((240:1)/256))}.
##' @param spoldfpar A list of arguments that allows the user to
##' modify SpatialPolygons display, using any argument available to
##' \code{plot}. Default is \code{list(col = "green")}.
##' @param spotdfpar A list of arguments that allows the user to
##' modify SpatialPoints display, using any argument available to
##' \code{plot}. Default is \code{list(pch = 3, col = "darkgreen")}.
##' @S3method plot ltraj
##' @seealso See \code{\link[adehabitatLT]{plot.ltraj}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##'
##' plot(puechcirc)
##' plot(puechcirc, ppar = list(pch = 2, cex = .5), lpar = list(lty = 2,
##'     col = grey(.5)))
##' plot(puechcirc, perani = FALSE)
##' \dontrun{
##' plot(puechcirc, perani = FALSE, mfrow = c(1, 2))}
##' plot(puechcirc, id = "JE93", perani = FALSE)
##'
##' data(puechabonsp)
##'
##' plot(puechcirc, perani = FALSE, spixdf = puechabonsp$map[,1])
##' plot(puechcirc, perani = FALSE, spixdf = puechabonsp$map[,1],
##'     ppar = list(pch = 2, cex = .5), lpar = list(lty = 2, col = "white"),
##'     spixdfpar = list(col = gray((1:240)/256)))
##'
##' cont <- getcontour(puechabonsp$map[,1])
##' plot(puechcirc, spoldf = cont)
##' plot(puechcirc, spoldf = cont, ppar = list(pch = 2, cex = .5),
##'     lpar = list(lty = 2, col = grey(.5)), spoldfpar = list(col = "cornsilk",
##'         border = grey(.5)))
plot.ltraj <- function(x, id = unique(unlist(lapply(x, attr,
    which = "id"))), burst = unlist(lapply(x, attr, which = "burst")),
    spixdf = NULL, spoldf = NULL, spotdf = NULL, xlim = NULL,
    ylim = NULL, center = FALSE, addpoints = TRUE, addlines = TRUE,
    perani = TRUE, final = TRUE, mfrow, ppar = list(pch = 21,
        col = "black", bg = "white"), lpar = list(), spixdfpar = list(col = gray((240:1)/256)),
    spoldfpar = list(col = "green"), spotdfpar = list(pch = 3,
        col = "darkgreen"), ...)
{
    ## Check the SpatialPoints(DataFrame)
    if (!is.null(spotdf)) {
        if (!inherits(spotdf, "SpatialPoints"))
            stop("spotdf should inherit the class SpatialPoints")
    }
    ## End of modification
    if (!is.null(spoldf)) {
        if (!inherits(spoldf, "SpatialPolygons"))
            stop("spoldf should inherit the class SpatialPolygons")
    }
    if (!is.null(spixdf)) {
        if (!inherits(spixdf, "SpatialPixelsDataFrame"))
            stop("spixdf should inherit the class SpatialPixelsDataFrame")
    }
    if (!inherits(x, "ltraj"))
        stop("x should be an object of class ltraj")
    x <- x[id = id]
    x <- x[burst = burst]
    typeII <- attr(x, "typeII")
    x <- lapply(x, function(i) {
        jj <- i[!is.na(i$x), ]
        attr(jj, "id") <- attr(i, "id")
        attr(jj, "burst") <- attr(i, "burst")
        return(jj)
    })
    class(x) <- c("ltraj", "list")
    attr(x, "typeII") <- typeII
    attr(x, "regular") <- is.regular(x)
    uu <- lapply(x, function(i) {
        i[, c("x", "y")]
    })
    if (!perani)
        idc <- "burst"
    else idc <- "id"
    id <- unique(unlist(lapply(x, function(i) {
        attr(i, idc)
    })))
    ## If more than one plot, force margins to the minimum, otherwise,
    ## use user specific margins
    ## if (length(id) > 1)
    ##     opar <- par(mar = c(0.1, 0.1, 2, 0.1), mfrow = n2mfrow(length(id)))
    if (length(id) > 1)
        mar = c(0.1, 0.1, 2, 0.1)
    else mar = par("mar")
    ## Allows user interaction with 'mfrow'
    if (missing(mfrow)) {
        if (length(id) > 12)
            mfrow <- c(3, 4)
        else mfrow <- n2mfrow(length(id))
    }
    ## If more plots than allowed in the graphic windows, use 'ask =
    ## TRUE'
    if (length(id) > prod(mfrow))
        ask <- TRUE
    else ask <- par("ask")
    opar <- par(mfrow = mfrow, mar = mar, ask = ask)
    on.exit(par(opar))
    ## End of modification
    if (is.null(xlim)) {
        if (perani) {
            ## Note the use of 'adehabitatLT::' to ensure the use of the
            ## 'id' function from this package, and avoid conflicts
            ## (e.g. with 'plyr::id')
            idtt <- unique(adehabitatLT::id(x))
            oo <- lapply(idtt, function(i) unlist(lapply(x[id = i],
                function(j) j$x)))
            ## If center, 'xlim' centered around the range of x
            ## maxxl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
            ## xlim <- lapply(oo, function(ki) c(min(ki), min(ki) +
            ##     maxxl))
            if (center) {
                maxxl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
                xlim <- lapply(oo, function(ki) c(min(ki) - (maxxl -
                    diff(range(ki)))/2, max(ki) + (maxxl -
                    diff(range(ki)))/2))
            }
            ## Else the same for all
            else {
                xrange <- range(unlist(oo))
                xlim <- lapply(oo, function(i) xrange)
            }
            ## End of modification
        }
        else {
            ## If center, 'xlim' centered around the range of x
            ## maxxl <- max(unlist(lapply(x, function(ki) diff(range(ki$x)))))
            ## xlim <- lapply(x, function(ki) c(min(ki$x), min(ki$x) +
            ##     maxxl))
            if (center) {
                maxxl <- max(unlist(lapply(x, function(ki) diff(range(ki$x)))))
                xlim <- lapply(x, function(ki) c(min(ki$x) - (maxxl -
                    diff(range(ki$x)))/2, max(ki$x) + (maxxl -
                    diff(range(ki$x)))/2))
            }
            ## Else the same for all
            else {
                xrange <- range(unlist(lapply(x, function(i) i$x)))
                xlim <- lapply(x, function(i) xrange)
            }
            ## End of modification
        }
    }
    else {
        xlim <- split(rep(xlim, length(id)), gl(length(id), 2))
    }
    if (is.null(ylim)) {
        if (perani) {
            ## Note the use of 'adehabitatLT::' to ensure the use of the
            ## 'id' function from this package, and avoid conflicts
            ## (e.g. with 'plyr::id')
            idtt <- unique(adehabitatLT::id(x))
            oo <- lapply(idtt, function(i) unlist(lapply(x[id = i],
                function(j) j$y)))
            ## If center, 'ylim' centered around the range of y
            ## maxyl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
            ## ylim <- lapply(oo, function(ki) c(min(ki), min(ki) +
            ##     maxyl))
            if (center) {
                maxyl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
                ylim <- lapply(oo, function(ki) c(min(ki) - (maxyl -
                    diff(range(ki)))/2, max(ki) + (maxyl -
                    diff(range(ki)))/2))
            }
            ## Else the same for all
            else  {
                yrange <- range(unlist(oo))
                ylim <- lapply(oo, function(i) yrange)
            }
            ## End of modification
        }
        else {
            ## If center, 'ylim' centered around the range of y
            ## maxyl <- max(unlist(lapply(x, function(ki) diff(range(ki$y)))))
            ## ylim <- lapply(x, function(ki) c(min(ki$y), min(ki$y) +
            ##     maxyl))
            if (center) {
                maxyl <- max(unlist(lapply(x, function(ki) diff(range(ki$y)))))
                ylim <- lapply(x, function(ki) c(min(ki$y) - (maxyl -
                    diff(range(ki$y)))/2, max(ki$y) + (maxyl -
                    diff(range(ki$y)))/2))
            }
            ## Else the same for all
            else {
                yrange <- range(unlist(lapply(x, function(i) i$y)))
                ylim <- lapply(x, function(i) yrange)
            }
            ## End of modification
        }
    }
    else {
        ylim <- split(rep(ylim, length(id)), gl(length(id), 2))
    }
    names(xlim) <- id
    names(ylim) <- id
    for (i in id) {
        if (!is.null(spixdf)) {
            ## In case of 'mfrow = c(1, 1)' (i.e. one plot per window)
            ## if (length(id) == 1) {
            if (isTRUE(all.equal(par("mfrow"), c(1, 1)))) {
            ## End of modification
                ## Allows spixdf's image modification
                ## image(spixdf, col = colspixdf, xlim = xlim[i][[1]],
                ##   ylim = ylim[i][[1]])
                do.call(image, c(list(spixdf, xlim = xlim[i][[1]],
                  ylim = ylim[i][[1]]), spixdfpar))
                ## End of modifictaion
            }
            else {
                ## Allows spixdf's image modification
                ## image(spixdf, col = colspixdf, xlim = xlim[i][[1]],
                ##   ylim = ylim[i][[1]])
                do.call(image, c(list(spixdf, xlim = xlim[i][[1]],
                  ylim = ylim[i][[1]]), spixdfpar))
                ## End of modification
                title(main = i)
            }
        }
        else {
            if (length(id) == 1) {
                plot(1, 1, type = "n", asp = 1, xlim = xlim[i][[1]],
                  ylim = ylim[i][[1]], ...)
            }
            else {
                ## In case of 'mfrow = c(1, 1)' (i.e. one plot per window)
                ## plot(1, 1, type = "n", asp = 1, xlim = xlim[i][[1]],
                ##   ylim = ylim[i][[1]], axes = (length(id) ==
                ##     1, main = i, ...)
                plot(1, 1, type = "n", asp = 1, xlim = xlim[i][[1]],
                  ylim = ylim[i][[1]], axes = isTRUE(all.equal(par("mfrow"),
                    c(1, 1))), main = i, ...)
                ## End of modification
            }
        }
        ## The box needs to be drawn after the polygon if one is
        ## added, otherwise the polygon covers it. Probably easier to
        ## draw one box in every case after all other graphical
        ## functions.
        ## if (length(id) > 1)
        ##     box()
        ## End of modification
        if (!is.null(spoldf)) {
            ## Allows spoldf's plot modification
            ## plot(spoldf, add = TRUE, col = colspoldf)
            do.call(plot, c(list(spoldf, add = TRUE), spoldfpar))
            ## End of modification
        }
        ## Allows to plot a SpatialPoints(DataFrame)
        if (!is.null(spotdf)) {
            do.call(points, c(list(spotdf), spotdfpar))
        }
        ## End of modification
        if (addlines) {
            if (idc == "burst") {
                ## Allows line modification
                ## lines(x[burst = i][[1]]$x, x[burst = i][[1]]$y)
                do.call(lines, c(list(x = x[burst = i][[1]]$x,
                  y = x[burst = i][[1]]$y), lpar))
                ## End of modification
            }
            else {
                xtmp <- x[id = i]
                for (j in 1:length(xtmp)) {
                  ## Allows line modification
                  ## lines(xtmp[[j]]$x, xtmp[[j]]$y)
                  do.call(lines, c(list(x = xtmp[[j]]$x, y = xtmp[[j]]$y),
                    lpar))
                  ## End of modification
                }
            }
        }
        if (addpoints) {
            if (idc == "burst") {
                ## Allows point modification
                ## points(x[burst = i][[1]]$x, x[burst = i][[1]]$y,
                ##   pch = 21, col = "black", bg = "white")
                do.call(points, c(list(x = x[burst = i][[1]]$x,
                  y = x[burst = i][[1]]$y), ppar))
                ## End of modification
            }
            else {
                xtmp <- x[id = i]
                for (j in 1:length(xtmp)) {
                  ## Allows point modification
                  ## points(xtmp[[j]]$x, xtmp[[j]]$y, pch = 21,
                  ##   col = "black", bg = "white")
                  do.call(points, c(list(x = xtmp[[j]]$x, y = xtmp[[j]]$y),
                    ppar))
                  ## End of modification
                }
            }
        }
        if (final) {
            if (idc == "burst") {
                points(x[burst = i][[1]]$x[c(1, length(x[burst = i][[1]]$x))],
                  x[burst = i][[1]]$y[c(1, length(x[burst = i][[1]]$y))],
                  pch = c(2, 14), col = c("blue", "red"), cex = 2,
                  lwd = 2)
            }
            else {
                xtmp <- x[id = i]
                for (j in 1:length(xtmp)) {
                  points(xtmp[[j]]$x[c(1, length(xtmp[[j]]$x))],
                    xtmp[[j]]$y[c(1, length(xtmp[[j]]$x))], pch = c(2,
                      14), col = c("blue", "red"), cex = 2, lwd = 2)
                }
            }
        }
        ## Creates a 'box' in every case
        box()
        ## End of modification
    }
    if (length(id) > 1)
        par(opar)
}


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
##' adehabitatLT:::plotltr(puechcirc, "cos(rel.angle)")
##' plotltr(puechcirc, "cos(rel.angle)")
##' \dontrun{
##' plotltr(puechcirc, "cos(rel.angle)", ppar = list(pch = 2, cex = 2),
##'     lpar = list(lty = 2, lwd = 2), mfrow = c(2, 1))}
##'
##' adehabitatLT:::plotltr(puechcirc, "dist")
##' plotltr(puechcirc, "dist")
##' plotltr(puechcirc, "dist", ppar = list(pch = 3, col = "blue"),
##'     lpar = list(lty = 3, col = "red"), perani = TRUE)
##'
##' adehabitatLT:::plotltr(puechcirc, "dx")
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


## QIC
##
##' @title QIC: Quasi-likelihood under Independence Criterion
##' @param mod A \code{coxph} or \code{clogit} model.
##' @param ... Optionally more fitted model objects.
##' @param details Logical, whether to provide detailed outputs
##' (turned automatically to \code{TRUE} when several models are
##' fitted).
##' @return If \code{details = FALSE}, simply returns the QIC. If
##' \code{details = TRUE}, returns a data frame presenting the QIC,
##' the quasi-likelihood, the number of observations, the number of
##' events, the number of paramaters and the trace. If several models
##' are provided, two additional columns present the delta QIC and the
##' model weights.
##' @author Mathieu Basille \email{basille@@ase-research.org} and
##' Thierry Duchesne
##' @export QIC
##' @S3method QIC coxph
QIC <- function(mod, ..., details = FALSE)
    UseMethod("QIC")
QIC.coxph <- function(mod, ..., details = FALSE)
{
    ## If several models are provided
    if (!missing(...)) {
        ## Check if all models are CoxPH or a CLogit
        if (any(!sapply(list(mod, ...), inherits, "coxph")))
            stop("Object of class 'coxph' or 'clogit' expected.")
        ## Check if robust variances were estimated for all models
        if (any(!sapply(list(mod, ...), function(x) exists("naive.var",
            x))))
            stop("QIC can be computed only if robust variances are estimated.")
        ## Compute the trace term for all models
        trace <- sapply(list(mod, ...), function(x) sum(diag(solve(x$naive.var) %*%
            x$var)))
        ## Extract the quasi-likelihood for all models
        quasi <- sapply(list(mod, ...), function(x) x$loglik[2])
        ## Compute the QIC
        QIC <- -2 * quasi + 2 * trace
        ## Extract the number of observations for all models
        n <- sapply(list(mod, ...), function(x) x$n)
        ## Extract the number of events for all models
        nevent <- sapply(list(mod, ...), function(x) x$nevent)
        ## Extract the number of parameters for all models
        K <- sapply(list(mod, ...), function(x) length(x$coefficients))
        ## Check the number of observations and events
        if (any(n != n[1L]) | any(nevent != nevent[1L]))
            warning("Models are not all fitted to the same number of observations or events.")
        ## Prepare a detailed output
        val <- data.frame(QIC = QIC, QuasiLL = quasi, n = n,
            nevent = nevent, K = K, Trace = trace, deltaQIC = QIC -
                min(QIC), weight = modWeights(mod, ..., criterion = "QIC",
                names = FALSE))
        ## Assign the names of the models as row names
        row.names(val) <- as.character(match.call()[-1L])
        ## Return the data frame
        return(val)
    }
    else {
        ## Check if model is CoxPH or a CLogit
        if (!inherits(mod, "coxph"))
            stop("Object of class 'coxph' or 'clogit' expected.")
        ## Check if robust variances were estimated
        if (!exists("naive.var", mod))
            stop("QIC can be computed only if robust variances are estimated.")
        ## Compute the trace term
        trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
        ## Extract the quasi-likelihood
        quasi <- mod$loglik[2]
        ## If 'details', return a detailed output
        if (details) {
            val <- data.frame(QIC = -2 * quasi + 2 * trace, QuasiLL = quasi,
                n = mod$n, nevent = mod$nevent, K = length(mod$coefficients),
                Trace = trace)
            ## Assign the name of the model as row names
            row.names(val) <- as.character(match.call()$mod)
            ## Return the data frame
            return(val)
        }
        ## Otherwise, just return the QIC
        else return(-2 * quasi + 2 * trace)
    }
}


## rec
##
##' Modified version of \code{\link[adehabitatLT]{rec}} that keeps the
##' original \code{row.names}.
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
        for (i in 1:length(x)) {
            if (!all(row.names(x[[i]]) %in% row.names(lif[[i]]))) {
                x[[i]] <- x[[i]][row.names(x[[i]]) %in% row.names(lif[[i]]),
                  ]
                attr(x[[i]], "infolocs") <- lif[[i]]
            }
            if (!all(row.names(lif[[i]]) %in% row.names(x[[i]]))) {
                lif[[i]] <- lif[[i]][row.names(lif[[i]]) %in%
                  row.names(x[[i]]), ]
                attr(x[[i]], "infolocs") <- lif[[i]]
            }
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


## rdsteps
##
## "infolocs<-" <- function(ltraj, value)
## {
##     if (!inherits(ltraj, "ltraj"))
##         stop("ltraj should be of class ltraj")
##     if (length(value)!=length(ltraj))
##         stop("the assignment should be a list of the same length as ltraj")
##     for (i in (1:length(ltraj))) {
##         df <- value[[i]]
##         if (!inherits(df, "data.frame"))
##             stop("value should be a list of data.frame")
##         if (nrow(df)!=nrow(ltraj[[i]]))
##             stop(paste("The burst", i,
##                        "does not have the same number of elements in ltraj and value"))
##         attr(ltraj[[i]], "infolocs") <- df
##     }
##     return(ltraj)
## }
##
### Mathieu Basille, basille@ase-research.org
### Last modified: 2012-11-19
### 'rdSteps' computes random steps. It takes a data frame as input,
### mandatory fields are 'id', 'x/ystart', 'rel./abs.angle', and
### 'dist'.
###   'distMax' allows to use a maximum distance in the empirical
###     distributions of distances used for random steps.
###   'simult' allows to draw simultaneously step lenght and turning
###     angle from the same step (instead of two independent
###     distributions)
###   'other' allows to draw random steps from distributions of all
###     other individuals.
### If they don't exist, two new columns 'xend' and 'yend' are created
### with the end of the random step coordinates computed based on
### 'dx/dy'.
### A new column 'case' is created (or overwritten) which takes the
### value 1 for observed steps and 0 for random ones.
### The 'protect' argument is to copy the value of the set of
### 'protected' variables from the observed step to the random
### ones. For example, use 'protect = c("Area", "Sex")' to copy the
### value of Area and Sex.
##
## Working version without 'simult'
## rdSteps <- function(df, emp = df, nr = 10, distMax = Inf, other = TRUE,
##     id = "id", xstart = "x", ystart = "y", date = "date", dx = "dx",
##     dy = "dy", dist = "dist", dt = "dt", abs.angle = "abs.angle",
##     rel.angle = "rel.angle", xend = "xend", yend = "yend", case = "case",
##     protect = NULL, reproducible = TRUE)
## {
##     if (!exists(xend, df))
##         df[, xend] <- df[, xstart] + df[, dx]
##     if (!exists(yend, df))
##         df[, yend] <- df[, ystart] + df[, dy]
##     df[, case] <- 1
##     angles <- na.omit(emp[, rel.angle])
##     idA <- emp[!is.na(emp[, rel.angle]), id]
##     dists <- na.omit(emp[emp[, dist] <= distMax, dist])
##     idD <- emp[!is.na(emp[, dist]) & emp[, dist] <= distMax,
##         id]
##     rdStepsId <- function(ldfk) {
##         idk <- as.character(ldfk[1, id])
##         if (other)
##             return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
##                 ], anglesk = angles[idA != idk], distsk = dists[idD !=
##                 idk], i = i))))
##         else return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
##             ], anglesk = angles, distsk = dists, i = i))))
##     }
##     rdStep <- function(pt, anglesk = anglesk, distsk = distsk,
##         i = i) {
##         if (is.na(pt[, xstart]) | is.na(pt[, rel.angle]))
##             return()
##         else {
##             if (reproducible)
##                 set.seed(i)
##             rhord <- sample(anglesk, nr, replace = TRUE)
##             alphard <- pt[, abs.angle] - pt[, rel.angle] + rhord
##             if (reproducible)
##                 set.seed(i)
##             distrd <- sample(distsk, nr, replace = TRUE)
##             rd <- pt
##             rd[1, ] <- NA
##             rd[, id] <- pt[1, id]
##             rd <- rd[rep(1, nr), ]
##             rd[, xstart] <- pt[1, xstart]
##             rd[, ystart] <- pt[1, ystart]
##             rd[, date] <- pt[1, date]
##             rd[, dx] <- cos(alphard * 180/pi) * distrd
##             rd[, dy] <- sin(alphard * 180/pi) * distrd
##             rd[, dist] <- distrd
##             rd[, dt] <- pt[1, dt]
##             rd[, abs.angle] <- alphard
##             rd[, rel.angle] <- rhord
##             rd[, xend] <- rd[, xstart] + rd[, dx]
##             rd[, yend] <- rd[, ystart] + rd[, dy]
##             rd[, case] <- 0
##             if (!is.null(protect))
##                 rd[, protect] <- pt[1, protect]
##             return(rbind(pt, rd))
##         }
##     }
##     ldf <- split(df, f = df[, id])
##     return(do.call("rbind", lapply(ldf, rdStepsId)))
## }
##
## rdSteps <- function(df, emp = df, nr = 10, distMax = Inf, simult = FALSE,
##     other = TRUE, id = "id", xstart = "x", ystart = "y", date = "date",
##     dx = "dx", dy = "dy", dist = "dist", dt = "dt", abs.angle = "abs.angle",
##     rel.angle = "rel.angle", xend = "xend", yend = "yend", case = "case",
##     protect = NULL, reproducible = TRUE) {
##     if (!exists(xend, df))
##         df[, xend] <- df[, xstart] + df[, dx]
##     if (!exists(yend, df))
##         df[, yend] <- df[, ystart] + df[, dy]
##     df[, id] <- factor(df[, id])
##     emp[, id] <- factor(emp[, id])
##     df[, case] <- 1
##     if (simult) {
##         tmp <- na.omit(emp[, c(dist, rel.angle, id)])
##         dists <- tmp[, dist]
##         angles <- tmp[, rel.angle]
##         idA <- idD <- tmp[, id]
##     }
##     else {
##         angles <- na.omit(emp[, rel.angle])
##         idA <- emp[!is.na(emp[, rel.angle]), id]
##         dists <- na.omit(emp[emp[, dist] <= distMax, dist])
##         idD <- emp[!is.na(emp[, dist]) & emp[, dist] <= distMax,
##             id]
##     }
##     rdStepsId <- function(ldfk) {
##         idk <- as.character(ldfk[1, id])
##         if (other)
##             return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
##                 ], anglesk = angles[idA != idk], distsk = dists[idD !=
##                 idk], i = i))))
##         else return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
##             ], anglesk = angles, distsk = dists, i = i))))
##     }
##     rdStep <- function(pt, anglesk = anglesk, distsk = distsk,
##         i = i) {
##         if (is.na(pt[, xstart]) | is.na(pt[, rel.angle]))
##             return()
##         else {
##             if (simult) {
##                 if (reproducible)
##                   set.seed(i)
##                 spl <- sample(1:length(anglesk), nr, replace = TRUE)
##                 rhord <- anglesk[spl]
##                 alphard <- pt[, abs.angle] - pt[, rel.angle] +
##                   rhord
##                 distrd <- distsk[spl]
##             }
##             else {
##                 if (reproducible)
##                   set.seed(i)
##                 rhord <- sample(anglesk, nr, replace = TRUE)
##                 alphard <- pt[, abs.angle] - pt[, rel.angle] +
##                   rhord
##                 if (reproducible)
##                   set.seed(i)
##                 distrd <- sample(distsk, nr, replace = TRUE)
##             }
##             rd <- pt
##             rd[1, ] <- NA
##             rd[, id] <- pt[1, id]
##             rd <- rd[rep(1, nr), ]
##             rd[, xstart] <- pt[1, xstart]
##             rd[, ystart] <- pt[1, ystart]
##             rd[, date] <- pt[1, date]
##             rd[, dx] <- cos(alphard * 180/pi) * distrd
##             rd[, dy] <- sin(alphard * 180/pi) * distrd
##             rd[, dist] <- distrd
##             rd[, dt] <- pt[1, dt]
##             rd[, abs.angle] <- alphard
##             rd[, rel.angle] <- rhord
##             rd[, xend] <- rd[, xstart] + rd[, dx]
##             rd[, yend] <- rd[, ystart] + rd[, dy]
##             rd[, case] <- 0
##             if (!is.null(protect))
##                 rd[, protect] <- pt[1, protect]
##             return(rbind(pt, rd))
##         }
##     }
##     ldf <- split(df, f = df[, id])
##     return(do.call("rbind", lapply(ldf, rdStepsId)))
## }
##
## data(puechcirc)
## puechcirc ## class ltraj
## uu <- ld(puechcirc)
## bli <- rdSteps(uu, simul = TRUE)
##
## load("~/Travail/Data/Québec/Côte-Nord/Loup/Outputs/LoupLocsCN.RData")
## head(LoupLocsCN)
## bla <- LoupLocsCN[1:100, ]
## bla$Id <- factor(bla$Id)
## bli <- rdSteps(bla, id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Sector", "Pack", "Id", "Sex"))
## blu <- rdSteps(bla, emp = LoupLocsCN, id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Sector", "Pack", "Id", "Sex"))
## rdSteps(bla[1, ], id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend")
## rdSteps(bla[1, ], emp = LoupLocsCN, id = "Id", xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Schedule", "Burst"))
## rdSteps(bla[2, ], id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend")
## rdSteps(bla[2, ], emp = LoupLocsCN, id = "Id", xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Schedule", "Burst"))
## debug(rdSteps)
## debug(rdStepsId)
## debug(rdStep)
## bli <- rdSteps(LoupLocsCN, id = "Id", xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Sector", "Pack", "Id", "Sex"))


## scatter.ENFA
##
## scatter.enfa <- function(x, xax = 1, yax = 2, pts = FALSE,
##     nc = TRUE, grid = TRUE, percent = 95, clabel = 1, side = c("top",
##         "bottom", "axes", "none"), posieig = c("none", "top",
##         "bottom"), Adensity, Udensity, Aangle, Uangle, Aborder,
##     Uborder, Acol, Ucol, Alty, Ulty, Abg, Ubg, Ainch, Uinch,
##     ...)
## {
##     side <- match.arg(side)
##     if (!inherits(x, "enfa"))
##         stop("Object of class 'enfa' expected")
##     old.par <- par(no.readonly = TRUE)
##     on.exit(par(old.par))
##     par(mar = c(0.1, 0.1, 0.1, 0.1))
##     x1 <- x$li[, xax]
##     x1 <- c(x1 - diff(range(x1)/50), x1 + diff(range(x1))/50)
##     xlim <- range(x1)
##     y1 <- x$li[, yax]
##     y1 <- c(y1 - diff(range(y1)/50), y1 + diff(range(y1))/50)
##     ylim <- range(y1)
##     pmar <- t(x$mar * x$cw) %*% as.matrix(x$co[, c(xax, yax)])
##     scatterutil.base(dfxy = x$li[, c(xax, yax)], xax = 1, yax = 2,
##         xlim = xlim, ylim = ylim, grid = grid, addaxes = FALSE,
##         cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "",
##         csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL,
##         area = NULL, add.plot = FALSE)
##     if (pts) {
##         if (missing(Acol))
##             Acol <- gray(0.8)
##         if (missing(Ucol))
##             Ucol <- "black"
##         if (missing(Abg))
##             Abg <- gray(0.8)
##         if (missing(Ubg))
##             Ubg <- "black"
##         if (missing(Ainch))
##             Ainch <- 0.03
##         if (missing(Uinch))
##             Uinch <- Ainch * max(x$pr)
##         symbols(x$li[, c(xax, yax)], circles = rep(1, length(x$pr)),
##             fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
##         symbols(x$li[x$pr > 0, c(xax, yax)], circles = x$pr[x$pr >
##             0], fg = Ucol, bg = Ubg, inches = Uinch, add = TRUE)
##         abline(v = 0)
##         abline(h = 0)
##         if (nc)
##             symbols(pmar, circles = 1, fg = "black", bg = "white",
##                 inches = Ainch * 2, add = TRUE)
##     }
##     else {
##         if (missing(Adensity))
##             Adensity <- NULL
##         if (missing(Udensity))
##             Udensity <- NULL
##         if (missing(Aangle))
##             Aangle <- 45
##         if (missing(Uangle))
##             Uangle <- 45
##         if (missing(Aborder))
##             Aborder <- NULL
##         if (missing(Uborder))
##             Uborder <- NULL
##         if (missing(Acol))
##             Acol <- gray(0.95)
##         if (missing(Ucol))
##             Ucol <- gray(0.6)
##         if (missing(Alty))
##             Alty <- NULL
##         if (missing(Ulty))
##             Ulty <- NULL
##         pcff <- function(xy) {
##             mo <- apply(xy, 2, mean)
##             dis <- apply(xy, 1, function(x) sum((x - mo)^2))
##             xy <- xy[dis < quantile(dis, percent/100), ]
##             return(xy[chull(xy[, 1], xy[, 2]), ])
##         }
##         mcpA <- pcff(x$li[, c(xax, yax)])
##         mcpU <- pcff(x$li[rep(1:length(x$pr), x$pr), c(xax, yax)])
##         polygon(mcpA, density = Adensity, angle = Aangle, border = Aborder,
##             col = Acol, lty = Alty)
##         polygon(mcpU, density = Udensity, angle = Uangle, border = Uborder,
##             col = Ucol, lty = Ulty)
##         abline(v = 0)
##         abline(h = 0)
##         if (nc)
##             points(pmar, pch = 21, bg = "white", cex = 1.5)
##     }
##     dfarr <- x$co[, c(xax, yax)]
##     born <- par("usr")
##     k1 <- min(dfarr[, 1])/born[1]
##     k2 <- max(dfarr[, 1])/born[2]
##     k3 <- min(dfarr[, 2])/born[3]
##     k4 <- max(dfarr[, 2])/born[4]
##     k <- c(k1, k2, k3, k4)
##     dfarr <- 0.75 * dfarr/max(k)
##     s.arrow(dfarr, clabel = clabel, addaxes = FALSE, add.plot = TRUE)
##     if (side != "none") {
##         if (xax == 1)
##             xleg <- "mar"
##         else xleg <- paste("sp", xax - 1)
##         if (yax == 1)
##             yleg <- "mar"
##         else yleg <- paste("sp", yax - 1)
##         xl <- par("usr")[1]
##         xr <- par("usr")[2]
##         yd <- par("usr")[3]
##         yu <- par("usr")[4]
##         if (side == "top") {
##             tra <- paste(" xax =", xleg, "\n yax =", yleg)
##             wd <- strwidth(tra, cex = 1)
##             ht <- strheight(tra, cex = 1) * 1.5
##             rect(xl, yu - ht, xl + wd, yu, col = "white", border = 0)
##             text(xl + wd/2, yu - ht/2, tra, cex = 1)
##         }
##         if (side == "bottom") {
##             tra <- paste(" xax =", xleg, "\n yax =", yleg)
##             wd <- strwidth(tra, cex = 1)
##             ht <- strheight(tra, cex = 1) * 1.5
##             rect(xl, yd + ht, xl + wd, yd, col = "white", border = 0)
##             text(xl + wd/2, yd + ht/2, tra, cex = 1)
##         }
##         if (side == "axes") {
##             trax <- paste(" xax =", xleg)
##             wdx <- strwidth(trax, cex = 1)
##             htx <- strheight(trax, cex = 1) * 1.5
##             rect(xr - 1.05 * wdx, 0 + 0.05 * htx, xr - 0.05 *
##                 wdx, 0 + 1.05 * htx, col = "white", border = 0)
##             text(xr - 1.05 * wdx/2, 0 + htx/2, trax, cex = 1)
##             tray <- paste(" yax =", yleg)
##             wdy <- strwidth(tray, cex = 1)
##             hty <- strheight(tray, cex = 1) * 1.5
##             rect(0 + 0.05 * wdy, yu - 1.05 * hty, 0 + 1.05 *
##                 wdy, yu - 0.05 * hty, col = "white", border = 0)
##             text(0 + wdy/2, yu - 1.05 * hty/2, tray, cex = 1)
##         }
##     }
##     add.scatter.eig(x$s, x$nf, xax = xax - 1, yax = yax - 1,
##         posi = posieig, sub = "Sp. eigenvalues")
##     box()
## }
##
## library(adehabitatHS)
## data(lynxjura)
## map <- lynxjura$map
## locs <- lynxjura$locs
## locs <- locs[slot(locs, "data")[,2]!="D",]
## slot(map,"data")[,4] <- sqrt(slot(map,"data")[,4])
## tab <- slot(map, "data")
## pr <- slot(count.points(locs, map), "data")[,1]
## pc <- dudi.pca(tab, scannf = FALSE)
## (enfa1 <- enfa(pc, pr, scannf = FALSE, nf = 3))
##
## scatter(enfa1)
## scatter(enfa1, grid = FALSE, side = "axes", posieig = "bottom")
## scatter(enfa1, grid = FALSE, side = "axes", posieig = "top", xax = 2, yax = 3)


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


## trajdyn
##
##' Modified version of \code{\link[adehabitatLT]{trajdyn}}, which
##' allows to 1) increment by \code{by} points/steps, 2) only display
##' \code{k} previous points/steps, 3) show the current step as a
##' vector, 4) to modify point and line display, 5) show the current
##' burst/loc/(infolocs) and 6) to update interactively a new or
##' existing variable in \code{infolocs}.
##'
##' \code{y} selects the number of previous points/steps to increment
##' at each step. It defaults to \code{1}, i.e. an increment of 1
##' point/step. Choosing anything different than a positive number
##' sets it back to \code{1}.
##'
##' \code{k} selects the number of previous points/steps to
##' display. It defaults to \code{Inf}, i.e. all
##' points/steps. Choosing anything different than a positive number
##' sets it back to \code{Inf}.
##'
##' \code{v} shows (or removes) the current step as a vector. Note
##' that the vector connects the current location to the next
##' available location (even if there are NAs in the data set).
##'
##' \code{s} shows the current buffer and localisation numbers,
##' together with the associated infolocs (if it exists). If a
##' \code{nvar} is requested, only this variable is shown.
##'
##' The argument \code{nvar} allows to work on a given variable: if a
##' variable of that name already exists in \code{infolocs(x)}, the
##' values of the variable are retrieved from there; otherwise, a
##' variable filled with NAs is used. If a \code{nvar} is requested,
##' \code{d} "deletes" the value and resets it to \code{NA}; every
##' other letter not in use in any option is saved in \code{nvar}. The
##' letters available are: "c", "e", "f", "h", "j", "m", "t", "u",
##' "w", "x".
##'
##' On a QWERTY keyboard:
##'
##' we t u
##'   f hj
##'  xcm
##'
##' On a AZERTY keyboard:
##'
##'  e t u
##'   f hj m
##' wxc
##'
##' @title Interactive Display of Objects of Class \code{ltraj}
##' @param by The number of previous points/steps to increment at each
##' step. Default is an increment of 1 point/step.
##' @param only The number of previous points/steps to
##' display. Default is \code{Inf}, i.e. all points/steps.
##' @param ppar A list of arguments that allows the user to modify
##' point display, using any argument available to
##' \code{points}. Default is \code{list(pch = 16)}.
##' @param lpar A list of arguments that allows the user to modify
##' line display, using any argument available to
##' \code{lines}. Default is \code{list(lwd = 2)}.
##' @param nvar A character string giving the name of a variable.
##' @return If a \code{nvar} is provided, return the original ltraj
##' with updated values in \code{infolocs(nvar)}.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' \dontrun{
##' data(puechcirc)
##'
##' ## Use of `by` and `only` to select the previous k points/steps:
##' trajdyn(puechcirc, by = 10, only = 20)
##'
##' ## Use of `nvar` to dynamically fill in new data:
##' (newtraj <- trajdyn(puechcirc, nvar = "Var"))
##' }
trajdyn <- function (x, burst = attr(x[[1]], "burst"), hscale = 1, vscale = 1,
    by = 1, only = Inf, recycle = TRUE, ppar = list(pch = 16), lpar = list(lwd = 2),
    nvar = NULL, display = c("guess", "windows", "tk"), ...)
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class 'ltraj'")
    e1 <- new.env(parent = baseenv())
    ## Prepare 'info': get the whole infolocs of the ltraj
    info <- infolocs(x)
    ## If a nvar is requested
    if (!is.null(nvar)) {
        ## Save the ltraj for backup and prepare a list in e1
        ltr_bkp <- x
        ## If info is NULL, set 'status_info' to "noinfo" (defaults to
        ## 'nvar_exists')
        status_info <- "nvar_exists"
        if (is.null(info))
            status_info <- "noinfo"
        ## Check that 'nvar' is a character
        if (!inherits(nvar, "character"))
            stop("nvar should be a character string")
        ## Get the column nvar in infolocs
        info <- infolocs(x, nvar)
        ## If info is NULL, create an empty info list and set
        ## 'status_info' to "novar"
        if (is.null(info)) {
            status_info <- "novar"
            info <- lapply(x, function(le) setNames(data.frame(rep(NA,
                nrow(le)), row.names = row.names(le)), nvar))
        }
    }
    ## Store info in e1
    assign("info", info, envir = e1)
    typeII <- attr(x, "typeII")
    x <- lapply(x, function(i) {
        jj <- i[!is.na(i$x), ]
        attr(jj, "id") <- attr(i, "id")
        attr(jj, "burst") <- attr(i, "burst")
        return(jj)
    })
    class(x) <- c("ltraj", "list")
    attr(x, "typeII") <- typeII
    attr(x, "regular") <- is.regular(x)
    u <- x
    ## With 'addvec', the ltraj parameters needs to be recalculated
    ## (otherwise dx/dt can be NAs)
    x <- rec(x)
    assign("x", x[burst = burst], envir = e1)
    assign("v", x[burst = burst], envir = e1)
    assign("ajouli", FALSE, envir = e1)
    assign("ajoupo", FALSE, envir = e1)
    assign("ajoubu", FALSE, envir = e1)
    assign("addpoints", TRUE, envir = e1)
    assign("addlines", TRUE, envir = e1)
    ## Prepare the step vector status (FALSE by default)
    assign("addvec", FALSE, envir = e1)
    assign("lim", TRUE, envir = e1)
    assign("buadd", burst, envir = e1)
    assign("K", 1, envir = e1)
    assign("N", nrow(get("x", envir = e1)[[1]]), envir = e1)
    assign("cusr", rep(0 + NA, 4), envir = e1)
    assign("cplt", rep(0 + NA, 4), envir = e1)
    ## Assign 'by' in e1
    assign("by", by, envir = e1)
    ## Assign 'only' in e1
    assign("only", only, envir = e1)
    opt <- options(warn = -1)
    on.exit(options(opt))
    dsp <- substring(match.arg(display), 1, 1)
    if (dsp == "g")
        dsp <- switch(.Platform$OS.type, windows = "w", "t")
    if (dsp == "t" && !require(tkrplot))
        stop("'tkrplot' package needed\n")
    if (dsp == "t")
        assign("hoho", 1, envir = e1)
    replot <- function() {
        opar <- par(mar = c(0, 0, 0, 0), bg = "white")
        tmptmp <- get("x", envir = e1)
        attr(tmptmp[[1]], "id") <- " "
        assign("x", tmptmp, envir = e1)
        if (get("lim", envir = e1)) {
            assign("xlim", range(get("x", envir = e1)[[1]]$x),
                envir = e1)
            assign("ylim", range(get("x", envir = e1)[[1]]$y),
                envir = e1)
        }
        plot(get("x", envir = e1), id = attr(get("x", envir = e1)[[1]],
            "id"), addlines = FALSE, addp = FALSE, final = FALSE,
            xlim = get("xlim", envir = e1), ylim = get("ylim",
                envir = e1), ...)
        assign("cusr", par("usr"), envir = e1)
        assign("cplt", par("plt"), envir = e1)
        scatterutil.sub(as.character(get("x", envir = e1)[[1]]$date[get("K",
            envir = e1)]), 1, "topleft")
        if (get("ajoubu", envir = e1)) {
            lapply(u[burst = get("buadd", envir = e1)], function(zz) {
                if (get("addpoints", envir = e1))
                  points(zz[, c("x", "y")], pch = 16, col = "grey")
                if (get("addlines", envir = e1))
                  lines(zz[, c("x", "y")], pch = 16, col = "grey")
            })
        }
        ## addlines before addpoints
        if (get("addlines", envir = e1))
            if (get("K", envir = e1) > 1)
                ## Plot only the last 'only' steps, and allows for line
                ## modification:
                ## lines(get("x", envir = e1)[[1]][1:get("K", envir = e1),
                ##   c("x", "y")], lwd = 2)
                do.call(lines, c(get("x", envir = e1)[[1]][tail(1:get("K",
                  envir = e1), ifelse(is.null(get("only", envir = e1)),
                  get("K", envir = e1), get("only", envir = e1))),
                  c("x", "y")], lpar))
        if (get("addpoints", envir = e1))
            ## Plot only the last 'only' points, and allows for point
            ## modification:
            ## points(get("x", envir = e1)[[1]][1:get("K", envir = e1),
            ##     c("x", "y")], pch = 16)
            do.call(points, c(get("x", envir = e1)[[1]][tail(1:get("K",
                envir = e1), ifelse(is.null(get("only", envir = e1)),
                get("K", envir = e1), get("only", envir = e1))),
                c("x", "y")], ppar))
        if (get("ajouli", envir = e1))
            lines(c(get("a1", envir = e1)[1], get("a2", envir = e1)[1]),
                c(get("a1", envir = e1)[2], get("a2", envir = e1)[2]),
                lwd = 2, col = "red")
        if (get("ajoupo", envir = e1))
            points(get("a5", envir = e1)[1], get("a5", envir = e1)[2],
                pch = 16, col = "red", cex = 1.7)
        iti <- unlist(get("x", envir = e1)[[1]][get("K", envir = e1),
            c("x", "y")])
        points(iti[1], iti[2], col = "blue", pch = 16, cex = 1.4)
        ## Add the step vector
        if (get("addvec", envir = e1)) {
            xx1 <- get("x", envir = e1)[[1]][get("K", envir = e1), "x"]
            yy1 <- get("x", envir = e1)[[1]][get("K", envir = e1), "y"]
            xx2 <- xx1 + get("x", envir = e1)[[1]][get("K", envir = e1), "dx"]
            yy2 <- yy1 + get("x", envir = e1)[[1]][get("K", envir = e1), "dy"]
            arrows(xx1, yy1, xx2, yy2, lwd = 3, length = .1, col = "blue")
        }
        par(opar)
    }
    ## Remove the final \n
    help.txt <- paste("\n-------- to obtain this help, type 'h' ------------------",
        "n/p            -- Next/Previous relocation", "a              -- show All relocations",
        "y              -- browse the trajectory bY n steps",
        "g              -- Go to...", "0-9            -- show a given part of the path",
        "k              -- display only the K previous steps",
        "s              -- Show burst/relocation/(infolocs)",
        "b              -- change Burst", "i              -- add/remove other bursts on the graph",
        "z/o            -- Zoom in/Out", "Left-Click     -- measure the distance between two points",
        "Right-Click    -- identify a relocation",
        "r/l/v          -- add or remove points/Lines/Vector",
        "q              -- Quit", "---------------------------------------------------------\n", sep = "\n")
    ## If a nvar is requested, complete the help text
    if (!is.null(nvar))
        help.txt <- paste0(help.txt,
            "d              -- Delete nvar for this relocation\n",
            "other letters  -- Set nvar to the letter value\n",
            "---------------------------------------------------------\n\n")
    assign("D", 0, envir = e1)
    assign("a1", 0, envir = e1)
    assign("a2", 0, envir = e1)
    if (dsp == "t") {
        tt <- tcltk::tktoplevel()
        tcltk::tkwm.title(tt, "Exploration of Animal Movements")
        img <- tkrplot::tkrplot(tt, replot, hscale = hscale,
            vscale = vscale)
        txt <- tcltk::tktext(tt, bg = "white", font = "courier 10")
        scr <- tcltk::tkscrollbar(tt, repeatinterval = 5, command = function(...) tcltk::tkyview(txt,
            ...))
        tcltk::tkconfigure(txt, yscrollcommand = function(...) tcltk::tkset(scr,
            ...))
        tcltk::tkpack(img, side = "top")
        tcltk::tkpack(txt, side = "left", fill = "both", expand = TRUE)
        tcltk::tkpack(scr, side = "right", fill = "y")
        iw <- as.numeric(tcltk::tcl("image", "width", tcltk::tkcget(img,
            "-image")))
        ih <- as.numeric(tcltk::tcl("image", "height", tcltk::tkcget(img,
            "-image")))
    }
    showz <- function() switch(dsp, w = replot(), t = {
        tkrplot::tkrreplot(img)
    })
    type <- function(s) switch(dsp, w = cat(s), t = {
        tcltk::tkinsert(txt, "end", s)
        tcltk::tksee(txt, "end")
    })
    type(help.txt)
    cc <- function(x, y) {
        if (dsp == "t") {
            x <- (as.double(x) - 1)/iw
            y <- 1 - (as.double(y) - 1)/ih
        }
        px <- (x - get("cplt", envir = e1)[1])/(get("cplt", envir = e1)[2] -
            get("cplt", envir = e1)[1])
        py <- (y - get("cplt", envir = e1)[3])/(get("cplt", envir = e1)[4] -
            get("cplt", envir = e1)[3])
        ux <- px * (get("cusr", envir = e1)[2] - get("cusr",
            envir = e1)[1]) + get("cusr", envir = e1)[1]
        uy <- py * (get("cusr", envir = e1)[4] - get("cusr",
            envir = e1)[3]) + get("cusr", envir = e1)[3]
        c(ux, uy)
    }
    mm.w <- function(buttons, x, y) {
        if (buttons == 0) {
            i <- get("D", envir = e1)
            if (i == 0) {
                assign("a1", cc(x, y), envir = e1)
                assign("D", 1, envir = e1)
            }
            if (i == 1) {
                assign("a2", cc(x, y), envir = e1)
                assign("D", 0, envir = e1)
                di <- sqrt(sum((get("a2", envir = e1) - get("a1",
                  envir = e1))^2))
                cat(paste("distance:", round(di, 6), "\n"))
                lines(c(get("a1", envir = e1)[1], get("a2", envir = e1)[1]),
                  c(get("a1", envir = e1)[2], get("a2", envir = e1)[2]),
                  lwd = 2, col = "red")
            }
            return()
        }
        if (buttons == 2) {
            w <- get("v", envir = e1)[[1]][1:get("K", envir = e1),
                ]
            assign("a3", cc(x, y), envir = e1)
            di <- sqrt((w$x - get("a3", envir = e1)[1])^2 + (w$y -
                get("a3", envir = e1)[2])^2)
            print(w[which.min(di), ])
            cat("\n")
            points(w[which.min(di), c("x", "y")], pch = 16, col = "red",
                cex = 1.7)
            return()
        }
    }
    mm.t <- function(x, y) {
        i <- get("D", envir = e1)
        if (i == 0) {
            assign("a1", cc(x, y), envir = e1)
            assign("D", 1, envir = e1)
        }
        if (i == 1) {
            assign("a2", cc(x, y), envir = e1)
            assign("D", 0, envir = e1)
            di <- sqrt(sum((get("a2", envir = e1) - get("a1",
                envir = e1))^2))
            type(paste("distance:", di, "\n"))
            assign("ajouli", TRUE, envir = e1)
            showz()
            assign("ajouli", FALSE, envir = e1)
        }
        return()
    }
    mm.t2 <- function(x, y) {
        w <- get("v", envir = e1)[[1]][1:get("K", envir = e1),
            ]
        assign("a3", cc(x, y), envir = e1)
        di <- sqrt((w$x - get("a3", envir = e1)[1])^2 + (w$y -
            get("a3", envir = e1)[2])^2)
        assign("a5", unlist(w[which.min(di), c("x", "y")]), envir = e1)
        assign("ajoupo", TRUE, envir = e1)
        showz()
        assign("ajoupo", FALSE, envir = e1)
        tmp <- w[which.min(di), ]
        se <- unlist(lapply((max(nchar(names(tmp)) + nchar(sapply(tmp,
            as.character)) + 1) - nchar(names(tmp)) - nchar(sapply(tmp,
            as.character))), function(zz) paste(rep(" ", zz),
            collapse = "")))
        so <- unlist(lapply(1:length(tmp), function(i) paste(paste(names(tmp)[i],
            as.character(tmp[1, i]), sep = se[i]), "\n")))
        type(paste("Relocation", row.names(w)[which.min(di)],
            ":\n"))
        sapply(so, type)
        type("\n")
        return()
    }
    mm.mouse <- function(buttons, x, y) {
        assign("a8", cc(x, y), envir = e1)
        return()
    }
    mm.mouset <- function(x, y) {
        assign("a8", cc(x, y), envir = e1)
        return()
    }
    kb <- function(A) {
        key <- tolower(A)
        if (key == "q") {
            if (dsp == "t")
                tcltk::tkdestroy(tt)
            return("OK - Finished")
        }
        if (key %in% c(0:9)) {
            if (key > 0)
                assign("K", round(seq(1, get("N", envir = e1),
                  length = 11))[as.numeric(key) + 1], envir = e1)
            if (key == 0)
                assign("K", 1, envir = e1)
            showz()
        }
        if (key == "z") {
            assign("tmppx", (get("cusr", envir = e1)[1:2] - get("cusr",
                envir = e1)[1])/2, envir = e1)
            assign("xlim", c((get("a8", envir = e1)[1] - (get("tmppx",
                envir = e1)[2] - get("tmppx", envir = e1)[1])/2),
                (get("a8", envir = e1)[1] + (get("tmppx", envir = e1)[2] -
                  get("tmppx", envir = e1)[1])/2)), envir = e1)
            assign("tmppy", (get("cusr", envir = e1)[3:4] - get("cusr",
                envir = e1)[3])/2, envir = e1)
            assign("ylim", c((get("a8", envir = e1)[2] - (get("tmppy",
                envir = e1)[2] - get("tmppy", envir = e1)[1])/2),
                (get("a8", envir = e1)[2] + (get("tmppy", envir = e1)[2] -
                  get("tmppy", envir = e1)[1])/2)), envir = e1)
            assign("lim", FALSE, envir = e1)
            showz()
        }
        if (key == "o") {
            assign("lim", TRUE, envir = e1)
            showz()
        }
        if (key == "n") {
            if (get("K", envir = e1) <= get("N", envir = e1))
                ## Browse by increments of 'by'
                ## assign("K", get("K", envir = e1) + 1, envir = e1)
                assign("K", get("K", envir = e1) + get("by",
                  envir = e1), envir = e1)
            if (get("K", envir = e1) > get("N", envir = e1)) {
                if (recycle)
                  assign("K", 1, envir = e1)
                if (!recycle) {
                  assign("K", get("N", envir = e1), envir = e1)
                  cat("End of burst !\n")
                }
            }
            showz()
        }
        if (key == "l") {
            assign("addlines", !get("addlines", envir = e1),
                envir = e1)
            showz()
        }
        if (key == "g") {
            if (dsp == "w") {
                recom <- TRUE
                while (recom) {
                  rr <- readline("Enter a relocation number: ")
                  recom <- FALSE
                  if (!(rr %in% row.names(get("x", envir = e1)[[1]]))) {
                    cat("invalid number\n")
                    recom <- TRUE
                  }
                }
                assign("K", which(row.names(get("x", envir = e1)[[1]]) ==
                  as.numeric(rr)), envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(row.names(get("x", envir = e1)[[1]])[1])
                tu <- tcltk::tktoplevel(tt, width = 500, height = 50)
                tcltk::tkwm.title(tu, "Enter a relocation number")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable = lv, width = 50)
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    rr <- tcltk::tclvalue(lv)
                    if (!(rr %in% row.names(get("x", envir = e1)[[1]]))) {
                      tcltk::tkmessageBox(message = "invalid number",
                        type = "ok")
                    }
                    else {
                      assign("K", which(row.names(get("x", envir = e1)[[1]]) ==
                        as.numeric(rr)), envir = e1)
                      showz()
                      tcltk::tkdestroy(tu)
                    }
                  })
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }
        if (key == "r") {
            assign("addpoints", !get("addpoints", envir = e1),
                envir = e1)
            showz()
        }
        ## Display the step vector
        if (key == "v") {
            assign("addvec", !get("addvec", envir = e1),
                envir = e1)
            showz()
        }
        if (key == "b") {
            assign("K", 1, envir = e1)
            if (dsp == "w") {
                assign("hoho", select.list(unlist(lapply(u, function(y) attr(y,
                  "burst")))), envir = e1)
                type(paste("Choice of the burst:", get("hoho",
                  envir = e1), "\n\n"))
                assign("x", u[burst = get("hoho", envir = e1)],
                  envir = e1)
                assign("v", u[burst = get("hoho", envir = e1)],
                  envir = e1)
                assign("N", nrow(get("x", envir = e1)[[1]]),
                  envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(unlist(lapply(u, function(y) attr(y,
                  "burst"))))
                bubu <- unlist(lapply(u, function(y) attr(y,
                  "burst")))
                tu <- tcltk::tktoplevel(tt)
                tcltk::tkwm.title(tu, "Choose a burst of relocations")
                tcltk::tkwm.resizable(tu, 0, 0)
                tfr <- tcltk::tkframe(tu)
                tli <- tcltk::tklistbox(tfr, bg = "white", font = "courier 12",
                  listvariable = lv)
                scr2 <- tcltk::tkscrollbar(tfr, repeatinterval = 5,
                  command = function(...) tcltk::tkyview(tli,
                    ...))
                tcltk::tkconfigure(tli, yscrollcommand = function(...) tcltk::tkset(scr2,
                  ...))
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    assign("hoho", ifelse(nchar(tcltk::tclvalue(tcltk::tkcurselection(tli))) ==
                      0, 1, as.numeric(tcltk::tclvalue(tcltk::tkcurselection(tli))) +
                      1), envir = e1)
                    type(paste("Choice of the burst:", bubu[get("hoho",
                      envir = e1)], "\n\n"))
                    tcltk::tkdestroy(tu)
                  })
                tcltk::tkpack(tli, side = "left", fill = "both",
                  expand = TRUE)
                tcltk::tkpack(scr2, side = "right", fill = "y")
                tcltk::tkpack(tfr, side = "right", fill = "y")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
                assign("x", u[burst = bubu[get("hoho", envir = e1)]],
                  envir = e1)
                assign("v", u[burst = bubu[get("hoho", envir = e1)]],
                  envir = e1)
                assign("N", nrow(get("x", envir = e1)[[1]]),
                  envir = e1)
                showz()
            }
        }
        if (key == "i") {
            if (get("ajoubu", envir = e1)) {
                assign("ajoubu", FALSE, envir = e1)
                showz()
            }
            else {
                if (dsp == "w") {
                  assign("buadd", select.list(unlist(lapply(u,
                    function(y) attr(y, "burst"))), multiple = TRUE),
                    envir = e1)
                  if (length(get("buadd", envir = e1) > 0)) {
                    type(paste("show bursts:", paste(get("buadd",
                      envir = e1), collapse = " "), "\n\n"))
                    assign("ajoubu", TRUE, envir = e1)
                    showz()
                  }
                }
                if (dsp == "t") {
                  lv <- tcltk::tclVar(unlist(lapply(u, function(y) attr(y,
                    "burst"))))
                  bubu <- unlist(lapply(u, function(y) attr(y,
                    "burst")))
                  tu <- tcltk::tktoplevel(tt)
                  tcltk::tkwm.title(tu, "Choose one or several bursts")
                  tcltk::tkwm.resizable(tu, 0, 0)
                  tfr <- tcltk::tkframe(tu)
                  tli <- tcltk::tklistbox(tfr, bg = "white",
                    font = "courier 12", listvariable = lv, selectmode = "multiple")
                  scr2 <- tcltk::tkscrollbar(tfr, repeatinterval = 5,
                    command = function(...) tcltk::tkyview(tli,
                      ...))
                  tcltk::tkconfigure(tli, yscrollcommand = function(...) tcltk::tkset(scr2,
                    ...))
                  submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                    command = function() {
                      argg <- ifelse(nchar(tcltk::tclvalue(tcltk::tkcurselection(tli))) ==
                        0, 1, 0)
                      if (argg == 0) {
                        assign("ajoubu", TRUE, envir = e1)
                        assign("buadd", bubu[as.numeric(unlist(strsplit(tcltk::tclvalue(tcltk::tkcurselection(tli)),
                          " "))) + 1], envir = e1)
                        type(paste("show bursts:", paste(get("buadd",
                          envir = e1), collapse = " "), "\n\n"))
                        showz()
                        tcltk::tkdestroy(tu)
                      }
                    })
                  tcltk::tkpack(tli, side = "left", fill = "both",
                    expand = TRUE)
                  tcltk::tkpack(scr2, side = "right", fill = "y")
                  tcltk::tkpack(tfr, side = "right", fill = "y")
                  tcltk::tkpack(submit.but, side = "bottom")
                  tcltk::tkwait.window(tu)
                  assign("x", u[burst = bubu[get("hoho", envir = e1)]],
                    envir = e1)
                  assign("v", u[burst = bubu[get("hoho", envir = e1)]],
                    envir = e1)
                  assign("N", nrow(get("x", envir = e1)[[1]]),
                    envir = e1)
                  showz()
                }
            }
        }
        if (key == "p") {
            if (get("K", envir = e1) > 1)
                ## Browse by increments of 'by'
                ## assign("K", get("K", envir = e1) - 1, envir = e1)
                assign("K", get("K", envir = e1) - get("by",
                  envir = e1), envir = e1)
            ## Condition becomes K <= 1
            ## if (get("K", envir = e1) == 1) {
            if (get("K", envir = e1) <= 1) {
                if (recycle)
                  assign("K", get("N", envir = e1), envir = e1)
                if (!recycle) {
                  assign("K", 1, envir = e1)
                  cat("Beginning of burst!\n")
                }
            }
            showz()
        }
        if (key == "a") {
            assign("K", get("N", envir = e1), envir = e1)
            showz()
        }
        if (key == "h")
            type(help.txt)
        ## If 'y', change 'by'
        if (key == "y") {
            if (dsp == "w") {
                rr <- as.numeric(readline("Increment by n steps: "))
                assign("by", ifelse(is.na(rr) | rr < 0, 1, rr),
                  envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(ifelse(is.null(get("by",
                  envir = e1)), "NULL", as.character(get("by",
                  envir = e1))))
                tu <- tcltk::tktoplevel(tt, width = 500, height = 50)
                tcltk::tkwm.title(tu, "Increment by n steps")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable = lv, width = 50)
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    rr <- as.numeric(tcltk::tclvalue(lv))
                    assign("by", ifelse(is.na(rr) | rr < 0, 1,
                      rr), envir = e1)
                    showz()
                    tcltk::tkdestroy(tu)
                  })
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }
        ## If 'k', change 'only'
        if (key == "k") {
            if (dsp == "w") {
                rr <- as.numeric(readline("Display only k steps: "))
                assign("only", ifelse(is.na(rr) | rr < 0, Inf,
                    rr), envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(ifelse(is.null(get("only",
                  envir = e1)), "NULL", as.character(get("only",
                  envir = e1))))
                tu <- tcltk::tktoplevel(tt, width = 500, height = 50)
                tcltk::tkwm.title(tu, "Display only k steps")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable = lv, width = 50)
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    rr <- as.numeric(tcltk::tclvalue(lv))
                    assign("only", ifelse(is.na(rr) | rr < 0,
                      Inf, rr), envir = e1)
                    showz()
                    tcltk::tkdestroy(tu)
                  })
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }
        ## If 's', show burst/loc/infolocs
        ## Function 'showloc', which will be reused later...
        showloc <- function() {
            bubu <- unlist(lapply(u, function(y) attr(y, "burst")))
            text <- paste0("Burst: ", bubu[get("hoho", envir = e1)], "; Loc: ",
                row.names(get("x", envir = e1)[[1]])[get("K", envir = e1)])
            ## If there is no infolocs, just show burst/loc
            if (is.null(info))
                type(paste0(text, "\n"))
            ## Else show info
            else {
                ## If no nvar is requested, all variables shown
                if (is.null(nvar)) {
                  infk <- get("info", envir = e1)[[get("hoho",
                    envir = e1)]][which(row.names(get("info",
                    envir = e1)[[get("hoho", envir = e1)]]) ==
                    row.names(get("x", envir = e1)[[1]])[get("K", envir = e1)]), ]
                  infk <- apply(infk, 2, function(x) if(is.numeric(x)) round(x, 3))
                  inftext <- paste(names(infk), as.character(infk), sep = ": ",
                    collapse = "; ")
                }
                ## If a nvar is requested, only show nvar
                else {
                  infv <- get("info", envir = e1)[[get("hoho",
                    envir = e1)]][which(row.names(get("info",
                    envir = e1)[[get("hoho", envir = e1)]]) ==
                    row.names(get("x", envir = e1)[[1]])[get("K", envir = e1)]),
                    nvar]
                  infv <- ifelse(is.numeric(infv), round(infv, 3), infv)
                  inftext <- paste0(nvar, ": ", infv)
                }
                type(paste0(text, "; ", inftext, "\n"))
            }
        }
        if (key == "s") {
            showloc()
        }
        ## If a nvar is requested, all letters not attributed can be
        ## used to save in the new variable
        if (!is.null(nvar) & key %in% c("c", "d", "e", "f", "h",
            "j", "m", "t", "u", "w", "x")) {
            ## Get back the list prepared
            infonew <- get("info", envir = e1)
            ## If D, delete the attribute and set it back to NA (note
            ## that the nvar take initial NAs into account)
            if (key == "d")
                infonew[[get("hoho", envir = e1)]][which(row.names(infonew[[get("hoho",
                  envir = e1)]]) == row.names(get("x", envir = e1)[[1]])[get("K",
                  envir = e1)]), nvar] <- NA
            ## Otherwise, set it to the value of the key
            else infonew[[get("hoho", envir = e1)]][which(row.names(infonew[[get("hoho",
                envir = e1)]]) == row.names(get("x", envir = e1)[[1]])[get("K",
                envir = e1)]), nvar] <- key
            ## Update lvar in e1
            assign("info", infonew, envir = e1)
            ## Displays the burst/loc/nvar
            showloc()
        }
        return()
    }
    showz()
    toto <- switch(dsp, w = getGraphicsEvent("", onKeybd = kb,
        onMouseDown = mm.w, onMouseMove = mm.mouse), t = {
        tcltk::tkbind(tt, "<Key>", kb)
        tcltk::tkbind(img, "<Button-1>", mm.t)
        tcltk::tkbind(img, "<Motion>", mm.mouset)
        tcltk::tkbind(img, "<Button-3>", mm.t2)
        tcltk::tkwait.window(tt)
    })
    ## If a nvar is requested, return the ltraj with modified infolocs
    if (!is.null(nvar)) {
        ## If 'noinfo', use 'info' directly
        if (status_info == "noinfo")
            infolocs(ltr_bkp) <- get("info", envir = e1)
        ## If 'novar', simply merge 'info' with the existing
        ## 'infolocs'
        else if (status_info == "novar")
            infolocs(ltr_bkp) <- mapply(cbind, infolocs(ltr_bkp),
                get("info", envir = e1), SIMPLIFY = FALSE)
        ## If 'nvar' already exists in 'infolocs', replace the column
        ## 'nvar' by the new one in 'info'
        else infolocs(ltr_bkp) <- mapply(function(x, y) {
            x[[nvar]] <- y[[nvar]]
            return(x)
        }, infolocs(ltr_bkp), get("info", envir = e1), SIMPLIFY = FALSE)
        return(ltr_bkp)
    }
}
