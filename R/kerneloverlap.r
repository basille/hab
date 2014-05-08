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


## summary.kerneloverlap
##
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
