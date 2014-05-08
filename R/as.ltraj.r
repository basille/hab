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
