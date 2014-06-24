## plot.ltraj
##
##' New arguments to allow a better control by the user, plus a
##' computation of xlim/ylim that ensures that the bounding box of the
##' maps are the same for each burst (see parameter
##' \code{center}). Also enables the plot of a
##' \code{SpatialPoints(DataFrame)} in the background.
##'
##' It is possible to use point and line parameters globally for every
##' trajectory displayed. In this case, \code{ppar} and \code{lpar} need
##' just be a list of graphical parameters, such as \code{list(pch = 21,
##' col = "black", bg = "white")}. It is also possible to use parameters
##' for single steps, using as graphical parameter a list of vectors of
##' length equal to each trajectory. Such information can be based on
##' \code{infolocs}, see Example.
##'
##' @title Graphical Display of an Object of Class "ltraj"
##' @param na.rm Logical, whether to remove missing locations.
##' @param spotdf An object of class \code{SpatialPoints}.
##' @param center Logical, whether to center each plot around the
##' current burst (default = \code{FALSE}).
##' @param mfrow A vector of the form \code{c(nr, nc)}, which allows
##' the user to define the numbers of rows (\code{nr}) and columns
##' (\code{nc}) in the device (the default uses
##' \code{n2mfrow(length(id))} if \code{length(id) <= 12}, and
##' \code{mfrow = c(3, 4)} otherwhise).
##' @param ppar A list of arguments that allows the user to modify point
##' display, using any argument available to \code{points}. Default is
##' \code{list(pch = 21, col = "black", bg = "white")}. See Details.
##' @param lpar A list of arguments that allows the user to modify line
##' display, using any argument available to \code{lines}. Default is
##' \code{list()}, i.e. an empty list. See Details.
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
##' @seealso See \code{\link[adehabitatLT]{plot.ltraj}} for further
##' details on the function and all available arguments.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##' ##'
##' ## Point and line parameters
##' plot(puechcirc)
##' plot(puechcirc, ppar = list(pch = 2, cex = .5), lpar = list(lty = 2,
##'     col = grey(.5)))
##' ##'
##' ## id/perani and mfrow
##' plot(puechcirc, perani = FALSE)
##' \dontrun{
##' plot(puechcirc, perani = FALSE, mfrow = c(1, 2))}
##' plot(puechcirc, id = "JE93", perani = FALSE)
##' ##'
##' ## Using parameters for single steps
##' info <- list(data.frame(col = sample(c("red",
##'          "grey"), 80, rep = TRUE), stringsAsFactors = FALSE),
##'          data.frame(col = sample(c("blue", "darkred"),
##'              69, rep = TRUE), stringsAsFactors = FALSE),
##'          data.frame(col = sample(c("darkgreen", "purple"),
##'              66, rep = TRUE), stringsAsFactors = FALSE))
##' info <- mapply(function(x, y) {
##'     row.names(x) <- row.names(y)
##'     return(x)
##' }, info, puechcirc, SIMPLIFY = FALSE)
##' infolocs(puechcirc) <- info
##' ##'
##' ## Per burst (default)
##' plot(puechcirc, ppar = list(pch = 19, col = infolocs(puechcirc,
##'     "col", simplify = TRUE)), lpar = list(col = infolocs(puechcirc,
##'     "col", simplify = TRUE)), na.rm = FALSE)
##' ##'
##' ## Per animal
##' plot(puechcirc, ppar = list(pch = 19, col = infolocs(puechcirc,
##'     "col", simplify = TRUE, perani = TRUE)), lpar = list(col = infolocs(puechcirc,
##'     "col", simplify = TRUE, perani = TRUE)), na.rm = FALSE, perani = TRUE)
##' ##'
##' ## Using a SpatialPixelsDataFrame
##' data(puechabonsp)
##' ##'
##' plot(puechcirc, perani = FALSE, spixdf = puechabonsp$map[,1])
##' plot(puechcirc, perani = FALSE, spixdf = puechabonsp$map[,1],
##'     ppar = list(pch = 2, cex = .5), lpar = list(lty = 2, col = "white"),
##'     spixdfpar = list(col = gray((1:240)/256)))
##' ##'
##' ## Using a SpatialPolygonsDataFrame
##' cont <- getcontour(puechabonsp$map[,1])
##' plot(puechcirc, spoldf = cont)
##' plot(puechcirc, spoldf = cont, ppar = list(pch = 2, cex = .5),
##'     lpar = list(lty = 2, col = grey(.5)), spoldfpar = list(col = "cornsilk",
##'         border = grey(.5)))
plot.ltraj <- function(x, id = unique(adehabitatLT::id(x)), burst =
    adehabitatLT::burst(x), na.rm = TRUE, spixdf = NULL, spoldf = NULL,
    spotdf = NULL, xlim = NULL, ylim = NULL, center = FALSE,
    addpoints = TRUE, addlines = TRUE, perani = FALSE, final = TRUE,
    mfrow, ppar = list(pch = 21, col = "black", bg = "white"),
    lpar = list(), spixdfpar = list(col = gray((240:1)/256)),
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
    ## Are there any lists in point/line parameters?
    plist <- sapply(ppar, is.list)
    llist <- sapply(lpar, is.list)
    ## Check the length of list parameters
    if (any(plist)) {
        if (perani) {
            for (k in (1:length(ppar))[plist])
                if (!isTRUE(all.equal(unlist(lapply(ppar[[k]], length)), unlist(lapply(bindltraj(x), nrow)), check.attributes = FALSE)))
                    stop("Point parameters for individual locations must have the same length as the corresponding burst")
        }
        else {
            for (k in (1:length(ppar))[plist])
                if (!isTRUE(all.equal(unlist(lapply(ppar[[k]], length)), unlist(lapply(x, nrow)), check.attributes = FALSE)))
                    stop("Point parameters for individual locations must have the same length as the corresponding burst")
        }
    }
    if (any(llist)) {
        if (perani)
            {
            for (k in (1:length(lpar))[llist])
                if (!isTRUE(all.equal(unlist(lapply(lpar[[k]], length)), unlist(lapply(bindltraj(x), nrow)), check.attributes = FALSE)))
                    stop("Line parameters for individual steps must have the same length as the corresponding burst")
        }
        else {
            for (k in (1:length(lpar))[llist])
                if (!isTRUE(all.equal(unlist(lapply(lpar[[k]], length)), unlist(lapply(x, nrow)), check.attributes = FALSE)))
                    stop("Line parameters for individual steps must have the same length as the corresponding burst")
        }
    }
    ## End of modification
    ## Parameter na.rm = TRUE/FALSE
    if (na.rm) {
        ## Remove NAs from individual point/line parameters
        nas <- lapply(x, function(i) !is.na(i$x))
        ## Note the use of 'adehabitatLT::' to ensure the use of the 'id'
        ## function from this package, and avoid conflicts (e.g. with
        ## 'plyr::id')
        names(nas) <- adehabitatLT::id(x)
        ## Only if the list of parameter is of length > 0
        if (length(ppar) > 0 & any(plist))
            for (k in (1:length(ppar))[plist])
                ppar[[k]] <- mapply(function(x, y) {
                  x[y]
                }, ppar[[k]], nas, SIMPLIFY = FALSE)
        if (length(lpar) > 0 & any(llist))
            for (k in (1:length(lpar))[llist])
                lpar[[k]] <- mapply(function(x, y) {
                  x[y]
                }, lpar[[k]], nas, SIMPLIFY = FALSE)
        ## Remove NAs from trajectory and rebuild it
        x <- na.omit(x)
        #typeII <- attr(x, "typeII")
        #x <- lapply(x, function(i) {
        #    jj <- i[!is.na(i$x), ]
        #    attr(jj, "id") <- attr(i, "id")
        #    attr(jj, "burst") <- attr(i, "burst")
        #    return(jj)
        #})
        #class(x) <- c("ltraj", "list")
        #attr(x, "typeII") <- typeII
        #attr(x, "regular") <- is.regular(x)
    }
    ## End of modification
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
                ## Allows for NAs
                maxxl <- max(unlist(lapply(oo, function(kk) diff(range(kk,
                  na.rm = TRUE)))))
                xlim <- lapply(oo, function(ki) c(min(ki, na.rm = TRUE) -
                  (maxxl - diff(range(ki, na.rm = TRUE)))/2,
                  max(ki, na.rm = TRUE) + (maxxl - diff(range(ki,
                    na.rm = TRUE)))/2))
            }
            ## Else the same for all
            else {
                ## Allows for NAs
                xrange <- range(unlist(oo), na.rm = TRUE)
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
                ## Allows for NAs
                maxxl <- max(unlist(lapply(x, function(ki) diff(range(ki$x,
                  na.rm = TRUE)))))
                xlim <- lapply(x, function(ki) c(min(ki$x, na.rm = TRUE) -
                  (maxxl - diff(range(ki$x, na.rm = TRUE)))/2,
                  max(ki$x, na.rm = TRUE) + (maxxl - diff(range(ki$x,
                    na.rm = TRUE)))/2))
            }
            ## Else the same for all
            else {
                ## Allows for NAs
                xrange <- range(unlist(lapply(x, function(i) i$x)),
                  na.rm = TRUE)
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
                ## Allows for NAs
                maxyl <- max(unlist(lapply(oo, function(kk) diff(range(kk,
                  na.rm = TRUE)))))
                ylim <- lapply(oo, function(ki) c(min(ki, na.rm = TRUE) -
                  (maxyl - diff(range(ki, na.rm = TRUE)))/2,
                  max(ki, na.rm = TRUE) + (maxyl - diff(range(ki,
                    na.rm = TRUE)))/2))
            }
            ## Else the same for all
            else  {
                ## Allows for NAs
                yrange <- range(unlist(oo), na.rm = TRUE)
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
                ## Allows for NAs
                maxyl <- max(unlist(lapply(x, function(ki) diff(range(ki$y,
                  na.rm = TRUE)))))
                ylim <- lapply(x, function(ki) c(min(ki$y, na.rm = TRUE) -
                  (maxyl - diff(range(ki$y, na.rm = TRUE)))/2,
                  max(ki$y, na.rm = TRUE) + (maxyl - diff(range(ki$y,
                    na.rm = TRUE)))/2))
            }
            ## Else the same for all
            else {
                ## Allows for NAs
                yrange <- range(unlist(lapply(x, function(i) i$y)),
                  na.rm = TRUE)
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
        ## ppark for individual point display parameters
        ppark <- ppar
        if (any(plist)) {
            for (k in (1:length(ppark))[plist])
                ppark[k] <- ppark[[k]][i]
        }
        ## End of modification
        ## lpark for individual line display parameters
        lpark <- lpar
        if (any(llist))
          for (k in (1:length(lpark))[llist])
            lpark[k] <- lpark[[k]][i]
        ## End of modification
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
                #do.call(lines, c(list(x = x[burst = i][[1]]$x,
                #  y = x[burst = i][[1]]$y), lpar))
                do.call(segments, c(list(x0 = x[burst = i][[1]]$x,
                  y0 = x[burst = i][[1]]$y, x1 = x[burst = i][[1]]$x +
                    x[burst = i][[1]]$dx, y1 = x[burst = i][[1]]$y +
                    x[burst = i][[1]]$dy), lpark))
                ## End of modification
            }
            else {
                ## Bind all identical ids together
                ## xtmp <- x[id = i]
                xtmp <- do.call(rbind, x[id = i])
                ## End of modification
                ## Allows line modification
                ## for (j in 1:length(xtmp)) {
                ##   lines(xtmp[[j]]$x, xtmp[[j]]$y)
                ## }
                do.call(segments, c(list(x0 = xtmp$x, y0 = xtmp$y,
                  x1 = xtmp$x + xtmp$dx, y1 = xtmp$y + xtmp$dy),
                  lpark))
                ## End of modification
            }
        }
        if (addpoints) {
            if (idc == "burst") {
                ## Allows point modification
                ## points(x[burst = i][[1]]$x, x[burst = i][[1]]$y,
                ##   pch = 21, col = "black", bg = "white")
                do.call(points, c(list(x = x[burst = i][[1]]$x,
                  y = x[burst = i][[1]]$y), ppark))
                ## End of modification
            }
            else {
                ## Bind all identical ids together
                ## xtmp <- x[id = i]
                xtmp <- do.call(rbind, x[id = i])
                ## End of modification
                ## Allows line modification
                ## for (j in 1:length(xtmp)) {
                ##   points(xtmp[[j]]$x, xtmp[[j]]$y, pch = 21,
                ##     col = "black", bg = "white")
                ## }
                do.call(points, c(list(x = xtmp$x, y = xtmp$y),
                  ppark))
                ## End of modification
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
