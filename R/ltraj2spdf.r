## ltraj2sldf
##
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


## ltraj2spdf
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
##' summary(adehabitatLT::ltraj2spdf(puechcirc))
##' summary(ltraj2spdf(puechcirc, strict = FALSE))
##' summary(ltraj2spdf(puechcirc, strict = FALSE, proj4string = CRS("+init=epsg:27572")))
##'
##' ## Conversion to SLDF:
##' summary(adehabitatLT::ltraj2sldf(puechcirc, byid = TRUE))
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
