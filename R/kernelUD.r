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
