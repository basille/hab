## na.omit.ltraj
##
##' \code{na.omit} removes missing locations from a \code{ltraj} object.
##' @title Handle Missing Values in Objects of Class 'ltraj'
##' @param x An object of class \code{ltraj}.
##' @param rec Logical, whether to recompute descriptive parameters of the
##' trajectory (in particular dx, dy, and angles). Use \code{FALSE} with
##' care.
##' @return A ltraj object, without missing locations.
##' @S3method na.omit ltraj
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##' puechcirc
##' na.omit(puechcirc)
na.omit.ltraj <- function(x, rec = TRUE)
{
  ## Check ltraj
  if (!inherits(x, "ltraj"))
    stop("x should be an object of class ltraj")
  ## Get typeII attribute
  typeII <- attr(x, "typeII")
  ## Get the position of non NAs,
  nas <- lapply(x, function(i) !is.na(i$x))
  #    names(nas) <- id(x)
  ## Get infolocs and remove it from the ltraj
  info <- infolocs(x)
  x <- removeinfo(x)
  ## If there is infolocs, remove lines with NAs
  if (!is.null(info))
    info <- mapply(function(x, y) {
      x[y, , drop = FALSE]
    }, info, nas, SIMPLIFY = FALSE)
  ## Remove NAs from the ltraj
  x <- lapply(x, function(i) {
    jj <- i[!is.na(i$x), ]
    attr(jj, "id") <- attr(i, "id")
    attr(jj, "burst") <- attr(i, "burst")
    return(jj)
  })
  ## Set back class and ltraj attributes
  class(x) <- c("ltraj", "list")
  attr(x, "typeII") <- typeII
  attr(x, "regular") <- is.regular(x)
  ## Recompute ltraj parameters
  if (rec)
      x <- rec(x)
  ## Associate infolocs without NAs
  if (!is.null(info))
      infolocs(x) <- info
  return(x)
}
