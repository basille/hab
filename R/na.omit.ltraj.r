## na.omit.ltraj
##
##' \code{na.omit} removes missing locations from a \code{ltraj} object.
##' @title Handle Missing Values in Objects of Class 'ltraj'
##' @param object An object of class \code{ltraj}.
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
na.omit.ltraj <- function(object, rec = TRUE)
{
  ## Check ltraj
  if (!inherits(object, "ltraj"))
    stop("object should be an object of class ltraj")
  ## Get typeII attribute
  typeII <- attr(object, "typeII")
  ## Get the position of non NAs,
  nas <- lapply(object, function(i) !is.na(i$x))
  ## Get infolocs and remove it from the ltraj
  info <- infolocs(object)
  object <- removeinfo(object)
  ## If there is infolocs, remove lines with NAs
  if (!is.null(info))
    info <- mapply(function(x, y) {
      x[y, , drop = FALSE]
    }, info, nas, SIMPLIFY = FALSE)
  ## Remove NAs from the ltraj
  object <- lapply(object, function(i) {
    jj <- i[!is.na(i$x), ]
    attr(jj, "id") <- attr(i, "id")
    attr(jj, "burst") <- attr(i, "burst")
    return(jj)
  })
  ## Set back class and ltraj attributes
  class(object) <- c("ltraj", "list")
  attr(object, "typeII") <- typeII
  attr(object, "regular") <- is.regular(object)
  ## Recompute ltraj parameters
  if (rec)
      object <- rec(object)
  ## Associate infolocs without NAs
  if (!is.null(info))
      infolocs(object) <- info
  return(object)
}
