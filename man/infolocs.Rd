\name{infolocs}
\alias{infolocs}
\title{Infolocs of an Object of Class ltraj}
\usage{
infolocs(ltraj, which)
}
\description{
Modified version of \code{\link[adehabitatLT]{infolocs}}
that returns \code{NULL} if \code{infolocs} exists, but
\code{which} is not a colum of it.
}
\examples{
## Load puechcirc data, and add some random infolocs:
data(puechcirc)
info <- list(data.frame(A = rnorm(80), B = runif(80), C = rpois(80,
    1)), data.frame(A = rnorm(69, 10), B = runif(69, 10, 11),
    C = rpois(69, 10)), data.frame(A = rnorm(66, -10), B = runif(66,
    -10, -9), C = rpois(66, 100)))
## Watch the row names:
info <- mapply(function(x, y) {
    row.names(x) <- row.names(y)
    return(x)
}, info, puechcirc, SIMPLIFY = FALSE)
infolocs(puechcirc) <- info

## Try to retrieve the column `toto`:
infolocs(puechcirc, "toto")
}
\author{
Modified by Mathieu Basille
\email{basille@ase-research.org}
}
\seealso{
See \code{\link[adehabitatLT]{infolocs}} for further
details on the function and all available arguments.
}

