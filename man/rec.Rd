\name{rec}
\alias{rec}
\title{Recalculates the descriptive parameters of a ltraj}
\usage{
rec(x, slsp = c("remove", "missing"))
}
\description{
Modified version of \code{\link[adehabitatLT]{rec}} that
keeps the original \code{row.names}.
}
\examples{
data(puechcirc)
(bla <- rec(puechcirc))
head(puechcirc[[2]])
head(bla[[2]])
}
\author{
Modified by Mathieu Basille
\email{basille@ase-research.org}
}
\seealso{
See \code{\link[adehabitatLT]{rec}} for further details on
the function and all available arguments.
}

