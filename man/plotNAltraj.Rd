\name{plotNAltraj}
\alias{plotNAltraj}
\title{Highlighting the Patterns in Missing Values in Trajects}
\usage{
plotNAltraj(x, perani = FALSE, addlines = TRUE, mfrow, ppar = list(pch =
  16, cex = 0.3), lpar = list(), ...)
}
\arguments{
  \item{perani}{logical.  If \code{FALSE} (default), one
  plot is drawn for each value of \code{burst}. If
  \code{TRUE}, one plot is drawn for each value of
  \code{id}, and the several bursts are superposed on the
  same plot for a given animal.}

  \item{addlines}{Logical.  Indicates whether lines should
  be added to the plot.}

  \item{mfrow}{A vector of the form \code{c(nr, nc)}, which
  allows the user to define the numbers of rows (\code{nr})
  and columns (\code{nc}) in the device (the default uses
  \code{n2mfrow(length(id))} if \code{length(id) <= 12},
  and \code{mfrow = c(3, 4)} otherwhise).}

  \item{ppar}{A list of arguments that allows the user to
  modify point display, using any argument available to
  \code{points}. Default is \code{list(pch = 21, col =
  "black", bg = "white")}.}

  \item{lpar}{A list of arguments that allows the user to
  modify line display, using any argument available to
  \code{lines}. Default is \code{list()}, i.e. an empty
  list.}
}
\description{
Five new arguments to allow a better control by the user.
}
\examples{
data(puechcirc)
adehabitatLT:::plotNAltraj(puechcirc)
plotNAltraj(puechcirc, perani = TRUE, addlines = FALSE, mfrow = c(1,
    2), ppar = list(pch = 15, cex = 0.5))
}
\author{
Modified by Mathieu Basille
\email{basille@ase-research.org}
}
\seealso{
See \code{\link[adehabitatLT]{plotNAltraj}} for further
details on the function and all available arguments.
}

