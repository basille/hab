\name{plot.ltraj}
\alias{plot.ltraj}
\title{Graphical Display of an Object of Class "ltraj"}
\usage{
plot.ltraj(x, id = unique(adehabitatLT::id(x)),
  burst = adehabitatLT::burst(x), na.rm = TRUE, spixdf = NULL,
  spoldf = NULL, spotdf = NULL, xlim = NULL, ylim = NULL,
  center = FALSE, addpoints = TRUE, addlines = TRUE, perani = FALSE,
  final = TRUE, mfrow, ppar = list(pch = 21, col = "black", bg = "white"),
  lpar = list(), spixdfpar = list(col = gray((240:1)/256)),
  spoldfpar = list(col = "green"), spotdfpar = list(pch = 3, col =
  "darkgreen"), ...)
}
\arguments{
  \item{na.rm}{Logical, whether to remove missing
  locations.}

  \item{spotdf}{An object of class \code{SpatialPoints}.}

  \item{center}{Logical, whether to center each plot around
  the current burst (default = \code{FALSE}).}

  \item{mfrow}{A vector of the form \code{c(nr, nc)}, which
  allows the user to define the numbers of rows (\code{nr})
  and columns (\code{nc}) in the device (the default uses
  \code{n2mfrow(length(id))} if \code{length(id) <= 12},
  and \code{mfrow = c(3, 4)} otherwhise).}

  \item{ppar}{A list of arguments that allows the user to
  modify point display, using any argument available to
  \code{points}. Default is \code{list(pch = 21, col =
  "black", bg = "white")}. See Details.}

  \item{lpar}{A list of arguments that allows the user to
  modify line display, using any argument available to
  \code{lines}. Default is \code{list()}, i.e. an empty
  list. See Details.}

  \item{spixdfpar}{A list of arguments that allows the user
  to modify SpatialPixelsDataFrame display, using any
  argument available to \code{image}. Default is
  \code{list(col = gray((240:1)/256))}.}

  \item{spoldfpar}{A list of arguments that allows the user
  to modify SpatialPolygons display, using any argument
  available to \code{plot}. Default is \code{list(col =
  "green")}.}

  \item{spotdfpar}{A list of arguments that allows the user
  to modify SpatialPoints display, using any argument
  available to \code{plot}. Default is \code{list(pch = 3,
  col = "darkgreen")}.}
}
\description{
New arguments to allow a better control by the user, plus a
computation of xlim/ylim that ensures that the bounding box
of the maps are the same for each burst (see parameter
\code{center}). Also enables the plot of a
\code{SpatialPoints(DataFrame)} in the background.
}
\details{
It is possible to use point and line parameters globally
for every trajectory displayed. In this case, \code{ppar}
and \code{lpar} need just be a list of graphical
parameters, such as \code{list(pch = 21, col = "black", bg
= "white")}. It is also possible to use parameters for
single steps, using as graphical parameter a list of
vectors of length equal to each trajectory. Such
information can be based on \code{infolocs}, see Example.
}
\examples{
data(puechcirc)

## Point and line parameters
plot(puechcirc)
plot(puechcirc, ppar = list(pch = 2, cex = .5), lpar = list(lty = 2,
    col = grey(.5)))

## id/perani and mfrow
plot(puechcirc, perani = FALSE)
\dontrun{
plot(puechcirc, perani = FALSE, mfrow = c(1, 2))}
plot(puechcirc, id = "JE93", perani = FALSE)

## Using parameters for single steps
infolocs(puechcirc) <- list(data.frame(col = sample(c("red",
    "grey"), 80, rep = TRUE), stringsAsFactors = FALSE),
    data.frame(col = sample(c("blue", "darkred"),
        69, rep = TRUE), stringsAsFactors = FALSE),
    data.frame(col = sample(c("darkgreen", "purple"),
        66, rep = TRUE), stringsAsFactors = FALSE))

## Per burst (default)
plot(puechcirc, ppar = list(pch = 19, col = infolocs(puechcirc,
    "col", simplify = TRUE)), lpar = list(col = infolocs(puechcirc,
    "col", simplify = TRUE)), na.rm = FALSE)

## Per animal
plot(puechcirc, ppar = list(pch = 19, col = infolocs(puechcirc,
    "col", simplify = TRUE, perani = TRUE)), lpar = list(col = infolocs(puechcirc,
    "col", simplify = TRUE, perani = TRUE)), na.rm = FALSE, perani = TRUE)

## Using a SpatialPixelsDataFrame
data(puechabonsp)

plot(puechcirc, perani = FALSE, spixdf = puechabonsp$map[,1])
plot(puechcirc, perani = FALSE, spixdf = puechabonsp$map[,1],
    ppar = list(pch = 2, cex = .5), lpar = list(lty = 2, col = "white"),
    spixdfpar = list(col = gray((1:240)/256)))

## Using a SpatialPolygonsDataFrame
cont <- getcontour(puechabonsp$map[,1])
plot(puechcirc, spoldf = cont)
plot(puechcirc, spoldf = cont, ppar = list(pch = 2, cex = .5),
    lpar = list(lty = 2, col = grey(.5)), spoldfpar = list(col = "cornsilk",
        border = grey(.5)))
}
\author{
Modified by Mathieu Basille
\email{basille@ase-research.org}
}
\seealso{
See \code{\link[adehabitatLT]{plot.ltraj}} for further
details on the function and all available arguments.
}

