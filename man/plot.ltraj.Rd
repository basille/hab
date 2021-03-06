% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plot.ltraj.r
\name{plot.ltraj}
\alias{plot.ltraj}
\title{Graphical Display of an Object of Class "ltraj"}
\usage{
plot.ltraj(x, id = unique(adehabitatLT::id(x)),
  burst = adehabitatLT::burst(x), by = c("burst", "id", "none"),
  na.rm = TRUE, spixdf = NULL, spoldf = NULL, spotdf = NULL,
  xlim = NULL, ylim = NULL, center = FALSE, addpoints = TRUE,
  addlines = TRUE, box = FALSE, final = TRUE, mfrow, ppar = list(pch =
  21, col = "black", bg = "white"), lpar = list(), spixdfpar = list(col =
  gray((240:1)/256)), spoldfpar = list(col = "green"), spotdfpar = list(pch
  = 3, col = "darkgreen"), ...)
}
\arguments{
\item{by}{Character, replaces \code{perani}. Either \code{"burst"}
(identical to \code{perani = FALSE}), \code{"id"} (identical to
\code{perani = TRUE}), or \code{"none"} to plot all bursts at
once.}

\item{na.rm}{Logical, whether to remove missing locations.}

\item{spotdf}{An object of class \code{SpatialPoints}.}

\item{center}{Logical, whether to center each plot around the
current burst (default = \code{FALSE}).}

\item{box}{Logical, whether to add a \code{box} after plotting
each individual plot (useful when a map overlaps the plot
borders, default = \code{FALSE}).}

\item{mfrow}{A vector of the form \code{c(nr, nc)}, which allows
the user to define the numbers of rows (\code{nr}) and columns
(\code{nc}) in the device (the default uses
\code{n2mfrow(length(id))} if \code{length(id) <= 12}, and
\code{mfrow = c(3, 4)} otherwhise).}

\item{ppar}{A list of arguments that allows the user to modify point
display, using any argument available to \code{points}. Default is
\code{list(pch = 21, col = "black", bg = "white")}. See Details.}

\item{lpar}{A list of arguments that allows the user to modify line
display, using any argument available to \code{lines}. Default is
\code{list()}, i.e. an empty list. See Details.}

\item{spixdfpar}{A list of arguments that allows the user to
modify SpatialPixelsDataFrame display, using any argument
available to \code{image}. Default is \code{list(col =
gray((240:1)/256))}.}

\item{spoldfpar}{A list of arguments that allows the user to
modify SpatialPolygons display, using any argument available to
\code{plot}. Default is \code{list(col = "green")}.}

\item{spotdfpar}{A list of arguments that allows the user to
modify SpatialPoints display, using any argument available to
\code{plot}. Default is \code{list(pch = 3, col = "darkgreen")}.}
}
\description{
New arguments to allow a better control by the user, plus a
computation of xlim/ylim that ensures that the bounding box of the
maps are the same for each burst (see parameter
\code{center}). Also enables the plot of a
\code{SpatialPoints(DataFrame)} in the background.
}
\details{
It is possible to use point and line parameters globally for every
trajectory displayed. In this case, \code{ppar} and \code{lpar}
need just be a list of graphical parameters, such as
\code{list(pch = 21, col = "black", bg = "white")}. It is also
possible to use parameters for single bursts or individual, using
atomic vectors with the same numbers of named elements than
bursts/id, such as \code{list(col = c(b1 = "blue", b2 =
"red"))}. Finally, it is also possible to use parameters for
single steps, using as graphical parameter a list of vectors of
length equal to each trajectory. Such information can be based on
\code{infolocs}, see Examples.
}
\examples{
data(puechcirc)

## Point and line parameters
plot(puechcirc)
plot(puechcirc, ppar = list(pch = 2, cex = .5), lpar = list(lty = 2,
    col = grey(.5)))
## By burst
plot(puechcirc, ppar = list(col = c(CH930824 = "blue", CH930827 = "green",
    JE930827 = "red"), pch = 20), lpar = list(col = c(CH930824 = "blue",
    CH930827 = "green", JE930827 = "red")))

## Parameter 'by', 'mfrow' and 'id'
plot(puechcirc, by = "id")
plot(puechcirc, by = "id", mfrow = c(1, 2))
plot(puechcirc, id = "JE93")
plot(puechcirc, by = "none", ppar = list(col = c(CH930824 = "blue",
    CH930827 = "green", JE930827 = "red"), pch = 20), lpar = list(col =
    c(CH930824 = "blue", CH930827 = "green", JE930827 = "red")))

## Using parameters for single steps
info <- list(data.frame(col = sample(c("red",
         "grey"), 80, rep = TRUE), stringsAsFactors = FALSE),
         data.frame(col = sample(c("blue", "darkred"),
             69, rep = TRUE), stringsAsFactors = FALSE),
         data.frame(col = sample(c("darkgreen", "purple"),
             66, rep = TRUE), stringsAsFactors = FALSE))
info <- mapply(function(x, y) {
    row.names(x) <- row.names(y)
    return(x)
}, info, puechcirc, SIMPLIFY = FALSE)
infolocs(puechcirc) <- info

## By burst (default)
plot(puechcirc, ppar = list(pch = 19, col = infolocs(puechcirc,
    "col", simplify = TRUE)), lpar = list(col = infolocs(puechcirc,
    "col", simplify = TRUE)), na.rm = FALSE)

## By animal
plot(puechcirc, by = "id",
    ppar = list(pch = 19, col = infolocs(puechcirc, "col", by = "id",
        simplify = TRUE)),
    lpar = list(col = infolocs(puechcirc, "col", by = "id",
        simplify = TRUE)), na.rm = FALSE)

## Using a SpatialPixelsDataFrame
data(puechabonsp)

plot(puechcirc, spixdf = puechabonsp$map[,1])
plot(puechcirc, spixdf = puechabonsp$map[,1],
    ppar = list(pch = 2, cex = .5), lpar = list(lty = 2, col = "white"),
    spixdfpar = list(col = gray((1:240)/256)))

## Using a SpatialPolygonsDataFrame
cont <- adehabitatMA::getcontour(puechabonsp$map[,1])
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

