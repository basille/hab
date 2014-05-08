\name{kernelUD}
\alias{kernelUD}
\title{Estimation of Kernel Home-Range}
\usage{
kernelUD(xy, h = "href", grid = 60, same4all = FALSE, hlim = c(0.1,
  1.5), kern = c("bivnorm", "epa"), extent = 0.5, boundary = NULL)
}
\description{
Modified version of \code{\link[adehabitatHR]{kernelUD}},
which silences the use of \code{same4all} when working on a
\code{SpatialPixels} as a \code{grid}.
}
\examples{
## Prepare the data:
data(puechabonsp)
loc <- puechabonsp$relocs
elev <- puechabonsp$map

## Compute the kernel estimation
ker1 <- kernelUD(puechabonsp$relocs[,1], grid = elev, same4all = TRUE)

## Summarizing the overlaps
summary(kerneloverlap(ker1, conditional = TRUE))
}
\author{
Modified by Mathieu Basille
\email{basille@ase-research.org}
}
\seealso{
See \code{\link[adehabitatHR]{kernelUD}} for further
details on the function and all available arguments.
}

