\name{trajdyn}
\alias{trajdyn}
\title{Interactive Display of Objects of Class \code{ltraj}}
\usage{
trajdyn(x, burst = attr(x[[1]], "burst"), hscale = 1, vscale = 1,
  by = 1, only = Inf, recycle = TRUE, ppar = list(pch = 16),
  lpar = list(lwd = 2), nvar = NULL, display = c("guess", "windows",
  "tk"), ...)
}
\arguments{
  \item{by}{The number of previous points/steps to
  increment at each step. Default is an increment of 1
  point/step.}

  \item{only}{The number of previous points/steps to
  display. Default is \code{Inf}, i.e. all points/steps.}

  \item{ppar}{A list of arguments that allows the user to
  modify point display, using any argument available to
  \code{points}. Default is \code{list(pch = 16)}.}

  \item{lpar}{A list of arguments that allows the user to
  modify line display, using any argument available to
  \code{lines}. Default is \code{list(lwd = 2)}.}

  \item{nvar}{A character string giving the name of a
  variable.}
}
\value{
If a \code{nvar} is provided, return the original ltraj
with updated values in \code{infolocs(nvar)}.
}
\description{
Modified version of \code{\link[adehabitatLT]{trajdyn}},
which allows to 1) increment by \code{by} points/steps, 2)
only display \code{k} previous points/steps, 3) show the
current step as a vector, 4) to modify point and line
display, 5) show the current burst/loc/(infolocs) and 6) to
update interactively a new or existing variable in
\code{infolocs}.
}
\details{
\code{y} selects the number of previous points/steps to
increment at each step. It defaults to \code{1}, i.e. an
increment of 1 point/step. Choosing anything different than
a positive number sets it back to \code{1}.

\code{k} selects the number of previous points/steps to
display. It defaults to \code{Inf}, i.e. all points/steps.
Choosing anything different than a positive number sets it
back to \code{Inf}.

\code{v} shows (or removes) the current step as a vector.
Note that the vector connects the current location to the
next available location (even if there are NAs in the data
set).

\code{s} shows the current buffer and localisation numbers,
together with the associated infolocs (if it exists). If a
\code{nvar} is requested, only this variable is shown.

The argument \code{nvar} allows to work on a given
variable: if a variable of that name already exists in
\code{infolocs(x)}, the values of the variable are
retrieved from there; otherwise, a variable filled with NAs
is used. If a \code{nvar} is requested, \code{d} "deletes"
the value and resets it to \code{NA}; every other letter
not in use in any option is saved in \code{nvar}. The
letters available are: "c", "e", "f", "h", "j", "m", "t",
"u", "w", "x".

On a QWERTY keyboard:

we t u f hj xcm

On a AZERTY keyboard:

e t u f hj m wxc
}
\examples{
\dontrun{
data(puechcirc)

## Use of `by` and `only` to select the previous k points/steps:
trajdyn(puechcirc, by = 10, only = 20)

## Use of `nvar` to dynamically fill in new data:
(newtraj <- trajdyn(puechcirc, nvar = "Var"))
}
}
\author{
Modified by Mathieu Basille
\email{basille@ase-research.org}
}

