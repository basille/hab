\name{modWeights}
\alias{modWeights}
\title{Model weights using Information Theory}
\usage{
modWeights(..., criterion = QIC, names = TRUE)
}
\arguments{
  \item{...}{All fitted model objects.}

  \item{criterion}{The function to be used as the
  Information Theory criterion. Only tested with \code{QIC}
  (default).}

  \item{names}{Logical, whether to assign names to the
  result.}
}
\value{
A numeric vector with the corresponding weights.
}
\description{
Compute model weights using Information Theory criteria.
}
\author{
Mathieu Basille \email{basille@ase-research.org}
}

