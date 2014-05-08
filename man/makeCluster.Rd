\name{makeCluster}
\alias{makeCluster}
\title{Create independent clusters}
\usage{
makeCluster(seq, strata, case, ndrop, nclust = NULL, nstep = NULL)
}
\arguments{
  \item{seq}{A numeric, indicating the continuous sequences
  of steps. See Details.}

  \item{strata}{A numeric, indicating the strata of each
  observed and random steps.}

  \item{case}{A numeric, with 1 for observed steps and 0
  for random steps.}

  \item{ndrop}{The number of steps to drop between each
  cluster. They will be NAs in the resulting vector.}

  \item{nclust}{The ideal number of clusters to obtain. See
  Details.}

  \item{nstep}{The number of steps consisting in a single
  cluster. See Details.}
}
\value{
A vector of the same length as \code{seq}, \code{strata}
and \code{case}, with the cluster number for each step, and
three attributes giving the number of clusters
(\code{nclust}), the number of successive steps kept per
cluster (\code{nstep}) and the number of successive steps
dropped between each cluster (\code{ndrop}).
}
\description{
Create independant clusters from a sequence of (observed +
random) steps.
}
\details{
There must be one and only one case per strata. In
addition, the sequences must be orderd in ascending order
(there is no check on this), and the case and strata should
correspond to that order.

It is only necessary to provide \code{nclust} to get at
least that many clusters (if possible). The argument
\code{nloc} can be set instead, if one wants exactly a
number of successive steps in a cluster.
}
\examples{
## case (1 observed + 9 random)
case <- rep(rep(1:0, c(1, 9)), 100)
## 100 stratas of 10 steps (1 observed + 9 random)
strata <- rep(1:100, each = 10)
## 5 sequences of 20*10 steps
seq <- rep(1:5, each = 200)
head(data.frame(case, strata, seq), 22)
## Make at least 15 clusters by dropping 2 steps between each of them
makeCluster(seq, strata, case, ndrop = 2, nclust = 15)
## Same result with 'nstep = 7' direcly
makeCluster(seq, strata, case, ndrop = 2, nstep = 7)
}
\author{
Mathieu Basille \email{basille@ase-research.org}
}

