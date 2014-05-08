\name{QIC}
\alias{QIC}
\title{QIC: Quasi-likelihood under Independence Criterion}
\usage{
QIC(mod, ..., details = FALSE)
}
\arguments{
  \item{mod}{A \code{coxph} or \code{clogit} model.}

  \item{...}{Optionally more fitted model objects.}

  \item{details}{Logical, whether to provide detailed
  outputs (turned automatically to \code{TRUE} when several
  models are fitted).}
}
\value{
If \code{details = FALSE}, simply returns the QIC. If
\code{details = TRUE}, returns a data frame presenting the
QIC, the quasi-likelihood, the number of observations, the
number of events, the number of paramaters and the trace.
If several models are provided, two additional columns
present the delta QIC and the model weights.
}
\description{
QIC: Quasi-likelihood under Independence Criterion
}
\author{
Mathieu Basille \email{basille@ase-research.org} and
Thierry Duchesne
}

