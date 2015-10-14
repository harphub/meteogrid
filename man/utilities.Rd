\name{utilities}
\alias{i2a}
\alias{matsum}
\title{Simple utilities}
\description{Some simple functions that may be useful. Just simple shortcuts.}
\usage{
i2a(x,n)
matsum(x, ...)
}
\arguments{
  \item{x,n}{integers}
  \item{x}{A matrix or array }
  \item{...}{More arguments to \code{summary}}.
}
\value{\code{i2a} formats x as an integer of the required length with zero padding. \code{matsum} is simply the summary of all the matrix columns at once, i.e. x treated as a vector.
}
\keyword{file}
