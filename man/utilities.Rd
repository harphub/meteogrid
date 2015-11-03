\name{utilities}
\alias{i2a}
\alias{trimk}
\alias{matsum}
\title{Simple utilities}
\description{Some simple functions that may be useful. Just simple shortcuts.}
\usage{
i2a(i,n)
trim(x)
matsum(x, ...)
}
\arguments{
  \item{i,n}{integers}
  \item{x}{A matrix or array (matsum), a character string(trim)}
  \item{...}{More arguments to \code{summary}}.
}
\value{\itemize{
\item{\code{i2a} formats i as an integer of the required length with zero padding.}
\item{\code{matsum} is simply the summary of all the matrix columns at once, i.e. x treated as a vector.}
\item{\code{trim} strip sleading and trailing blanks from a string.}
}}
\keyword{file}
