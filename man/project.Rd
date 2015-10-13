\name{project}
\alias{project}
\title{Projections}
\description{A wrapper for calling proj4 code for standard projections.}
\usage{
project(x, y, proj=.Last.Domain()$projection, inv=FALSE)
}
\arguments{
  \item{x,y}{Longitude and latitude values (inv=FALSE) or projected values (inv=TRUE). The two vectors must have equal length. NA's are allowed. If y is missing, x is expected to be a vector of length 2 or to have two columns.}
  \item{proj}{A proj4 projection description, given as a list. In the future, a single projection string may also be supported. Default value is the projection of the last geodomain to be plotted.}
  \item{inv}{TRUE indicates inverse projection}
}
\value{A data.frame with columns x and y.
}
\details{ This routine calls the proj.4 library of geographical projections. Longitude and latitude values are always in (decimal) degrees, not radians.
}
\keyword{file}
