\name{gridpoints}
\alias{gridpoints}
\title{Indicate grid points on the current image.}
\description{
  This little function will return the co-ordinates and put a marker on given grid points. It is
  essentially the same ase \code{points}, but it accepts grid indices in
  stead of geographical co-ordinates.
}
\usage{
gridpoints(x,y=NULL,plot=TRUE,...)
}
\arguments{
  \item{x,y}{Grid indices (column and row) of the points to plot. Either
    2 vectors, or a 2-column matrix. }
  \item{...}{Graphical parameters for the generic \code{points} function.}
}
\value{A list with following members:
  \item{data}{The value of the field in the chosen grid point.}
  \item{x,y}{The actual coordinates of the chosen grid point.}
  \item{i,j}{The matrix indices of the chosen grid point.}
}
\seealso{
  \code{\link{DomainPoints}},\link{points}
  }
