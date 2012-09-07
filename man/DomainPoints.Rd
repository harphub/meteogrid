\name{DomainPoints}
\alias{DomainPoints}
\alias{DomainPoints.geofield}
\alias{DomainPoints.default}
\title{Find Lat-Lon (or X-Y) coordinates of all grid points.}
\description{This function returns the co-ordinates of all the grid
  points of a given \code{geofield} domain.}
\usage{DomainLatLon(x,type='lalo')}
\arguments{
  \item{x}{May be either a geofield object or only its domain
    descrtiption (\code{attributes(x,"domain")}.}
  \item{type}{Default is "lalo": the Lat-Lon co-ordinates of all grid
    points. Alternative is "xy" which may be useful for Lambert
    projected grids. It returns the X-Y co-ordinates of the grid points
    as they appear in a graph.}
}
\value{A list containing 2 real-valued matrices \code{lon} and
  \code{lat} (or \code{x} and \code{y})..}
\seealso{\code{DomainExtent}
}
\keyword{file}

