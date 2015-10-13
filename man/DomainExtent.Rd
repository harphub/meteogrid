\name{DomainExtent}
\alias{DomainExtent}
\title{Finding the (Lat-Lon) extent of a domain.}
\description{
A function that returns the basic dimensions of a domain.
}
\usage{
DomainExtent(geo)
}
\arguments{
  \item{geo}{May be either a geofield object or only its domain
    descrtiption (\code{attributes(x,"domain")}.}
  }
  \value{
    The result is a list with elements:
    \item{lonlim,latlim}{The extreme longitude and latitude values of he
      domain.}
    \item{clonlat}{The (lon,lat) of the center of the domain.}
    \item{x0,x1,y0,y1}{The SW and NE points of the domain in projection
      co-ordinates.}
    \item{dx,dy}{The interval between grid points in meters.}
    \item{nx,ny}{The number of points in X and Y direction.}
}
\seealso{\code{DomainPoints}
}
\keyword{file}
