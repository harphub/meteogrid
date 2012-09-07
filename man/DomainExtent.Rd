\name{DomainExtent}
\alias{DomainExtent.default}
\alias{DomainExtent.geofield}
\title{Finding the (Lat-Lon) extent of a domain.}
\description{}
\usage{
DomainExtent(domain)
}
\arguments{
  \item{domain}{May be either a geofield object or only its domain
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
\details{
}
\examples{
}
\seealso{\code{DomainPoints}
}
\keyword{file}
