\name{lalopoint}
\alias{lalopoint}
\title{Find the grid point that is closest to given co-ordinates.}
\description{
  In the given grid, the point closest to the (Lon,Lat) co-ordinates is
  found. 
}
\usage{
lalopoint(geofield,lon,lat,minimise='lalo')
}
\arguments{
  \item{geofield}{A \code{geofiled} record.}
  \item{lon,lat}{longitude and latitude.}
  \item{minimise}{Which "distance" to minimise: can be 'lalo' (minimise lat-lon) or 'proj' (minimise geographical distance in the projection).  }
}
\value{A list with following members:
  \item{data}{The value of the field in the chosen grid point.}
  \item{x,y}{The actual coordinates of the chosen grid point.}
  \item{i,j}{The matrix indices of the chosen grid point.}
}
\seealso{
  \code{\link{DomainPoints}}
  }
