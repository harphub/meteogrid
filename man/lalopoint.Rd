\name{lalopoint}
\alias{lalopoint}
\title{Find the grid point that is closest to given co-ordinates.}
\description{
  In the given grid, the point closest to the (Lon,Lat) co-ordinates is
  found. OBSOLETE! Use \code{point.interp}.
}
\usage{
lalopoint(geofield,lon,lat,minimise='lalo')
}
\arguments{
  \item{geofield}{A \code{geofield, FAfile, geodomain, FAframe} object.}
  \item{lon,lat}{longitude and latitude.}
  \item{minimise}{Which "distance" to minimise: can be 'lalo' (minimise lat-lon) or 'proj' (minimise geographical distance in the projection).  }
}
\value{A list with following members:
  \item{data}{The value of the field in the chosen grid point. If \code{data} is not a \code{geofield}, this is set to NA.}
  \item{lonlat}{The actual coordinates of the chosen grid point.}
  \item{index}{The matrix indices of the chosen grid point.}
}
\seealso{
  \code{point.interp}
}
