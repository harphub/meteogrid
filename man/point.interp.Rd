\name{point.interp}
\alias{point.interp}
\alias{point.closest}
\title{Approximate field values at given Lat-Lon co-ordinates.}
\description{Approximate field values for given Lat-Lon co-ordinates. Either by interpolation or by taking the closest grid point value.}
\usage{
point.interp(lon,lat,infield)
point.closest(lon,lat,infield,method="lalo")
}
\arguments{
  \item{infield}{A \code{geofield} object.}
  \item{lon}{A vector of longitues.}
  \item{lat}{A vector of latitudes (same length as lon).}
  \item{method}{Either "lalo" (minimisation in LatLon) or "proj".}

}
\value{A vector containg the interpolated values for the given co-ordinates.
}
\details{ \code{point.interp} uses bilinear interpolation from the \code{fields} package.
}
\examples{
}
\seealso{\code{\link[geogrid]{regrid}}, \code{\link[geogrid]{lalopoint}}
}
\keyword{file}
