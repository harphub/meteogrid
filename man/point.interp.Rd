\name{point.interp}
\alias{point.interp}
\alias{point.closest}
\alias{point.bilin}
\alias{point.bicubic}
\alias{point.index}
\title{Approximate field values at given Lat-Lon co-ordinates.}
\description{Approximate field values for given Lat-Lon co-ordinates. Either by interpolation or by taking the closest grid point value.}
\usage{
point.interp(lon,lat,method="bilin",...)
point.closest(lon,lat,infield,mask=NULL,init=FALSE,weights=NULL)
point.index(lon,lat,domain,clip=FALSE)
}
\arguments{
  \item{infield}{A \code{geofield} object. Or a \code{geodomain} if you are only interested in initialising the weights.}
  \item{lon}{A vector of longitues.}
  \item{lat}{A vector of latitudes (same length as lon).}
  \item{method}{Can be "bilin", "bicubic" or "closest".}
  \item{init}{If TRUE, the interpolation weights (or closest indices) are returned. These can be used for later calls.}
  \item{weights}{If given, the weights are not re-calculated. For interpolation of many fields, it is often faster to initialise once and use these weights.}
  \item{mask}{A matrix of the same domain, indicating which points should be masked from the interpolation. Usually a land/sea mask.}
  \item{clip}{If TRUE, fractional indices that are lower than 0.5 or higher than gridsize+0.5 are set to NA.}

}
\value{A vector containg the interpolated values for the given co-ordinates. \code{point.index} returns the (fractional) matrix indices, which may be of interest.
}
\details{
\code{point.bilin} etc. are called by the main routine \code{point.interp}, depending on the chosen method. You don't really have to call them explicitely.
}
\examples{
}
\seealso{\code{\link[geogrid]{regrid}}, }
}
\keyword{file}
