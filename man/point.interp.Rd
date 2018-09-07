\name{point.interp}
\alias{point.interp}
\alias{point.interp.init}
\alias{point.closest}
\alias{point.closest.init}
\alias{point.bilin}
\alias{point.bilin.init}
\alias{point.bicubic}
\alias{point.bicubic.init}
\alias{point.index}
\title{Approximate field values at given Lat-Lon co-ordinates.}
\description{Approximate field values for given Lat-Lon co-ordinates. Either by interpolation or by taking the closest grid point value. The \code{.init} routines allow to calculate the interpolation weights only once. The routines for different methods don't have to be called explicitely. In principle, you should only call \code{point.interp} and \code{point.interp.init} and specify a method.}
\usage{
point.interp(lon, lat, infield, method="bilin", mask=NULL, pointmask=NULL, force=FALSE, weights=NULL)
point.interp.init(lon, lat, domain=.Last.domain(), method="bilin", mask=NULL, pointmask=NULL, force=FALSE)
point.index(lon, lat, domain=.Last.domain(), clip=TRUE)
point.closest(lon, lat, infield, mask=NULL, pointmask=NULL, force=FALSE, weights=NULL)
point.closest.init(lon, lat, domain=.Last.domain(), mask=NULL, pointmask=NULL, force=FALSE)
point.bilin(lon, lat, infield, mask=NULL, pointmask=NULL, force=FALSE, weights=NULL)
point.bilin.init(lon, lat, domain=.Last.domain(), mask=NULL, pointmask=NULL, force=FALSE)
point.bicubic(lon, lat, infield, weights=NULL)
point.bicubic.init(lon, lat, domain=.Last.domain())
}
\arguments{
  \item{infield}{A \code{geofield} object. Or a \code{geodomain} if you are only interested in initialising the weights.}
  \item{domain}{A geodomain or geofield object. The \code{.init} routines may be called with only a domain specification, as the weights do not depend on the field values. If a geofield is provided, the domain attribute is taken.}
  \item{lon}{A vector of longitues.}
  \item{lat}{A vector of latitudes (same length as lon).}
  \item{method}{Can be "bilin", "bicubic" or "closest".}
  \item{weights}{These weights are the result of the \code{point.[method].init} routines. If provided, the weights are not re-calculated. For interpolation of many fields, it is often faster to initialise once and use these weights.}
  \item{mask}{A matrix of the same domain, indicating which points should be masked from the interpolation. Usually a land/sea mask. Maks is not supported for bicubic interpolation.}
  \item{pointmask}{An optional vector that can be used in combination with the mask if some of the points are actually meant to be interpolated from masked points (e.g. temperature at sea). Default is a vector with constant value 1 for all points.}
  \item{force}{When using a mask, it is possible that all 4 points are masked out. If \code{force=TRUE}, the output is set to \code{NA}. The default is to avoid this and revert to using all 4 points for the interpolation.}
  \item{clip}{If TRUE, fractional indices that are lower than 0.5 or higher than gridsize+0.5 are set to NA. If FALSE, points outside the domain will still get an interpolated value.}

}
\value{\code{point.interp()} returns a vector containing the interpolated values for the given co-ordinates. \code{point.index} returns the (fractional) matrix indices, which may be of interest. \code{point.interp.init()} returns a list with interpolation weights that can be used to speed up future calls.
}
\details{
\code{point.bilin} etc. are called by the main routine \code{point.interp}, depending on the chosen method. You don't really have to call them explicitely.
}
\examples{
\dontrun{
w1 <- point.interp.init(lon, lat, domain=.Last.domain(), method="bilin")
f1 <- point.interp(lon, lat, infield, weights=w1, method="bilin")
}
}
\seealso{\code{\link[meteogrid]{regrid}} }
\keyword{file}
