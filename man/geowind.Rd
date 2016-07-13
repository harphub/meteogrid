\name{geowind}
\alias{geowind}
\alias{geowind.init}
\title{Convert wind components from model axes to geographic E/N axes.}
\description{
Convert model wind from projection to Lat/Lon axes (or back) using the map factor. This is mainly important for calculating the wind direction. Only Lambert, Rotated LatLon, rotated Mercator and Polar Stereogrphic projections are supported.
}

\usage{
geowind(u, v, inv=FALSE, init=NULL)
geowind.init(domain)
}

\arguments{
\item{u, v}{u and v components of the model wind. These must be \code{geofield} objects or vectors of identical length. }
\item{inv}{if inv=TRUE, the vector is rotated from N/E back to the grid axes.
}
\item{init}{Either NULL or a list with components \code{mapfactor} and \code{angle} (e.g. the output from \code{geowind.init}). If \code{u, v} are not geofields, init can not be NULL (because angle can no longer be calculated). Both components must be the same dimension as \code{u, v}. If init=NULL, the local rotation angle and map factor are calculated prior to calculating the rotated wind vectors. }
\item{domain}{A geofield or geodomain object.}
}
\value{
A list with components \code{U,V} that have the same type as the input. The wind vector has been rotated to the geographical N/E axes, or inverse. \code{geowind.init} returns a list with components \code{angle,mapfactor} that can serve as \code{init} value for later calls to \code{geowind}.}
\details{
\code{geowind.init()} is merely meant as a means to save computation time if the same correction must be done many times.
}
\seealso{\code{\link{wind.dirspeed}}
}
\examples{
\dontrun{
  u <- FAdec('ICMSHABOF+0001',"S020WIND.U")
  v <- FAdec('ICMSHABOF+0001',"S020WIND.V")
  uvg <- geowind(u, v) ### a list containing 2 geofield objects
}
}
\keyword{file}
