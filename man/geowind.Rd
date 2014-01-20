\name{geowind}
\alias{geowind}
\title{Convert wind components from model axes to geographic E/N axes.}
\description{
Convert model wind from projection to Lat/Lon axes (or back) using the map factor. This is mainly important for calculating the wind direction. Only Lambert, Rotated LatLon, rotated Mercator and Polar Stereogrphic projections are supported.
}

\usage{
geowind(U, V, inv=FALSE)
}

\arguments{
\item{U,V}{ u and v components of the model wind. These must be \code{geofield} objects.
\item{inv} { if inv=TRUE, the vector is rotated from N/E back to the grid axes.
   }
\value{
A list with components \code{U,V} that have the same type as the input. The wind vectot has been rotated to the geographical N/E axes, or inverse. }
\details{
}
\seealso{\code{\link{FAdec}},\code{\link{wind.direction}}
}
\examples{
  u <- FAdec('ICMSHABOF+0001',"S020WIND.U")
  v <- FAdec('ICMSHABOF+0001',"S020WIND.V")
  uvg <- geowind(u,v) ### a list containing 2 geofield objects
}
\keyword{file}
