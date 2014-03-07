\name{geowind}
\alias{geowind}
\title{Convert wind components from model axes to geographic E/N axes.}
\description{
Convert model wind from projection to Lat/Lon axes (or back) using the map factor. This is mainly important for calculating the wind direction. Only Lambert, Rotated LatLon, rotated Mercator and Polar Stereogrphic projections are supported.
}

\usage{
geowind(u, v, inv=FALSE,init=FALSE,angle=NULL)
}

\arguments{
\item{u,v}{ u and v components of the model wind. These must be \code{geofield} objects. If \code{init} is TRUE, u may be a \code{geodomain} object.}
\item{inv}{ if inv=TRUE, the vector is rotated from N/E back to the grid axes.
}
\item{init}{ if init=TRUE, the local rotation angle and map factor are calculated and returned. The result can serve as \code{angle} input to speed up later calls.}
\item{angle}{ If given, this should be the result of an initialisation run of geowind for the current domain. If it is NULL (default), the angles and map factor are calculated. }
}
\value{
A list with components \code{U,V} that have the same type as the input. The wind vector has been rotated to the geographical N/E axes, or inverse. If init=TRUE, th eresult is a list with components \code{angle,mapfactor} that can serve as input for later calls.}
\details{
}
\seealso{\code{\link{wind.dirspeed}}
}
\examples{
  u <- FAdec('ICMSHABOF+0001',"S020WIND.U")
  v <- FAdec('ICMSHABOF+0001',"S020WIND.V")
  uvg <- geowind(u,v) ### a list containing 2 geofield objects
}
\keyword{file}
