\name{wind.dirspeed}
\alias{wind.dirspeed}
\alias{wind.uv}
\title{Calculate wind direction and speed from u and v components.}
\description{
\code{wind.dispeed} is a simple calculation of the wind direction for given u and v. Notice that these should already be converted to natural directions along E/N axes. \code{wind.uv} is the inverse.
}

\usage{
wind.dirspeed(u, v,fieldname=c("Wind direction","Wind speed"),rad=FALSE)
wind.uv(wspeed,wdir,fieldname=c("U","V"),rad=FALSE)
}

\arguments{
\item{u,v}{ u and v components of the model wind. These may be either vectors or \code{geofield} objects. They must be rotated to true N/E axis for the wind direction to be correct! u may also be a list with components (u,v), in which case v should be missing. 
   }
\item{wspeed,wdir}{Wind speed and direction direction (radians or degrees).}
\item{fieldname}{A character string giving the names for the result. Only necessary if u and v are \code{geofield} objects.}
\item{rad}{If TRUE, the wind direction is in radians. If FALSE, it is in (decimal) degrees.}
}
\value{
A list with elements \code{wdir} and \code{wspeed}.}
\details{
Wind direction is given in degrees by default. If the wind speed is lower than 10E-6, the result is set to NA.
}
\seealso{
\code{\link{geowind}}
}
\keyword{file}
