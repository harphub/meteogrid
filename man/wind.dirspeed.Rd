\name{wind.dirspeed}
\alias{wind.dirspeed}
\title{Calculate wind direction and speed from u and v components.}
\description{
A simple calculation of the wind direction for given u and v. Notice that these should already be converted to natural directions along E/N axes.
}

\usage{
wind.dirspeed(u, v,fieldname=c("Wind direction","Wind speed"))
}

\arguments{
\item{u,v}{ u and v components of the model wind. These may be either vectors or \code{geofield} objects. They must be rotated to true N/E axis for the wind direction to be correct! 
   }
\item{fieldname}{A character string giving the names for the result. Only necessary if u and v are \code{geofield} objects.}
}
\value{
A list with elements \code{wdir} and \code{wspeed}.}
\details{
Wind direction is given in degrees. If the wind speed is lower than 10E-6, the result is set to NA.
}
\seealso{
\code{\link{geowind}}
}
\examples{
  u <- FAdec('ICMSHABOF+0001',"S020WIND.U")
  v <- FAdec('ICMSHABOF+0001',"S020WIND.V")
  UV <- geowind(u,v)
  DS <- wind.dirspeed(UV$U,UV$V) ### a list containg 2 geofield objects

}
\keyword{file}
