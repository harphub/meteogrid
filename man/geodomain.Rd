\name{geodomain}
\alias{geodomain}
\alias{print.geodomain}
\alias{is.geodomain}
\alias{compare.geodomain}
\title{Geodomain data class}
\description{The main object class for geographical domains (grids), described by a projection and a grid.}
\usage{
\method{print}{geodomain}(x, ...)
is.geodomain(x)
compare.geodomain(x,y,eps=1e-10)
}
\arguments{
\item{x,y}{\code{geodomain} class objects.}
\item{eps}{Maximal numerical de=ifference allowed for 'identical' values.}
\item{...}{Not used.}
}
\details{The \code{geodomain} class consists of a list describe a regular gridded domain. An object is list with following elements:
  \itemize{
    \item{projection}{Details of the grid's projection. A list of proj.4 parameters.}
    \item{nx,ny}{Number of gridpoints}
    \item{SW,NE}{(Longitude,Latitude) pairs describing two outer cornes of the grid.}
    \item{dx,dy}{Grid resolution. Either in degrees (LatLon grid) or in meter (projected grid). If missing, these can in fact be calculated from the other components.}
  }
\code{compare.geodomain} compares two domains numerically and returns TRUE if they are equal. This is NOT very reliable. \code{is.geodomain} is a simple shortcut.
}
\seealso{\code{\link{geofield}}, \code{\link{plot.geodomain}},
         \code{\link{DomainExtent}}, \code{\link{DomainPoints}}
}
\keyword{file}
