\name{subgrid}
\alias{subgrid}
\alias{zoomgrid}
\title{Defining subdomains}
\description{These functions may be used to define subdomains (i.e. a
  rectangular subset of grid points) of any \code{geofield} domain.}
\usage{
subgrid(geo, x1, x2, y1, y2, reso = 1) 
zoomgrid(geo, x, y, zoom = 50)
}
\arguments{
  \item{geo}{A \code{geofield} or \code{geodomain} object.}
  \item{x1,y1}{Matrix index of the SW corner of the subdomain.}
  \item{x2,y2}{Matrix index of the NE corner of the subdomain.}
  \item{reso}{Resolution of the subdomain. It may be any positive
    integer. If reso is an integer larger than 1, the resolution is
    dropped by taking only a subset of points with interval reso. Default is 1 (all gridpoints
    are kept).}
  \item{x,y}{Matrix index of the central point of a 'zoomed' grid.}
  \item{zoom}{The number of matrix points at either side of the central
    point of the subgrid.}
 }
\value{A \code{geofield} object with a regular subset of the grid. If the input is a \code{geodomain}, then the result is too.
}
\details{\code{zoomgrid} is in fact just a wrapper for a particular type of call to \code{subgrid}, which is a bit more general. Note that there is no interpolation of any kind involved in these routines. It is simply a subselection of grid boxes. This could also be interpreted as a version of \code{regrid(...,method="closest")}.
}
\seealso{\code{\link{regrid}}
}
\keyword{file}
