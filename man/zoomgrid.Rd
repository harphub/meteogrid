\name{subgrid}
\alias{subgrid}
\alias{zoomgrid}
\title{Defining subdomains}
\description{These functions may be used to define subdomains (i.e. a
  rectangular subset of grid points) of any \code{geofield} domain.}
\usage{
subgrid(infield, x1 = 1, x2 = attr(infield, "domain")$nx, y1 = 1, 
    y2 = attr(infield, "domain")$ny, reso = 1) 
zoomgrid(infield, x, y, zoom = 50)
}
\arguments{
  \item{infield}{A \code{geofield} object.}
  \item{x1,y1}{Matrix index of the SW corner of the subdomain.}
  \item{x2,y2}{Matrix index of the NE corner of the subdomain.}
  \item{reso}{Resolution of the subdomain. It may be any positive
    integer. If reso is an integer larger than 1, the resolution is
    dropped by taking only a subset of points with interval reso. Default is 1 (all gridpoints
    are kept).}
  \item{x,y}{Matrix index of the central point of a "zoomed" grid.}
  \item{zoom}{The number of matrix points at either side of the central
    point of the subgrid.}
 }
\value{A \code{geofield} object.
}
\details{
}
\examples{
}
\seealso{
}
\keyword{file}
