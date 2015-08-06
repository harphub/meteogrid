\name{Make.domain}
\alias{Make.domain}
\title{Create a geographic grid.}
\description{This function returns a \code{geofield} domain with the
  given parameters.}
\usage{Make.domain(projtype="lambert",clonlat,nxny,dxdy,reflat,reflon,tilt)}
\arguments{
  \item{projtype}{The projection type. Possible values are
    \code{"lambert"},\code{"mercator"}, \code{"omerc"} (oblique
    mercator) and \code{"lalo"}.}
  \item{clonlat}{ The c(lon,lat) co-ordinates of the center point of the domain.}
  \item{nxny}{The c(nx,ny) number of grid points, i.e. the dimensions of
    the required grid.}
  \item{dxdy}{The resolution in meter. For lat-lon projection, the unit
    is (decimal) degrees.}
  \item{reflat}{The reference latitude of the projection.}
  \item{reflon}{The reference longitude of the projection.}
}
\value{A \code{geodomain} domain.}
\seealso{}
\examples{
## Not run:
\dontrun{
newdom <- Make.domain("lambert",clonlat=c(50.5,3.5),nxny=c(128,128),dxdy=c(8500,8500),reflon=50.5,reflat=3.5)
plot(newdom,add=FALSE)
}
## End(Not run)
}
\keyword{file}

