\name{Make.domain}
\alias{Make.domain}
\alias{Make.domain.RLL}
\title{Create a geographic grid.}
\description{This function returns a \code{geodomain} with the
  given parameters.}
\usage{
Make.domain(projtype="lambert", clonlat, nxny, dxdy, exey=NULL,
            reflat=clonlat[2],reflon=clonlat[1],tilt=0, earth=list(R = 6371229))
Make.domain.RLL(Lon1, Lat1, SPlon, SPlat, SPangle=0, nxny, dxdy)
}
\arguments{
  \item{projtype}{The projection type. Possible values include
    \code{"lambert"},\code{"mercator"}, \code{"omerc"} (oblique
    mercator) and \code{"lonlat"}.}
  \item{clonlat}{ The c(lon,lat) co-ordinates of the center point of the domain.}
  \item{nxny}{The c(nx,ny) number of grid points, i.e. the dimensions of
    the required grid.}
  \item{dxdy}{The resolution in meter. For lat-lon projection, the unit
    is (decimal) degrees.}
  \item{exey}{Extension zone (ref. ALADIN domains). This is only added as an attribute and does not influence the domain itself. But if the domain is later saved in FA format, this is added to the grid dimensions.}
  \item{reflat}{The reference latitude of the projection.}
  \item{reflon}{The reference longitude of the projection.}
  \item{tilt}{Only for rotated mercator.}
  \item{Lon1,Lat1,SPlon,SPlat,SPangle}{HIRLAM-style rotated LatLon grid: rotation of South Pole.}
  \item{earth}{The earth shape description (in proj4 format). Default is a sperical earth .}
}
\value{A \code{geodomain} domain. \code{Make.RLL} is a special function for Rotated LatLon using the HIRLAM conventions, which require location of the rotated South Pole.}
\examples{
## Not run:
\dontrun{
newdom <- Make.domain("lambert",clonlat=c(50.5,3.5),nxny=c(128,128),dxdy=c(8500,8500),reflon=50.5,reflat=3.5)
plot(newdom,add=FALSE)
}
## End(Not run)
}
\keyword{file}

