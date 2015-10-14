\name{regrid}
\alias{regrid}
\alias{regrid.init}
\title{Interpolating to another domain}
\description{Using spatial interpolation this
  function transforms a geofield to another geometry. }
\usage{
regrid(infield, newdomain=.Last.domain(),method="bilin",mask=NULL,weights=NULL)
regrid.init(olddomain, newdomain=.Last.domain(),method="bilin", mask=NULL)
}
\arguments{
  \item{infield}{A \code{geofield} object.}
  \item{olddomain}{A \code{geodomain} object. It may aso be a \code{geofield}, but only the domain attribute is used for initialising the interpolation.}
  \item{newdomain}{A \code{geodomain} or a \code{geofield}. In the second case, the new domain attribute is taken from this second field. The default value is the \code{.Last.domain}, the last domain that has been plotted. }
#  \item{...}{Options for \code{interp}.}
  \item{method}{The interpolation method can be bilinear ("bilin"), bicubic ("bicub") or nearest neighbour ("n" or "c").}
  \item{mask}{A mask (grid of same dimension as infield), points with value 0 (or FALSE) are not used for interpolation.}
  \item{weights}{Output from \code{regrid.init}. By providing such weights, you can save a lot of time if many fields need interpolation.}
 }
\value{A \code{geofield} object.
}
\details{ This routine uses interpolation to transform the data
  to a new domain (for instance, different grid resolution, different projection...).
}
\examples{
\dontrun{
regrid(field1,field2) # this will output the data of field1 interpolated to the grid of field 2

iview(field1)
cview(regrid(field2),add=TRUE) # add contours of field2, interpolated to the last plotted grid, which is field1.
}
}
\seealso{\code{\link[geogrid]{zoomgrid}},\code{\link[geogrid]{subgrid}},\code{\link[geogrid]{Make.domain}},\code{\link[geogrid]{point.interp}}
}
\keyword{file}
