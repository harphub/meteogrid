\name{regrid}
\alias{regrid}
\title{Interpolating to another domain}
\description{Using bilinear spatial interpolation from the fields library, this
  function transforms a geofield to another geometry. }
\usage{
regrid(infield, newdomain=.Last.domain)
}
\arguments{
  \item{infield}{A \code{geofield} object.}
  \item{newdomain}{A \code{geodomain} or a \code{geofield}. In the second case, the new domain attribute is taken from this second field. The default value is the \code{.Last.domain}, the last domain that has been plotted. }
#  \item{...}{Options for \code{interp}.}

 }
\value{A \code{geofield} object.
}
\details{ This routine uses bilinear interpolation to transform the data
  to a new domain (for instance, different grid resolution, different projection...).
}
\examples{
regrid(field1,field2) # this will output the data of field1 interpolated to the grid of field 2

iview(field1)
cview(regrid(field2),add=TRUE) # add contours of field2, interpolated to the last plotted grid, which is field1.
}
\seealso{\code{\link[geogrid]{zoomgrid}},\code{\link[geogrid]{subgrid}},\code{\link[geogrid]{Make.domain}},\code{\link[geogrid]{point.interp}}
}
\keyword{file}
