\name{upscale}
\alias{upscale}
\alias{upscale_factor}
\alias{upscale_regrid}
\alias{upscale_regrid_init}
\title{Aggregating to another, coarser domain}
\description{Change resolution/geometry of a geofield by averaging. }
\usage{
upscale(infield, factor=NULL, newdomain=NULL, method="mean", weights=0, ...)
upscale_factor(infield, factor, method="mean", ...)
upscale_regrid(infield, newdomain, method="mean", weights=NULL, ...)
upscale_regrid_init(olddomain, newdomain)
}
\arguments{
  \item{infield}{A \code{geofield} object.}
  \item{factor}{If given, the new grid is simply a coarser version of the original. \code{factor} should be a vector of 2 integers (if it is scalar, the value is repeated). It dscribes how many gridd cells in the X and Y direction are aggregated to one new grid box.}
  \item{newdomain, olddomain}{A \code{geodomain} or a \code{geofield}. In the second case only the domain attribute is used.}
  \item{method}{The aggregation method. When aggregating by a given factor, this could be a function like \code{mean}, \code{median}, \code{max}... For ergridding, only 'mean' is supported.}
  \item{weights}{Initialised weights for faster computation, output from \code{upscale_regrid_init}. Only used for regridding (factor=NULL).}
  \item{...}{Options for the method (e.g. \code{na.rm=TRUE}).}
 }
\value{A \code{geofield} object.
}
\details{ These routines aggregate all grid boxes of the original domain into the new domain. When a factor is used, the result is exactly as expected. When regridding to a different domain, the average is taken of all original cells whose center falls inside the new grid cell. There is no compensation (yet?) for partial grid cells.
}
\seealso{\code{\link[geogrid]{zoomgrid}},\code{\link[geogrid]{subgrid}},\code{\link[geogrid]{regrid}}
}
\keyword{file}
