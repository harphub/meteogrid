\name{domainbox}
\alias{domainbox}
\title{Draw contour of a geodomain}
\description{
Function for adding the outline of a geographical domain to the current plot.
}
\usage{
domainbox(domain, add.dx=TRUE, ...)
}
\arguments{
\item{domain}{A \code{geodomain} object. Or any object (like \code{geofield}) thta has a domain attribute.}
\item{add.dx}{Usually, the co-ordinates are the centers of the grid boxes, not the corners. So we add dx/2 to get the full plot window. Set to FALSE in other cases.}
}
\value{
The outline of the domain is added to a pre-existing plot. If the projection is the same, this may be a square, but otherwise it is re-projected.
}
\seealso{\code{\link[geogrid]{iview}}
}
\keyword{file}
