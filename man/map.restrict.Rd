\name{map.restrict}
\alias{map.restrict}
\alias{periodicity}
\title{Restrict map to plot domain}
\description{
This function clips a set of map boundaries exactly to the plot domain.
}
\usage{
map.restrict(bxy,xlim,ylim,xperiod=NA_real_,xfrac=0.5,yperiod=NA_real_,yfrac=NA_real_)
periodicity(domain=.Last.domain())
}
\arguments{
  \item{bxy}{Number of points to be clicked.}
  \item{xlim,ylim}{Min/max values of the domain.}
  \item{xperiod, xfrac, yperiod, yfrac}{Uually, only the x axis may be periodic for certain projections. \code{xperiod} is the numerical value of the total map projection (the domain may be smaller). A jump larger than xfrac * xperiod is assumed to be passing over the boundaries and is thus clipped.}
  \item{domain}{A geodomain.}
}
\details{
Usually, a user will never use these functions directly, but they are available just in case. \code{periodicity} is a simple function that guesses the periodicity of a map, while \code{map.restrict} clips a map to the correct boundaries.
}
\examples{
\dontrun{
lalo <- laloclick(2)   # lalo now contains the co-ordinates
points(proj4(lalo),pch='x')   # indicate the points on the present map.
}
}
\seealso{
\code{\link{locator}}
}
