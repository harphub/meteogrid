\name{iview}
\alias{iview}
\alias{cview}
\alias{fcview}
\alias{vview}
\title{Visualisation shortcuts}
\description{
Routines for visualisation of geogrid objects. They use the geogrid
methods of image, contour, fcontour and vecplot.
}
\usage{
iview(x,nlevels=15,col=irainbow(nlevels),title=paste(x$info$name,"\n",x$time),...)
cview(x,nlevels=15,title=paste(x$info$name,"\n",x$time),...)
fcview(x,nlevels=15,color.palette=irainbow,title=paste(x$info$name,"\n",x$time),...)
vview(U,V,,...)
}
\arguments{
\item{record}{A \code{geofield} object as returned by \code{FAdec}. It
  may in some cases also be a plain data matrix, if it can be
  interpreted by the user-supplied function \code{as.geofield}.}
\item{nlevels}{Number of levels.}
\item{palette}{Colours for the levels.}
\item{geocol}{Colour of geogaphical boundaries}
}
\value{
  An image (\code{iview}), contour (\code{cview}) or filled contour
  (\code{fcview}) plot of the data, overlaid by geographical boundaries.
}
\details{
  These functions are in fact simple "wrappers" for the more general
  class functions \code{image.geofield} etc. The (...) indicates any
  parameter that is accepted by these functions.
}
\seealso{\code{\link[geogrid]{limage}},\code{\link[geogrid]{limage.geofield}},\code{\link[geogrid]{vecplot}}
}
\keyword{file}
