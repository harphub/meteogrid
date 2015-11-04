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
iview(x,nlevels=15,color.palette=irainbow,
      title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
      legend=FALSE,mask=NULL,
      drawmap=TRUE, maplwd=.5, mapcol='black', map.database='world', ...)
cview(x,nlevels=15,
      title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
      mask=NULL, add=FALSE,
      drawmap=!add, maplwd=.5, mapcol="black", map.database="world", ...)
fcview(x,nlevels=15,color.palette=irainbow,
       title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
       legend=TRUE,mask=NULL,
       drawmap=TRUE, maplwd=.5, mapcol='black', map.database='world', ...)
vview(U, V, add=FALSE, aspcorrect=TRUE,
      drawmap=TRUE, maplwd=.5, mapcol="black", map.database='world', ...)
}
\arguments{
\item{x}{A \code{geofield} object.}
\item{U,V}{\code{geofield} objects interpreted as U and V components of a wind field.}
\item{nlevels}{Number of (colour) levels. You may also give the levels explicitely.}
\item{color.palette}{Colours for the levels.}
\item{title}{Title for the plot.}
\item{add}{If TRUE, the plot is superposed on the current image,}
\item{mask}{Defines a mask such that all grid points with mask==FALSE are set to NA. The mask may be a matrix with the same dimension as x, an expression involving x, or even a string that is parsed and which uses 'x' to represent the data.}
\item{legend}{If TRUE, a color legend is drawn at the right of the image.}
\item{drawmap}{If TRUE, a world map is superposed on the image.}
\item{maplwd, mapcol}{Line width and colour for the world map.}
\item{map.database}{Which map data to use. Default is 'world' from the \code{maps} package.}
\item{aspcorrect}{In Lat/Lon co-ordinates, the aspect ratio is not 1:1. If this argument is TRUE, a correction is attempted. For other projections, this argument is ignored.}
\item{...}{Further options for \code{limage}, \code{vecplot} or \code{contour}.}
}
\value{
  An image (\code{iview}), contour (\code{cview}) or filled contour
  (\code{fcview}) plot of the data, overlaid by geographical boundaries. \code{vview} plots a vector field (e.g. wind).
}
\details{
  These functions are in fact simple "wrappers" for the more general
  functions \code{limage} etc. The (...) indicates any
  parameter that is accepted by these functions.
}
\examples{
\dontrun{
# mask as a string that is evaluated literally (x is the input data)
iview(x=field,mask="x>0") # all negative values are left blank
# mask as an expression, assuming lsm is a known geofield
iview(x=field,mask=lsm) # all points where lsm=FALSE are left blank
}}
\seealso{\code{\link{limage}},\code{\link{vecplot}},\code{\link[graphics]{contour}}
}
\keyword{file}
