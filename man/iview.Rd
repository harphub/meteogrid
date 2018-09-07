\name{iview}
\alias{iview}
\alias{cview}
\alias{fcview}
\alias{vview}
\title{Visualisation shortcuts}
\description{
Routines for visualisation of meteogrid objects. They use the meteogrid
methods of image, contour, fcontour and vecplot.
}
\usage{
iview(x,nlevels=15,color.palette=irainbow,
      levels=pretty(x[1:gdomain$nx,1:gdomain$ny], nlevels),
      col=color.palette(length(levels)-1),
      title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
      legend=FALSE, mask=NULL, na.col=par("bg"),
      drawmap=TRUE, maplwd=.5, mapcol='black', map.database='world',
      interior=TRUE, fill=FALSE, ...)
cview(x,nlevels=15,
      title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
      mask=NULL, add=FALSE,
      drawmap=!add, maplwd=.5, mapcol="black", map.database="world",
      interior=TRUE, fill=FALSE, ...)
fcview(x,nlevels=15,color.palette=irainbow,
       levels=pretty(x[1:gdomain$nx,1:gdomain$ny], nlevels),
       col=color.palette(length(levels)-1),
       title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
       legend=TRUE,mask=NULL,
       drawmap=TRUE, maplwd=.5, mapcol='black', map.database='world',
      interior=TRUE, fill=FALSE, ...)
vview(U, V, add=FALSE, aspcorrect=TRUE,
      drawmap=TRUE, maplwd=.5, mapcol="black", map.database='world',
      interior=TRUE, fill=FALSE, ...)
}
\arguments{
\item{x}{A \code{geofield} object.}
\item{U,V}{\code{geofield} objects interpreted as U and V components of a wind field.}
\item{nlevels}{Number of (colour) levels. You may also give the levels explicitely.}
\item{color.palette}{Colours for the levels.}
\item{levels}{A fixed set of breaks. If \code{levels} is defined, \code{nlevels} is ignored. If missing, the value range is devided in approximately \code{nlevels} breaks.}
\item{col}{A fixed vector of colours. If provided, the length must be exactly equal to \code{length(levels)-1}. So it is only practicle to define this if you also have a fixed set of levels.}
\item{title}{Title for the plot.}
\item{add}{If TRUE, the plot is superposed on the current image,}
\item{mask}{Defines a mask such that all grid points with mask==FALSE are set to NA. The mask may be a matrix with the same dimension as x, an expression involving x, or even a string that is parsed and which uses 'x' to represent the data.}
\item{na.col}{By default, masked grid points (and any other NA values) are transparent. With \code{na.col} you can give them another colour.}
\item{legend}{If TRUE, a color legend is drawn at the right of the image.}
\item{drawmap}{If TRUE, a world map is superposed on the image.}
\item{fill}{If TRUE, the map is drawn as filled polygons, not just lines. This is most interesting for contours and vectors, but may in rare cases be used with images, if the colours are transparent.}
\item{interior}{If FALSE, the interior political borders are not drawn.}
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
