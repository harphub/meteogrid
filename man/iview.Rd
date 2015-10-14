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
      legend=TRUE,mask=NULL,...)
cview(x,nlevels=15,
      title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
      mask=NULL,mapcol="black",add=FALSE,drawmap=!add,maplwd=.5,
      map.database="world",...)
fcview(x,nlevels=15,color.palette=irainbow,
       title=paste(attr(x,"info")$name,"\n",attr(x,"time")),
       legend=TRUE,mask=NULL,...)
vview(U,V,,...)
}
\arguments{
\item{x}{A \code{geofield} object.}
\item{nlevels}{Number of (colour) levels. You may also give the levels explicitely.}
\item{color.palette}{Colours for the levels.}
\item{title}{Title for the plot.}
\item{mask}{A masking function, may be an expression or a string that will be parsed.}
\item{legend}{If TRUE, a color legend is drawn at the right of the image.}
\item{...}{Further options for \code{limage}, \code{vecplot} or \code{contour}.}
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
\examples{
\dontrun{
# mask as a string that is evaluated literally (x is the input data)
iview(x=field,mask="x>0") # all positive values are left blank
# mask as an expression, assuming lsm is a known geofield
iview(x=field,mask= !lsm) # all points where lsm=FALSE are left blank
}}
\seealso{\code{\link[geogrid]{limage}},\code{\link[geogrid]{limage.geofield}},\code{\link[geogrid]{vecplot}}
}
\keyword{file}
