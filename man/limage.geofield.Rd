\name{limage.geofield}
\alias{limage.geofield}
\title{Visualisation of gridded geographic data.}
\description{These are class functions for visualising \code{geofield}
  objects. They are derived from the default functions, but add geographical information from the map database. \code{iview} and \code{fcview} are just simple shortcuts for calling \code{limage} with some special settings. }
\usage{
\method{limage}{geofield}(x,smooth=FALSE,drawmap = TRUE, maplwd = 0.5,
      mapcol = "black", 
      map.database="world", ...)
}
\arguments{
\item{x}{A field in \code{geofield} format.}
\item{mapcol}{Colour of the map contours. }
\item{map.database}{A character string for the map database. The default is "world".}
\item{drawmap}{If TRUE, the map contours are superposed on the data plot.}
\item{maplwd}{Line width of the map.}
\item{smooth}{If FALSE, the image consists of the actual grid boxes (cfr. image), if TRUE the values are smoothed (cfr. filled.contour). Note that in the first case, the domain size is extended by half a grid interval in all directions!}

\item{...}{Other graphical parameters for limage. These function don't allow passing extra graphical parameters for the map boundaries. For these, you must call plot.geodomain directly. }
}
\details{There are also shortcuts \code{iview, cview, fcview, vview} that go a little further by interpreting any matrix as a geofield (using the \code{as.geofield} function) and have different colour defaults. 
}
\examples{
\dontrun{
fcview(ttt,legend=T,legend.title=title(main="meter",cex.main=.6,line=1),legend.width=1/7)
}
}
\seealso{
  \code{\link[graphics]{image}},\link{iview},\code{\link[graphics]{filled.contour}},\code{\link{limage}},\code{\link{plot.geodomain}}
}
\keyword{file}
