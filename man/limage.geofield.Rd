\name{limage.geofield}
\alias{limage.geofield}
\title{Visualisation of gridded geographic data.}
\description{These are class functions for visualising \code{geofield}
  objects. They are derived from the default functions, but add geographical information from the map database. \code{iview} and \code{fcview} are just simple shortcuts for calling \code{limage} with some special settings. }
\usage{
limage.geofield(field,smooth=FALSE,drawmap = TRUE, maplwd = 0.5,
      mapcol = "black", mapreso = 0.02,
      map.database="worldNew", map.interior=TRUE, ...)
}
\arguments{
\item{field}{A field in \code{geofield} format.}
\item{mapcol}{Colour of the map contours. }
\item{mapreso}{Map resolution factor. The full map database is "thinned" by this factor. 0 gives highest resolution. For large domains, higher value (e.g. .1 or more) will improve plotting speed. For backward compatibility, mapreso may also be a character string, in which case it is interpreted as the map database (usually 'world' or 'worldHires'), but this should now be set using map.database .}
\item{map.database}{A character string for the map database. The default is "mapNew" which requires the library \code{mapNew} to be installed. This is an updated world map (e.g. European borders are updated) with pretty good resolution. The \code{map} package also allows "world" (very coarse resolution for local domains) and "worldHires". Both of these are howerver pretty old (1990). When using a high resolution set, mapreso may have to be increased.}
\item{drawmap}{If TRUE, the map contours are superposed on the data plot.}
\item{maplwd}{Line width of the map.}
\item{smooth}{If FALSE, the image consists of the actual grid boxes (cfr. image), if TRUE the values are smoothed (cfr. filled.contour). Note that in the first case, the domain size is extended by half a grid interval in all directions!}

\item{...}{Other graphical parameters for limage. These function don't allow passing extra graphical parameters for the map boundaries. For these, you must call plot.geodomain directly. }
}
\value{
}
\details{There are also shortcuts \code{iview, cview, fcview, vview} that go a little further by interpreting any matrix as a geofield (using the \code{as.geofield} function) and have different colour defaults. 
}
\examples{
fcview(ttt,legend=T,legend.title=title(main="meter",cex.main=.6,line=1),legend.width=1/7)
}
\seealso{
  \code{\link[base]{image}},\link{iview},\code{\link[base]{filled.contour}},\code{\link{limage}},\code{\link{plot.geodomain}}
}
\keyword{file}
