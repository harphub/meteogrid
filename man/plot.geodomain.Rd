\name{plot.geodomain}
\alias{plot.geodomain}
\alias{plot.geofield}
\title{Plotting a map for a given domain.}
\description{This function plots country borders for a given domain.  }
\usage{
plot.geodomain=function(x=.Last.domain,add=T,
             maplwd=1,mapcol='light grey',
             add.dx=T,drawmap=!add, box=drawmap,
             mapreso=0.02,map.database='worldHires',...)
}
\arguments{
\item{domain}{A \code{geodomain} format. Default is the last domain plotted.}
\item{mapcol}{Colour of the map contours. }
\item{mapreso}{Map resolution factor. The full map database is "thinned" by this factor. 0 gives highest resolution. For large domains, higher value (e.g. .1 or more) will improve plotting speed. For backward compatibility, mapreso may also be a character string, in which case it is interpreted as the map database (usually 'world' or 'worldHires'), but this should now be set using map.database .}
\item{map.database}{A character string for the map database. The map package allows "world" (very coarse resolution when for small local domains) and "worldHires". When using a high resolution set, mapreso may have to be increased.}
\item{drawmap}{If TRUE, the map contours are superposed on the data plot.}
\item{maplwd}{Line width of the map.}
\item{add.dx}{If TRUE, the domain size is extended by half a grid interval in all directions, i.e. the co-ordinates are treated as the centers of grid boxes.}
\item{box}{If TRUE, a rectangular box is drawn around the domain.}
\item{...}{Other standard graphical parameters.}
}
\value{
}
\details{There are also shortcuts \code{iview, cview, fcview, vview} that go a little further by interpreting any matrix as a geofield (using the \code{as.geofield} function) and have different colour defaults. 
}
\examples{
ttt=FAdec(...)
### to plot the domain of ttt (not the data itself):
plot(ttt,mapreso=.01)
}
\seealso{
  \code{\link{geofield}},\code{\link{plot}}
}
\keyword{file}
