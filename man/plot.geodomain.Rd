\name{plot.geodomain}
\alias{plot.geodomain}
\alias{plot.geofield}
\title{Plotting a map for a given domain.}
\description{This function plots country borders for a given domain.  }
\usage{
\method{plot}{geodomain}(x=.Last.domain(), add=!is.null(.Last.domain()),
             maplwd=1, mapcol='black',
             add.dx=TRUE, drawmap=!add, box=drawmap,
             map.database='world',...)
\method{plot}{geofield}(x, ...)
}
\arguments{
\item{x}{A an object of \code{geofield} or \code{geodomain} class. Default is the last domain plotted.}
\item{mapcol}{Colour of the map contours. }
\item{map.database}{A character string for the map database. Default is 'world' from the \code{maps} package.}
\item{drawmap}{If TRUE, the map contours are superposed on the data plot.}
\item{maplwd}{Line width of the map.}
\item{add.dx}{If TRUE, the domain size is extended by half a grid interval in all directions, i.e. the co-ordinates are treated as the centers of grid boxes.}
\item{box}{If TRUE, a rectangular box is drawn around the domain.}
\item{add}{If TRUE, the plot is superposed on the current plot.}
\item{...}{Other standard graphical parameters.}
}
\details{These plot methods only plot a map. For \code{geofield} objects there are various other functions that plot the actual data.
}
\examples{
\dontrun{
ttt <- FAdec(...)
### to plot the domain of ttt (not the data itself):
plot(ttt, mapreso=.01, add=FALSE)
}
}
\seealso{
  \code{\link{geofield}},\code{\link{geodomain}}, \code{\link{plot}}, \code{\link{iview}}, \code{\link{domainbox}}
}
\keyword{file}
