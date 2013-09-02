\name{seamask}
\alias{seamask}
\alias{as.ocean}
\title{Mask out the ocean}
\description{
These routines allow you to plot the ocean as a polygon (in stead of the land surface), which may be useful for some graphs. This routine is rather experimental and will not work under all circumstances! }

\usage{
seamask(col='white',add.dx=TRUE,drawmap=TRUE,map.database='worldcoast')
as.ocean(data,XYLIM=par('usr'))
}
\arguments{
\item{col}{What colour to give to the ocean. Typically white or blue.}
\item{add.dx}{Does the existing figure have 'dx' added? True for 'iview', FALSE for fcview.}
\item{data}{as.ocean expects a list or data frame with components 'x' and 'y'. These must consist of closed polygons, separated by NA.}
\item{XYLIM}{The limits opf the 'ocean'. This may usually be set to the size of the current plot region (even if the map polygons go beyond).}
}
\value{
  An image (\code{iview}), contour (\code{cview}) or filled contour
  (\code{fcview}) plot of the data, overlaid by geographical boundaries.
}
\details{
  seamask() redefines the plot region. If you resize the plot window afterwards, the result will be completely destroyed! Also, it may be necessary to close and re-open the plot window after using this function.

}
\seealso{\code{\link[geogrid]{limage.geofield}},\code{\link[geogrid]{vecplot}}
}
\examples{
}
\keyword{file}
