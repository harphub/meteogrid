\name{getmask}
\alias{getmask}
\alias{getmap}
\alias{getbox}
\title{Domain map and boundaries.}
\description{
These routines return the map (polygons) for drawing over a plot. Or the inversei.e. the sea becomes the inside of the polygons, e.g. to erase values over the sea.
}
\usage{
getmap(domain=.Last.domain(), interior=TRUE, 
                   fill=FALSE, map.database="world")
getbox(domain=.Last.domain())
getmask(domain=.Last.domain(), map.database="world")
}
\arguments{
\item{domain}{A geodomain}
\item{interior}{Draw interior (political) country borders.}
\item{fill}{Return the map as polygons, not simple lines.}
\item{map.database}{The map database to be used.}
}
\details{
 \code{getmask} depends on the \code{sf} package.
}
\examples{
  \dontrun{
  iview(myfield, fill=TRUE)
  plot(getmask(), col="white", border=1, add=TRUE)
  }
}
\seealso{\code{\link[geogrid]{limage}},\code{\link[geogrid]{vecplot}}
}
\keyword{file}
