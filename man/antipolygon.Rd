\name{antipolygon}
\alias{antipolygon}
\title{Colour everything outside a given polygon.}
\description{
This routine does the inverse of 'polygon': it colours everything outside a given polygnon. This may be useful for e.g. blanking out a plot outside certain borders. This routine is quite simple and may not work under all circumstances.}
\usage{
antipolygon(data,col='white',xylim=DomainExtent(.Last.domain()))
}
\arguments{
\item{data}{A polygon (x and y columns).}
\item{col}{Colout to use.}
\item{xylim}{Boundaries for the 'antipolygon'. Usually the boundaries of the presentt plot.}
}
\details{
 antipolygon() works by adding the 4 corner points to the polygon, thus 'inverting' it.
}
\seealso{\code{\link[geogrid]{limage}},\code{\link[geogrid]{vecplot}}
}
\keyword{file}
