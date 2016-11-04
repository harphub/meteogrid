\name{antipolygon}
\alias{antipolygon}
\title{Colour everything outside a given polygon.}
\description{
This routine does the inverse of 'polygon': it colours everything outside a given polygnon. This may be useful for e.g. blanking out a plot outside certain borders. This routine is quite simple and may not work under all circumstances.}
\usage{
antipolygon(data,xylim=DomainExtent(.Last.domain()), bg='white')
}
\arguments{
\item{data}{A polygon (x and y columns).}
\item{bg}{Colour to use.}
\item{xylim}{Boundaries for the 'antipolygon'. Usually the boundaries of the presentt plot.}
}
\details{
 antipolygon() works by adding the 4 corner points to the polygon, thus 'inverting' it.
}
\examples{ \dontrun{
mask <- project(map("world","Belgium",fill=0,plot=FALSE))
xylim <- DomainExtent(datafield)
antipolygon(mask,xylim=xylim)
} }
\seealso{\code{\link[geogrid]{limage}},\code{\link[geogrid]{vecplot}}
}
\keyword{file}
