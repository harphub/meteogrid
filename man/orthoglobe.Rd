\name{orthoglobe}
\alias{orthoglobe}
\title{Draw a global map in orthorombic projection.}
\description{
This routine simply drwas a world map that may be useful as a background for other data.
}
\usage{
orthoglobe(reflon=0,reflat=90,map.database='world',...)
}
\arguments{
\item{reflon, reflat}{Pole location for the projection.}
\item{map.database}{Map database.}
\item{...}{Options for 'plot' (colour, line type etc.).}
}
\seealso{\code{\link[meteogrid]{limage}},\code{\link[meteogrid]{vecplot}}
}
\keyword{file}
