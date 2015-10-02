\name{DrawLatLon}
\alias{DrawLatLon}
\title{Add a grid of meridians and latitude lines.}
\description{}
\usage{
DrawLatLon(nx=9, ny=9, labels = TRUE, pretty = TRUE, cex = 1, 
            col = "blue", lty = 2, font = 2, ...)
}
\arguments{
  \item{nx,ny}{The (approximate) number of lines to be drawn.}
  \item{labels}{If TRUE, the latitude and longitude values are printed at the bottm and left of the map.}
  \item{pretty}{If TRUE, "nice" values are used.}
  \item{...}{Graphical parameters (cex=text size, lty=line type) .}
}
\value{Lines of constant latitude and longitude are added to the current map. The lines are at regular intervals with appropriate values to get aproximately the wanted number of lines.
}
\seealso{\code{DomainPoints}
}
\keyword{file}
