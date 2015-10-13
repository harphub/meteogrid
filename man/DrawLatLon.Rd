\name{DrawLatLon}
\alias{DrawLatLon}
\title{Add a grid of meridians and latitude lines.}
\description{This function adds a grid of lat/lon lines to an existing plot.
}
\usage{
DrawLatLon(nx=9, ny=9, lines = TRUE, labels = TRUE, cex = 1, 
            col = "grey", lty = 2, font = 2, lab.x = 2, lab.y = 2, ...)
}
\arguments{
  \item{nx,ny}{The (approximate) number of lines to be drawn.}
  \item{lines, labels}{If labels=TRUE, the latitude and longitude values are printed at the bottom and left of the map. If lines=TRUE, a grid of lines is drawn on the plot. If both are FALSE nothing happens, so what's the point...}
  \item{lty, cex, font, lab.x, lab.y}{Graphical parameters.}
  \item{...}{More graphical parameters.}
}
\value{Lines of constant latitude and longitude are added to the current map. The lines are at regular intervals with appropriate values to get aproximately the wanted number of lines.
}
\seealso{\code{DomainPoints}
}
\keyword{file}
