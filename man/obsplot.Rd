\name{obsplot}
\alias{obsplot}
\title{Plot point values on a map with a colour function.}
\description{
Add some coloured markers for observations on an existing map.
}
\usage{
obsplot(x,y,z,breaks=5,pretty=TRUE,legend.pos=NULL,
        add=TRUE,domain=.Last.domain(),col=irainbow,...)
}
\arguments{
  \item{x,y}{Vectors with co-ordinates. If y is missing, x is suposed to have two columns x and y.}
  \item{z}{Vector with data (same length as x and y).}
  \item{breaks}{The number of breaks in the colour code.}
  \item{pretty}{If TRUE, the breaks are chosen at 'nice' values.}
  \item{legend.pos}{Legend position.}
  \item{add}{If TRUE, add to current plot. Otherwise, start a new plot with dmain.}
  \item{domain}{By default, the last plotted domain.}
  \item{col}{Colour scheme.}
  \item{...}{Options for 'points'.}
}
\value{
  A set of points is plotted at the x,y points, with colour taken from the chosen colour scheme. The return value (hidden) is a list wit hthe projected x,y points, the data values, levels and colours.
}
\seealso{
  \code{\link{DomainPoints}},\code{\link{lalopoint}}
  }
