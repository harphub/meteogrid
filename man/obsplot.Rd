\name{obsplot}
\alias{obsplot}
\title{Plot point values on a map with a colour function.}
\description{
Add some coloured markers for observations on an existing map.
}
\usage{
obsplot(x,y,data,breaks=5,pretty=TRUE,legend.pos=NULL,
        add=TRUE,domain=.Last.domain,col=irainbow,...)
}
\arguments{
  \item{x,y}{Vectors with co-ordinates.}
  \item{data}{Vector with data (same length as x and y).}
}
\value{
}
\seealso{
  \code{\link{DomainPoints}},\code{\link{lalopoint}}
  }
