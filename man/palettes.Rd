\name{palettes}
\alias{col.oro}
\alias{col.precip}
\alias{col.temp}
\alias{irainbow}
\alias{col.cloud}
\alias{grainbow}
\title{Colour paletes for images.}
\description{Some simple colour palettes that may be used for images.}
\usage{
col.oro(n)
}
\arguments{
  \item{n}{Integer: number of colours. irainbow and grainbow accept negative values, in which case the order is reversed.}
}
\value{
All these functions return a vector of n colours that span the chosen colour palette. \code{col.cloud} is a bit special, as it actually returns transparancy values. \code{grainbow} is kind of a stupid name for grey-scaling. \code{irainbow} and \code{grainbow} accept negative values for n, in which case the order is reversed. For best results, some of these functions require RColorBrewer. If this is not available, a similar choice is made from standard colours.}
\keyword{file}
