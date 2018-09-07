\name{proj4.str2list}
\alias{proj4.str2list}
\alias{proj4.list2str}
\title{Change proj4 strings to list and back.}
\description{
These routines reformat a proj4 projection string to a list of arguments or vice versa.}
\usage{
proj4.str2list(pp)
proj4.list2str(pp) 
}
\arguments{
\item{pp}{A proj4 projection, either a single character string or a list of named arguments.}
}
\details{
proj.4 strings are almost ubiquitous as definitions of projections. Internally, \code{meteogrid} works with a list of arguments rather than a single character strings. These functions provide simple conversion between the two formats.
}
\seealso{\code{\link[meteogrid]{project}}}
\keyword{file}
