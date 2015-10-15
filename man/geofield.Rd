\name{geofield}
\alias{geofield}
\alias{as.geofield}
\alias{is.geofield}
\alias{print.geofield}
\title{Geofield class format}
\description{The main object class for gridded data. Consists of a matrix with a number of attributes.}
\usage{
as.geofield(x=NA, domain, time = "", info = list()) 
is.geofield(x)
\method{print}{geofield}(x, ...)
}
\arguments{
  \item{x}{a data matrix. If a vecor (or single value), this is cast to a matrix with dimensions required by the domain.}
  \item{domain}{A geodomain that should describe a domain with exactly the dimensions of x.}
  \item{time, info}{Time and a list of information that are added as attributes to the resulting geofield object.}
  \item{...}{Not used.}
}
\value{
\item{is.geofield}{returns TRUE if x is a geofield object.}    
\item{as.geofield}{returns a geofield object. If x is a matrix and a domain is specified, this is added as an attribute.}
}
\details{The geofield class consists a data matrix and some attributes that describe the domain and the field value. An object is a data matrix with the
  following attributes:
  \itemize{
  \item{domain}{The domain description, a \code{geodomain} class object. This is itself a list with elements:
    \code{projection}: Details of the grid's projection.\\
    \code{nx,ny}: Number of gridpoints
    \code{SW,NE}: (Latitude,Longitude) pairs describing the outer cornes of the grid.
  }
  \item{time}{Date and time information of the data field.}
  \item{info}{Other information.}
  }
}
\seealso{\code{\link{geodomain}}
}
\keyword{file}
