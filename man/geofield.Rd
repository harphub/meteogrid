\name{geofield}
\alias{as.geofield}
\alias{is.geofield}
\alias{geofield}
\title{Geofield data format}
\description{}
\usage{
is.geofield(x)
as.geofield(x, grid, time = "", info = "", domain) 
}
\arguments{
  \item{x}{a data matrix}
}
\value{
\item{is.geofield}{returns TRUE if x is a geofield object.}    
\item{as.geofield}{returns a geofield object. If x is a matrix and a domain is specified, this is added as an attribute. If no domain is specified, the function makes an educated guess (e.g. a 300x300 matrix defaults to the Aladin-France domain). This function can be edited for local applications.
}
}
\details{The geofield format consists a data matrix and some attributes that describe the domain and the field value. An object is a data matrix with the
  following attributes:
  \itemize{
  \item{domain}{The domain description. This is itself a list with elements:
    \code{projection}: Details of the grid's projection.\\
    \code{nx,ny}: Number of gridpoints
    \code{SW,NE}: (Latitude,Longitude) pairs describing the outer cornes of the grid.
  }
  \item{time}{Date and time information of the data field.}
  \item{info}{Other information.}
  }
}
\examples{
}
\seealso{
}
\keyword{file}
