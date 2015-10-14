\name{laloclick}
\alias{laloclick}
\title{Click on the map and get the Lat-Lon co-ordinates.}
\description{
  With this command you caqn click n times on a map, and get the Lat-Lon co-ordinates of the points. 
}
\usage{
laloclick(n=1, ...)
}
\arguments{
  \item{n}{Number of points to be clicked.}
  \item{...}{Other arguments to 'locator'.}
}
\value{An matrix dim(n,2) containing longitude and latitude.
}
\examples{
\dontrun{
lalo <- laloclick(2)   # lalo now contains the co-ordinates
points(proj4(lalo),pch='x')   # indicate the points on the present map.
}
}
\seealso{
\code{\link{locator}}
}
