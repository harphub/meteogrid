\name{vecplot}
\alias{vecplot}
\title{Plotting vector fields}
\description{A simple and primitive function for plotting vector
  fields.}
\usage{
vecplot(U,V,x=1:dim(U)[1],y=1:dim(U)[2],thinx=10,thiny=thinx,
                aspect, add=FALSE, arrowcolor=1, arrowsize=.03,
                maxscale=.9, lwd=.5, rescale, xlab="", ylab="", ...)
}
\arguments{
  \item{U,V}{2 matrices containing the X and Y components of the field.}
  \item{x,y}{2 vectors with the x and y coordinates of the grid points.}
  \item{thinx,thiny}{Thinning of the grid points. Vectors are only drawn
                      every thinx,thiny points.}
  \item{aspect}{Aspect ratio. By default it is set to 1, i.e. the X and Y
    axes have the same units. For e.g. Lat-Lon maps, this is wrong!}
  \item{arrowcolor}{Color of the vectors.}
  \item{arrowsize}{Size of the actual "arrow" part of the vector.}
  \item{maxscale}{Scaling of vector lengths vs grid interval.}
  \item{add}{If TRUE, the vectors are added to an existing plot. }
  \item{lwd}{line width}
  \item{rescale}{Another way to rescale the vector lengths. }
  \item{xlab,ylab}{Labels for the X and Y axes. Default is to leave them clear.}
  \item{...}{Options for 'arrows'.}
 }
 \value{ A vector field plot.
}
\details{
   The vector field is plotted using the "arrow" function. This is quite
   promitive. Also, it means that the vectors are bound to the
   co-ordinate system of the underlying plot. In other words: the X and
   Y axis should have the same units. If this is not the case, the
   aspect ratio (asp) should be adapted manually.
}
\keyword{file}
