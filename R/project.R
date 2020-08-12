#' Projections
#' 
#' A wrapper for calling PROJ code for standard projections.
#' 
#' This routine calls the proj.4 (or later) library of geographical projections. Longitude
#' and latitude values are always in (decimal) degrees, not radians.
#' 
#' @param x,y Longitude and latitude values (inv=FALSE) or projected values
#' (inv=TRUE). The two vectors must have equal length. NA's are allowed. If y
#' is missing, x is expected to be a vector of length 2 or to have two columns.
#' @param proj A proj4 projection description, given as a list or a string. Default value is the
#'   projection of the last geodomain to be plotted.
#' @param inv TRUE indicates inverse projection
#' @return A data.frame with columns x and y.
#' @keywords file
#' @export project

project <- function(x, y, proj=.Last.domain()$projection, inv=FALSE)
{
  if (missing(y)) {
# apparantly, is.list gives TRUE for data.frames, but let's be careful:
    if (is.list(x) | is.data.frame(x)) {
      if (all(is.element(c("x","y"), names(x)))) {
        y <- x$y
        x <- x$x
    } else if (all(is.element(c("lon","lat"), tolower(substring(names(x), 1, 3))))) {
        y <- x$lat
        x <- x$lon
    } else {
### we assume the first two list elements are longitude & latitude
### whatever the actual names
        y <- as.numeric(x[[2]])
        x <- as.numeric(x[[1]])
      }
      if (is.null(x) | is.null(y)) stop("No valid x and y co-ordinates found.")
    }
    else if (is.vector(x)) {
      y <- x[2]
      x <- x[1]
    } else {
      y <- x[,2]
      x <- x[,1]
    }
  }
  if (length(x) != length(y)) stop("x and y don't have the same length.")
  if (!is.numeric(x) | !is.numeric(y)) stop("Co-ordinates must be numeric values.")
  if (missing(proj)) {
    if(!is.null(.Last.domain())) proj <- .Last.domain()$projection
    else return("No projection.")
  }

  if (proj$proj == "latlong") return(data.frame(x=x,y=y))
### longitude should probably be in the interval [-180,180[
### unless e.g. if my global data is on a globe [0,360[
### we assume that the meridian = MinLon + 180
### so meridian-180 must not be transported.
## BUG: this may introduce wrapping problems
#    meridian <- checknull(proj$lon0, default=0)
#    x[x <  (meridian-180)] <- x[x <  (meridian-180)] + 360
#    x[x >= (meridian+180)] <- x[x >= (meridian+180)] - 360
 
  if (inv && is.list(proj) && proj$proj=='omerc') {
### to circumvent some bugs [ot]merc inverse in PROJ.4 (versions 4.7 - 4.9) 
### If they ever solve this bug, I'll have to change this!
    if (proj$alpha<0 ) x <- -x
    else y <- -y
  }

  if (is.list(proj)) {
    proj_string <- proj4.list2string(proj, join=TRUE)
  } else if (length(proj) > 1) {
    proj_string <- paste("+", proj, collapse=" ")
  } else {
    proj_string <- proj
  }
  result <- mg_project(x=x, y=y, proj_string, inv=inv)

### again the same proj.4 bug:
  if (!inv) {
      if (is.list(proj) && proj$proj=='omerc') {
        if (proj$alpha<0 ) result$x <- -result$x
        else result$y <- -result$y
      }
  } else {
    if (is.list(proj) && proj$proj=="tmerc") {
      result$y[y<0] <- -result$y[y<0]
    }
  }

  return(result)

}

