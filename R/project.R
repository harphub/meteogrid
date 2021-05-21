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

project <- function(x, y=NULL, proj=.Last.domain()$projection,
                    inv=FALSE, proj4fix=proj_version()[0]==4)
{
  if (is.null(y) && is.numeric(x) && length(x) == 2) {
    xy <- list(x=x[1], y=x[2])
  } else {
    xy <- xy.coords(x, y)
  }
  if (is.null(proj)) {
    stop("No projection.")
  }
  # projection may be a proj4 string or a list of named (proj4) options
  if (is.list(proj)) {
    proj_string <- proj4.list2str(proj, join=TRUE)
  } else if (length(proj) > 1) {
    proj_string <- paste("+", proj, collapse=" ")
  } else {
    proj_string <- proj
  }
  # FIXME:
  # make sure it's not a *rotated* latlong !!!
  if (grepl("+proj=latlong", proj_string, fixed=TRUE)) {
    cat("Nothing to project!\n")
    return(data.frame(x=xy$x, y=xy$y))
  }
 
  if (proj4fix && inv && grepl("+proj=omerc", proj_string, fixed=TRUE)) {
### to circumvent some bugs [ot]merc inverse in PROJ.4 (versions 4.7 - 4.9) 
### DO NOT RUN THIS WITH LATER VERSION (>= v5.0)!
    if (proj4.str2list(proj_string)$alpha<0 ) xy$x <- -xy$x
    else xy$y <- -xy$y
  }
  result <- mg_project(x=xy$x, y=xy$y, proj_string, inv=inv)

### again the same proj.4 bug:
  if (proj4fix && !inv && grepl("+proj=omerc", proj_string, fixed=TRUE)) {
    if (proj4.str2list(proj_string)$alpha<0 ) result$x <- -result$x
    else result$y <- -result$y
  } else if (proj4fix && inv && grepl("+proj=tmerc", proj_string, fixed=TRUE)) {
    result$y[y<0] <- -result$y[y<0]
  }

  return(result)

}

