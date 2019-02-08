#-------------------------------------------#
# Part of R-package meteogrid                 #
# Copyright (c) 2003-2017 Alex Deckmyn      #
#   Royal Meteorological Institute, Belgium #
# alex.deckmyn@meteo.be                     #
# Released under GPL-3 license              #
#-------------------------------------------#

#########################
### SUBGRID           ###
#########################

subgrid <- function(geo, x1, x2, y1, y2, reso=1) {
  if (inherits(geo,"geofield")) gdomain <- attr(geo,"domain")
  else if (inherits(geo,"geodomain")) gdomain <- geo
  else stop("subgrid requires a geofield or geodomain as input.")

  xsub <- seq(x1, x2, by=reso)
  ysub <- seq(y1, y2, by=reso)
  nx <- length(xsub)
  ny <- length(ysub)

  gxy <- DomainPoints(geo, "xy")
  clonlat <- as.numeric(project(x = (gxy$x[xsub[1], ysub[1]] + gxy$x[xsub[nx], ysub[ny]])/2,
				y = (gxy$y[xsub[1], ysub[1]] + gxy$y[xsub[nx], ysub[ny]])/2,
				proj = gdomain$projection,
				inv = TRUE))
  if (is.null(gdomain$dx)) {
    dx <- (gxy$x[xsub[nx], ysub[ny]] - gxy$x[xsub[1], ysub[1]]) / (nx-1)
  } else {
    dx <- gdomain$dx * reso
  }
  if (is.null(gdomain$dy)) {
    dy <- (gxy$y[xsub[nx], ysub[ny]] - gxy$y[xsub[1], ysub[1]]) / (ny-1)
  } else {
    dy <- gdomain$dy * reso
  }

  newdomain <- list(projection=gdomain$projection, 
                    nx = nx, ny = ny,
		    dx = dx, dy = dy,
		    clonlat = clonlat)
  class(newdomain) <- "geodomain"

  # use as.numeric()' to coerce a data.frame to a vector:
  if (inherits(geo, "geofield")) {
    as.geofield(geo[xsub,ysub], domain=newdomain, time=attr(geo,"time"),
                info=c(attr(geo,"info"), extra="SUBFIELD"))
  } else {
    newdomain
  }
}

################

zoomgrid <- function(geo, x, y, zoom=50){
  if (inherits(geo,"geofield")) gdomain <- attr(geo,"domain")
  else if (inherits(geo,"geodomain")) gdomain <- geo
  else stop("zoomgrid requires a geofield or geodomain as input.")

  if (x-zoom < 1) xmin <- 1
  else if (x+zoom > gdomain$nx) xmin <- gdomain$nx-2*zoom
  else xmin <- x-zoom

  if (y-zoom < 1) ymin <- 1
  else if (y+zoom > gdomain$ny) ymin <- gdomain$ny-2*zoom
  else ymin <- y-zoom
  xmax <- xmin + 2*zoom
  ymax <- ymin + 2*zoom

  subgrid(geo, xmin, xmax, ymin, ymax, reso=1)
}


