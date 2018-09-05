#-------------------------------------------#
# Part of R-package geogrid                 #
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

  xsub <- seq(x1,x2,by=reso)
  ysub <- seq(y1,y2,by=reso)
  subnx <- length(xsub)
  subny <- length(ysub)

  newlalo <- DomainPoints(geo,"lalo")
  newlalo$lon <- newlalo$lon[xsub,ysub]
  newlalo$lat <- newlalo$lat[xsub,ysub]

  gdomain$SW <- c(newlalo$lon[1,1],newlalo$lat[1,1])
  gdomain$NE <- c(newlalo$lon[subnx,subny],newlalo$lat[subnx,subny])

  gdomain$nx <- length(xsub)
  gdomain$ny <- length(ysub)
  gdomain$dx <- gdomain$dx * reso
  gdomain$dy <- gdomain$dy * reso
  cxy <- (project(gdomain$NE, proj=gdomain$projection) + project(gdomain$SW, proj=gdomain$projection)) / 2
  # use as.numeric()' to coerce a data.frame to a vector:
  gdomain$clonlat <- as.numeric(project(cxy, proj=gdomain$projection, inv=TRUE))
  if (inherits(geo,"geofield")) {
    as.geofield(geo[xsub,ysub], domain=gdomain, time=attr(geo,"time"),
             info=c(attr(geo,"info"), extra="SUBFIELD"))
  } else gdomain
}

################

zoomgrid <- function(geo, x, y, zoom=50){
  if (inherits(geo,"geofield")) gdomain <- attr(geo,"domain")
  else if (inherits(geo,"geodomain")) gdomain <- geo
  else stop("subgrid requires a geofield or geodomain as input.")

  if (x-zoom < 1) xmin <- 1
  else if (x+zoom > gdomain$nx) xmin <- gdomain$nx-2*zoom
  else xmin <- x-zoom

  if (y-zoom < 1) ymin <- 1
  else if (y+zoom > gdomain$ny) ymin <- gdomain$ny-2*zoom
  else ymin <- y-zoom
  xmax <- xmin+2*zoom
  ymax <- ymin+2*zoom

  subgrid(geo,xmin,xmax,ymin,ymax,reso=1)
}


