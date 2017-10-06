#-------------------------------------------#
# Part of R-package geogrid                 #
# Copyright (c) 2003-2016 Alex Deckmyn      #
#   Royal Meteorological Institute, Belgium #
# alex.deckmyn@meteo.be                     #
# Released under GPL-3 license              #
#-------------------------------------------#

######################
### geogrid
######################

is.geofield <- function(x){
  inherits(x,"geofield")
}

is.geodomain <- function(x){
  inherits(x,"geodomain")
}

print.geofield <- function(x, ...){
  cat(paste(attr(x,"info")$origin,":",attr(x,"info")$name),"\n")
  cat("Time:\n")
  cat(attr(x,"time"),"\n")
  cat("Domain summary:\n")
  print(attr(x,"domain"))
  cat("Data summary:\n")
  cat(summary(as.vector(x)),"\n")
}

print.geodomain = function(x, ...){
  cat(x$nx,"x",x$ny,"domain\n")
  cat("Projection summary:\n")
  cat("proj=",x$projection$proj,"\n")
  if (!is.null(x$NE)) cat("NE = (",x$NE[1],",",x$NE[2],")\n")
  if (!is.null(x$SW)) cat("SW = (",x$SW[1],",",x$SW[2],")\n")
  if (!is.null(x$clonlat)) cat("centre = (",x$clonlat[1],",",x$clonlat[2],")\n")
#  print.noquote(domain$projection)
}

as.geofield <- function (x=NA, domain, time = "", info = list()) {
  if (is.geofield(x)) return(x)
  if (is.vector(x)) x <- array(x,dim=c(domain$nx,domain$ny))
  attr(x,"domain") <- domain
  attr(x,"time") <- time
  attr(x,"info") <- info
  class(x) <- c("geofield",class(x))
  return(x)
}

compare.geodomain <- function(x, y ,eps=1e-10){
### TRUE if they are equal, FALSE if they are not
### this is not exhaustive, but I am in a hurry!
  OK <- (x$nx == y$nx) &  (x$ny == y$ny) &
  (max(abs(x$NE-y$NE)) < eps) & (max(abs(x$SW-y$SW)) < eps)
  if (!OK) return(FALSE)
### compare the projection: may have a different order!
  n1 <- names(x$projection)
  n2 <- names(y$projection)
  ndif <- setdiff(n1,n2)
  if (length(ndif) > 0) return(FALSE)
  for(nn in n1){
    if (is.character(x[[nn]])) OK <- (x[[nn]]==y[[nn]])
    else OK <- (abs(x[[nn]]-y[[nn]])<eps)
    if (!OK) return (FALSE)
  }
  return(TRUE)
}

#####################################

DomainExtent <- function(geo){
### We look for the extreme LatLon values of a domain domain
### these are useful to select the right part from the world map
### first find the SW and NE point, then draw a rectangular box
### and project it back to LatLon
### You could use DomainPoints in stead, but then you'd be projecting
### all domain points, which is a bit of an overkill...

  if (!inherits(geo,"geodomain")) geo <- attr(geo,"domain")

  if (!is.null(geo$clonlat) && !is.null(geo$dx) && !is.null(geo$dy)) {
    dx <- geo$dx
    dy <- geo$dy
    clonlat <- geo$clonlat
    cxy <- project(x=geo$clonlat[1], y=geo$clonlat[2], proj=geo$projection)
    x0  <- cxy$x - dx*(geo$nx -1)/2
    x1  <- cxy$x + dx*(geo$nx -1)/2
    y0  <- cxy$y - dy*(geo$ny -1)/2
    y1  <- cxy$y + dy*(geo$ny -1)/2
  } else {
    xy <- project(list(x=c(geo$SW[1], geo$NE[1]),y=c(geo$SW[2], geo$NE[2])),
                proj=geo$projection)
    x0 <- xy$x[1]
    y0 <- xy$y[1]
    x1 <- xy$x[2]
    y1 <- xy$y[2]
    dx <- (x1-x0)/(geo$nx-1)
    dy <- (y1-y0)/(geo$ny-1)
    xc <- (x0+x1)/2
    yc <- (y0+y1)/2
    cll <- project(list(x=xc,y=yc),proj=geo$projection,inv=TRUE)
    clonlat <- c(cll$x, cll$y)
  }
  borders <- project(list(x=c(seq(x0,x1,length=geo$nx),rep(x1,geo$ny),
                              seq(x0,x1,length=geo$nx),rep(x0,geo$ny)),
                          y=c(rep(y0,geo$nx),seq(y0,y1,length=geo$ny),
                              rep(y1,geo$nx),seq(y0,y1,length=geo$ny))),
                      proj=geo$projection,inv=TRUE)
### ATTENTION: if the map crosses or comes close to the date line (meridian 180/-180) this will not
### work correctly. In the map command we must then set lonlim=NULL !!!

  lonlim <- range(borders$x, na.rm=TRUE)
  llr <- lonlim[2]-lonlim[1]
  if (llr <= 1 | llr > 300) lonlim <- NULL

  list(lonlim = lonlim , latlim = range(borders$y, na.rm=TRUE),
       clonlat=clonlat,
       x0=x0, y0=y0, x1=x1, y1=y1, dx=dx, dy=dy, nx=geo$nx, ny=geo$ny)
}

##########################################################
DomainPoints <- function (geo,type="lalo"){
### return lat's and lon's of all domain points (or leave in projection if type = "xy")
  
  if (!inherits(geo,"geodomain")) geo <- attr(geo,"domain")

  if (!is.null(geo$clonlat) && !is.null(geo$dx) && !is.null(geo$dy)) {
    cxy <- project(x=geo$clonlat[1], y=geo$clonlat[2], proj=geo$projection)
    xy <- data.frame(x= cxy$x + c(-1, +1)*dx*(geo$nx-1)/2 ,
                     y= cxy$y + c(-1, +1)*dy*(geo$yx-1)/2) 
  } else {
    lalo <- list(x=c(geo$SW[1],geo$NE[1]),y=c(geo$SW[2],geo$NE[2]))
    xy <- project(lalo, proj = geo$projection)
  }
  xydomain <- expand.grid(x = seq(xy$x[1], xy$x[2], length = geo$nx),
                          y = seq(xy$y[1], xy$y[2], length = geo$ny))
  if (type == "lalo") {
    lalolist <- project(xydomain, proj = geo$projection, inv = TRUE)
    list(lon = matrix(lalolist$x, ncol = geo$ny, nrow = geo$nx),
         lat = matrix(lalolist$y, ncol = geo$ny, nrow = geo$nx))
  }
  else if (type == "xy") {
    list(x=matrix(xydomain$x, ncol = geo$ny, nrow = geo$nx),
         y=matrix(xydomain$y, ncol = geo$ny, nrow = geo$nx))
  } else {
    print("Unknown type.")
  }
}

gridpoints <- function(x ,y=NULL, ...){
### indicate particular grid points on a map (like points, but with grid indices)
  AllCoords <- DomainPoints(.Last.domain(),type="xy")
  xy <- xy.coords(x,y)

  xlist <- AllCoords$x[cbind(xy$x,xy$y)]
  ylist <- AllCoords$y[cbind(xy$x,xy$y)]

  points(xlist,ylist,...)
}

###########################################################



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

  if (inherits(geo,"geofield")) {
    as.geofield(geo[xsub,ysub], domain=gdomain,time=attr(geo,"time"),
             info=c(attr(geo,"info"),extra="SUBFIELD"))
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


###########################
### CREATE NEW DOMAINS  ###
###########################

Make.domain <- function(projtype="lambert", clonlat, nxny, dxdy, exey=NULL,
                        reflat=clonlat[2], reflon=clonlat[1], tilt=0,
                        earth=list(R=6371229)){
#                        earth=list(a=6371229.0, es=0.0)){
  if (length(dxdy)==1) dxdy <- rep(dxdy,2)
  if (projtype %in% c("lcc", "lambert")) {
### Lambert (as used in ALADIN: only 1 reference latitude)
    projection <- list(proj="lcc", lon_0=reflon, lat_1=reflat, lat_2=reflat)
  } else if (projtype %in% c("merc", "omerc", "tmerc", "mercator")){
### Rotated & tilted Mercator
### as used in ALADIN: reflon=clon, reflat=clat !!!
### we simply *ignore* reflon and reflat (coming from FA file they are 0 !)
#    if (reflat != clonlat[2] || reflon != clonlat[1]) warning('This domain is not ALADIN-compatible!')
    rlon <- clonlat[1]
    rlat <- clonlat[2]
    if (abs(rlat) < 0.01 && abs(tilt) < 0.01) {
# not necessary: this is OK in "omerc"
      projection <- list(proj="merc",lon_0=rlon)
    } else if (abs(abs(tilt)-90) < 1.0E-5) {
# this would crash in "omerc" (alpha ~ 0)
      projection <- c(proj = "tmerc", lon_0 = rlon, lat_0 = rlat)
    } else if (abs(tilt) < 1.E-5) {
      projection <- c(proj = "somerc", lon_0 = rlon, lat_0 = rlat)
    } else if (tilt>0) {
      projection <- list(proj = "omerc", lonc = rlon,
                            lat_0 = rlat, alpha = -90 + tilt, no_rot=NA)
    } else {
      projection <- list(proj = "omerc", lonc = rlon,
                            lat_0 = rlat, alpha = 90 + tilt, no_rot=NA)
    }
  } else if (projtype %in% c("latlong")) {
    projection=list(proj="latlong")
  } else if (projtype %in% c("ob_tran", "RotLatLon")){
    projection <- list(proj="ob_tran","o_proj"="latlong",
                       "o_lat_p"=-reflat,"o_lon_p"=0,"lon_0"=reflon)
  } else if (projtype %in% c("stere", "stereographic")) {
    projection <- list(proj="stere",lon_0=reflon,lat_0=reflat)
  } else {
    stop("Unknown projection.")
  }
  projection <- c(projection, earth)
### project the centre point

#  cxy <- project(list(x=clonlat[1], y=clonlat[2]), proj=projection)

#  SW0 <- c(cxy$x,cxy$y) - dxdy*(nxny - 1)/2
#  NE0 <- SW0 + dxdy*(nxny - 1)

#  lims <- project(list(x=c(SW0[1],NE0[1]),y=c(SW0[2],NE0[2])), proj=projection, inv=TRUE)
#  SW <- c(lims$x[1],lims$y[1])
#  NE <- c(lims$x[2],lims$y[2])
### and the output is...
  result <- list(projection=projection, nx=nxny[1], ny=nxny[2], dx=dxdy[1], dy=dxdy[2],
                 clonlat=clonlat)
  if (!is.null(exey)) result <- c(result, ex=exey[1], ey=exey[2])
  class(result) <- "geodomain"
  result
}

Make.domain.RLL <- function(Lon1,Lat1,SPlon,SPlat,SPangle=0,nxny,dxdy){
### This is for Rotated LatLon as used by Hirlam: central meridian is vertical.
### In the future, this should be merged with Make.domain
### but that is not trivial: here, the centre point is of no real consequence
### and you have to define the rotated South Pole
  if (length(dxdy)==1) dxdy <- rep(dxdy,2)
  Lon2 <- Lon1 + (nxny[1]-1)*dxdy[1]
  Lat2 <- Lat1 + (nxny[2]-1)*dxdy[2]
  projection <- list(proj="ob_tran","o_proj"="latlong",
                       "o_lat_p"=-SPlat,"o_lon_p"=0,"lon_0"=SPlon)

  RR <- project(list(x=c(Lon1,Lon2)*pi/180,y=c(Lat1,Lat2)*pi/180),proj=projection,inv=TRUE)
  SW <- c(RR$x[1],RR$y[1])
  NE <- c(RR$x[2],RR$y[2])

  result <- list(projection=projection,nx=nxny[1],ny=nxny[2],SW=SW,NE=NE)
  class(result) <- "geodomain"
  result
}

