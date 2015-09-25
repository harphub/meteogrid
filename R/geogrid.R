#--------------------------------------#
# Part of R-package geogrid            #
# Copyright Alex Deckmyn               #
# Released under GPL-3 license         #
#--------------------------------------#

######################
### geogrid
######################

is.geofield <- function(x){
  inherits(x,"geofield")
}

is.geodomain <- function(x){
  inherits(x,"geodomain")
}

print.geofield <- function(x,...){
  cat(paste(attr(x,"info")$origin,":",attr(x,"info")$name),"\n")
  cat("Time:\n")
  cat(attr(x,"time"),"\n")
  cat("Domain summary:\n")
  print(attr(x,"domain"))
  cat("Data summary:\n")
  cat(summary(as.vector(x)),"\n")
}

print.geodomain = function(x,...){
  cat(x$nx,"x",x$ny,"domain\n")
  cat("Projection summary:\n")
  cat("proj=",x$projection$proj,"\n")
  cat("NE = (",x$NE[1],",",x$NE[2],")\n")
  cat("SW = (",x$SW[1],",",x$SW[2],")\n")
#  print.noquote(domain$projection)
}

compare.geodomain <- function(domain1,domain2,eps=1e-10){
### TRUE if they are equal, FALSE if they are not
### this is not exhaustive, but I am in a hurry!
  (domain1$nx == domain2$nx) &
  (domain1$ny == domain2$ny) &
  (domain1$projection$proj == domain2$projection$proj) &
  (max(abs(domain1$NE-domain2$NE)) < eps) &
  (max(abs(domain1$SW-domain2$SW)) < eps)
}

#####################################

DomainExtent <- function(geo,...){
  UseMethod("DomainExtent")
}

DomainExtent.geofield <- function(geo,...){
  DomainExtent(attr(geo,"domain"))
}

DomainExtent.geodomain <- function(geo,...){
### We look for the extreme LatLon values of a domain domain
### these are useful to select the right part from the world map
### first find the SW and NE point, then draw a rectangular box
### and project it back to LatLon
### You could use DomainPoints in stead, but then you'd be projecting
### all domain points, which is a bit of an overkill...
###  if(is.geofield(domain)) domain <- attr(domain,"domain")
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

  clonlat <- project(list(x=xc,y=yc),proj=geo$projection,inv=TRUE)
  borders <- project(list(x=c(seq(x0,x1,length=geo$nx),rep(x1,geo$ny),
                              seq(x0,x1,length=geo$nx),rep(x0,geo$ny)),
                          y=c(rep(y0,geo$nx),seq(y0,y1,length=geo$ny),
                              rep(y1,geo$nx),seq(y0,y1,length=geo$ny))),
                      proj=geo$projection,inv=TRUE)
### ATTENTION: if the map crosses or comes close to the date line (meridian 180/-180) this will not
### work correctly. In the map command we must then set lonlim=NULL !!!

  lonlim <- range(borders$x,na.rm=TRUE)
  llr <- lonlim[2]-lonlim[1]
  if (llr <= 1 | llr > 300) lonlim <- NULL

  list(lonlim = lonlim , latlim = range(borders$y,na.rm=TRUE),
       clonlat=c(clonlat$x,clonlat$y),
       x0=x0,y0=y0,x1=x1,y1=y1,dx=dx,dy=dy,nx=geo$nx,ny=geo$ny)
}

##########################################################
DomainPoints <- function(geo,...){
  UseMethod("DomainPoints")
}

DomainPoints.geofield <- function(geo,...){
  DomainPoints(attr(geo,"domain"),...)
}

DomainPoints.geodomain <- function (geo,type="lalo",...){
### return lat's and lon's of all domain points (or leave in projection if type = "xy")

  lalo <- list(x=c(geo$SW[1],geo$NE[1]),y=c(geo$SW[2],geo$NE[2]))
  xy <- project(lalo, proj = geo$projection)
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

gridpoints <- function(x,y=NULL,...){
### indicate particular grid points on a map (like points, but with grid indices)
  AllCoords <- DomainPoints(.Last.domain(),type="xy")
  xy <- xy.coords(x,y)

  xlist <- AllCoords$x[cbind(xy$x,xy$y)]
  ylist <- AllCoords$y[cbind(xy$x,xy$y)]

  points(xlist,ylist,...)
}

###########################################################

domain2lalo <- function(infield){
### return the LatLon co-ordinates of all grid points (X,Y), and the field values as Z
### all in 1 data frame (I think that's better than a list)
### You could use this for spatial interpolations...
  lalodomain <- DomainPoints(infield,"lalo")
  data.frame(x=as.vector(lalodomain$lon),y=as.vector(lalodomain$lat),
             z = as.vector(infield[1:attr(infield,'domain')$nx, 1:attr(infield,'domain')$ny]))
}


#########################
### SUBGRID           ###
#########################

subgrid <- function(geo,...){
  UseMethod("subgrid")
}

subgrid.geodomain <- function(geo,x1=1,x2=geo$nx,
                             y1=1,y2=geo$ny,reso=1,...) {
  xsub <- seq(x1,x2,by=reso)
  ysub <- seq(y1,y2,by=reso)
  subnx <- length(xsub)
  subny <- length(ysub)

  newlalo <- DomainPoints(geo,"lalo")
  newlalo$lon <- newlalo$lon[xsub,ysub]
  newlalo$lat <- newlalo$lat[xsub,ysub]
  geo$SW <- c(newlalo$lon[1,1],newlalo$lat[1,1])
  geo$NE <- c(newlalo$lon[subnx,subny],newlalo$lat[subnx,subny])

  geo$nx <- length(xsub)
  geo$ny <- length(ysub)
  geo
}

subgrid.geofield <- function(geo,x1=1,x2=attr(geo,"domain")$nx,
                             y1=1,y2=attr(geo,"domain")$ny,reso=1,...){
  xsub <- seq(x1,x2,by=reso)
  ysub <- seq(y1,y2,by=reso)

  subfield  <- geo[xsub,ysub]
  subdomain <- subgrid(attr(geo,"domain"),x1=x1,x2=x2,y1=y1,y2=y2,reso=reso)

  as.geofield(subfield,domain=subdomain,time=attr(geo,"time"),
             info=c(attr(geo,"info"),extra="SUBFIELD"))
}

################

zoomgrid <- function(infield,x,y,zoom=50){
  infield <- as.geofield(infield)
  if (x-zoom < 1) xmin <- 1
  else if (x+zoom > attr(infield,"domain")$nx) xmin <- attr(infield,"domain")$nx-2*zoom
  else xmin <- x-zoom

  if (y-zoom < 1) ymin <- 1
  else if (y+zoom > attr(infield,"domain")$ny) ymin <- attr(infield,"domain")$ny-2*zoom
  else ymin <- y-zoom
  xmax <- xmin+2*zoom
  ymax <- ymin+2*zoom

  subgrid(infield,xmin,xmax,ymin,ymax,reso=1)
}


#############################################################
### Some utilities for working with maps and projections  ###
#############################################################

DrawLatLon <- function(nx=9,ny=9,labels=TRUE,pretty=TRUE,
                     cex=1,col="blue",lty=2,font=2,lab.x=2, lab.y=2, ...) {
### add lattitude and longitude lines to  a map

### adapted from "map.grid" in mapproj library
### uses map.wrap from maps library
  pretty.range <- function(lim,...) {
    # like pretty but ensures that the range is identical:
    # range(pretty.range(x)) == range(x)
    x <- pretty(lim,...)
    if(abs(x[1]-lim[1]) > abs(x[2]-lim[1])) x = x[-1]
    n <- length(x)
    if(abs(x[n]-lim[2]) > abs(x[n-1]-lim[2])) x = x[-n]
#    x[1] <- lim[1]; x[length(x)] = lim[2]
    x
  }
  auto.format <- function(x) {
    # from minka library
    # use the minimal number of digits to make x's unique
    # similar to abbrev
    for(digits in 0:6) {
      s <- formatC(x,digits=digits,format="f")
      if(all(duplicated(s) == duplicated(x))) break
    }
    s
  }
  # by default, use limits of last map
  if(is.null(.Last.domain())) stop("Sorry, no projection has been defined.")
  glimits <- DomainExtent(.Last.domain())
  xlim <- glimits$lonlim
  ylim <- glimits$latlim

  if(pretty) {
    x <- pretty.range(xlim,n=nx)
    y <- pretty.range(ylim,n=ny)
  } else {
    x <- seq(xlim[1],xlim[2],len=nx)
    y <- seq(ylim[1],ylim[2],len=ny)
  }
  p <- project(expand.grid(x=c(seq(xlim[1],xlim[2],len=100),NA),y=y),
               proj=.Last.domain()$projection)
  p <- map.wrap(p)

  p <- map.restrict(p , c(glimits$x0,glimits$x1) , c(glimits$y0,glimits$y1),
                    xperiod=periodicity()$xper,yperiod=periodicity()$yper)
  lines(p,col=col,lty=lty,...)

  p2 <- expand.grid(y = c(seq(ylim[1], ylim[2], len = 100), NA),x=x)
### we must put the x and y components in the right order, because "project"
### just takes "p[[1]]", not p$x.
  p2 <- list(x=p2$x,y=p2$y)
  p2 <- project(p2, proj = .Last.domain()$projection)
  p2 <- map.restrict(p2, c(glimits$x0, glimits$x1), c(glimits$y0,glimits$y1),
                    xperiod=periodicity(.Last.domain()),yperiod=NA)
  lines(p2, col=col,lty=lty,...)

  if(labels) {
    if(lab.x<1 | lab.x>length(x)) lab.x=min(2,length(x))
    if(lab.y<1 | lab.y>=length(y)) lab.y=min(2,length(y)-1)
    tx <- x[lab.x]
    xinc <- median(diff(x))
    ty <- y[length(y)-lab.y]
    yinc <- median(diff(y))
    text.coord1 <- project(expand.grid(x = x + xinc * 0.05, y = ty + yinc * 0.5),
                        proj = .Last.domain()$projection)
### don't write labels outside the domain
### we should include a little buffer zone

    xbuf <- .02*(glimits$x1-glimits$x0)
    ybuf <- .02*(glimits$y1-glimits$y0)

    text.coord1$x[(text.coord1$x > (glimits$x1-xbuf))]=NA
    text.coord1$x[(text.coord1$y > (glimits$y1-ybuf))]=NA

    text(text.coord1, labels = auto.format(x),
         cex = cex, adj = c(0, 0), col = col, font = font,
         ...)

    text.coord2 <- project(expand.grid(x = tx + xinc * 0.5, y = y +
                             yinc * 0.05), proj = .Last.domain()$projection)
    text.coord2$x[(text.coord2$x > (glimits$x1-xbuf))]=NA
    text.coord2$x[(text.coord2$y > (glimits$y1-ybuf))]=NA

    text(text.coord2, labels = auto.format(y),
         cex = cex, adj = c(0, 0), col = col, font = font,
         ...)
  }
}

###########################
### CREATE NEW DOMAINS  ###
###########################

Make.domain <- function(projtype="lambert",clonlat,nxny,dxdy,reflat=clonlat[2],reflon=clonlat[1],tilt){
### Lambert (as used in ALADIN: only 1 reference latitude)
  if(length(dxdy)==1) dxdy <- rep(dxdy,2)
  if(projtype=="lambert")
    projection <- list(proj="lcc",lon_0=reflon,lat_1=reflat,lat_2=reflat,a=6371229.0,es=0.0)

### Rotated & tilted Mercator
### as used in ALADIN: reflon=clon, reflat=clat !!!
  else if(projtype=="mercator"){
    if(reflat!=clonlat[2] | reflon!=clonlat[1]) warning('This domain is not ALADIN-compatible!')
    if(abs(reflat)<.01)
      projection <- list(proj="merc",lon_0=reflon,a=6371229.0,es=0.0)
    else{
      if(abs(tilt)<1.0E-7) projection <- list(proj = "somerc", lonc = reflon,
                            lat_0 = reflat, a = 6371229, es = 0)
      else if (abs(abs(tilt)-90) < 1.0E-7)  projection <- list(proj = "tmerc", lonc = reflon,
                            lat_0 = reflat, a = 6371229, es = 0)
      else {
        if(tilt>0) projection <- list(proj = "omerc", lonc = reflon,
                            lat_0 = reflat, alpha = -90 + tilt, a = 6371229,
                            es = 0,no_rot=NA)
        else projection <- list(proj = "omerc", lonc = reflon,
                            lat_0 = reflat, alpha = 90 + tilt, a = 6371229,
                            es = 0,no_rot=NA)
      }
    }
  }
  else if(projtype=="latlong")
    projection=list(proj="latlong")
  else if(projtype=="RotLatLon"){
      projection <- list(proj="ob_tran","o_proj"="latlong",
                       "o_lat_p"=-reflat,"o_lon_p"=0,"lon_0"=reflon)
  }
  else stop("Unknown projection.")

### project the center point

  cxy <- project(list(x=clonlat[1],y=clonlat[2]),proj=projection)

  SW0 <- c(cxy$x,cxy$y) - dxdy*(nxny - 1)/2
  NE0 <- SW0 + dxdy*(nxny - 1)

  lims <- project(list(x=c(SW0[1],NE0[1]),y=c(SW0[2],NE0[2])),proj=projection,inv=TRUE)
  SW <- c(lims$x[1],lims$y[1])
  NE <- c(lims$x[2],lims$y[2])
### and the output is...
  result <- list(projection=projection,nx=nxny[1],ny=nxny[2],SW=SW,NE=NE,center=clonlat)
  class(result) <- "geodomain"
  result
}

MakeRLL <- function(Lon1,Lat1,SPlon,SPlat,SPangle=0,nxny,dxdy){
### This is for Rotated LatLon as used by Hirlam: central meridian is vertical.
### In the future, this should be merged with Make.domain
### but that is not trivial: here, the center point is of no real consequence
### and you have to define the rotated South Pole
  if(length(dxdy)==1) dxdy <- rep(dxdy,2)
  Lon2 <- Lon1+nxny[1]*dxdy[1]
  Lat2 <- Lat1+nxny[2]*dxdy[2]
  projection <- list(proj="ob_tran","o_proj"="latlong",
                       "o_lat_p"=-SPlat,"o_lon_p"=0,"lon_0"=SPlon)

  RR <- project(list(x=c(Lon1,Lon2)*pi/180,y=c(Lat1,Lat2)*pi/180),proj=projection,inv=TRUE)
  SW <- c(RR$x[1],RR$y[1])
  NE <- c(RR$x[2],RR$y[2])

  result <- list(projection=projection,nx=nxny[1],ny=nxny[2],SW=SW,NE=NE)
  class(result) <- "geodomain"
  result
}

##################################
### Interface to PROJ4 library ###
##################################

project <- function(x,y,proj=.Last.domain()$projection,inv=FALSE) {

  if (missing(y)) {
# apparantly, is.list gives TRUE for data.frames, but lets be careful:
    if(is.list(x) | is.data.frame(x)) {y <- x$y;x <- x$x}
    else if(is.vector(x)) {
      y <- x[2]
      x <- x[1]
    } else {
      y <- x[,2]
      x <- x[,1]
    }
  }

  if (missing(proj)) {
    if(!is.null(.Last.domain())) proj <- .Last.domain()$projection
    else return("No projection.")
  }

  if (proj$proj == "latlong") {
### longitude should be in the interval [-180,180[
### unless e.g. if my global data is on a globe [0,360[
### we assume that the meridian = MinLon + 180
### so meridian-180 must not be transported.
    meridian <- if(is.null(proj$lon0)) 0 else proj$lon0
#    x <- ifelse(x < meridian-180,x+360,x)
#    x <- ifelse(x >= meridian+180,x-360,x)
## much faster (!):
    x[which(x <  (meridian-180))] <- x[which(x <  (meridian-180))] + 360
    x[which(x >= (meridian+180))] <- x[which(x >= (meridian+180))] - 360
    
    data.frame(x=x,y=y)
  } else  {
    npoints <- as.integer(length(x))
    npar <- as.integer(length(proj))
    par <- paste(names(proj),lapply(proj,function(x) if(is.na(x)) "" else paste("=",x,sep="")),sep="")
### SIMPLER: par=paste(names(proj),'=',proj,sep='')
### but this doesn't allow options without =x value

    if (!inv) {
      x <- x/180*pi
      y <- y/180*pi
    }
### to fix what *I think* is a bug in PROJ4 (never had a reply)
### If they ever solve this bug, I'll have to change this!
    if (proj$proj=='omerc' & inv ){
      if (proj$alpha<0 ) x <- -x
      else y <- -y
    }
    result <- .C("Rproj4",x=x,y=y,npoints=npoints,par=par,
                 npar=npar,inv=as.integer(inv),NAOK=TRUE,PACKAGE="geogrid")
### again the same proj.4 bug:
    if (proj$proj=='omerc' & !inv) {
      if (proj$alpha<0 ) result$x <- -result$x
      else result$y <- -result$y
    }
    if (!inv) data.frame(x=result$x,y=result$y)
    else data.frame(x=result$x*180/pi,y=result$y*180/pi)
  }

}

periodicity <- function(domain=.Last.domain()){
### function that decides whether a domain/projection is periodic.
### Note that e.g. in any LatLon domain, the X co-ordinate is periodic
### even if the domain does cover the globe.
### Only cylindrical projections are periodic, I suppose.
### TO DO: there are more, but I am not interested in them.
### even
  xper <- switch(domain$projection$proj,
            "latlong"   = 360,
            "ob_tran"   = 2*pi, ### BUG: this could be any oblique projection
            "merc"     = 2 * pi * domain$projection$a,
            "omerc"    = 2 * pi * domain$projection$a,
            "tmerc"    = 2 * pi * domain$projection$a,
            "somerc"   = 2 * pi * domain$projection$a,
            NA_real_)
  yper <- NA_real_
  list(xper=xper,yper=yper)
}

