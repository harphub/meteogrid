######################
### geogrid
######################
### Alex Deckmyn
### Royal Meteorological Institute of Belgium
### alex.deckmyn@oma.be
######################
### routines for viewing gridded geographical/meteorological  data
### uses Proj.4 library . It requires the proj4 package.
###
### 03/03/2004 : adapted image.geofield to be closer to image.default
###              changed proj4 format to resemble mapproject
###              dropped Rutils dependence
### 23/06/2004 : adapted DomainLatLon, subdomain, lalopoint
### 01/07/2004 : changed class name to geofield
### 22/02/2005 : turn geofield class into a matrix with attributes
### 05/01/2006 : small changes, DomainPoints function replaces DomainLatLon
### 08/08/2006 : added easy creation and plotting of domains
###              (Make.lccdomain, plot.geodomain)
### 16/08/2006 : bugfix in DomainPoints (Luc Gerard)
### 14/09/2006 : little fix for titles in cview and fcview
### 15/09/2006 : assume aspect=1 in vecplot (i.e. the axes have the same scale)
###              this should correct the direction of the arrows
###              For Lat-Lon: local correction for aspect ratio!
### 26/09/2006 : Bugfix in domain2lalo (Luc Gerard)
### 04/10/2006 : Bugfix in plot.geodomain (Luc Gerard)
###              + a little bit of cleaning
###              + renamed DrawLatLon
### 17/01/2007 : New function Make.domain : includes different (ALADIN) projections
### 18/01/2007 : Add "orthoglobe" for a nice plot of the globe.
### 19/01/2007 : Port map (-> "map4")
###              add little function "gridpoint"
### 23/02/2007 : Added "regrid": interpolation from one domain to another.
### 08/05/2007 : Added interp.point (interpolation to 1 lat-lon point)
### 23/01/2008 : Added Rotated LatLon projection (for Hirlam output)
###              replace akima by bilinear interpolation (using "fields" library)-> FASTER
### 24/01/2008 : adapt to CRAN-package "proj4" (so my routine drops out!)
###              projections are now defined as a list (but backward compatible)
### xx/xx/2008 : Little improvements to domain plotting: by using "image(col=0)"
###              you get the exact borders (no more boundary lines beyond the frame)
### 14/10/2008 : New function "limage" for legend support in image and filled.contour.
###              Rewrite of iview and fcview for decent legend support.
###              Code cleaning.
### 17/11/2008 : Some bug fixes (Jure Cedilnik)
###              new function domainbox
### 02/02/2009 : bugfix domainbox
### ../06/2009 : several fixes for Rotated Mercator.
### 21/10/2009 : bugfix in lalopoint (Tomasz Kułakowski)
### 08/03/2010 : bugfix for legend.title in limage.default
### 09/06/2010 : bugfix in project - latlong was not processed correctly
### 13/12/2010 : bug fix for contour.geofield(drawmap=FALSE) (Luc Gerard)
### 21/01/2011 : bugfixes for global maps (reported by Piet Termonia)
### 03/05/2011 : bugfix as.geofield
### 16/09/2011 : add compare.geodomain
### 29/05/2012 : add useRaster=TRUE as an option to iview.
###              Check that version>=2.13 before calling image.
### 29/06/2012 : Add "mask" option to lalopoint (e.g. to find closest LAND point)
### 21/03/2013 : New regrid implementation:
###              Bilinear, Nearest neighbour and bicubic spline interpolations.
### 11/07/2013 : Replace "mapNew" by a simple lat/lon table, so yo udon't need "maps".
##########################################

is.geofield <- function(x){
  inherits(x,"geofield")
}

is.geodomain <- function(x){
  inherits(x,"geodomain")
}

"print.geofield" <- function(field){
  print.noquote(paste(attr(field,"info")$origin,":",attr(field,"info")$name))
  print.noquote("Time:")
  print.noquote(attr(field,"time"))
  print.noquote("Domain summary:")
  print(attr(field,"domain"))
  print.noquote("Data summary:")
  print.noquote(summary(as.vector(field)))
}

print.geodomain = function(domain){
  print.noquote(paste(domain$nx,"x",domain$ny,"domain"))
  print.noquote("Projection summary:")
  print.noquote(paste("proj=",domain$projection$proj))
  print.noquote(paste("NE = (",domain$NE[1],",",domain$NE[2],")"))
  print.noquote(paste("SW = (",domain$SW[1],",",domain$SW[2],")"))
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

DomainExtent <- function(x,...){
  UseMethod("DomainExtent")
}

DomainExtent.geofield <- function(field){
  DomainExtent(attr(field,"domain"))
}

DomainExtent.geodomain <- function(domain){
### We look for the extreme LatLon values of a domain domain
### these are useful to select the right part from the world map
### first find the SW and NE point, then draw a rectangular box
### and project it back to LatLon
### You could use DomainPoints in stead, but then you'd be projecting
### all domain points, which is a bit of an overkill...
###  if(is.geofield(domain)) domain <- attr(domain,"domain")
  xy <- project(list(x=c(domain$SW[1], domain$NE[1]),y=c(domain$SW[2], domain$NE[2])),
                proj=domain$projection)
  x0 <- xy$x[1]
  y0 <- xy$y[1]
  x1 <- xy$x[2]
  y1 <- xy$y[2]
  dx <- (x1-x0)/(domain$nx-1)
  dy <- (y1-y0)/(domain$ny-1)
  xc <- (x0+x1)/2
  yc <- (y0+y1)/2

  clonlat <- project(list(x=xc,y=yc),proj=domain$projection,inv=TRUE)
  borders <- project(list(x=c(seq(x0,x1,length=domain$nx),rep(x1,domain$ny),
                              seq(x0,x1,length=domain$nx),rep(x0,domain$ny)),
                          y=c(rep(y0,domain$nx),seq(y0,y1,length=domain$ny),
                              rep(y1,domain$nx),seq(y0,y1,length=domain$ny))),
                      proj=domain$projection,inv=TRUE)
### ATTENTION: if the map crosses or comes close to the date line (meridian 180/-180) this will not
### work correctly. In the map command we must then set lonlim=NULL !!!

  lonlim <- range(borders$x,na.rm=TRUE)
  llr <- lonlim[2]-lonlim[1]
  if ( llr <= 1 | llr > 300 ) lonlim <- NULL

  list(lonlim = lonlim , latlim = range(borders$y,na.rm=TRUE),
       clonlat=c(clonlat$x,clonlat$y),
       x0=x0,y0=y0,x1=x1,y1=y1,dx=dx,dy=dy,nx=domain$nx,ny=domain$ny)
}

##########################################################
DomainPoints <- function(x,...){
  UseMethod("DomainPoints")
}

DomainPoints.geofield <- function(field,...){
  DomainPoints(attr(field,"domain"),...)
}

DomainPoints.geodomain <- function (domain,type="lalo"){
### return lat's and lon's of all domain points (or leave in projection if type = "xy")

  lalo <- list(x=c(domain$SW[1],domain$NE[1]),y=c(domain$SW[2],domain$NE[2]))
  xy <- project(lalo, proj = domain$projection)
  xydomain <- expand.grid(x = seq(xy$x[1], xy$x[2], length = domain$nx),
                          y = seq(xy$y[1], xy$y[2], length = domain$ny))
  if(type=="lalo") {
    lalolist <- project(xydomain, proj = domain$projection, inv = TRUE)
    list(lon = matrix(lalolist$x, ncol = domain$ny, nrow = domain$nx),
         lat = matrix(lalolist$y, ncol = domain$ny, nrow = domain$nx))
  }
  else if (type=="xy") {
    list(x=matrix(xydomain$x, ncol = domain$ny, nrow = domain$nx),
         y=matrix(xydomain$y, ncol = domain$ny, nrow = domain$nx))
  }
  else print("Unknown type.")
}

gridpoints <- function(x,y=NULL,...){
### indicate particular grid points on a map (like points, but with grid indices)
  AllCoords <- DomainPoints(.Last.domain,type="xy")
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

subgrid <- function(x,...){
  UseMethod("subgrid")
}

subgrid.geodomain <- function(domain,x1=1,x2=domain$nx,
                             y1=1,y2=domain$ny,reso=1) {
  xsub <- seq(x1,x2,by=reso)
  ysub <- seq(y1,y2,by=reso)
  subnx <- length(xsub)
  subny <- length(ysub)

  newlalo <- DomainPoints(domain,"lalo")
  newlalo$lon <- newlalo$lon[xsub,ysub]
  newlalo$lat <- newlalo$lat[xsub,ysub]
  domain$SW <- c(newlalo$lon[1,1],newlalo$lat[1,1])
  domain$NE <- c(newlalo$lon[subnx,subny],newlalo$lat[subnx,subny])

  domain$nx <- length(xsub)
  domain$ny <- length(ysub)
  domain
}

subgrid.geofield <- function(infield,x1=1,x2=attr(infield,"domain")$nx,
                             y1=1,y2=attr(infield,"domain")$ny,reso=1){
  xsub <- seq(x1,x2,by=reso)
  ysub <- seq(y1,y2,by=reso)

  subfield  <- infield[xsub,ysub]
  subdomain <- subgrid(attr(infield,"domain"),x1=x1,x2=x2,y1=y1,y2=y2,reso=reso)

  as.geofield(subfield,domain=subdomain,time=attr(infield,"time"),
             info=c(attr(infield,"info"),extra="SUBFIELD"))
}

################

zoomgrid <- function(infield,x,y,zoom=50){
  infield <- as.geofield(infield)
  if(x-zoom < 1) xmin <- 1
  else if (x+zoom > attr(infield,"domain")$nx) xmin <- attr(infield,"domain")$nx-2*zoom
  else xmin <- x-zoom

  if(y-zoom < 1) ymin <- 1
  else if (y+zoom > attr(infield,"domain")$ny) ymin <- attr(infield,"domain")$ny-2*zoom
  else ymin <- y-zoom
  xmax <- xmin+2*zoom
  ymax <- ymin+2*zoom

  subgrid(infield,xmin,xmax,ymin,ymax,reso=1)
}



########## OLD INTERPOLATION ROUTINES :

regrid.old <- function (infield, newdomain=.Last.domain)
{
### Bilinear interpolation to different grid
  require(fields)
  if(is.geofield(newdomain)) newdomain=attr(newdomain,"domain")
  newpoints <- DomainPoints(newdomain)
  projpoints <- project(list(x=as.vector(newpoints$lon), y=as.vector(newpoints$lat)),
                      proj = attr(infield,"domain")$projection)
  glimits <- DomainExtent(attr(infield, "domain"))
  x <- seq(glimits$x0, glimits$x1, length = attr(infield, "domain")$nx)
  y <- seq(glimits$y0, glimits$y1, length = attr(infield, "domain")$ny)

  result <- interp.surface(list(x=x, y=y, z=infield[1:glimits$nx, 1:glimits$ny]),cbind(projpoints$x,projpoints$y))

  as.geofield(matrix(result, ncol = newdomain$ny, nrow = newdomain$nx),
              domain = newdomain,time=attr(infield,"time"),info=attr(infield,"info"))
}

point.interp.old <- function(lon,lat,infield){
### interpolate to a given (set of) point(s):
  require(fields)
  gdomain <- attr(infield, "domain")
  glimits <- DomainExtent(gdomain)
  x <- seq(glimits$x0, glimits$x1, length = gdomain$nx)
  y <- seq(glimits$y0, glimits$y1, length = gdomain$ny)

  projpoints <- project(list(x=lon,y=lat),proj=gdomain$projection)

  result <- interp.surface(list(x=x, y=y, z=infield[1:glimits$nx, 1:glimits$ny]),
               cbind(projpoints$x,projpoints$y))
  result
}

point.closest.old <- function(lon,lat,infield,...){
### simply return value of the closest grid points
### as calculated by lalopoint
  n <- length(lon)
  x <- rep(NA,n)
  for(i in 1:n) x[i] <- lalopoint(infield,lon[i],lat[i],...)$data
  x
}

lalopoint <- function(data,lon,lat,minimise='lalo',mask=NULL){
### find the closest domain point to the given co-ordinates
### by minimising distance in either LonLat or projected co-ordinates
### This is in fact not exactly the same as minimising geographical distance!
  ldata <- is.geofield(data)
  if (ldata) domain <- attr(data, "domain")
  else domain <- data

  lalodomain <- DomainPoints(domain,"lalo")
  nx <- domain$nx
  ny <- domain$ny

  if(minimise=='lalo') dist <- sqrt((lalodomain$lon-lon)^2+(lalodomain$lat-lat)^2)
  else if (minimise=='proj') {
    xydomain <- DomainPoints(domain, "xy")
    ppp <- project(list(x=lon,y=lat), proj = domain$projection)
    dist <- sqrt((xydomain$x - ppp$x)^2 + (xydomain$y - ppp$y)^2)
  }
  else stop("unknown minimisation.")

  if (!is.null(mask)) {
      if (is.character(mask))
          mask <- eval(parse(text = mask))
      else if (is.expression(mask))
          mask <- eval(mask)
      dist[eval(expression(mask))] <- NA
  }
# We wouldn't expect two points to be exactly as close (numerically)
# But you never know. So in this case we chose at random.
# Default method is to set the index of both to 1.5, which screws up everything.

#  dist.rank[] <- rank(dist,ties.method="random",na.last=TRUE)
#  ij <- which(dist.rank==1,arr.ind=TRUE)

  ij <- which(dist==min(dist,na.rm=TRUE),arr.ind=TRUE)
  i <- ij[1,1]
  j <- ij[1,2]

  il <- lalodomain$lon[i,j]
  jl <- lalodomain$lat[i,j]

  if(ldata)
    list(data = data[i, j], lonlat = c(il, jl), index = c(i,j))
  else
    list(data = NA, lonlat = c(il, jl), index = c(i,j) )
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
  if(!exists(".Last.domain")) stop("Sorry, no projection has been defined.")
  glimits <- DomainExtent(.Last.domain)
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
               proj=.Last.domain$projection)
  p <- map.wrap(p)

  p <- map.restrict(p , c(glimits$x0,glimits$x1) , c(glimits$y0,glimits$y1),
                    xper=periodicity()$xper,yper=periodicity()$yper)
  lines(p,col=col,lty=lty,...)

  p2 <- expand.grid(y = c(seq(ylim[1], ylim[2], len = 100), NA),x=x)
### we must put the x and y components in the right order, because "project"
### just takes "p[[1]]", not p$x.
  p2 <- list(x=p2$x,y=p2$y)
  p2 <- project(p2, proj = .Last.domain$projection)
  p2 <- map.restrict(p2, c(glimits$x0, glimits$x1), c(glimits$y0,glimits$y1),
                    xperiod=periodicity(.Last.domain),yperiod=NA)
  lines(p2, col=col,lty=lty,...)

  if(labels) {
    if(lab.x<1 | lab.x>length(x)) lab.x=min(2,length(x))
    if(lab.y<1 | lab.y>=length(y)) lab.y=min(2,length(y)-1)
    tx <- x[lab.x]
    xinc <- median(diff(x))
    ty <- y[length(y)-lab.y]
    yinc <- median(diff(y))
    text.coord1 <- project(expand.grid(x = x + xinc * 0.05, y = ty + yinc * 0.5),
                        proj = .Last.domain$projection)
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
                             yinc * 0.05), proj = .Last.domain$projection)
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

project <- function(x,y,proj=.Last.domain$projection,inv=FALSE)
{

  if(missing(y)){
# apparantly, is.list gives TRUE for data.frames, but lets be careful:
    if(is.list(x) | is.data.frame(x)) {y <- x$y;x <- x$x}
    else if(is.vector(x)) {y <- x[2];x <- x[1]}
    else {y <- x[,2];x <- x[,1]}
  }

  if(missing(proj)){
    if(exists(".Last.domain")) proj <- .Last.domain$projection
    else return("No projection.")
  }

  if(proj$proj == "latlong")  {
### longitude should be in the interval [-180,180[
### unless e.g. if my global data is on a globe [0,360[
### we assume that the meridian = MinLon + 180
### so meridian-180 must not be transported.
    meridian <- ifelse(is.null(proj$lon0),0,proj$lon0)
    x <- ifelse(x < meridian-180,x+360,x)
    x <- ifelse(x >= meridian+180,x-360,x)
    data.frame(x=x,y=y)
  }
  else  {
    npoints <- as.integer(length(x))
    npar <- as.integer(length(proj))
    par <- paste(names(proj),lapply(proj,function(x) if(is.na(x)) "" else paste("=",x,sep="")),sep="")
### SIMPLER: par=paste(names(proj),'=',proj,sep='')
### but this doesn't allow options without =x value

    if(!inv) {x <- x/180*pi ; y <- y/180*pi}
### to fix what *I think* is a bug in PROJ4 (never had a reply)
### If they ever solve this bug, I'll have to change this!
    if(proj$proj=='omerc' & inv ){
      if (proj$alpha<0 ) x <- -x
      else y <- -y
    }
    result <- .C("Rproj4",x=x,y=y,npoints=npoints,par=par,
                 npar=npar,inv=as.integer(inv),NAOK=TRUE)
### again the same proj.4 bug:
    if(proj$proj=='omerc' & !inv){
      if (proj$alpha<0 ) result$x <- -result$x
      else result$y <- -result$y
    }
    if (!inv) data.frame(x=result$x,y=result$y)
    else  data.frame(x=result$x*180/pi,y=result$y*180/pi)
  }

}

periodicity <- function(domain=.Last.domain){
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



