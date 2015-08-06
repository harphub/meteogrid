#--------------------------------------#
# Part of R-package geogrid            #
# Copyright Alex Deckmyn               #
# RMI Belgium, alex.deckmyn@meteo.be   #
# Released under GPL-3 license         #
#--------------------------------------#

# geowind: rotation of wind fields
# not easy to do in a generic way, so we only implement for a few projections:
# rotated LatLon and Lambert.

# TO DO: add more projections? (rotated) mercator, polar-sterographic
#        check mapfactor
#        output as speed/direction?

### rotate model wind to geographic wind
### the conversion is done for the whole domain

###########################################################
### 1. very basic functions for (u,v) <-> (wdir,wspeed) ###
###########################################################

wind.dirspeed <- function(u,v,fieldname=c("Wind direction","Wind speed"),rad=FALSE){
  if(missing(v) & is.list(u)) {
    v <- u[[2]]
    u <- u[[1]]
  }
  MINSPEED <- 10E-6

  wspeed <- sqrt(u^2 + v^2)
  wdir <- ifelse(abs(u)>MINSPEED,
                      ( -180 - atan(v/u) * 180/pi + sign(u)*90 ) %% 360,
                      ifelse(abs(v)<MINSPEED,NA, ifelse(v>0,180,0) ) )
  if(rad) wdir <- wdir*pi/180.
  if(is.geofield(u)) {
    attributes(wdir) <- attributes(u)
    attributes(wdir)$info$name <- fieldname[1]
    attributes(wspeed) <- attributes(u)
    attributes(wspeed)$info$name <- fieldname[2]
  }
  if(is.vector(u)) return(data.frame(wdir=wdir,wspeed=wspeed))
  else return(list(wdir=wdir,wspeed=wspeed))
}

wind.uv <- function(wspeed,wdir,fieldname=c("U","V"),rad=FALSE){
  if(missing(wdir) & is.list(wspeed)){
    wdir <- wspeed$wdir
    wspeed <- wspeed$wspeed
  }
  if(!rad) wdir <- wdir*pi/180.
  u <- wspeed*cos(wdir)
  v <- wspeed*sin(wdir)
  if(is.geofield(wdir)) {
    attributes(u)$info$name <- fieldname[1]
    attributes(v)$info$name <- fieldname[2]
  }
  if(is.vector(wspeed)) return(data.frame(U=u,V=v))
  else return(list(U=u,V=v))
}


###########################################################
### 2. main routine for rotation grid axes <-> N/E axes ###
###########################################################

geowind <- function(u,v,inv=FALSE,init=NULL){
  if( is.null(init) ){
    domain <- attributes(u)$domain
    ww <- geowind.init(domain)
  }
  else ww <- init

  if(inv) {
    ww$angle <- -ww$angle
    ww$mapfactor <- 1/ww$mapfactor
  }
  U <- ( cos(ww$angle) * u + sin(ww$angle) * v) * ww$mapfactor
  V <- (-sin(ww$angle) * u + cos(ww$angle) * v) * ww$mapfactor
  if(!inv){
    attributes(U)$info$name <- paste(attributes(u)$info$name,
            "Rotated to N/E axes.")
    attributes(V)$info$name <- paste(attributes(v)$info$name,
            "Rotated to N/E axes.")
  }
  else {
    attributes(U)$info$name <- paste(attributes(u)$info$name,
            "grid axes.")
    attributes(V)$info$name <- paste(attributes(v)$info$name,
            "grid axes.")
  }
  if(is.vector(u)) data.frame(U=U,V=V)
  else list(U=U,V=V)
}

geowind.init <- function(domain){
### ATTENTION: in the output
### $angle is the angle that must be added to the wind direction (=wind origin!) from grid
### to get the actual geographic value
### checked for LCC,RLL
  if (is.geofield(domain)) domain <- attributes(domain)$domain
  ww <- switch(domain$projection$proj,
          "ob_tran" = geowind.RLL(domain),
          "lcc" = geowind.LCC(domain),
          "stere" = geowind.PS(domain),
          "omerc" = geowind.RM(domain),
          stop(paste("unimplemented projection: ",domain$projection$proj))
        )
  ww
}

######################################################################
### 3. functions that calculate local rotation angle and mapfactor ###
######################################################################

## rotated Lat/Lon
geowind.RLL <- function(domain){
  rad <- pi/180.
  SPlat <- -domain$projection$o_lat_p
  SPlon <- domain$projection$lon_0
  sinP <- sin( (SPlat + 90) * rad )
  cosP <- cos( (SPlat + 90) * rad )

  lalo <- DomainPoints(domain,"lalo")
  Rlalo <- DomainPoints(domain,"xy")

  sinRX <- sin(Rlalo$x)
  cosRX <- cos(Rlalo$x)
  cosRY <- cos(Rlalo$y)

  lon0 <- rad*(lalo$lon - SPlon)
  sinX0 <- sin(lon0)
  cosX0 <- cos(lon0)

  sinA <- sinP*sinX0/cosRY
  cosA <- cosP*sinX0*sinRX + cosX0*cosRX
# fix numerical issues (acos will crash on 1.+eps):
  cosA[which(cosA > 1)] <- 1
  cosA[which(cosA < -1)] <- -1
  angle <- acos(cosA) * (2*(sinA>=0)-1)
# the result has domain ]-PI,+PI[
## ATTENTION: angle calculated here is for the wind origin
  list(angle = angle,mapfactor = 1)
}

## Lambert conformal conical
geowind.LCC <- function(domain){
### version in Rfa (based on Luc Gerard?)
### 12/2014: change sign of angle for consistency
  rad <- pi/180.

  reflat <- domain$projection$lat_1
  reflon <- domain$projection$lon_0
  refsin <- sin(reflat * rad)
  refcos <- cos(reflat * rad)

  lalo <- DomainPoints(domain,"lalo")

  mapfactor <- (refcos/cos(lalo$lat * rad))^(1 - refcos) * ((1 + refsin)/(1 + sin(lalo$lat * rad)))^refsin
  angle <- refsin * (lalo$lon - reflon) * rad
  list(angle=angle,mapfactor=mapfactor)
}

geowind.LCC2 <- function(domain){
### as coded in GL
### no mapfactor??? for the rest, identical result if lat1=lat2
  rad <- pi/180.
  lalo <- DomainPoints(domain,"lalo")

  lat1 <- domain$projection$lat_1
  lat2 <- domain$projection$lat_2
  lon0 <- domain$projection$lon_0

  if (abs(lat1-lat2)<1.E-6) {n <- sin(lat1 * rad)
  }else n <- log(cos(lat1*rad)/cos(lat2*rad)) / log(tan((lat2/2+45)*rad)/tan((lat1/2+45)*rad))

  diff <- lalo$lon - lon0
  diff[diff < -180] <- diff[diff < -180]+360
  diff[diff >  180] <- diff[diff >  180]-360
  alpha <- diff * n * rad * sign(lat1)
  list(angle= -alpha,mapfactor=1)
}

## Polar Stereographic
geowind.PS <- function(domain){
  warning("Polar Stereographic wind rotation: unvalidated!!!")
  rad <- pi/180.
  lalo <- DomainPoints(domain,"lalo")

  reflon <- domain$projection$lon_0
  reflat <- domain$projection$lat_0

  diff <- (lalo$lon - reflon)%%360
  diff[diff > 180] <- diff[diff > 180]-360

  angle <- diff * rad * (2*(reflat>=0)-1)
  list(angle=angle,mapfactor=1)
}

## Rotated Mercator
geowind.RM <- function(domain){
  stop("Rotated Mercator wind rotation: unfinished!!!")
  rad <- pi/180.
  lalo <- DomainPoints(domain,"lalo")
# TO DO
  cosA[which(cosA > 1)] <- 1
  cosA[which(cosA < -1)] <- -1
  angle <- acos(cosA) * (2*(sinA>=0)-1)
  list(angle=angle,mapfactor=1)
}


