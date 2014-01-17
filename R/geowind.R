# geowind: rotation of wind fields
# not easy to do in a generic way, so we only implement for a few projections:
# rotated LatLon and Lambert.

# TO DO: add more projections? (rotated) mercator, polar-sterographic
#        check mapfactor
#        output as speed/direction?
#        RLL can be simplified A LOT: the GL code is insanely complex

### write functions that calculate anngle and mapfactor
geowind <- function(U,V,inv=FALSE){
  domain <- attributes(U)$domain
  ppp <- domain$projection
  if (ppp$proj=="ob_tran") ww <- geowind.RLL(domain)
  else if(ppp$proj=="lcc") ww <- geowind.LCC(domain)
  else if(ppp$proj=="stere") ww <- geowind.PS(domain)
  else if(ppp$proj=="omerc") ww <- geowind.RM(domain) # TODO: may also be tmerc,somerc,merc
  else error("unimplemented projection")

  if(inv) {
    ww$angle <- -ww$angle
    ww$mapfactor <- 1/ww$mapfactor
  }
  U <- (cos(ww$angle) * U - sin(ww$angle) * V) * ww$mapfactor
  V <- (sin(ww$angle) * U + cos(ww$angle) * V) * ww$mapfactor
  if(!inv){
    attributes(U)$info$name <- paste(attributes(U)$info$name, 
            "Rotated to N/E axes.")
    attributes(V)$info$name <- paste(attributes(V)$info$name, 
            "Rotated to N/E axes.")
  }
  else {
    attributes(U)$info$name <- paste(attributes(U)$info$name, 
            "grid axes.")
    attributes(V)$info$name <- paste(attributes(V)$info$name, 
            "grid axes.")
  }
  list(U=U,V=V)
}

geowind.RLL <- function(domain){
  rad <- pi/180.
  SPlat <- -domain$projection$o_lat_p
  SPlon <- domain$projection$lon_0
  sinP <- sin( (SPlat + 90) * rad )
  cosP <- cos( (SPlat + 90) * rad )

  lalo <- DomainPoints(domain,"lalo")

# sinRY <- cosP*sinY + sinP*cosY*cosX

  sinA <- sinP*sinX0/cosRY
  cosA <- cosP*sinX0*sinRX + cosX0*cosRX
# fix numerical issues (acos will crash on 1.+eps):
  cosA[which(cosA > 1)] <- 1
  cosA[which(cosA < -1)] <- -1
  angle <- acos(cosA) * (2*(sinA>=0)-1)
  list(angle=angle,mapfactor=1)
}

geowind.LCC <- function(domain){
### version in Rfa (based on Luc Gerard?)
  rad <- pi/180.

  reflat <- domain$projection$lat_1
  reflon <- domain$projection$lon_0
  refsin <- sin(reflat * rad)
  refcos <- cos(reflat * rad)

  lalo <- DomainPoints(domain,"lalo")

  mapfactor <- (refcos/cos(lalo$lat * rad))^(1 - refcos) * ((1 + refsin)/(1 + sin(lalo$lat * rad)))^refsin
  angle <- -refsin * (lalo$lon - reflon) * rad
  list(angle=angle,mapfactor=mapfactor)
  list(U=U,V=V)
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

# Polar Stereographic
geowind.PS <- function(domain){
  rad <- pi/180.
  lalo <- DomainPoints(domain,"lalo")

  reflon <- domain$projection$lon_0
  reflat <- domain$projection$lat_0

  diff <- (lalo$lon - reflon)%%360
  diff[diff > 180] <- diff[diff > 180]-360

  angle <- diff * rad * (2*(reflat>=0)-1) 
  list(angle=angle,mapfactor=1)
}

# Rotated Mercator
geowind.RM <- function(domain){
  rad <- pi/180.
  lalo <- DomainPoints(domain,"lalo")
# TO DO
  cosA[which(cosA > 1)] <- 1
  cosA[which(cosA < -1)] <- -1
  angle <- acos(cosA) * (2*(sinA>=0)-1)
  list(angle=angle,mapfactor=1)
}


geowind.RLL0 <- function(U,V,inv=FALSE) {
### transform wind field from Rotated LatLon to normal (or back)
### R version of the code in GL: a bit complex...
   rad <- pi/180.
  domain <- attributes(U)$domain
  SPlat <- -domain$projection$o_lat_p
  SPlon <- domain$projection$lon_0
  sinP <- sin( (SPlat + 90) * rad )
  cosP <- cos( (SPlat + 90) * rad )

  lalo <- DomainPoints(U,"lalo")
  Rlalo <- DomainPoints(U,"xy")
  Rlalo$lat <- Rlalo$y * 180/pi
  Rlalo$lon <- Rlalo$x * 180/pi

#  sinX <- sin(lalo$lon * rad)
#  cosX <- cos(lalo$lon * rad)
  sinY <- sin(lalo$lat * rad)
  cosY <- cos(lalo$lat * rad)

  sinRX <- sin(Rlalo$lon * rad)
  cosRX <- cos(Rlalo$lon * rad)
  sinRY <- sin(Rlalo$lat * rad)
  cosRY <- cos(Rlalo$lat * rad)

  lon0 <- rad*(lalo$lon - SPlon)
  sinX0 <- sin(lon0)
  cosX0 <- cos(lon0)

### these should be competely inverse. Are they? RX <-> X etc?
### This doesn't look like a rotation at all...
### IN FACT: A==D,C==-B , and for (inv) it's just a matter of C -> -C
### so we can simplify a lot!!!
### just find rotation angle: angle=acos(A) * sign2(B) |sign2(0)=1 sign2(x)= 2*(x>=0)-1
  if(inv){ # from latlon to RotLatLon
    A <- cosP*sinX0*sinRX + cosX0*cosRX
    B <- cosP*cosX0*sinY*sinRX - sinP*cosY*sinRX - sinX0*sinY*cosRX
    C <- sinP*sinX0/cosRY
    D <- (sinP*cosX0*sinY + cosP*cosY)/cosRY
  }
  else { # from RotLatLon to LatLon
    A <- cosP*sinX0*sinRX + cosX0*cosRX
    B <- cosP*sinX0*cosRX*sinRY + sinP*sinX0*cosRY - cosX0*sinRX*sinRY
    C <- -sinP*sinRX/cosY
    D <- (cosP*cosRY - sinP*cosRX*sinRY)/cosY
  }
  result <- list(U= A*U + B*V, V=C*U+D*V)
  result
}
