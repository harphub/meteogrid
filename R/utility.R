#-------------------------------------------#
# Part of R-package geogrid                 #
# Copyright (c) 2003-2016 Alex Deckmyn      #
#   Royal Meteorological Institute, Belgium #
# Released under GPL-3 license              #
#-------------------------------------------#

### format an integer as a string (simplified)
i2a <- function(i,n=ceiling(log(i)/log(10)) )
             formatC(i,width=n,flag=0,format="d")

### trim whitespace from character strings
trim <- function(x) sub(pattern=" +$",replacement="",
                        x=sub(pattern="| +",replacement="",x))

### To have the summary of a matrix or array
matsum <- function(x,...) summary(as.vector(x), ...)

### give a default value for NULL
checknull <- function(x, default=0) if (is.null(x)) default else x

### click on the map and get the Lat-Lon coordinates
laloclick <- function(n=1, ...){
  xy <- locator(n, ...)
  project(xy,proj=.Last.domain()$projection,inv=TRUE)
}

### indicate particular grid points on a map (like points, but with grid indices)
gridpoints <- function(x ,y=NULL, ...){
  AllCoords <- DomainPoints(.Last.domain(),type="xy")
  xy <- xy.coords(x,y)

  xlist <- AllCoords$x[cbind(xy$x,xy$y)]
  ylist <- AllCoords$y[cbind(xy$x,xy$y)]

  points(xlist,ylist,...)
}


### For fun: plot an orthographic projection of the world.
orthoglobe <- function(reflon=0,reflat=90,map.database="world",...){
  projection <- list(proj="ortho",lon_0=reflon,lat_0=reflat,a=6371229.0,es=0.0)
  boundaries <- maps::map(database=map.database, plot=FALSE)
  geo <- project(boundaries, proj =projection,inv = FALSE)

  plot(geo,type="l",xlab="",ylab="",axes=FALSE,...)
  box()
  domain <- list(projection=projection)
  .Last.domain(domain)
}

antipolygon <- function(data, xylim=DomainExtent(.Last.domain()), bg="white"){
### inspired by the version in the R-package GEOmap
### this colours everything outside the polygon given by "data"
### up to the limits in xylim.
### typical usage:
### mask=project(map("world","Belgium",fill=0,plot=FALSE))
### xylim=DomainExtent(datafield)
### antipolygon(mask,xylim=xylim)
    x <- c(data$x, data$x[1],xylim$x0,xylim$x0,xylim$x1,xylim$x1,xylim$x0)
    y <- c(data$y, data$y[1],xylim$y0,xylim$y1,xylim$y1,xylim$y0,xylim$y0)
    polygon(x, y, border = bg, col = bg, xpd = TRUE)
}


########################
### COLOUR PALETTES  ###
########################
### for nice colour palettes, install RColorBrewer
col.oro <- function(n){
  GYR <- colors()[c(615,81,258,384,38,504,634)]
  colorRampPalette(GYR)(n)
}

col.precip <- function(n){
  WYB <- colors()[c(0,513,105,125,29,477)]
  colorRampPalette(WYB)(n)
}

col.temp <- function(n){
  BYR <- colors()[c(477,619,615,383,38,573,504,35,645)]
  colorRampPalette(BYR)(n)
}

irainbow=function(n) {
  if(requireNamespace("RColorBrewer",quietly=TRUE)){
    Spectral <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
  } else Spectral <- colors()[c(107,619,12,105,384,383,38,573,504,35,645)]

  if (n>0) colorRampPalette(Spectral)(n)
  else rev(colorRampPalette(Spectral)(-n))
}

### transparent white (probably requires Cairo)
col.cloud=function(n) rgb(1,1,1,(0:(n-1))/(n-1))

#"irainbow" <- function (n) {
#  if(n>0)  rev(rainbow(n, end = 0.65))
#  else rainbow(-n, end = 0.65)
#}
### grey scale rainbow (black=low,white=high)
### useful for B&W graphs in journals that charge extra for colour print :-)
"grainbow" <- function (n) {
  if (n<0) grey(seq(1,0,length=-n))
  else grey(seq(0,1,length=n))
}
