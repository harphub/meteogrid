#--------------------------------------#
# Part of R-package geogrid            #
# © Alex Deckmyn                       #
# Released under GPL-3 license         #
#--------------------------------------#

### format an integer as a string (simplified)
i2a=function(x,n=ceiling(log(x)/log(10)) )
             formatC(x,width=n,flag=0,format="d")

### To have the summary of a matrix or array
matsum <- function(x,...) summary(as.vector(x))

### click on the map and get the Lat-Lon coordinates
laloclick = function(n=1,type='n',...){
  xy=locator(n,type=type,...)
  project(xy,proj=.Last.domain$projection,inv=TRUE)
}

### For fun: plot an orthographic projection of the world.
orthoglobe = function(reflon=0,reflat=90,map.database="worldmap",...){
  projection = list(proj="ortho",lon_0=reflon,lat_0=reflat,a=6371229.0,es=0.0)
  if(!exists(map.database)) data(list=map.database)
  boundaries <- worldmap
  geo = project(boundaries, proj =projection,inv = FALSE)

  plot(geo,type="l",xlab="",ylab="",axes=FALSE,...)
  box()
  domain = list(projection=projection)
  assign(".Last.domain",domain,envir=globalenv())
}

as.ocean=function(data, XYLIM=par('usr')){
### convert land boundaries into an 'ocean boundary'
###  data is a data.frame or list of polygons
###  with coordinates 'x' and 'y'
###  polygons are seperated by NA

### outer border of the 'ocean':
  bx=XYLIM[c(1,2,2,1,1)]
  by=XYLIM[c(3,3,4,4,3)]

  ee=(1:length(data$x))[is.na(data$x)]-1
  if (!is.na(data$x)[length(data$x)]) ee=c(ee,length(data$x))

  px=c(data$x[1:ee[1]],XYLIM[1])
  py=c(data$y[1:ee[1]],XYLIM[3])

  if( length(ee)>1) for(i in 1:(length(ee)-1)){
     px=c(px,data$x[ (ee[i]+2):ee[(i+1)] ],XYLIM[1])
     py=c(py,data$y[ (ee[i]+2):ee[(i+1)] ],XYLIM[3])
  }
  data.frame(x=c(bx,px),y=c(by,py))
}

box.resize=function(glim){
  usr=par('usr')
  plt=par('plt')
#  fin=par('fin')
#  mai=par('mai')

  f1=plt[1]+(glim[1]-usr[1])*(plt[2]-plt[1])/(usr[2]-usr[1])
  f2=plt[2]-(usr[2]-glim[2])*(plt[2]-plt[1])/(usr[2]-usr[1])
  f3=plt[3]+(glim[3]-usr[3])*(plt[4]-plt[3])/(usr[4]-usr[3])
  f4=plt[4]-(usr[4]-glim[4])*(plt[4]-plt[3])/(usr[4]-usr[3])
#  print(c(f1,f2,f3,f4))

#  pin=c(fin[1]*(f2-f1),fin[2]*(f4-f3))
  par(plt=c(f1,f2,f3,f4),usr=glim)
### the following 'box' command is necessary for the par change to take effect
### I don't really know why
  box(lty=0)
}

seamask=function(col='white',add.dx=TRUE,
                 drawmap=TRUE,map.database='worldcoast'){
# This doesn't work very well!!!
  domain=.Last.domain
  opar=par('plt','usr')
  glimits=DomainExtent(domain)
  LXY=c(glimits$x0,glimits$x1,glimits$y0,glimits$y1)
  if(add.dx) LXY=LXY+c(-glimits$dx,glimits$dx,-glimits$dy,glimits$dy)/2

  if(!exists(map.database)) data(list=map.database)
  map <- eval(parse(text=map.database))
  boundaries <- map[map$x >= glimits$lonlim[1] &
                    map$x <= glimits$lonlim[2] &
                    map$y >= glimits$latlim[1] &
                    map$y <= glimits$latlim[2] ,]
  borders=project(boundaries,proj=domain$projection)
  box.resize(LXY)
  invborders=as.ocean(borders)
  polygon(invborders,col=col,border=FALSE)
  par(opar)
  box(lty=0)
  plot(domain,box=TRUE,drawmap=drawmap,add.dx=add.dx)
}

antipolygon <- function(data,col="white",xylim){
### inspired by the version in GEOmap
### this colours everything outside the polygon given by "data"
### up to the limits in xylim.
### typical usage:
### mask=project(map("worldHires","Belgium",fill=0,plot=FALSE))
### xylim=DomainExtent(datafield)
### antipolygon(mask,xylim=xylim)
    x <- c(data$x, data$x[1],xylim$x0,xylim$x0,xylim$x1,xylim$x1,xylim$x0)
    y <- c(data$y, data$y[1],xylim$y0,xylim$y1,xylim$y1,xylim$y0,xylim$y0)
    polygon(x, y, border = col, col = col, xpd = TRUE)
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
  if(is.element("package:RColorBrewer",search()))
    Spectral <- rev(brewer.pal(11,"Spectral"))
  else Spectral <- colors()[c(107,619,12,105,384,383,38,573,504,35,645)]

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
