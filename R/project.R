#-------------------------------------------#
# Part of R-package geogrid                 #
# Copyright (c) 2003-2016 Alex Deckmyn      #
#   Royal Meteorological Institute, Belgium #
# Released under GPL-3 license              #
#-------------------------------------------#

##################################
### Interface to PROJ4 library ###
##################################

project <- function(x, y, proj=.Last.domain()$projection, inv=FALSE)
{
  if (missing(y)) {
# apparantly, is.list gives TRUE for data.frames, but let's be careful:
    if (is.list(x) | is.data.frame(x)) {
#      y <- x$y
#      x <- x$x
### we assume the first two list elements are longitude & latitude
### whatever the actual names (could be x/y, lon/lat, Lon/Lat...)
      y <- as.numeric(x[[2]])
      x <- as.numeric(x[[1]])
      if (is.null(x) | is.null(y)) stop("No valid x and y co-ordinates found.")
    }
    else if (is.vector(x)) {
      y <- x[2]
      x <- x[1]
    } else {
      y <- x[,2]
      x <- x[,1]
    }
  }
  if (length(x) != length(y)) stop("x and y don't have the same length.")
  if (!is.numeric(x) | !is.numeric(y)) stop("Co-ordinates must be numeric values.")
  if (missing(proj)) {
    if(!is.null(.Last.domain())) proj <- .Last.domain()$projection
    else return("No projection.")
  }

  if (proj$proj == "latlong") {
### longitude should be in the interval [-180,180[
### unless e.g. if my global data is on a globe [0,360[
### we assume that the meridian = MinLon + 180
### so meridian-180 must not be transported.
    meridian <- if (is.null(proj$lon0)) 0 else proj$lon0
#    x <- ifelse(x < meridian-180,x+360,x)
#    x <- ifelse(x >= meridian+180,x-360,x)
## much faster (!):
    x[which(x <  (meridian-180))] <- x[which(x <  (meridian-180))] + 360
    x[which(x >= (meridian+180))] <- x[which(x >= (meridian+180))] - 360
    
    data.frame(x=x,y=y)
  } else  {
    npoints <- as.integer(length(x))
    npar <- as.integer(length(proj))
    par <- proj4.list2str(proj)
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

#################
### utilities ###
#################

### clip a map exactly to the domain borders:
.map.restrict1 <- function(bx,by,xlim,xperiod=NA_real_,xfrac=0.5){
  lx <- length(bx)

  result <- .C("maprestrict1",bx=as.double(bx),by=as.double(by),lx=as.integer(lx),
            x0=as.double(xlim[1]),x1=as.double(xlim[2]),
            nx=double(2*lx),ny=double(2*lx),
            newlen=integer(1),xperiod=as.numeric(xperiod),
            xfrac=as.numeric(xfrac),NAOK=TRUE)
  data.frame(x=result$nx[2:result$newlen],y=result$ny[2:result$newlen])
}

map.restrict <- function(bxy,xlim,ylim,xperiod=NA_real_,xfrac=0.5,yperiod=NA_real_,yfrac=NA_real_){
  bxy <- .map.restrict1(bxy$x,bxy$y,xlim,xperiod,xfrac)
  byx <- .map.restrict1(bxy$y, bxy$x, ylim, yperiod, yfrac)
  list(x = byx$y, y = byx$x)
}
######################
### periodicity depends on the projection
### may be 360, 2*pi or ~40.000 in Mercator !
### a non-periodic LatLon domain still has a "period" of 360 !

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

### turn a proj4 string into a list
proj4.str2list <- function(pp){
  tryNum <- function(x) if ( !is.na(suppressWarnings(as.numeric(x))) ) as.numeric(x) else x
  p1 <- strsplit(pp,"+",fixed=TRUE)[[1]][-1]  # vector of "x=n" strings
  p2 <- gsub(' ','',p1) # remove all blanks
  p3 <- strsplit(p2,"=")
  parNam <- vapply(p3, function(x) x[1], FUN.VALUE="a")
  prj <- lapply(p3, function(x) if (length(x)==1) NA else tryNum(x[2]))
  names(prj) <- parNam
  prj
}

proj4.list2str <- function(pp){
  paste(names(pp), lapply(pp, function(x) if (is.na(x)) "" else paste("=",x,sep="")),sep="")
}


