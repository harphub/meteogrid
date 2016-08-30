#--------------------------------------#
# Part of R-package geogrid            #
# Copyright Alex Deckmyn               #
# Released under GPL-3 license         #
#--------------------------------------#

###################
### UPSCALING   ###
###################
### routines for going to a coarser resolution
### either staying in the same grid or regridding
### in these routines: no interpolation, but aggregation (mean or median)
### of the smaller (original) grid boxes to the coarser resolution

### upscaling can be either by a factor (e.g. by two grid cells)
### or to a new grid
upscale <- function(infield, factor=NULL, newdomain=NULL, method="mean", ... ) {
  if (!is.null(factor)) {
    result <- upscale_factor(infield, factor=factor, method=method, ... )
  }
  else {
    result <- upscale_regrid(infield, newdomain=newdomain, method=method, ... )
  }
  result
}

upscale_factor <- function(infield, factor, method="mean", ... ){
  if (length(factor)==1) factor <- c(factor, factor)
  ### define the upscaled domain 
  olddomain <- attr(infield, "domain")
  newdomain <- olddomain
  newdomain$dx <- olddomain$dx * factor[1]
  newdomain$dy <- olddomain$dy * factor[2]
# take the "floor" for the new nx,ny 
  newdomain$nx <- floor(olddomain$nx / factor[1])
  newdomain$ny <- floor(olddomain$ny / factor[2])
# TO DO: currently we stick to the lower left border
# this could be an option...
  sw0 <- project(olddomain$SW, proj=olddomain$projection, inv=FALSE)
  sw1 <- c(sw0$x, sw0$y) + c(olddomain$dx, olddomain$dy) * (factor - 1)/2
  newdomain$SW <- as.numeric(project(sw1[1], sw1[2], proj=newdomain$projection, inv=TRUE))

  ne1 <- sw1 + c( (newdomain$nx-1) * newdomain$dx, (newdomain$ny-1) * newdomain$dy)
  newdomain$NE <- as.numeric(project(ne1[1], ne1[2], proj=newdomain$projection, inv=TRUE))

### the actual data
### from Daan: we do this by reshaping the matrix to a 4d array
### notice that we drop a few points in some cases
### IT IS MUCH FASTER TO DO SUM and devide, than take the MEAN !!!
### in the future, dedicated C code may be even faster (certainly for median...)
  zz <- array(infield[1:(factor[1] * newdomain$nx),1:(factor[2] * newdomain$ny)],
                    c(factor[1], newdomain$nx, factor[2], newdomain$ny))
  if (method=="mean") result <- apply(zz, c(2,4), sum, ...)/(factor[1]*factor[2])
  else result <- apply(zz, c(2,4), eval(parse(text=method)), ...) # this can be SLOW!!!

  as.geofield(result, domain=newdomain)
}

### take the mean value of all cells whose center falls in the new grid cell
upscale_regrid <- function(infield, newdomain, method="mean", ... ) {
  if (!is.geodomain(newdomain)) {
    if (inherits(newdomain, "geofield") || inherits(newdomain, "FAfile")) newdomain <- attributes(newdomain)$domain
    else if (inherits(newdomain, "FAframe") && requireNamespace("Rfa")) newdomain <- Rfa::FAdomain(newdomain)
    else stop("new domain not well defined!")
  }


  gnx <- newdomain$nx
  gny <- newdomain$ny
  if (method != "mean") stop("Only mean is available for upscale regridding.") 
  opoints <- DomainPoints(infield)
#  opoints$value <- as.vector(infield)
  pind <- point.index(as.vector(opoints$lon), as.vector(opoints$lat), domain=newdomain, clip=FALSE)

  result <- .C("upscale_by_mean", npoints=as.integer(prod(dim(infield))),
                                  px=as.integer(round(pind$i)), py=as.integer(round(pind$j)), 
                                  pval=as.numeric(infield),
                                  gnx=as.integer(gnx), gny=as.integer(gny),
                                  gcount=integer(gnx * gny),
                                  gval=numeric(gnx * gny),
                                  NAOK=TRUE, PACKAGE="geogrid")
# we don't really use gcount, but it is available if necessary...
  as.geofield(matrix(result$gval, nrow=gnx), domain=newdomain)
}



