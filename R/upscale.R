###################
### UPSCALING   ###
###################
### routines for going to a coarser resolution
### either staying in the same grid or regridding
### in these routines: no interpolation, but aggregation (mean or median)
### of the smaller (original) grid boxes to the coarser resolution

### upscaling can be either by a factor (e.g. by two grid cells)
### or to a new grid
upscale <- function(infield, factor=NULL, newdomain=NULL, method="mean", weights=0, ... ) {
  if (!is.null(factor)) {
    result <- upscale_factor(infield, factor=factor, method=method, ... )
  }
  else {
    result <- upscale_regrid(infield, newdomain=newdomain, method=method, weights=weights, ... )
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
  newdomain$SW <- as.numeric(project(sw1, proj=newdomain$projection, inv=TRUE))

  ne1 <- sw1 + c( (newdomain$nx-1) * newdomain$dx, (newdomain$ny-1) * newdomain$dy)
  newdomain$NE <- as.numeric(project(ne1, proj=newdomain$projection, inv=TRUE))

## domains defined by their centre should also be supported...
  if (!is.null(olddomain$clonlat)) newdomain$clonlat <- as.numeric(project( (ne1+sw1)/2, proj=newdomain$projection, inv=TRUE))

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

### take the mean value of all cells whose centre falls in the new grid cell
upscale_regrid <- function(infield, newdomain, method="mean", weights=NULL, ... ) {
  if (method != "mean") stop("upscale_regrid only supports method=mean.")
  newdomain <- as.geodomain(newdomain)

  if (is.null(weights)) {
    weights <- upscale_regrid_init(infield, newdomain)
  }

  result <- upscale_by_mean( px=weights$px, py=weights$py, 
                             pval=as.numeric(infield),
                             gnx=newdomain$nx, gny=newdomain$ny)
  print(dim(result))

  as.geofield(result, domain=newdomain)
}

# initialising upscale_regrid: pre-calculate the grid indices
upscale_regrid_init <- function(olddomain, newdomain) {
  olddomain <- as.geodomain(olddomain)
  newdomain <- as.geodomain(newdomain)

  opoints <- DomainPoints(olddomain)
  npoints <- olddomain$nx * olddomain$ny
  pind <- point.index(domain=newdomain, lon=as.vector(opoints$lon), lat=as.vector(opoints$lat), clip=FALSE)

  return(data.frame(px = round(pind$i), py=round(pind$j)))
}

