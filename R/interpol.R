##############################
### INTERPOLATION          ###
##############################
### routines for nearest neighbour, bilinear and bicubic interpolation
### regrid is a function for interpolating data from one grid to another
### It simply considers the new grid as a set of lat/lon points for interpolation.
### "init" allows you to return the indices and weights for the interpolation
### "weights" : if this is provided, they are not re-calculated.
### "mask" lets you define a land/sea mask. only points where mask==TRUE are used.
###             (FIX ME: not implemented for bicubic)

### regrid with a single L/S mask will give nonsense (sea points get NA), needs a newmask too 
### and even then, you MUST avoid NA's caused by "errors" in the LSM of "newdomain"
### That means an extra option "force=FALSE".
regrid <- function (infield, newdomain=.Last.domain(), method="bilin",
                    mask=NULL, newmask=NULL, weights=NULL)
{
### regridding: bilinear, bi-cubic or nearest neighbour, and now also upscaling by mean, subgrid
  if (!is.null(weights)) {
    if (!is.null(attr(weights, "newdomain")))  {
      newdomain <- attr(weights, "newdomain")
    } else {
      newdomain <- as.geodomain(newdomain)
    }
    if (!is.null(attr(weights, "method")))  {
      method <- attr(weights, "method")
    }
  }

  if (method %in% c("mean", "median")) {
    return(upscale_regrid(infield=infield, newdomain=newdomain, method=method, weights=weights))
  } else if (method == "subgrid") {
    if (is.null(weights)) weights <- regrid.init(olddomain=infield, newdomain=newdomain, method="subgrid")
    if (length(dim(infield)) > 2) {
      edim(infield)[-(1:2)]
      result <- apply(infield, 3:length(dim(infield)),
                      function(x) x[weights$i0:weights$i1, weights$j0:weights$j1])
      dim(result) <- c(newdomain$nx, newdomain$ny, dim(infield)[-(1:2)])
    } else {
      result <- infield[weights$i0:weights$i1, weights$j0:weights$j1]
      edim <- NULL
    }
    return(as.geofield(result,
                       domain = newdomain,
                       info   = attr(infield, "info"),
                       extra_dim = edim))

  } else {
# speed-up: if you already have the weights, no need to calculate lon/lat of the new domain!
    if (is.null(weights)) weights <- regrid.init(olddomain=infield, newdomain=newdomain, method=method, 
                                                 mask=mask, newmask=newmask)

    result <- point.interp(infield=infield, method=method, weights=weights)
    # 3d data
#    cat("dimensions", paste(dim(infield),sep="x"), paste(dim(result), sep="x"), "\n")
    edim <- if (length(dim(infield)) > 2) dim(infield)[-(1:2)] else NULL
    return(as.geofield(as.numeric(result),
                       domain = newdomain,
                       info   = attr(infield, "info"),
                       extra_dim = edim))
  }
}

regrid.init <- function (olddomain, newdomain=.Last.domain(),
                         method="bilin", mask=NULL, newmask=NULL)
{
### regridding: either bilinear of nearest neighbour
### olddomain and newdomain may be geofields 
  olddomain <- as.geodomain(olddomain)
  newdomain <- as.geodomain(newdomain)
  if (method %in% c("mean")) {
    result <- upscale_regrid_init(olddomain, newdomain)
  } else if (method %in% c("median")) {
    stop("Method", method, "not yet supported. Sorry.")
  } else if (method == "subgrid") {
    # identify how newdomain is a subgrid of olddomain
    if (newdomain$dx != olddomain$dx || newdomain$dy != olddomain$dy) {
      stop("Grids do not conform. Clearly not a subgrid.")
    }
    odp <- DomainPoints(olddomain, type="xy")
    ndp <- DomainPoints(newdomain, type="xy") # you only need SW and NE point...
    i0 <- which.min(abs(odp$x[,1] - ndp$x[1,1]))
    i1 <- which.min(abs(odp$x[,1] - ndp$x[newdomain$nx, newdomain$ny]))
    j0 <- which.min(abs(odp$y[1,] - ndp$y[1,1]))
    j1 <- which.min(abs(odp$y[1,] - ndp$y[newdomain$nx, newdomain$ny]))
    if (i1 - i0 + 1 != newdomain$nx || j1 - j0 + 1 != newdomain$ny) {
      stop("Domains do not correspond. Not a subgrid.")
    }
    result <- list(i0 = i0, i1 = i1, j0 = j0, j1 = j1)
    attr(result, "method") <- "subgrid"
  } else {
    newpoints <- DomainPoints(newdomain)
    if (is.null(mask) != is.null(newmask)) stop("When using Land/Sea masks, you *must* provide both domains!")
    result <- point.interp.init(domain=olddomain,
                                lon=as.vector(newpoints$lon),lat=as.vector(newpoints$lat),
                                method=method, mask=mask,
                                pointmask=as.vector(newmask), force=FALSE)
  }
  attr(result, "olddomain") <- olddomain
  attr(result, "newdomain") <- newdomain
  result
}


############ METHODS #############

### fractional indices of points whithin a grid
### clip means: only consider points that are less than half a grid box out of the domain
point.index <- function(domain=.Last.domain(), lon, lat, clip=TRUE){
  domain <- as.geodomain(domain)

  if (missing(lat)){
    if (is.matrix(lon)){
      lat <- lon[,2]
      lon <- lon[,1]
    }
    else if (is.list(lon)){
      lat <- lon[[2]]
      lon <- lon[[1]]
    }
    else stop("lat is missing!")
  }
  glimits <- DomainExtent(domain)
  projpoints <- project(list(x=lon, y=lat), proj = domain$projection)
  i <- (projpoints$x - glimits$x0)/glimits$dx + 1
  j <- (projpoints$y - glimits$y0)/glimits$dy + 1

  if (clip) {
    i[i<0.5 | i>glimits$nx+1/2] <- NA
    j[j<0.5 | j>glimits$ny+1/2] <- NA
  }
  data.frame(i=i, j=j)
}


point.interp <- function(infield, lon=NULL, lat=NULL, method="bilin",
                         mask=NULL, pointmask=NULL, force=FALSE, weights=NULL){
  if (!is.null(weights) && !is.null(attr(weights, "method"))) method <- attr(weights, "method")
  if (substring(method, 1, 3)=="bil") {
    point.bilin(infield=infield, lon=lon, lat=lat, mask=mask,
                pointmask=pointmask, force=force, weights=weights)
  } else if (substring(method,1,3)=="bic") {
    if (!is.null(mask) | !is.null(pointmask) | force) {
      warning("Mask is not supported for bicubic interpolation.")
    }
    point.bicubic(infield=infield, lon=lon, lat=lat, weights=weights)
  } else if (is.element(substring(method,1,1),c("n","c"))) {
    point.closest(infield=infield, lon=lon, lat=lat, mask=mask, pointmask=pointmask,
                  force=force, weights=weights)
  } else stop(paste("Unknown interpolation method",method))
}

point.interp.init <- function(domain=.Last.domain(), lon, lat, 
                              method="bilin", mask=NULL, pointmask=NULL, 
                              force=FALSE){
  if (substring(method,1,3)=="bil") {
    point.bilin.init(domain=domain, lon=lon, lat=lat, mask=mask, pointmask=pointmask, force=force)
  } else if (substring(method,1,3)=="bic") {
    if (!is.null(mask) | !is.null(pointmask) | force) warning("Mask is not supported for bicubic interpolation.")
    point.bicubic.init(domain=domain, lon=lon, lat=lat) 
  } else if (is.element(substring(method,1,1),c("n","c"))) {
    point.closest.init(domain=domain, lon=lon, lat=lat, mask=mask, pointmask=pointmask, force=force)
  } else stop(paste("Unknown interpolation method",method))
}


### bilinear interpolation
point.bilin.init <- function(domain=.Last.domain(), lon, lat,
                             mask=NULL, pointmask=NULL, force=FALSE){
  domain <- as.geodomain(domain)

  nx <- domain$nx
  ny <- domain$ny
  index <- point.index(domain=domain, lon=lon, lat=lat)
  fi <- floor(index$i)
  fj <- floor(index$j)
  fi[fi<1 | fi>(nx-1)] <- NA
  fj[fj<1 | fj>(ny-1)] <- NA
  ci <- fi + 1
  cj <- fj + 1
#  ci[ci<2 | ci>nx] <- NA
#  cj[cj<2 | cj>ny] <- NA
  not_in_domain <- (is.na(ci) | is.na(cj) | is.na(fi) | is.na(fj))
# interpolation weights (sum is always 1 or NA)
  di <- index$i - fi
  dj <- index$j - fj
  w <- matrix(0, ncol=4, nrow=length(lon))
  w[,1] <- (1-di)*(1-dj)  # w00
  w[,2] <- dj*(1-di)      # w01
  w[,3] <- di*(1-dj)      # w10
  w[,4] <- di*dj          # w11

  if (any(not_in_domain)) {
    message("Warning: ", sum(not_in_domain), " points are outside of the domain.")
  }

  if (!is.null(mask)) {
    # if some of the 4 points are masked: set weight to zero, and renormalise the remaining weights to sum=1
    if (is.null(pointmask)) pointmask <- rep(1, length(fi))
    w[mask[cbind(fi,fj)] != pointmask, 1] <- 0
    w[mask[cbind(fi,cj)] != pointmask, 2] <- 0
    w[mask[cbind(ci,fj)] != pointmask, 3] <- 0
    w[mask[cbind(ci,cj)] != pointmask, 4] <- 0
    wsum <- rowSums(w) # w00 + w01 + w10 + w11
    # if some indices are NA, this should not crash, so use which()!
    # it should be OK to use 0 and not some delta=1.E-9 :
    wnn <- which(wsum > 0 & wsum < 1. - 1.E-9)
    w[wnn,1] <- w[wnn,1] / wsum[wnn]
    w[wnn,2] <- w[wnn,2] / wsum[wnn]
    w[wnn,3] <- w[wnn,3] / wsum[wnn]
    w[wnn,4] <- w[wnn,4] / wsum[wnn]
# if all 4 weights are 0, make sure the result is NA, not 0:
# unless force==FALSE: then take the original weights
    if (any(wsum == 0, na.rm=TRUE)) {
      w0 <- which(wsum == 0)
      message("Warning: ", length(w0), "points have no source points with similar L/S mask!")
      if (force) {
        message("Set to NA.")
        w[w0,] <- NA
      } else { # restore original
        message("Set to original (4 point) interpolation.")
        ddi=di[w0]
        ddj=dj[w0]
        w[w0, 1] <- (1-ddi)*(1-ddj)
        w[w0, 2] <- ddj*(1-ddi)
        w[w0, 3] <- ddi*(1-ddj)
        w[w0, 4] <- ddi*ddj
      }
    }
  }
#  list(w00=w00, w10=w10, w01=w01, w11=w11,
#       F00=cbind(fi,fj), F01=cbind(fi,cj),
#       F10=cbind(ci,fj), F11=cbind(ci,cj))
  result <- data.frame(w00=w[,1], w01=w[,2], w10=w[,3], w11=w[,4],
      F00=I(cbind(fi,fj)), F01=I(cbind(fi,cj)),
      F10=I(cbind(ci,fj)), F11=I(cbind(ci,cj)))
  attr(result, "method") <- "bilin"
  attr(result, "mask") <- is.null(mask)
  result
}


point.bilin <- function(infield, lon=NULL, lat=NULL, mask=NULL, 
                        pointmask=NULL, force=FALSE, weights=NULL)
{
  if (is.null(weights)){
## for gaussian grid: call different init function!
#    if(inherits(infield,"gaussian")) weights <- point.bilin.gaussian.init(lon,lat,infield)
#    else
    weights <- point.bilin.init(infield, lon=lon, lat=lat,
                                mask=mask, pointmask=pointmask, force=force)
  }
  ndim <- length(dim(infield))
  if (ndim==2) {
    weights$w00*infield[weights$F00] + weights$w01*infield[weights$F01] +
              weights$w10*infield[weights$F10] + weights$w11*infield[weights$F11]
  } else {
    result <- apply(infield, 3:ndim, function(x)  
                    weights$w00*x[weights$F00] + weights$w01*x[weights$F01] +
                    weights$w10*x[weights$F10] + weights$w11*x[weights$F11])
    result
  }
}

### nearest neighbour (closest point)
point.closest.init <- function(domain=.Last.domain(), lon, lat,
                               mask=NULL, pointmask=NULL, force=FALSE) {
  domain <- as.geodomain(domain)

  nx <- domain$nx
  ny <- domain$ny
  index <- point.index(domain=domain, lon=lon, lat=lat)
# closest point is just a matter of rounding the index to closest integer
  i <- round(index$i)
  j <- round(index$j)
# boundary issues:
  i[(i<1)|(i>nx)] <- NA
  j[(j<1)|(j>ny)] <- NA

  if (!is.null(mask)) {
    if (is.null(pointmask)) pointmask=rep(1, length(i))
# with a mask, we need a bit more work
# we do this only for those points where the closest point is masked...
    ismasked <- which(mask[cbind(i,j)] != pointmask)
    nmasked <- length(ismasked)
    if (nmasked>0) {
      fi <- floor(index$i[ismasked])
      fj <- floor(index$j[ismasked])
      ci <- fi + 1
      cj <- fj + 1
      di <- index$i[ismasked] - fi
      dj <- index$j[ismasked] - fj
      ci[ ci<1 | ci>nx ] <- NA
      cj[ cj<1 | cj>ny ] <- NA
      fi[ fi<1 | fi>nx ] <- NA
      fj[ fj<1 | fj>ny ] <- NA
      dist <- data.frame(d00=ifelse(mask[cbind(fi,fj)] == pointmask[ismasked], di^2+dj^2,         NA),
                         d01=ifelse(mask[cbind(fi,cj)] == pointmask[ismasked], di^2+(1-dj)^2,     NA),
                         d10=ifelse(mask[cbind(ci,fj)] == pointmask[ismasked], (1-di^2)+dj^2,     NA),
                         d11=ifelse(mask[cbind(ci,cj)] == pointmask[ismasked], (1-di)^2+(1-dj)^2, NA))
      novalue <- is.na(dist$d00) &  is.na(dist$d01) & is.na(dist$d10) & is.na(dist$d11)
# if we don't force NA (=> never NA inside the domain), just leave the original closest point
      if (force && any(novalue)) {
        i[ismasked][novalue] <- NA
        j[ismasked][novalue] <- NA
      }
      if (any(!novalue)) {
        closest <- apply(dist[!novalue, ], 1, which.min)
        mi <- ifelse(closest <= 2, fi, ci)
        mj <- ifelse(closest%%2 == 1, fj, cj)
        i[ismasked][!novalue] <- mi
        j[ismasked][!novalue] <- mj
      }
    }
  }
#  list(index=cbind(i,j))
  result <- data.frame(i=i, j=j)
  attr(result, "method") <- "closest"
  attr(result, "mask") <- !is.null(mask)
  result
}

point.closest <- function(infield, lon=NULL, lat=NULL, mask=NULL, 
                          pointmask=NULL, force=FALSE, weights=NULL){
### where mask != pointmask (default 1): take next closest point
### but only from the four closest !
  if (is.null(weights)) weights <- point.closest.init(infield, lon=lon, lat=lat, mask, pointmask, force)
  index <- cbind(weights$i, weights$j)
  ndim <- length(dim(infield))
  if (ndim==2) infield[index]
  else apply(infield, 3:ndim, function(x) x[index])
}


### bicubic spline interpolation

### 1D routine (for testing purposes)
.interp.cubic.1D <- function(i, data){
  nx <- length(data)
  fi <- floor(i)
  di <- i - fi
  ci <- ceiling(i)
# if we are exactly on a point, floor and ceiling are the same
  ci[di==0] <- ci[di==0] + 1
  ffi <- fi - 1
  cci <- ci + 1

# FIX ME: simplest border handling
  cci[cci<1] <- 1
  ffi[ffi<1] <- 1
  cci[cci>nx] <- nx
  ffi[ffi>nx] <- nx

  result <- di*(-di^2/2+di-1/2)*data[ffi] + di^2*(di-1)/2*data[cci] +
             (1+3/2*di^3 -5/2*di^2) * data[fi] + di*(-3/2*di^2+2*di+1/2)*data[ci]

  result
}

point.bicubic.init <- function(domain=.Last.domain(), lon, lat){
### there's actually not much initialisation that you can do
### "weights" in fact only contains the index of neighbouring points
### because the weights are computed from the field values.
### But it may make a little difference.
  domain <- as.geodomain(domain)

  nx <- domain$nx
  ny <- domain$ny
  index <- point.index(domain=domain, lon=lon, lat=lat)

  fi <- floor(index$i)
  di <- index$i - fi
  ci <- fi + 1
  ffi <- fi - 1
  cci <- ci + 1
# FIX ME: simplest border handling
  cci[cci<1] <- 1
  ffi[ffi<1] <- 1
  cci[cci>nx] <- nx
  ffi[ffi>nx] <- nx

  fj <- floor(index$j)
  dj <- index$j - fj
  cj <- fj + 1
  ffj <- fj - 1
  ccj <- cj + 1

# FIX ME: simplest border handling
  ccj[ccj<1] <- 1
  ffj[ffj<1] <- 1
  ccj[ccj>ny] <- ny
  ffj[ffj>ny] <- ny
  result <- data.frame(di=di, ffi=ffi, fi=fi, ci=ci, cci=cci,
                       dj=dj, ffj=ffj, fj=fj, cj=cj, ccj=ccj)
  attr(result, "method") <- "bicubic"
  result
}


point.bicubic <- function(infield, lon=NULL, lat=NULL, weights=NULL){
  ndim <- length(dim(infield))
  if (ndim != 2) stop("bicubic interpolation only for 2d fields.")
  if (is.null(weights)) weights <- point.bicubic.init(infield, lon=lon, lat=lat)
  di <- weights$di
  dj <- weights$dj
  FFi <-  dj*(-dj^2/2+dj-1/2)*infield[cbind(weights$ffi,weights$ffj)] +
          dj^2*(dj-1)/2*infield[cbind(weights$ffi,weights$ccj)] +
          (1+3/2*dj^3 -5/2*dj^2) * infield[cbind(weights$ffi,weights$fj)] +
          dj*(-3/2*dj^2+2*dj+1/2)*infield[cbind(weights$ffi,weights$cj)]

  Fi <-  dj*(-dj^2/2+dj-1/2)*infield[cbind(weights$fi,weights$ffj)] +
         dj^2*(dj-1)/2*infield[cbind(weights$fi,weights$ccj)] +
         (1+3/2*dj^3 -5/2*dj^2) * infield[cbind(weights$fi,weights$fj)] +
         dj*(-3/2*dj^2+2*dj+1/2)*infield[cbind(weights$fi,weights$cj)]

  Ci <-  dj*(-dj^2/2+dj-1/2)*infield[cbind(weights$ci,weights$ffj)] +
         dj^2*(dj-1)/2*infield[cbind(weights$ci,weights$ccj)] +
         (1+3/2*dj^3 -5/2*dj^2) * infield[cbind(weights$ci,weights$fj)] +
         dj*(-3/2*dj^2+2*dj+1/2)*infield[cbind(weights$ci,weights$cj)]

  CCi <- dj*(-dj^2/2+dj-1/2)*infield[cbind(weights$cci,weights$ffj)] +
         dj^2*(dj-1)/2*infield[cbind(weights$cci,weights$ccj)] +
         (1+3/2*dj^3 -5/2*dj^2) * infield[cbind(weights$cci,weights$fj)] +
         dj*(-3/2*dj^2+2*dj+1/2)*infield[cbind(weights$cci,weights$cj)]

  result <- di*(-di^2/2+di-1/2)*FFi + di^2*(di-1)/2*CCi +
             (1+3/2*di^3 -5/2*di^2) * Fi + di*(-3/2*di^2+2*di+1/2)*Ci
  result
}

