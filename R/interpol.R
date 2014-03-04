##############################
### INTERPOLATION          ###
##############################
### routines for nearest neighbour, bilinear and bicubic interpolation
### regrid is a function for interpolating data from one grid to another
### It simply considers the new grid as a set of lat/lon points for interpolation.
### "init" allows you to return the indices and weights for the interpolation
### "weights" : if this is profided, they are not re-calculated. 
### "mask" lets you define a land/sea mask. only points where mask==TRUE are used. 
###             (FIX ME: not implemented for bicubic)
regrid <- function (infield, newdomain=.Last.domain,method="bilinear",
                    mask=NULL,init=FALSE,weights=NULL)
{
### regridding: either bilinear of nearest neighbour
  olddomain <- attr(infield,"domain")
  if(is.geofield(newdomain)) newdomain <- attr(newdomain,"domain")
  newpoints <- DomainPoints(newdomain)

### Bilinear interpolation
  if(substring(method,1,3)=="bil"){
    result <- point.bilin(as.vector(newpoints$lon),
                             as.vector(newpoints$lat),
                             infield,mask=mask,init=init,weights=weights)
  }
### Nearest neighbour interpolation
  else if (substring(method,1,1)=="n"){
    result <- point.closest(as.vector(newpoints$lon),
                             as.vector(newpoints$lat),
                             infield,mask=mask,init=init,weights=weights)
  }
### Bicubic spline
  else if (substring(method,1,3)=="bic"){
    result <- point.bicubic(as.vector(newpoints$lon),
                             as.vector(newpoints$lat),
                             infield,init=init,weights=weights)
  }
  else stop(paste("Unknown interpolation method",method))

  if(init) return(result)
  else return(
    as.geofield(matrix(result, ncol = newdomain$ny, nrow = newdomain$nx),
                domain = newdomain,time=attr(infield,"time"),
                info=attr(infield,"info")))
}

############ METHODS #############

### fractional indices of points whithin a grid
### clip 
point.index <- function(lon,lat,domain,clip=FALSE){
  if(is.geofield(domain)) domain <- attributes(domain)$domain

  glimits <- DomainExtent(domain)
  projpoints <- project(list(x=lon,y=lat),
                        proj = domain$projection)
  i <- (projpoints$x - glimits$x0)/glimits$dx + 1
  j <- (projpoints$y - glimits$y0)/glimits$dy + 1

  if(clip) { 
    i[i<0.5 | i>glimits$nx+1/2] <- NA
    j[j<0.5 | j>glimits$ny+1/2] <- NA
  }
  list(i=i,j=j)
}


point.interp <- function(lon,lat,method="bilin",...){
  if(substring(method,1,3)=="bil") point.bilin(lon,lat,...)
  else if (substring(method,1,3)=="bic") point.bicubic(lon,lat,...)
  else if (substring(method,1,1)=="n" | substring(method,1,1)=="c") point.closest(lon,lat,...)
  else stop(paste("Unknown interpolation method",method))
}
  
### bilinear interpolation
point.bilin.init <- function(lon,lat,domain,mask=NULL){
  if(is.geofield(domain)) domain <- attributes(domain)$domain
  nx <- domain$nx
  ny <- domain$ny
  index <- point.index(lon,lat,domain)
  fi <- floor(index$i)
  fj <- floor(index$j)
  ci <- fi + 1
  cj <- fj + 1
  di <- index$i - fi
  dj <- index$j - fj
  ci[ci<1|ci>nx] <- NA
  cj[cj<1|cj>ny] <- NA
  fi[fi<1|fi>nx] <- NA
  fj[fj<1|fj>ny] <- NA

# interpolation weights (sum is always 1)
  w00 <- (1-di)*(1-dj)
  w01 <- dj*(1-di)
  w10 <- di*(1-dj)
  w11 <- di*dj
 
  if(!is.null(mask)){
# if some of the 4 points are masked: set weight to zero, and renormalise the remaining weights to sum=1
    w00[!mask[cbind(fi,fj)]] <- 0
    w01[!mask[cbind(fi,cj)]] <- 0
    w10[!mask[cbind(ci,fj)]] <- 0
    w11[!mask[cbind(ci,cj)]] <- 0
    wsum <- w00 + w01 + w10 + w11
    wnn <- (wsum>1.E-9)
    w00[wnn] <- w00[wnn] / wsum[wnn] 
    w01[wnn] <- w01[wnn] / wsum[wnn] 
    w10[wnn] <- w10[wnn] / wsum[wnn] 
    w11[wnn] <- w11[wnn] / wsum[wnn] 
# if all weights are 0, make sure the result is NA, not 0:
    w00[!wnn] <- NA
  }
  list(w00=w00,w10=w10,w01=w01,w11=w11,
       F00=cbind(fi,fj),F01=cbind(fi,cj),
       F10=cbind(ci,fj),F11=cbind(ci,cj))
}


point.bilin <- function(lon,lat,infield,mask=NULL,init=FALSE,weights=NULL)
{
### How to introduce a L/S mask? Must adapt weights.
### That would be slow in R.
  if(is.null(weights)){ 
    weights <- point.bilin.init(lon,lat,infield,mask=mask)
    if(init) return(weights)
  }

  return(weights$w00*infield[weights$F00] + weights$w01*infield[weights$F01] +
         weights$w10*infield[weights$F10] + weights$w11*infield[weights$F11])
}

### nearest neighbour (closest point)
point.closest.init <- function(lon,lat,domain,mask=NULL) {
  if(is.geofield(domain)) domain <- attributes(domain)$domain
  nx <- domain$nx
  ny <- domain$ny
  index <- point.index(lon,lat,domain)
# closest point is just a matter of rounding the index to closest integer
  i <- round(index$i)
  j <- round(index$j)
# boundary issues: 
  i[(i<1)|(i>nx)] <- NA
  j[(j<1)|(j>ny)] <- NA

  if(!is.null(mask)){
# with a mask, we need a bit more work
# we do this only for those points where the closest point is masked...
    ismasked <- !mask[cbind(i,j)]
    nmasked <- sum(ismasked)
    if(nmasked>0){
      fi <- floor(index$i[ismasked])
      fj <- floor(index$j[ismasked])
      ci <- fi + 1
      cj <- fj + 1
      di <- index$i[ismasked] - fi
      dj <- index$j[ismasked] - fj
      ci[ci<1|ci>nx] <- NA
      cj[cj<1|cj>ny] <- NA
      fi[fi<1|fi>nx] <- NA
      fj[fj<1|fj>ny] <- NA
      dist <- data.frame(d00=ifelse(mask[cbind(fi,fj)],(di^2+dj^2),NA),
                         d01=ifelse(mask[cbind(fi,cj)],(di^2+(1-dj)^2),NA),
                         d10=ifelse(mask[cbind(ci,fj)],((1-di^2)+dj^2),NA),
                         d11=ifelse(mask[cbind(ci,cj)],((1-di)^2+(1-dj)^2),NA))
      novalue <- is.na(dist$d00) &  is.na(dist$d01) & is.na(dist$d10) & is.na(dist$d11)
      if(sum(novalue)==nmasked) result[ismasked] <- NA
      else {
        closest <- rep(NA,nmasked)
        closest[!novalue] <- apply(dist[!novalue,],1,which.min)
        mi <- ifelse(closest <= 2,fi,ci)
        mj <- ifelse(closest%%2 == 1,fj,cj)
        i[ismasked] <- mi
        j[ismasked] <- mj
      }
    }
  }
  return(list(index=cbind(i,j)))
}

point.closest <- function(lon,lat,infield,mask=NULL,init=FALSE,weights=NULL){
### where mask=FALSE: take next closest point
### but only from the four closest !
  if(is.null(weights)){
    weights <- point.closest.init(lon,lat,infield,mask)
    if(init) return(weights)
  }
  return(infield[weights$index])
}


### bicubic spline interpolation

### 1D routine (for testing purposes)
interp.cubic <- function(i,data){
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

  return(result)
}

point.bicubic.init <- function(lon,lat,domain,mask=NULL){
### there's actually not much initialisation that you can do
### "weights" in fact only contains the index of neighbouring points
### because the weights are computed from the field values.
### But it may make a little difference.
  if(!is.null(mask)) stop("Mask not (yet) supported in bicubic interpolation.")
  if(is.geofield(domain)) domain <- attributes(domain)$domain

  nx <- domain$nx
  ny <- domain$ny
  index <- point.index(lon,lat,domain)

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
  return(list(ffi=ffi,fi=fi,ci=ci,cci=cci,ffj=ffj,fj=fj,cj=cj,ccj=ccj))
}


point.bicubic <- function(lon,lat,infield,init=FALSE,weights=NULL,mask=NULL){
  if(is.null(weights)){
    weights <- point.bicubic.init(lon,lat,infield,mask)
    if(init) return(weights)
  }
    
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

# For mask, we must calculate the 16 weights, set some to 0, normalise weights to sum=1
# But how much sence does it make to do masked bicubic?

  result <- di*(-di^2/2+di-1/2)*FFi + di^2*(di-1)/2*CCi + 
             (1+3/2*di^3 -5/2*di^2) * Fi + di*(-3/2*di^2+2*di+1/2)*Ci

  return(result)
}
##########################################3
