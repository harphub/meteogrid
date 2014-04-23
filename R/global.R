# basic support for (reduced) gaussian grids
# main purpose is to have a "regrid" function
# EXPERIMENTAL!!!
# the domain is specified by a list of latitudes & NLON
# lon=0 is always one of the points in the reduced grid. CHECK THIS!

# first initialise the weights
point.index.gaussian <- function(lon,lat,domain){
### stretching? Pole rotation?
### after that, you can call C code.
    if(is.geofield(domain)) domain <- attributes(domain)$domain
    NX=length(lat)
    lon <- lon %% 360
    INDEX=.C("gridindex",inlat=lat,inlon=lon,nx=as.integer(NX),
             gridlat=domain$latlist,nlon=domain$nlon,nlat=as.integer(length(domain$latlist)),
             iout1=numeric(NX),iout2=numeric(NX),jout=numeric(NX))
    list(iout1=INDEX$iout1,iout2=INDEX$iout2,jout=INDEX$jout)
}

point.bilin.gaussian.init <- function(lon,lat,domain){
    if(is.geofield(domain)) domain <- attributes(domain)$domain
    index=point.index.gaussian(lon,lat,domain)
## how to treat poles: just use the closest latitude
    N <- length(domain$latlist)
    index$jout[which(index$jout==0)] <- 1
    index$jout[which(index$jout==N+1)] <- N


    fi1 <- floor(index$iout1)
    ci1 <- fi1+1
    di1 <- index$iout1 - fi1
    fi2 <- floor(index$iout2)
    ci2 <- fi2+1
    di2 <- index$iout2 - fi2
    fj  <- floor(index$jout)
    cj  <- fj+1
    dj  <- index$jout - fj
### EXPERIMENTAL!!! CHECK THIS!!!
    w00 <- (1-di1)*(1-dj)
    w01 <- (1-di2)*dj
    w10 <- di1 * (1-dj)
    w11 <- di2 * dj
    list(w00 = w00, w10 = w10, w01 = w01, w11 = w11, F00 = cbind(fi1, 
        fj), F01 = cbind(fi2, cj), F10 = cbind(ci1, fj), F11 = cbind(ci2, 
        cj))
}
### I don't need this function: use point.bilin !!!  
### just add a check for gaussian: the init function is different.
point.bilin.gaussian <- function(lon,lat,infield,weights=NULL){
  domain <- attributes(infield)$domain
  if (is.null(weights)) weights <- point.bilin.gaussian.init(lat,lon,domain)
  return(weights$w00 * infield[weights$F00] + weights$w01 * 
        infield[weights$F01] + weights$w10 * infield[weights$F10] + 
        weights$w11 * infield[weights$F11])
}


### make a default global domain
### mercator or latlon, resolution...
global.lalo <- function(dx=0.5, dy=dx){
  NX<- round(360/dx)
  NY <- round(180/dy)
  Make.domain(projtype="latlong",nxny=c(NX,NY),clonlat=c(0,0),dxdy=c(dx,dy))
}


####################################################3
### SPECTRAL HARMONIC FIELDS

# initialse for 1 single latitude
legendre.init <- function(lat,nsmax){
  legendre <- numeric(((nsmax+1)*(nsmax+2))/2)
#  cat("length=",length(legendre),"\n")
  result <- .C("legendre_init",x=lat,nsmax=as.integer(nsmax),
                legendre=legendre)
  result$legendre
}

legendre.trans <- function(data,nsmax,LAT,reduced=FALSE){
# data has dimension NLAT x (nsmax+1)*(nsmax+2)/2
# for now, we assume the data is conveniently ordered
# in Daan's code: FM[2*m] ~ LF[n*(n+1)/2+m+1]*data[n*n + n +m-1]
#                 FM[2*m+1] ~ LF[n*(n+1)/2+m+1]*data[n*n + n -m-1] (imaginary part)
# so ideally, data is alread in complex form
  NLAT <- length(LAT)
  if(NLAT!=dim(data)[1]) stop("wrong dimension")
  for(ll in 1:NLAT){
    legendre <- legendre.init1(LAT[ll],nsmax)
    res1 <- legendre * data[ll,]
# we have to sum by "m"! Not trivial. INPUT: Data is ordered by n, then by m=-n:n,( then real/im)
# now we do FFT:
    result[ll,] <- Re(fft(res1))
  }
# reduced gaussian grid:
}

legendre.init2 <- function(lat,nsmax){
  NX=length(lat)
  legendre <- numeric(NX*(nsmax+1)*(nsmax+2)/2)
  cat("length=",length(legendre),"\n")
  result <- .C("legendre_init2",x=lat,nx=as.integer(NX),nsmax=as.integer(nsmax),
                legendre=legendre)
  matrix(result$legendre,nrow=NX,byrow=FALSE)
}

legendre.reorder <- function(indata,nsmax){
  NX=length(indata)
  if(NX != ((nsmax+1)*(nsmax+2))/2) stop("NSMAX inconsistent with data length.")
  .C("legendre_reorder",indata=indata,nsmax=as.integer(nsmax),
     outdata=numeric(NX))$outdata
}



