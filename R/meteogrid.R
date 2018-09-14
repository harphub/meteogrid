#-------------------------------------------#
# Part of R-package meteogrid               #
# Copyright (c) 2003-2018 Alex Deckmyn      #
#   Royal Meteorological Institute, Belgium #
# alex.deckmyn@meteo.be                     #
# Released under GPL-3 license              #
#-------------------------------------------#

######################
### meteogrid
######################

is.geofield <- function(x){
  inherits(x,"geofield")
}

is.geodomain <- function(x){
  inherits(x,"geodomain")
}

as.geodomain <- function(x) {
  if (inherits(x, "geodomain")) return(x)
  else if ("domain" %in% names(attributes(x))) return(attributes(x)$domain)
  else stop("domain not well defined!")
}

print.geofield <- function(x, ...){
  cat(paste(attr(x,"info")$origin,":",attr(x,"info")$name),"\n")
  if (length(dim(x) > 2)) {
    cat("Extra dimensions: ", paste(names(dim(x)), dim(x), sep="=", collapse=", "), "\n")
  }
  cat("Time:\n")
  cat(attr(x,"time"),"\n")
  cat("Domain summary:\n")
  print(attr(x,"domain"))
  cat("Data summary:\n")
  cat(summary(as.vector(x)),"\n")
}

print.geodomain = function(x, ...){
  cat(x$nx,"x",x$ny,"domain\n")
  cat("Projection summary:\n")
  cat("proj=",x$projection$proj,"\n")
  if (!is.null(x$NE)) cat("NE = (",x$NE[1],",",x$NE[2],")\n")
  if (!is.null(x$SW)) cat("SW = (",x$SW[1],",",x$SW[2],")\n")
  if (!is.null(x$clonlat)) cat("centre = (",x$clonlat[1],",",x$clonlat[2],")\n")
#  print.noquote(domain$projection)
}

as.geofield <- function (x=NA, domain, time = attr(domain, "time"),
                         info = attr(domain, "info"),
                         extra_dimensions=NULL) {
  mydomain <- as.geodomain(domain)
  if (is.geofield(x)) return(x)
  if (is.vector(x)) {
    x <- array(x, dim=c(mydomain$nx, mydomain$ny, extra_dimensions))
  } else {
    if (any(dim(x)[1:2] != c(mydomain$nx, mydomain$ny))) stop(
          "Wrong dimensions.", paste(dim(x), sep=" x "), " vs ", 
          mydomain$nx , " x ", mydomain$ny)
  }
  if (!is.null(extra_dimensions)) {
    names(dim(x))[-(1:2)] <- names(extra_dimensions)
  }

  attr(x, "domain") <- mydomain
  attr(x, "time")   <- time
  attr(x, "info")   <- if (is.null(info)) list() else info
  class(x) <- "geofield"
  return(x)
}

# a simple method for >2d geofields using x[[i,j...]]
# max 3 extra dimensions are allowed (e.g. time, level, members).
# NOT COMPLETELY CORRECT
# with 2 extra dimensions z[[n]] should give an error, not z[[n,]]
# but that is (almost?) impossible to solve with R code
# "ellipsis" can't deal with empty dimensions
# I guess that's exactly why "[" is primitive
.subset.geofield <- function(x, i, j, k) {
  dimx <- length(dim(x))
  if (dimx <= 2) stop("Subsetting a geofield requires extra dimensions")
  if (missing(i) && dimx >= 3) i <- 1:dim(x)[3]
  if (missing(j) && dimx >= 4) j <- 1:dim(x)[4]
  if (missing(k) && dimx >= 5) k <- 1:dim(x)[5]

  result <- switch(dimx - 2, x[,,i], x[,,i,j], x[,,i,j,k])
  as.geofield(result, x)
}

# ASSIGNMENT:
# MUCH HARDER: with missing dimensions, how do you find the data?
# only simple with 1 extra dimension
.subset.assign.geofield <- function(x, i, j, k, value) {
  stop("Please use standard subscripting [,,,] for geofield assignment.")
  dimx <- length(dim(x))
  if (dimx <= 2) stop("Subsetting a geofield requires extra dimensions")
  if (missing(i) && dimx >= 3) i <- 1:dim(x)[3]
  if (missing(j) && dimx >= 4) j <- 1:dim(x)[4]
  if (missing(k) && dimx >= 5) k <- 1:dim(x)[5]

  if (dimx==3) x[,,i] <- value
  x
}

# Much harder : using x[,,i]
# there is no "[.default", so you need to "unclass"
# which I suppose makes an extra copy of the whole array
# very messy and buggy, because you can not detect errors
# AVOID
.geosub1 <- function(x,i,j,k,l,m, ..., drop) {
  dimx=length(dim(x))
  if (missing(i) && dimx >= 1) i <- 1:dim(x)[1]
  if (missing(j) && dimx >= 2) j <- 1:dim(x)[2]
  if (missing(k) && dimx >= 3) k <- 1:dim(x)[3]
  if (missing(l) && dimx >= 4) l <- 1:dim(x)[4]
  if (missing(m) && dimx >= 5) m <- 1:dim(x)[5]

  result <- switch(dimx, unclass(x)[i], unclass(x)[i,j],
         unclass(x)[i,j,k], unclass(x)[i,j,k,l],
         unclass(x)[i,j,k,l,m])
  if (missing(i) && missing(j)) result <- as.geofield(result, x)
  result
}

 
