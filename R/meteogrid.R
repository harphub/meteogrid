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
  if (length(dim(x)) > 2) {
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
                         extra_dim=list()) {
  mydomain <- as.geodomain(domain)
  if (is.geofield(x)) return(x)

  # fix the dimensions  
  if (!is.list(extra_dim)) {
    dims <- c("x"=mydomain$nx, "y"=mydomain$ny, extra_dim)
    extra_dim <- lapply(extra_dim, function(x) 1:x)
    set_dimnames <- FALSE
  } else {
    dims <- c("x"=mydomain$nx, "y"=mydomain$ny, 
              vapply(extra_dim, FUN=length, FUN.VALUE=1))
    set_dimnames <- (length(dims) > 2)
  }

  if (is.vector(x)) {
    x <- array(x, dim=dims)
  } else {
    # extra_dimensions is probably not given in this case
    if (length(dim(x))>2 && length(extra_dim)==0) {
	    dims <- c(dims, dim(x)[-(1:2)])
    }
    # check all dimensions
    if (any(dim(x) != dims)) stop("Wrong dimensions.", 
				     paste(dim(x), collapse=" x "), " vs ",
				     paste(dims, collapse=" x "))
  }
  noname <- which(is.na(names(dims)) | names(dims)=="" )
  names(dims)[noname] <- paste0("D", noname)
  names(dim(x)) <- names(dims)
  # also set dimnames (not for x,y!)
  if (set_dimnames) dimnames(x) <- c(list(x=NULL, y=NULL), extra_dim)

  attr(x, "domain") <- mydomain
  attr(x, "time")   <- time
  attr(x, "info")   <- if (is.null(info)) list() else info
  class(x) <- "geofield"
  return(x)
}

# a simple method for >2d geofields using x[[i,j...]]
# max 3 extra dimensions are allowed (e.g. time, level, members).
# ATTENTION: with 2 extra dimensions z[[n]] gives z[[n,]]
#            it should throw an error
# but that is (almost?) impossible to solve with R code
# e.g. ellipsis (...) can't deal with empty dimensions
# I guess that's exactly why "[" is a primitive
.subset.geofield <- function(x, i, j, k) {
  dimx <- length(dim(x))
  if (dimx <= 2) stop("Subsetting a geofield requires extra dimensions")
  if (dimx > 3) warning("Working more than 3 dimensions in a geofield is dangerous.")
  if ( (dimx < 5 && !missing(k)) || (dimx < 4 && !missing(j))) stop("Bad dimensioning.")
  if (missing(i) && dimx >= 3) i <- 1:dim(x)[3]
  if (missing(j) && dimx >= 4) j <- 1:dim(x)[4]
  if (missing(k) && dimx >= 5) k <- 1:dim(x)[5]
  result <- as.geofield(switch(dimx - 2, x[,,i], x[,,i,j], x[,,i,j,k]),
			domain=x)
  # when dropping a dimension (e.g. length(i)=1), you fix a level/member/...
  # so that should be in the info field
  if (dimx >= 3 && length(i)==1) attributes(result)$info[[names(dim(x))[3] ]] <- dimnames(x)[[3]][i]
  if (dimx >= 4 && length(j)==1) attributes(result)$info[[names(dim(x))[4] ]] <- dimnames(x)[[4]][j]
  if (dimx == 5 && length(k)==1) attributes(result)$info[[names(dim(x))[5] ]] <- dimnames(x)[[5]][k]
  result

}

# ASSIGNMENT:
# MUCH HARDER: with missing dimensions, how do you find the data?
# only simple with 1 extra dimension
.subset.assign.geofield <- function(x, i, value) {
  dimx <- length(dim(x))
  if (dimx <= 2) stop("Subsetting a geofield requires extra dimensions")
  if (dimx > 3) stop("Please use standard subscripting [,,,] for geofield assignment.")
  x[,,i] <- value
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

 
