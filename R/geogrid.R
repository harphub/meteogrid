#-------------------------------------------#
# Part of R-package geogrid                 #
# Copyright (c) 2003-2017 Alex Deckmyn      #
#   Royal Meteorological Institute, Belgium #
# alex.deckmyn@meteo.be                     #
# Released under GPL-3 license              #
#-------------------------------------------#

######################
### geogrid
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

as.geofield <- function (x=NA, domain, time = "", info = list()) {
  if (is.geofield(x)) return(x)
  if (is.vector(x)) x <- array(x,dim=c(domain$nx,domain$ny))
  attr(x,"domain") <- domain
  attr(x,"time") <- time
  attr(x,"info") <- info
  class(x) <- c("geofield",class(x))
  return(x)
}

