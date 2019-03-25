######################
### meteogrid
######################
is.geodomain <- function(x){
  inherits(x,"geodomain")
}

as.geodomain <- function(x) {
  if (inherits(x, "geodomain")) return(x)
  else if (!is.null(attr(x, "domain"))) return(attr(x, "domain"))
  else stop("domain not well defined!")
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


