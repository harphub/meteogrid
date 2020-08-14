
#################
### utilities ###
#################

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
  radius <- domain$projection$R
  if (is.null(radius)) radius <- domain$projection$a
  if (is.null(radius)) radius <- NA_real_
  xper <- switch(domain$projection$proj,
            "latlong"   = 360,
            "ob_tran"   = 2*pi, ### BUG: this could be any oblique projection
            "merc"     = 2 * pi * radius,
            "omerc"    = 2 * pi * radius,
            "tmerc"    = 2 * pi * radius,
            "somerc"   = 2 * pi * radius,
            NA_real_)
  yper <- NA_real_
  list(xper=xper,yper=yper)
}

### turn a proj4 string into a list
proj4.str2list <- function(pp){
  tryNum <- function(x) if ( !is.na(suppressWarnings(as.numeric(x))) ) as.numeric(x) else x
  # split at every + preceded by white space or at start of string
  p1 <- strsplit(pp,"[ ]+[+]|^[+]")[[1]][-1] # vector of "x=n" strings
  p2 <- gsub(' ','',p1) # remove all blanks
  p3 <- strsplit(p2,"=")
  parNam <- vapply(p3, function(x) x[1], FUN.VALUE="a")
  prj <- lapply(p3, function(x) if (length(x)==1) NA else tryNum(x[2]))
  names(prj) <- parNam
  prj
}

#' export
proj4.list2str <- function(pp, join=TRUE){
  result <- paste0(names(pp), 
                   lapply(pp, function(x) if (is.na(x)) "" else paste0("=",x)))
  if (join) result <- paste0("+", result, collapse=" ")
  result
}


