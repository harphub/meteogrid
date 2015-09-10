# Some obsolete functions I would rather eliminate
# but which are requested by people too lazy to update their code
# BUT: they will not be documented!

lalopoint <- function(data,lon,lat,minimise='lalo',mask=NULL){
### find the closest domain point to the given co-ordinates
### by minimising distance in either LonLat or projected co-ordinates
### This is in fact not exactly the same as minimising geographical distance!
  ldata <- is.geofield(data)
  if (ldata) domain <- attr(data, "domain")
  else domain <- data

  lalodomain <- DomainPoints(domain,"lalo")
  nx <- domain$nx
  ny <- domain$ny

  if(minimise=='lalo') dist <- sqrt((lalodomain$lon-lon)^2+(lalodomain$lat-lat)^2)
  else if (minimise=='proj') {
    xydomain <- DomainPoints(domain, "xy")
    ppp <- project(list(x=lon,y=lat), proj = domain$projection)
    dist <- sqrt((xydomain$x - ppp$x)^2 + (xydomain$y - ppp$y)^2)
  }
  else stop("unknown minimisation.")

  if (!is.null(mask)) {
      if (is.character(mask))
          mask <- eval(parse(text = mask))
      else if (is.expression(mask))
          mask <- eval(mask)
      dist[eval(expression(mask))] <- NA
  }
# We wouldn't expect two points to be exactly as close (numerically)
# But you never know. So in this case we chose at random.
# Default method is to set the index of both to 1.5, which screws up everything.

#  dist.rank[] <- rank(dist,ties.method="random",na.last=TRUE)
#  ij <- which(dist.rank==1,arr.ind=TRUE)

  ij <- which(dist==min(dist,na.rm=TRUE),arr.ind=TRUE)
  i <- ij[1,1]
  j <- ij[1,2]

  il <- lalodomain$lon[i,j]
  jl <- lalodomain$lat[i,j]

  if(ldata)
    list(data = data[i, j], lonlat = c(il, jl), index = c(i,j))
  else
    list(data = NA, lonlat = c(il, jl), index = c(i,j) )
}
