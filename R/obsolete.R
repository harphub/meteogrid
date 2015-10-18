# Some obsolete functions I would rather eliminate
# but which are requested by people too lazy to update their code

lalopoint <- function(geo,lon,lat,minimise='proj',mask=NULL){
### find the closest domain point to the given co-ordinates
  if (minimise!='proj') return(.lalopoint0(geo,lon,lat,minimise,mask))
  ldata <- inherits(geo,"geofield")

  if (ldata ) domain <- attr(geo, "domain")
  else if (inherits(data,"geodomain")) domain <- geo
  else if ("domain" %in% names(attributes(geo))) domain <- attr(geo,"domain")
  else stop("geo is not a geographical object. Can not interprete.")

  lalodomain <- DomainPoints(domain,"lalo")

  index <- point.closest.init(lon=lon,lat=lat,domain=domain,mask=mask)$index

  if (ldata) dd <- geo[index] else dd <- NA
  ij <- index
  ll <- cbind(lalodomain$lon[index],lalodomain$lat[index])

  list(data = dd, lonlat = ll, index = ij)
} 

### The old legacy code for LatLon minimisation
### Only on life support because colleagues yell at me.
.lalopoint0 <- function(data,lon,lat,minimise='lalo',mask=NULL){
### find the closest domain point to the given co-ordinates
### by minimising distance in either LonLat or projected co-ordinates
### This is in fact not exactly the same as minimising geographical distance!
  ldata <- inherits(data,"geofield")
  if (ldata) domain <- attr(data, "domain")
## remove Rfa dependency...
#  else if (inherits(data,"FAfile")) domain <- attr(data, "domain")
#  else if (inherits(data,"FAframe")) domain <- FAdomain(data)
  else if (inherits(data,"geodomain")) domain <- data
  else stop("data is not a geographical object. Can not interprete.")
  if(length(lon)!=1) stop("lalopoint0 only accepts a single point.")

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


##########################

.domain2lalo <- function(infield){
### return the LatLon co-ordinates of all grid points (X,Y), and the field values as Z
### all in 1 data frame (I think that's better than a list)
### You could use this for spatial interpolations...
  lalodomain <- DomainPoints(infield,"lalo")
  data.frame(x=as.vector(lalodomain$lon),y=as.vector(lalodomain$lat),
             z = as.vector(infield[1:attr(infield,'domain')$nx, 1:attr(infield,'domain')$ny]))
}

