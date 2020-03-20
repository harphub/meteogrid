# people still yell at me for this:
lalopoint <- function(geo,lon,lat,minimise='proj',mask=NULL){
### find the closest domain point to the given co-ordinates
  if (minimise!='proj') stop("Only minimise='proj' is still supported.")
  ldata <- inherits(geo,"geofield")

  if (ldata ) domain <- attr(geo, "domain")
  else if (inherits(geo,"geodomain")) domain <- geo
  else if ("domain" %in% names(attributes(geo))) domain <- attr(geo,"domain")
  else stop("geo is not a geographical object. Can not interprete.")

  lalodomain <- DomainPoints(domain,"lalo")

  index <- point.closest.init(geo,lon=lon,lat=lat,mask=mask)$index

  if (ldata) dd <- geo[index] else dd <- NA
  ij <- index
  ll <- cbind(lalodomain$lon[index],lalodomain$lat[index])

  list(data = dd, lonlat = ll, index = ij)
} 


