# people still yell at me for this:
lalopoint <- function(geo, lon, lat, minimise='proj', mask=NULL){
### find the closest domain point to the given co-ordinates
  if (minimise!='proj') stop("Only minimise='proj' is still supported.")
  ldata <- inherits(geo, "geofield")

  domain <- as.geodomain(geo)

  lalodomain <- DomainPoints(domain, "lalo")

  index <- as.matrix(point.closest.init(domain, lon=lon, lat=lat, mask=mask))

  if (ldata) dd <- geo[index] else dd <- NA
  ij <- index
  ll <- cbind(lalodomain$lon[index],lalodomain$lat[index])

  list(data = dd, lonlat = ll, index = ij)
} 


