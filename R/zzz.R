.onLoad <- function(lib,pkg){
  if(is.element("RColorBrewer",.packages(all.available=TRUE))) {
    require(RColorBrewer)
  }
  else print("Warning: You don't have RColorBrewer installed! That package (available from CRAN) gives nicer colour palettes!")
}
