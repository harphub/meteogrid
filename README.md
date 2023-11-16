# meteogrid

_A toolbox for working with and viewing gridded (meteorological) data_

meteogrid (which used to be called geogrid until 2018) is a toolbox for gridded 
data, mainly focussed on meteorological (model) data. It offers tools for 
visualisation, interpolation etc. It does not contain routines for decoding of 
binary data formats. Therefore, meteogrid will typically be used together with 
[Rgrib2](https://github.com/harphub/Rgrib2) (package for GRIB decoding) or 
possibly [Rfa](https://github.com/harphub/Rfa) (internal to the ACCORD 
consortium).


## Installation

The installation requires that you first install the PROJ library for map 
projections. This may be obtained from http://proj.org

```r
# install.packages("remotes")
remotes::install_github("harphub/meteogrid")
```



