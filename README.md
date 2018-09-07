# meteogrid
##R package for geographical grids (mainly meteorology)
meteogrid`(which used to be called geogrid until 2018) is a toolbox for gridded data, mainly focussed on meteorological (model) data. It offers tools for visualisation, interpolation etc. It does not contain routines for decoding of binary data formats. Therefore, meteogrid will typically be used together with `Rgrib2` (package for GRIB decoding, available under GPLv3) or possibly `Rfa` (internal to the ALADIN & HIRLAM consortia).


##Installation
The installation requires that you first install the PROJ4 library for map projections. This may be obtained from http://proj4.org

##License
Copyright 2003-2016, Alex Deckmyn, Royal Meteorological Institute of Belgium
alex.deckmyn@meteo.be

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
inst/COPYING


