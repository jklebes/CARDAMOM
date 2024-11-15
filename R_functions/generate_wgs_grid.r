#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Function to create regular latitude / longitude grids for global WGS-84.
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

generate_wgs84_grid<-function(lat,long,resolution) {

    # Check whether we have two values for the resolution (longitude, latitude)
    if (length(resolution) == 1) {tmp = rep(resolution, length.out = 2)} else {tmp = resolution}

    # pt create raster of our desired resolution (in degrees) with the spatial extent set in the control file.
    # This function makes use of the WGS-84 coordinate system
    pt = rast(vals = 1, resolution = tmp, xmin = long[1], xmax = long[2], ymin = lat[1], ymax = lat[2], crs = "+init=epsg:4326")

    # extract the spatial information needed else where
    dims = dim(pt) ; lat_dim = dims[1] ; long_dim = dims[2]
    # extract the lat / long information needed
    long = crds(pt,df=TRUE, na.rm=FALSE)
    lat  = long$y ; long = long$x
    # reverse latitude vector to correct with R plotting orientation (i.e. make the world N/S rather than S/N)
    # It is very important to remember this has happened as the output from this fuction will not be consistent
    # with much of what is done with the calibration datasets.
    lat = lat[length(lat):1]

    # combine into output object
    output = list(cardamom_ext = pt, lat = lat, long = long, lat_dim = lat_dim, long_dim = long_dim)
    # tidy up
    rm(lat,long,lat_dim,long_dim,pt) ; gc(reset=TRUE,verbose=FALSE)

    # return back the solution
    return(output)

} # end function

## Use byte compile
generate_wgs84_grid<-cmpfun(generate_wgs84_grid)
