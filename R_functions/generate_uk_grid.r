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
# Function to create regular latitude / longitude grids for British National Grid
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

generate_uk_grid<-function(lat,long,resolution) {

    # Check whether we have two values for the resolution (longitude, latitude), (x,y)
    if (length(resolution) == 1) {tmp = rep(resolution, length.out = 2)} else {tmp = resolution}

    # Convert lat / long corners into OSGB36 format
    pt=data.frame(x=long,y=lat)
    # convert them into coordinates
    coordinates(pt)=~x+y
    # impose lat/long system (WGS-84) which the box is defined in
    proj4string(pt)=CRS("+init=epsg:4326")
    # transform to ORGB system so we can compare with the resolution
    # +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs
    pt=spTransform(pt,CRS("+init=epsg:27700"))

    # Create extent in GB grid
    pt = raster(vals = 1, resolution = tmp, ext = extent(pt), crs = "+init=epsg:27700")
    # Reproject onto the WGS84 grid
    pt = projectRaster(pt,crs = CRS("+init=epsg:4326"), method = "ngb")
    # extract the spatial information needed else where
    dims = dim(pt) ; lat_dim = dims[1] ; long_dim = dims[2]
    lat = coordinates(pt)[,2] ; long = coordinates(pt)[,1]
    # reverse latitude vector to correct with R plotting orientation (i.e. make the world N/S rather than S/N)
    lat = lat[length(lat):1]

    # combine into output object
    output = list(cardamom_exts = pt, lat = lat, long = long, lat_dim = lat_dim, long_dim = long_dim)
    # tidy up
    rm(lat,long,lat_dim,long_dim,pt) ; gc(reset=TRUE,verbose=FALSE)
    # And leave
    return(output)

} #end function

## Use byte compile
generate_uk_grid<-cmpfun(generate_uk_grid)
