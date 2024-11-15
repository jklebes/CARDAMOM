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
# Function determines all lat long coordinates needed for grid run mode
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

determine_lat_long_needed<- function(lat,long,resolution,grid_type,remove) {

    # check input data
    if (length(which(long > 180)) > 0) {stop("Long should be -180 to +180")}

    # generate UK or WGS-84 lat long grid
    if (grid_type == "UK") {
        output = generate_uk_grid(lat,long,resolution)
    } else if (grid_type=="wgs84") {
        output = generate_wgs84_grid(lat,long,resolution)
    } else {
        stop('have selected invalid grid type, the valid options are "UK" and "wgs84"')
    }
    # extract the latitude / longitude and extent/resolution information
    lat = output$lat ; long = output$long ; long_dim = output$long_dim ; lat_dim = output$lat_dim
    cardamom_ext = output$cardamom_ext

    # Create a grid specifically to be used for extracting the correct location from the gridded datasets
    obs_long_grid = array(long, dim=c(long_dim,lat_dim))
    obs_lat_grid = array(rev(lat), dim=c(long_dim,lat_dim)) # rev() accounts for the flipping of orientation conducted in the generate_*_grid()

    # remove the values we don't want
    if (length(remove) > 0) {lat = lat[-remove] ; long = long[-remove]}

    # output the result
    output = list(lat = lat, long = long, cardamom_ext = cardamom_ext,
                  obs_long_grid = obs_long_grid, obs_lat_grid = obs_lat_grid)
    rm(lat,long) ; gc(reset=TRUE, verbose=FALSE)
    return(output)

}
## Use byte compile
determine_lat_long_needed<-cmpfun(determine_lat_long_needed)
