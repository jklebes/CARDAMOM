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

    # Spatial grid
    output = generate_grid(grid_type,lat,long,resolution)
    # extract the latitude / longitude and extent/resolution information
    lat = output$lat ; long = output$long ; long_dim = output$long_dim ; lat_dim = output$lat_dim
    cardamom_ext = output$cardamom_ext
    rm(output)

    # Create a grid specifically to be used for extracting the correct location from the gridded datasets
    # i.e. these will be upside down from the human eye.
    # This reverses what was done in generate_grid()
    obs_long_grid = long[,lat_dim:1]
    obs_lat_grid = lat[,lat_dim:1] 

    # remove the values we don't want
    if (length(remove) > 0) {
        lat = lat[-remove] ; long = long[-remove]
    } else {
        lat = as.vector(lat) ; long = as.vector(long)
    }

    # output the result
    output = list(lat = lat, long = long, cardamom_ext = cardamom_ext,
                  obs_long_grid = obs_long_grid, obs_lat_grid = obs_lat_grid)
    rm(lat,long) ; gc(reset=TRUE, verbose=FALSE)
    return(output)

}
## Use byte compile
determine_lat_long_needed<-cmpfun(determine_lat_long_needed)
