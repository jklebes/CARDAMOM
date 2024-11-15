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
# Function extracts location specific information on the potential steady state value 
# of the wood stock pool from a loaded gridded dataset.
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_Cwood_potential<- function(i1,j1,timestep_days,spatial_type,resolution,
                                   grid_type,latlon_in,Cwood_potential_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("Cwood potential extracted for current location ",Sys.time(),sep=""))}

#   # find the nearest location
#   output = closest2d_2(1,Cwood_potential_all$lat,Cwood_potential_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract to local variable
   Cwood = Cwood_potential_all$biomass_gCm2[i1,j1]
   Cwood_unc = Cwood_potential_all$biomass_uncertainty_gCm2[i1,j1]

   # Convert any NaN to missing data flag -9999
   Cwood[which(is.na(Cwood))] = -9999 ; Cwood_unc[which(is.na(Cwood_unc))] = -9999

   # pass the information back
   return(list(Cwood_stock = Cwood, Cwood_stock_unc = Cwood_unc))

} # end function extract_Cwood_potential

## Use byte compile
extract_Cwood_potential<-cmpfun(extract_Cwood_potential)
