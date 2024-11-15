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
# Function extracts location specific information on the timeseries of 
# wood stock production drawn from an already loaded gridded dataset.
# 
# Author: T. Luke Smallman (08/12/2021)
#
#########################################################################################

extract_wood_productivity<- function(i1,j1,timestep_days,spatial_type,resolution,
                                     grid_type,latlon_in,Cwood_inc_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("Cwood productivity / increment extracted for current location ",Sys.time(),sep=""))}

#   # find the nearest location
#   output = closest2d_2(1,Cwood_inc_all$lat,Cwood_inc_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Create time series output variables
   Cwood_inc = rep(-9999, length(timestep_days))
   Cwood_inc_unc = rep(-9999, length(timestep_days))
   Cwood_inc_lag = rep(-9999, length(timestep_days))

   # Loop through each time step of the Cwood increment / production timeseries,
   # its associated uncertainty and period of effect
   for (t in seq(1, length(Cwood_inc_all$place_obs_in_step))) {
        # Prodictivity estimate
        Cwood_inc[Cwood_inc_all$place_obs_in_step[t]] = Cwood_inc_all$Cwood_increment_gCm2[i1,j1,t]
        # Its uncertainty
        tmp = min(Cwood_inc[Cwood_inc_all$place_obs_in_step[t]], Cwood_inc_all$Cwood_increment_uncertainty_gCm2[i1,j1,t])
        Cwood_inc_unc[Cwood_inc_all$place_obs_in_step[t]] = tmp
        # Its period of effect
        Cwood_inc_lag[Cwood_inc_all$place_obs_in_step[t]] = Cwood_inc_all$Cwood_increment_lag_step[i1,j1,t]
   }

   # Set any time series values with NaN to missing data flag (-9999)
   Cwood_inc[which(is.na(Cwood_inc))] = -9999
   Cwood_inc_unc[which(is.na(Cwood_inc_unc))] = -9999
   Cwood_inc_lag[which(is.na(Cwood_inc_lag))] = -9999

   # pass the information back
   return(list(Cwood_inc = Cwood_inc, Cwood_inc_unc = Cwood_inc_unc, Cwood_inc_lag = Cwood_inc_lag))

} # end function extract_wood_productivity

## Use byte compile
extract_wood_productivity<-cmpfun(extract_wood_productivity)
