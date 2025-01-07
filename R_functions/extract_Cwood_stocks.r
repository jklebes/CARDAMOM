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
# Function extracts location specific information on the timeseries of wood stock  
# information drawn from an already loaded gridded dataset.
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_Cwood_stocks<- function(i1,j1,timestep_days,spatial_type,resolution,grid_type,
                                latlon_in,Cwood_stock_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("Cwood stocks extracted for current location ",Sys.time(),sep=""))}

   # Create time series output variables
   Cwood_stock = rep(-9999, length(timestep_days))
   Cwood_stock_unc = rep(-9999, length(timestep_days))

   # Loop through each time step of the Cwood time series obs and
   # estimate average value
   for (t in seq(1, length(Cwood_stock_all$doy_obs))) {
        Cwood_stock[Cwood_stock_all$doy_obs[t]] = Cwood_stock_all$biomass_gCm2[i1,j1,t]
        tmp = min(Cwood_stock[Cwood_stock_all$doy_obs[t]], Cwood_stock_all$biomass_uncertainty_gCm2[i1,j1,t])
        Cwood_stock_unc[Cwood_stock_all$doy_obs[t]] = tmp
   }

   # Set any time series values with NaN to missing data flag (-9999)
   Cwood_stock[which(is.na(Cwood_stock))] = -9999
   Cwood_stock_unc[which(is.na(Cwood_stock_unc))] = -9999

   # pass the information back
   return(list(Cwood_stock = Cwood_stock, Cwood_stock_unc = Cwood_stock_unc))
   #return(list(Cwood_stock = Cwood_stock, Cwood_stock_unc = 250))

} # end function extract_Cwood_stocks

## Use byte compile
extract_Cwood_stocks<-cmpfun(extract_Cwood_stocks)
