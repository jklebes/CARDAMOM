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
# Function used to extract location specific information on leaf carbon per leaf area (gC/m2) 
# from already loaded gridded datasets
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_lca_prior<- function(i1,j1,spatial_type,resolution,grid_type,latlon_in,lca_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("LCA prior extracted for current location ",Sys.time(),sep=""))}

#   # find the nearest location
#   output = closest2d_2(1,lca_all$lat,lca_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract target location
   lca_gCm2 = lca_all$lca_gCm2[i1,j1]
   lca_unc_gCm2 = lca_all$lca_uncertainty_gCm2[i1,j1]

   # Convert any NaN to missing data flag -9999
   lca_gCm2[is.na(lca_gCm2)] = -9999 ; lca_unc_gCm2[is.na(lca_unc_gCm2)] = -9999

   # pass the information back
   return(list(lca_gCm2 = lca_gCm2, lca_unc_gCm2 = lca_unc_gCm2))

} # end function extract_lca_prior

## Use byte compile
extract_lca_prior<-cmpfun(extract_lca_prior)
