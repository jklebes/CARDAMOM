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
# Function extracts a location specific estimate of the initial Csom value 
# from the gridded dataset.
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_Csom_prior<-function(i1,j1,spatial_type,resolution,grid_type,
                             latlon_wanted,Csom_all) {

  # extract information on soil C content from soil grids database

  # Update the user
  if (use_parallel == FALSE) {print(paste("Csom prior data extracted for current location ",Sys.time(),sep=""))}

#  # find desired lat / long location within the soilgrid database
#  output = closest2d_2(1,Csom_all$lat,Csom_all$long,latlon_wanted[1],latlon_wanted[2])
#  i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

  # Extract target location
  Csom = Csom_all$Csom[i1,j1]
  Csom_unc = Csom_all$Csom_unc[i1,j1]

  # Convert any NaN to missing data flag -9999
  Csom[is.na(Csom) == TRUE] = -9999 ; Csom_unc[is.na(Csom_unc) == TRUE] = -9999

  # retun back to the user
  return(list(Csom_initial = Csom, Csom_initial_unc = Csom_unc))

} # end function extract_Csom_prior

## Use byte compile
extract_Csom_prior<-cmpfun(extract_Csom_prior)
