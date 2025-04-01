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
# Function to carry out stage 5 processes,
# i.e. dump of all outputs to netcdf files
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

cardamom_stage_5<-function(PROJECT) {

   # Inform the user
   print("Stage 5 write a netcdf dump of the CARDAMOM output")

   if (PROJECT$spatial_type == "grid") {
       # Create the netcdf file
       create_grid_output_nc(PROJECT)
   } else if (PROJECT$spatial_type == "site") {
       create_states_all_nc(PROJECT)
   } else {
       print("PROJECT$spatial_type does not have valid value")
   }

   # report to the user
   return(paste("CARDAMOM Report: 5 completed", sep=""))

} # end function cardamom_stage_5

## Use byte compile
cardamom_stage_5<-cmpfun(cardamom_stage_5)