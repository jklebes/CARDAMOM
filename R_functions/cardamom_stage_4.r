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
# Function to carry out stage 4 processes,
# i.e. creating generic figures of the analysis
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

cardamom_stage_4<-function(PROJECT) {

      print("Beginning stage 4: generating stardard outputs")

      # Generating site level plots or gridded
      if (PROJECT$spatial_type == "site" | grid_override) {
          # Generate figures of parameters and model values including
          # uncertainty information
          generate_uncertainty_figures(PROJECT)
      } else if (PROJECT$spatial_type == "grid") {
          # will generate spatial maps instead
          generate_parameter_maps(PROJECT)
          generate_simplified_stock_and_flux_maps(PROJECT)
      } else {
          stop('missing spatial_type definition (i.e. grid or site)')
      } # grid or site run

      # report to the user
      return(paste("CARDAMOM Report: 4 completed", sep=""))

} # end function cardamom_stage_4

## Use byte compile
cardamom_stage_4<-cmpfun(cardamom_stage_4)