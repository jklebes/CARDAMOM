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
# Function to load estimates of initial soil moisture from gridded dataset
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

load_soilwater_fields_for_extraction<-function(latlon_in,soilwater_source) {

  if (soilwater_source == "GLEAM") {

    # let the user know this might take some time
    print("Loading processed GLEAM soil water fraction for subsequent sub-setting ...")

    # open processed modis files
    input_file_1=paste(path_to_gleam,"/GLEAM_soil_moisture_prior_0.25.nc",sep="")
    data1=nc_open(input_file_1)

    # extract location variables
    lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
    # read the file
    soil_water=ncvar_get(data1, "GLEAM_soil_moisture_prior")
    soil_water_unc=ncvar_get(data1, "GLEAM_soil_moisture_prior_SD")

    # close files after use
    nc_close(data1)

    # output variables
    return(list(soil_water=soil_water,soil_water_unc = soil_water_unc,lat=lat,long=long))

  } else {
    # output variables
    return(list(soil_water=-9999,soil_water_unc = -9999,lat=-9999,long=-9999))
  }

} # function end load_soilwater_fields_for_extraction

## Use byte compile
load_soilwater_fields_for_extraction<-cmpfun(load_soilwater_fields_for_extraction)
