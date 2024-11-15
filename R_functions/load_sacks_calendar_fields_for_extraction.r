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
# Function to load prior estimates of arable crop sowing and harvest dates from gridded dataset
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

load_sacks_calendar_fields_for_extraction<-function(latlon_in,crop_management_source) {

  if (crop_management_source == "sacks_crop_calendar") {

    # let the user know this might take some time
    print("Loading processed 5 min Sacks crop calendar fields for subsequent sub-setting ...")

    # open processed modis files
    input_file_1=paste(path_to_crop_management,"/Wheat.Winter.crop.calendar.fill.nc",sep="")
    data1=nc_open(input_file_1)

    # extract location variables
    lat=ncvar_get(data1, "latitude") ; long=ncvar_get(data1, "longitude")
    # read the modis lai drivers
    plant=ncvar_get(data1, "plant")
    #	plant_range=ncvar_get(data1, "plant.range")
    harvest=ncvar_get(data1, "harvest")
    #	harvest_range=ncvar_get(data1, "harvest.range")

    if (length(dim(latlon_in)) > 1) {
      max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
      min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
    } else {
      max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
      min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
    }
    # determine which locations to keep form inputs
    keep_lat=which(lat > min_lat & lat < max_lat)
    keep_long=which(long > min_long & long < max_long)
    # filter the inputs
    plant=plant[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    plant_range=-9999 #plant_range[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    harvest=harvest[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    harvest_range=-9999 #harvest_range[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    # now filter the inputs spatial data
    lat=lat[keep_lat] ; long=long[keep_long]
    # close files after use
    nc_close(data1)

    # clean
    rm(keep_lat,keep_long,max_lat,max_long,min_lat,min_long) ; gc(reset=TRUE,verbose=FALSE)

    # output variables
    return(list(plant=plant,plant_range=plant_range,harvest=harvest,harvest_range=harvest_range,lat=lat,long=long))

  } else {
    # output variables
    return(list(plant=-9999,plant_range=-9999,harvest=-9999,harvest_range=-9999,lat=-9999,long=-9999))
  }

} # function end load_sacks_calendar_fields_for_extraction

## Use byte compile
load_sacks_calendar_fields_for_extraction<-cmpfun(load_sacks_calendar_fields_for_extraction)
