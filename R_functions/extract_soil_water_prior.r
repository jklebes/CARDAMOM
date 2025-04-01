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
# Function to extract a prior estimate of soil moisture content from gridded dataset
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_soilwater_initial<- function(spatial_type,resolution,grid_type,latlon_in,soilwater_all) {

  # Update the user
  if (use_parallel == FALSE) {print(paste("Initial soil water fraction extracted for current location ",Sys.time(),sep=""))}

  # convert input data long to conform to what we need
  check1=which(soilwater_all$long > 180) ; if (length(check1) > 0) { soilwater_all$long[check1] = soilwater_all$long[check1]-360 }

  # find the nearest location
  output=closest2d_2(1,soilwater_all$lat,soilwater_all$long,latlon_in[1],latlon_in[2])
  j1=unlist(output)[2];i1=unlist(output)[1]

  # return long to 0-360
  if (length(check1) > 0) { soilwater_all$long[check1] = soilwater_all$long[check1]+360 }
  # If resolution has been provides as single value then adjust this here
  if (length(resolution) == 1 & spatial_type == "grid") {tmp_res = resolution * c(1,1)} else {tmp_res = resolution}

  # work out number of pixels to average over
  if (spatial_type == "grid") {
      # resolution of the product
      product_res = c(abs(soilwater_all$long[2,1]-soilwater_all$long[1,1]),abs(soilwater_all$lat[1,2]-soilwater_all$lat[1,1]))
      if (grid_type == "wgs84") {
          # radius is ceiling of the ratio of the product vs analysis ratio
          radius = round(resolution / product_res, digits=0)
      } else if (grid_type == "UK") {
          # Estimate radius for UK grid assuming radius is determine by the longitude size
          # 6371e3 = mean earth radius (m)
          radius = round(rad2deg(sqrt((resolution / 6371e3**2))) / product_res, digits=0)
          #radius = max(0,floor(1*resolution*1e-3*0.5))
      } else {
          stop("have not specified the grid used in this analysis")
      }
  } else {
      radius = c(0,0)
  }

  answer=NA
  while (is.na(answer) == TRUE) {
     # work out average areas
     average_i = (i1-radius[1]):(i1+radius[1]) ; average_j = (j1-radius[2]):(j1+radius[2])
     average_i = max(1,(i1-radius[1])):min(dim(soilwater_all$soil_water)[1],(i1+radius[1]))
     average_j = max(1,(j1-radius[2])):min(dim(soilwater_all$soil_water)[2],(j1+radius[2]))
     # carry out averaging
     tmp=soilwater_all$soil_water[average_i,average_j] ; tmp[which(tmp == -9999)] = NA
     soilwater = mean(tmp, na.rm=TRUE)
     tmp = soilwater_all$soil_water_unc[average_i,average_j] ; tmp[which(tmp == -9999)] = NA
     soilwater_unc = mean(tmp, na.rm=TRUE)
     # error checking
     if (is.na(soilwater) | soilwater == 0) {radius = radius+1 ; answer = NA} else {answer = 0}
  }
  if (use_parallel == FALSE) {print(paste("NOTE: Initial soil water averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))}

  # restrict minimum uncertainty allowed to prevent errors
  # uncertainty provided is standard error averaged across 15 years of annual averages
  soilwater_unc = soilwater_unc * 3.872983 # sqrt(15) = 3.872983
  soilwater_unc = max(0.05,soilwater_unc)
  # pass the information back
  return(list(soil_water = soilwater,soil_water_unc = soilwater_unc))

} # end function extract_soilwater_initial

## Use byte compile
extract_soilwater_initial<-cmpfun(extract_soilwater_initial)
