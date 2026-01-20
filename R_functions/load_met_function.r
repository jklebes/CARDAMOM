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
# Companion function to load_met_fields_for_extraction. Supports use of parallel 
# processing.
# 
# Author: T. Luke Smallman (02/05/2024)
# Modification History:
# 1) T. Luke Smallman (11/03/2025) 
#    ERA5 workflow updated to read whole year files and to only read required spatial domain rather than whole domain
#
#########################################################################################

load_met_function<- function (year_to_do,varid,infile_varid,spatial_type,cardamom_ext,
                              path_to_met_source,met_source,wheat,xy_bounds = xy_bounds) {

    if (met_source == "ERA" | met_source == "isimip3a") {

        # Determine the next file we expect to read in
        input_file_1 = paste(path_to_met_source,varid[1],"_",year_to_do,".nc",sep="")        
        # open netcdf files
        data1 = nc_open(input_file_1)
        # read the met drivers
        var_in = ncvar_get(data1, infile_varid[1],start=c(xy_bounds[1],xy_bounds[3],1), 
                                                  count=c(xy_bounds[2],xy_bounds[4],-1))
        # keep count of time steps
        t_grid = dim(var_in)[3]
        # read in the location information
        lat = ncvar_get(data1, "lat", start=c(xy_bounds[1],xy_bounds[3]), count=c(xy_bounds[2],xy_bounds[4])) 
        long = ncvar_get(data1, "lon", start=c(xy_bounds[1],xy_bounds[3]), count=c(xy_bounds[2],xy_bounds[4]))

        # close files after use
        nc_close(data1)

        # Initialise the output variable
        var1_out = rep(NA, length(wheat)*dim(var_in)[3])

        # Loop each time step and aggregte before placing into to an out variables
        for (t in seq(1, dim(var_in)[3])) {

             # Inform the user
             #if (use_parallel == FALSE) {print(paste("...processing time step ",t," of ",dim(var_in)[3]," for ",year_to_do,sep=""))}

             # Convert to a raster, assuming standad WGS84 grid
             # This dependes on the lat / long / tmp1 spatially matching each other AND
             # latitude ranging -90/90 and longitude ranging -180/180 degrees
             var1 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(var_in[,,t]))
             var1 = rast(var1, crs = ("epsg:4326"), type="xyz")

             # Extract the epsg from the file
             epsg = crs(var1, describe = TRUE)$code
             # If we have an epsg then we want to know if it differs from the one desired by the analysis
             if (epsg != gsub("epsg:","",cardamom_grid_type)) {
                 # Ensure that the extent of the input object is consistent 
                 # with the possible extent of the selected epsg
                 var1 = crop(var1, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))
                 # If it does not match we need to reproject it
                 var1 = project(var1, cardamom_grid_type, method="near", align = FALSE) ; gc()
             }

             # Extend the extent of the overall grid to the analysis domain
             var1 = extend(var1,cardamom_ext)
             # Trim the extent of the overall grid to the analysis domain
             var1 = crop(var1,cardamom_ext)
             # Set any missing value flags to NA
             var1[which(as.vector(var1) == -9999)] = NA
             # Match resolutions
             if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                 # Resample to correct grid
                 var1 = resample(var1, cardamom_ext, method="average") ; gc() # "average"
             } # Aggrgeate to resolution

             # determine new (s)tart and (e)nd points for the output array
             s = (length(wheat)*(t-1)) + 1 ; e = length(wheat)*t
             # Load the current output into the output variable.
             # Don't use append as this requires copying of the whole array each time,
             # takes time.
             #var1_out = append(var1_out, values(var1)[wheat])
             var1_out[s:e] = values(var1)[wheat]

        } # t loop
        
        # clean up
        rm(var1,epsg,lat,long,var_in) ; gc(reset=TRUE,verbose=FALSE)

    } else if (met_source == "trendy") {

        # Check whether this is a leap year or not
        nos_days = nos_days_in_year(year_to_do)
        # Determine what the number of days per month are
        # Define days in month
        days_per_month = rep(31,12) ; days_per_month[c(9,4,6,11)] = 30
        if (nos_days == 365) {days_per_month[2] = 28} else {days_per_month[2] = 29}

        # open first file in the sequence
        input_file_1 = paste(path_to_met_source,"/",varid[1],"_",year_to_do,"_monthly.nc",sep="")
        # open netcdf files
        data1 = nc_open(input_file_1)
        # read the met drivers
        var_in = ncvar_get(data1, infile_varid[1])
        # read in the location information
        lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")

        # expand the one directional values here into 2 directional
        lat_dim = length(lat) ; long_dim = length(long)
        long = array(long,dim=c(long_dim,lat_dim))
        lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)
        # close files after use
        nc_close(data1)

        # Loop each time step and aggregte before placing into to an out variables
        var1_out = 0
        for (t in seq(1, dim(var_in)[3])) {

             # Inform the user
             #if (use_parallel == FALSE) {print(paste("...processing time step ",t," of ",dim(var_in)[3]," for ",year_to_do,sep=""))}

             # Convert to a raster, assuming standad WGS84 grid
             # This dependes on the lat / long / tmp1 spatially matching each other AND
             # latitude ranging -90/90 and longitude ranging -180/180 degrees
             var1 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(var_in[,,t]))
             var1 = rast(var1, crs = ("epsg:4326"), type="xyz")

             # Extract the epsg from the file
             epsg = crs(var1, describe = TRUE)$code
             # If we have an epsg then we want to know if it differs from the one desired by the analysis
             if (epsg != gsub("epsg:","",cardamom_grid_type)) {
                 # Ensure that the extent of the input object is consistent 
                 # with the possible extent of the selected epsg
                 var1 = crop(var1, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))
                 # If it does not match we need to reproject it
                 var1 = project(var1, cardamom_grid_type, method="near", align = FALSE) ; gc()
             }

             # Extend the extent of the overall grid to the analysis domain
             var1 = extend(var1,cardamom_ext)
             # Trim the extent of the overall grid to the analysis domain
             var1 = crop(var1,cardamom_ext)
             var1[which(as.vector(var1) == -9999)] = NA
             # Match resolutions of the datasets
             if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                 # Resample to correct grid
                 var1 = resample(var1, cardamom_ext, method="average") ; gc() 
             } # Aggrgeate to resolution

             # break out from the rasters so can manipulate
             var1_out = append(var1_out, rep(as.vector(unlist(var1))[wheat], times = days_per_month[t]))

        } # end of time loop

        # Remove initial value
        var1_out = var1_out[-1]

        # keep count of time steps
        t_grid = sum(days_per_month)

        # clean up
        rm(var_in,var1,t) ; gc(reset=TRUE,verbose=FALSE)

    } # end data source selection

    # return back to the user
    return(list(var_out=var1_out,t_grid=t_grid))

} # end function
## Use byte compile
load_met_function<-cmpfun(load_met_function)
