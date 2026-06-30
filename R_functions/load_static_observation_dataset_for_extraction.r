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
# Generic function to load observational datasets into CARDAMOM
# 
# Author: T. Luke Smallman (16/12/2024)
#
#########################################################################################

load_static_observation_dataset_for_extraction<-function(latlon_in,cardamom_ext,grid_type,
                                                         data_source,data_path,prefix,
                                                         est_var_name_in,unc_var_name_in,
                                                         est_var_name_out,unc_var_name_out) {

    # Select data_source option
    if (data_source == "Gridded_nc") {

        # let the user know this might take some time
        print(paste("Loading ",est_var_name_in," for subsequent sub-setting ...",sep=""))

        # check which file prefix we are using today
        # list all available files which we will then search
        input_files = list.files(data_path,full.names=TRUE,pattern="\\.nc$")
        #prefix = "MCD15A2H_LAI_(.)*" # (.)* wildcard characters for unix standard MCD15A2H_LAI_*
        #prefix = "net_biome_exchange_"

        ## Begin reading the files in

        # Check the expected file pattern is found in the available files
        this_year = input_files[grepl(paste("/",prefix,sep=""), input_files)]
        if (length(this_year) > 0) {

            # open the file
            data1 = nc_open(this_year)

            # Extract spatial information - latitude
            if (length(which(names(data1$var) == "lat")) > 0 | length(which(names(data1$dim) == "lat"))) {
                lat_in = ncvar_get(data1, "lat") 
            } else if (length(which(names(data1$var) == "latitude")) > 0 | length(which(names(data1$dim) == "latitude"))) {
                lat_in = ncvar_get(data1, "latitude") 
            } else {
                stop("......no variable or dimension called lat or latitude could be found")
            } # finding lat
            # Extract spatial information - longitude
            if (length(which(names(data1$var) == "lon")) > 0 | length(which(names(data1$dim) == "lon"))) {
                long_in = ncvar_get(data1, "lon") 
            } else if (length(which(names(data1$var) == "longitude")) > 0 | length(which(names(data1$dim) == "longitude"))) {
                long_in = ncvar_get(data1, "longitude") 
            } else {
                stop("......no variable or dimension called lat or latitude could be found")
            } # finding lat
            # Now check whether this is a 2D array or not
            if (length(dim(lat_in)) == 2 && length(dim(long_in)) == 2) {
                 # Do nothing, as this should be exactly what we want
             } else {
                 # We will assume that the arrays match the vectorisation order found in R.
                 # A warning will be issued, placing the onus on the user to make sure this is right
                 print("......The lat or lon information are provided as a vector and not as a 2D array, as specified in the dataset description documents")
                 print("......The code will construct the 2D array assuming the vectorisation order used in R and that the vectors represent the x and y coordinates.")
                 tmp1 = length(long_in) ; tmp2 = length(lat_in)
                 lat_in = t(array(lat_in, dim=c(tmp2,tmp1)))
                 long_in = array(long_in, dim=c(tmp1,tmp2))
                 rm(tmp1,tmp2) 
            } 
            # Extract the current global attributes
            global_attributes = ncatt_get(data1,0)
            # Check whether there is any information regarding the EPSG
            global_attributes = unlist(global_attributes)
            # Default assumption for netcdf files is for the epsg: 4326 i.e. the WGS-84 lat/long grid
            # But here we will search for any specific information
            epsg = 4326 ; aa = 1
            if (length(global_attributes) >= 1) {
                while (aa > 0) {
                   # Check whether the epsg is provided somewhere
                   if (grepl("epsg", global_attributes[aa], ignore.case = FALSE)) {
                       # We need to extract this information
                       epsg = gsub("[^0-9.-]", "", unlist(strsplit(global_attributes[aa],"epsg"))[2])
                       aa = -9999 # escape flag
                   } else if (grepl("EPSG", global_attributes[aa], ignore.case = FALSE)) {
                       # We need to extract this information
                       epsg = gsub("[^0-9.-]", "", unlist(strsplit(global_attributes[aa],"EPSG"))[2])
                       aa = -9999 # escape flag
                   } else { 
                       # Break the loop condition or increment
                       if (aa == length(global_attributes)) {aa = -9999} else { aa = aa + 1 }
                   }
                }
            } # assuming we have any global attributes

            # read the observation estimate 
            est_in = ncvar_get(data1, est_var_name_in) # Variable estimate
            # read error variable, if present
            if (length(which(grepl(unc_var_name_in,names(data1$var)) == TRUE)) > 0) {
                std_in = ncvar_get(data1, unc_var_name_in) # Variable standard deviation 
                std_present = TRUE
            } else {
                std_in = -9999 ; std_present = FALSE
            }                 
            # Close the current file
            nc_close(data1) 

            # Convert to a raster, assuming standard WGS84 grid
            est_out_tif = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(est_in))
            est_out_tif = rast(est_out_tif, crs = paste("epsg:",epsg,sep=""), type="xyz")
            rm(est_in)
            # Check that we have a uncertainty estimate
            if (std_present) {
                std_out_tif = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(std_in))
                std_out_tif = rast(std_out_tif, crs = paste("epsg:",epsg,sep=""), type="xyz")
                rm(std_in)
            } 

            # Extract the epsg from the file
            epsg = crs(est_out_tif, describe = TRUE)$code
            if (is.null(epsg) | epsg == "") { stop(paste("the use_lcm specification leads to a geotif which does not contain epsg information."))}
            # If we have an epsg then we want to know if it differs from the one desired by the analysis
            if (epsg != gsub("epsg:","",grid_type)) {
                # Ensure that the extent of the input object is consistent 
                # with the possible extent of the selected epsg
                est_out_tif = crop(est_out_tif, ext(unlist(crs(grid_type, describe=TRUE)$extent)))
                # If it does not match we need to reproject it
                est_out_tif = project(est_out_tif, grid_type, method="near", align = FALSE) ; gc()
                if (std_present) { 
                    # Ensure that the extent of the input object is consistent 
                    # with the possible extent of the selected epsg
                    std_out_tif = crop(std_out_tif, ext(unlist(crs(grid_type, describe=TRUE)$extent)))                                                  
                    std_out_tif = project(std_out_tif, grid_type, method="near", align = FALSE) ; gc() 
                }
            }
 
            # Use extend and crop to ensure we closely match the extent of the analysis
            est_out_tif = extend(est_out_tif,cardamom_ext) ; est_out_tif = crop(est_out_tif,cardamom_ext) 
            if (std_present) { std_out_tif = extend(std_out_tif,cardamom_ext) ; std_out_tif = crop(std_out_tif,cardamom_ext) }

            # Adjust spatial resolution of the datasets, this occurs in all cases
            if (res(est_out_tif)[1] != res(cardamom_ext)[1] | res(est_out_tif)[2] != res(cardamom_ext)[2]) {
                # Resample to correct grid.
                # Probably should be done via aggregate function to allow for correct error propogation
                est_out_tif = resample(est_out_tif, cardamom_ext, method="average") ; gc() 
                if (std_present) { std_out_tif = resample(std_out_tif, cardamom_ext, method="average") ; gc() }
            } # Aggrgeate to resolution

        } # is there information for the current year?

        # extract dimension information for the grid, 
        # note the axis switching between raster and actual array
        xdim = dim(est_out_tif)[2] ; ydim = dim(est_out_tif)[1]
        # extract the lat / long information needed
        long = crds(est_out_tif, df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))

        # Extract from the raster structure into arrays
        est_out = array(values(est_out_tif), dim=c(xdim,ydim))
        if (std_present) { std_out = array(values(std_out_tif), dim=c(xdim,ydim)) }

        # If the standard deviation exists, then we should ensure that 
        # both the estimate and uncertainty have common occurance of NaN
        if (std_present) {
            # Check which are NA in either the estimate or the uncertainty
            filter = which(is.na(est_out) | is.na(std_out))
            # and remove those which differ
            est_out[filter] = NA ; std_out[filter] = NA
        }
 
        # Set dummy value for the output uncertainty if required
        if (std_present == FALSE) {std_out = -9999} 

        # Create output object
        output_all = list(est_out, std_out, lat = lat, long = long) 
        # Update with the correct variable names
        names(output_all)[1:2]<-c(est_var_name_out,unc_var_name_out)

        # clean up variables
        rm(doy_in,est_out,std_out,lat,long) ; gc(reset=TRUE,verbose=FALSE)

        # Return to function
        return(output_all)

    } else if (data_source == "Gridded_tif") { 

        # let the user know this might take some time
        print(paste("Loading ",est_var_name_in," for subsequent sub-setting ...",sep=""))

        # list all available files which we will then search
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_files = list.files(data_path, full.names=TRUE, recursive = TRUE, pattern="\\.tif$")
        #prefix = "MCD15A2H_LAI_(.)*" # (.)* wildcard characters for unix standard MCD15A2H_LAI_*
        #prefix = "net_biome_exchange_"
   
        # Search files with the correct prefix only  
        input_files = input_files[grepl(paste("/",prefix,sep=""), input_files)]

        # Pull out those specifically for the variable of interest
        # based on them not being those with uncertainty
        est_input_files = input_files[grepl("uncertainty", input_files) == FALSE]
        # As there should only be a single static file, stop if not the case
        if (length(est_input_files) != 1) {stop(paste("Only a single file was expected for ",prefix," but a different number was found", sep=""))}
        # Read in the data file
        est_out = rast(est_input_files)
        # As we are using geotifs, we need to be sure we have matching estimate
        # and uncertainty files. Assuming we have asked for an uncertainty variable.
        if (unc_var_name_out != "") {
            unc_input_files = input_files[grepl("uncertainty", input_files)]
            if (length(unc_input_files) != 1) {stop(paste("Only a single file was expected for ",prefix," but a different number was found", sep=""))}
            if (length(est_input_files) != length(unc_input_files)) {
                stop(paste("The number of estimate files for ",prefix,
                           " and its uncertainty ",prefix,"_uncertainty do not match",sep=""))
            }
            # Read in the uncertainty file
            std_out = rast(unc_input_files)            
            std_present = TRUE
        } else {
            std_present = FALSE
        }

        # Extract the epsg from the file
        epsg = crs(est_out, describe = TRUE)$code
        if (is.null(epsg) | epsg == "") { stop(paste("the use_lcm specification leads to a geotif which does not contain epsg information."))}
        # If we have an epsg then we want to know if it differs from the one desired by the analysis
        if (epsg != gsub("epsg:","",grid_type)) {
            # Ensure that the extent of the input object is consistent 
            # with the possible extent of the selected epsg
            est_out = crop(est_out, ext(unlist(crs(grid_type, describe=TRUE)$extent)))
            # If it does not match we need to reproject it
            est_out = project(est_out, grid_type, method="near", align = FALSE) ; gc()
            if (std_present) { 
                # Ensure that the extent of the input object is consistent 
                # with the possible extent of the selected epsg
                std_out = crop(std_out, ext(unlist(crs(grid_type, describe=TRUE)$extent)))                                  
                std_out = project(std_out, grid_type, method="near", align = FALSE) ; gc() 
            }
        }

        # Use extend and crop to ensure we closely match the extent of the analysis
        est_out = extend(est_out,cardamom_ext) ; est_out = crop(est_out,cardamom_ext) 
        if (std_present) { std_out = extend(std_out,cardamom_ext) ; std_out = crop(std_out,cardamom_ext) }

        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(est_out)[1] != res(cardamom_ext)[1] | res(est_out)[2] != res(cardamom_ext)[2]) {
            # Resample to correct grid.
            # Probably should be done via aggregate function to allow for correct error propogation
            est_out = resample(est_out, cardamom_ext, method="average") ; gc() 
            if (std_present) { std_out = resample(std_out, cardamom_ext, method="average") ; gc() }
        } # Aggrgeate to resolution             

        # extract dimension information for the grid, 
        # note the axis switching between raster and actual array
        xdim = dim(est_out)[2] ; ydim = dim(est_out)[1]
        # extract the lat / long information needed
        long = crds(est_out, df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))

        # Extract from the raster structure into arrays
        est_out = array(values(est_out), dim=c(xdim,ydim))
        if (std_present) { std_out = array(values(std_out), dim=c(xdim,ydim)) }

        # If the standard deviation exists, then we should ensure that 
        # both the estimate and uncertainty have common occurance of NaN
        if (std_present) {
            # Check which are NA in either the estimate or the uncertainty
            filter = which(is.na(est_out) | is.na(std_out))
            # and remove those which differ
            est_out[filter] = NA ; std_out[filter] = NA
        }
 
        # Set dummy value for the output uncertainty if required
        if (std_present == FALSE) {std_out = -9999} 

        # Create output object
        output_all = list(est_out, std_out, lat = lat, long = long) 
        # Update with the correct variable names
        names(output_all)[1:2]<-c(est_var_name_out,unc_var_name_out)

        # clean up variables
        rm(est_out,std_out,lat,long) ; gc(reset=TRUE,verbose=FALSE)
        return(output_all)

    } else if (data_source == " " | data_source == "site_specific") {
    
        # Do nothing as this should be not needed or will be read later
    
    } else {

        # We have a problem as something unexpected has been specified
        stop(paste("The data_source = ",data_source," is not a recognised value (Gridded_nc, Gridded_tif, site_specific or blank)",sep=""))

    } # Gridded dataset or not

} # function end load_static_observation_dataset_for_extraction

## Use byte compile
load_static_observation_dataset_for_extraction<-cmpfun(load_static_observation_dataset_for_extraction)
