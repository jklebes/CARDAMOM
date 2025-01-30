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

load_observation_dataset_for_extraction<-function(latlon_in,cardamom_ext,grid_type,
                                                  data_source,data_path,prefix,
                                                  years_to_load,
                                                  est_var_name_in,unc_var_name_in,
                                                  lag_var_name_in,est_var_name_out,
                                                  unc_var_name_out,lag_var_name_out) {

    # Select data_source option
    if (data_source == "Gridded_nc") {

        # let the user know this might take some time
        print(paste("Loading ",est_var_name_in," for subsequent sub-setting ...",sep=""))

        # check which file prefix we are using today
        # list all available files which we will then search
        input_files = list.files(data_path,full.names=TRUE,pattern="\\.nc$")
        input_files = input_files[grepl(paste("/",prefix,sep=""),input_files)]
        #prefix = "MCD15A2H_LAI_(.)*" # (.)* wildcard characters for unix standard MCD15A2H_LAI_*
        #prefix = "net_biome_exchange_"

        # timing information on the number of day in a month
        month_days = rep(31,length.out=12)
        month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30 

        # Set default starting positions
        lat_done = FALSE ; missing_years = 0 ; years_with_obs = 0 ; keepers = 0 ; yrs = 1 ; doy_out = 0

        # Loop through all the analyses years and check whether files exist for it
        for (yr in seq(1, length(years_to_load))) {
             # create the prefix for the observation file files we will want for a given year
             # NOTE this assumes that the prefix is directly followed by the year of the file
             input_file_1 = paste(prefix,years_to_load[yr],sep="")
             # then check whether this pattern is found in the available files
             this_year = grepl(input_file_1, input_files) ; this_year = which(this_year == TRUE)
             # We should have a single file per year containing both the estiate and the standard deviation
             if (length(this_year) == 1) {
                 # Track the number of years we have information for
                 keepers = keepers+1
                 years_with_obs = append(years_with_obs,years_to_load[yr])
             } else {
                 missing_years = append(missing_years,years_to_load[yr])
             }
        } # loop through possible years
        missing_years = missing_years[-1]
        years_with_obs = years_with_obs[-1]

        # Warn the user if there are no data found
        if (length(missing_years) == length(years_to_load)) {
            print(paste("WARNINGS: ",est_var_name_in," have been requested but none found for the analysis time period",sep=""))
            # Create output object
            output_all = list(-9999, -9999, doy_obs = -9999, years = -9999, lat = -9999, long = -9999, missing_years = missing_years, data_available = FALSE) 
            # Update with the correct variable names
            names(output_all)[1:2]<-c(est_var_name_out,unc_var_name_out)
            # Return to function
            return(output_all)
        }

        # Flag to ensure we create output variables once
        lat_done = FALSE ; done_first_time = FALSE ; years_loaded = 0
        # Loop for year here
        for (yr in seq(1, length(years_to_load))) {

             ## Begin reading the files in now for real

             # Update the user as to our progress
             print(paste("...",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

             # Determine the unique file name pattern
             input_file_1 = paste(prefix,years_to_load[yr],sep="")

             # Then check whether this pattern is found in the available files
             this_year = input_files[grepl(input_file_1, input_files)]
             if (length(this_year) > 0) {

                 # Update year counter for output
                 years_loaded = append(years_loaded, years_to_load[yr])

                 # open the file
                 data1 = nc_open(this_year)

                 twodim = FALSE
                 # Get timing variable...
                 if (length(which(names(data1$var) == "doy")) > 0 | length(which(names(data1$dim) == "doy"))) {
                     doy_in = ncvar_get(data1, "doy") 
                 } else {
                     # We don't have the desired time variable
                     print(paste("......doy variable missing from ",est_var_name_in," Gridded_nc variable",sep=""))
                     if (data1$ndim == 2) {
                         print("......the code will assume that doy is the middle of the year, assuming only 2 dimensions are found in the file")
                         doy_in = 187 # middle day of the year
                         twodim = TRUE # flag to allow for correction to the dimension in the read variable
                     } else {
                         stop("......doy missing and there appears to be >2 dimension, i.e. more than x~y")
                     }
                 } # checking for doy variable
                 #...and accumulate for the overall vector
                 doy_out = append(doy_out,doy_in)
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
                 # Extract the current global attributes
                 global_attributes = ncatt_get(data1,0)
                 # Check whether there is any information regarding the EPSG
                 global_attributes = unlist(global_attributes)
                 # Default assumption for netcdf files is for the epsg: 4326 i.e. the WGS-84 lat/long grid
                 # But here we will search for any specific information
                 epsg = 4326 ; aa = 1
                 if (length(global_attributes) > 0) {
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
                 } # Assuming we have any global attributes

                 # read the observation estimate 
                 est_in = ncvar_get(data1, est_var_name_in) # Variable estimate
                 if (twodim) {est_in = array(est_in, dim=c(dim(est_in),1))}
                 # read error variable, if present
                 if (length(which(names(data1$var) == unc_var_name_in)) > 0) {
                     std_in = ncvar_get(data1, unc_var_name_in) # Variable standard deviation 
                     std_present = TRUE
                     if (twodim) {std_in = array(std_in, dim=c(dim(std_in),1))}
                 } else {
                     std_in = -9999 ; std_present = FALSE
                 }                 
                 # Close the current file
                 nc_close(data1) 

                 # Loop through the file to aggregate each time step in turn
                 for (t in seq(1, length(doy_in))) {

                      # Convert to a raster, assuming standard WGS84 grid
                      var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(est_in[,,t]))
                      var1 = rast(var1, crs = paste("epsg:",epsg,sep=""), type="xyz")
                      # Check that we have a uncertainty estimate
                      if (std_present) {
                          var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(std_in[,,t]))
                          var2 = rast(var2, crs = paste("epsg:",epsg,sep=""), type="xyz")
                      } 

                      # Extract the epsg from the file
                      epsg = crs(var1, describe = TRUE)$code
                      if (is.null(epsg) | epsg == "") { stop(paste("......the use_lcm specification leads to a geotif which does not contain epsg information."))}
                      # If we have an epsg then we want to know if it differs from the one desired by the analysis
                      if (epsg != gsub("epsg:","",grid_type)) {
                          # Ensure that the extent of the input object is consistent 
                          # with the possible extent of the selected epsg
                          var1 = crop(var1, ext(unlist(crs(grid_type, describe=TRUE)$extent)))                      
                          # If it does not match we need to reproject it
                          var1 = project(var1, grid_type, method="near", align = FALSE) ; gc()
                          if (std_present) { 
                              # Ensure that the extent of the input object is consistent 
                              # with the possible extent of the selected epsg
                              var2 = crop(var2, ext(unlist(crs(grid_type, describe=TRUE)$extent)))
                              var2 = project(var2, grid_type, method="near", align = FALSE) ; gc() 
                          }
                      }
 
                      # Use extend and crop to ensure we closely match the extent of the analysis
                      var1 = extend(var1,cardamom_ext) ; var1 = crop(var1,cardamom_ext) 
                      if (std_present) { var2 = extend(var2,cardamom_ext) ; var2 = crop(var2,cardamom_ext) }

                      # Adjust spatial resolution of the datasets, this occurs in all cases
                      if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                          # Resample to correct grid.
                          # Probably should be done via aggregate function to allow for correct error propogation
                          var1 = resample(var1, cardamom_ext, method="average") ; gc() 
                          if (std_present) { var2 = resample(var2, cardamom_ext, method="average") ; gc() }
                      } # Aggrgeate to resolution

                      # Combine estimate and uncertainty variables into a stacked raster
                      if (done_first_time == FALSE) {
                          # create output raster
                          est_out_tif = var1 ; rm(var1)
                          if (std_present) { std_out_tif = var2 ; rm(var2) }
                          done_first_time = TRUE                   
                      } else {                  
                          # add to existing output raster
                          add(est_out_tif) <- var1 ; rm(var1)
                          if (std_present) { add(std_out_tif) <- var2 ; rm(var2) }
                      }                     

               } # loop through available time steps in the current year

               # keep track of years actually ran
               yrs = yrs + 1
               # clean up allocated memeory
               gc()

           } # is there information for the current year?

      } # year loop

      # Extract spatial information just the once
      if (lat_done == FALSE) {
          # Set flag to true
          lat_done = TRUE
          # extract dimension information for the grid, 
          # note the axis switching between raster and actual array
          xdim = dim(est_out_tif)[2] ; ydim = dim(est_out_tif)[1]
          # extract the lat / long information needed
          long = crds(est_out_tif, df=TRUE, na.rm=FALSE)
          lat  = long$y ; long = long$x
          # restructure into correct orientation
          long = array(long, dim=c(xdim,ydim))
          lat = array(lat, dim=c(xdim,ydim))
      }

      # Correct for initialisation
      doy_out = doy_out[-1] ; years_loaded = years_loaded[-1]

      # Sanity check for dataset
      if (lat_done == FALSE) {stop(paste('No ',est_var_name_in,' information could be found...',sep=""))}

      # Extract from the raster structure into arrays
      est_out = array(NA, dim=c(xdim,ydim,length(doy_out)))
      for (d in seq(1, length(doy_out))) {
           est_out[,,d] = array(values(subset(est_out_tif, d)), dim=c(xdim,ydim))
      }
      if (std_present) { 
          # Extract from the raster structure into arrays
          std_out = array(NA, dim=c(xdim,ydim,length(doy_out)))      
          for (d in seq(1, length(doy_out))) {
               std_out[,,d] = array(values(subset(std_out_tif, d)), dim=c(xdim,ydim))
          }
      }

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
      if (length(missing_years) == 0) {missing_years = -9999}

      # Create output object
      output_all = list(est_out, std_out, doy_obs = doy_out, years = years_with_obs, lat = lat, long = long, missing_years = missing_years, data_available = TRUE) 
      # Update with the correct variable names
      names(output_all)[1:2]<-c(est_var_name_out,unc_var_name_out)

      # clean up variables
      rm(doy_in,est_out,std_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)

      # Return the output to the user
      return(output_all)

    } else if (data_source == "Gridded_tif") { 

        # let the user know this might take some time
        print(paste("Loading ",est_var_name_in," for subsequent sub-setting ...",sep=""))

        # list all available files which we will then search
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_files = list.files(data_path, full.names=TRUE, recursive = TRUE, pattern="\\.tif$")
        input_files = input_files[grepl(paste("/",prefix,sep=""),input_files)]
        #prefix = "MCD15A2H_LAI_(.)*" # (.)* wildcard characters for unix standard MCD15A2H_LAI_*
        #prefix = "net_biome_exchange_"
   
        # Search files with the correct prefix only  
        # Assumed file structure for estimate file
        # <variable_name>_<unit>_YYYY-DOY.tif
        # For uncertainty file
        # <variable_name>_<unit>_uncertainty_YYYY-DOY.tif
        input_files = input_files[grepl(prefix, input_files)]

        # Pull out those specifically for the variable of interest
        # based on them not being those with uncertainty
        est_input_files = input_files[grepl("uncertainty", input_files) == FALSE]
        # As we are using geotifs, we need to be sure we have matching estimate
        # and uncertainty files. Assuming we have asked for an uncertainty variable.
        if (unc_var_name_in != "") {
            unc_input_files = input_files[grepl("uncertainty", input_files)]
            if (length(est_input_files) != length(unc_input_files)) {
                stop(paste("The number of estimate files for ",prefix,
                           " and its uncertainty ",prefix,"_uncertainty do not match",sep=""))
            }
            std_present = TRUE
        }
            
        # Begin extraction of the year information found in the file name
        years_with_obs = gsub(data_path,"",est_input_files)
        years_with_obs = gsub(prefix,"",years_with_obs)
        years_with_obs = gsub("/","",years_with_obs)
        years_with_obs = gsub("\\.tif$","",years_with_obs)
        # Remove any further underscores, 
        # this should leave us with the time information alone
        years_with_obs = gsub("_","",years_with_obs)
        # Extract the first 4 characters as these should be YYYY
        # Check whether the file includes the option day of year (doy) 
        # information
        tmp1 = rep(NA, length(years_with_obs))
        doy_out = rep(NA, length(years_with_obs))
        for (y in seq(1, length(years_with_obs))) {
             # Extract the year information
             tmp1[y] = as.numeric(substr(years_with_obs[y],1,4))
             # Extract the day of year information
             doy_out[y] = as.numeric(substr(years_with_obs[y],6,8))
        } 
        years_with_obs = tmp1 ; rm(tmp1)
        # Therefore we can determine the years without data, 
        # i.e. the missing years
        missing_years = 0
        for (y in seq(1, length(years_to_load))) {
             if (length(which(years_with_obs == years_to_load[y])) > 0) {
                 # Do nothing
             } else {
                 missing_years = append(missing_years, as.numeric(years_to_load[y]))
             }
        }
        missing_years = missing_years[-1]

        # Warn the user if there are no data found
        if (length(missing_years) == length(years_to_load)) {
            print(paste("WARNINGS: ",est_var_name_in," have been requested but none found for the analysis time period",sep=""))
            # Create output object
            output_all = list(-9999, -9999, doy_obs = -9999, years = -9999, lat = -9999, long = -9999, missing_years = missing_years, data_available = FALSE) 
            # Update with the correct variable names
            names(output_all)[1:2]<-c(est_var_name_out,unc_var_name_out)
            # Return to function
            return(output_all)
        }

        # Loop through each year and extract if appropriate
        lat_done = FALSE ; done_first_time = FALSE 
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(years_to_load[1]) & years_with_obs[t] <= as.numeric(years_to_load[length(years_to_load)])) {

                 # Subset to the files found in the current year
                 est_files = est_input_files[grepl(years_with_obs[t],est_input_files)]
                 if (std_present) {
                     std_files = unc_input_files[grepl(years_with_obs[t],unc_input_files)]
                 }

                 # Loop through all available steps in the current year
                 for (tt in seq(1, length(est_files))) {

                      # Read in the estimate and uncertainty rasters
                      var1 = rast(est_files[tt])
                      if (std_present) { var2 = rast(std_files[tt]) }

                      # Extract the epsg from the file
                      epsg = crs(var1, describe = TRUE)$code
                      if (is.null(epsg) | epsg == "") { stop(paste("the use_lcm specification leads to a geotif which does not contain epsg information."))}
                      # If we have an epsg then we want to know if it differs from the one desired by the analysis
                      if (epsg != gsub("epsg:","",grid_type)) {
                          # Ensure that the extent of the input object is consistent 
                          # with the possible extent of the selected epsg
                          var1 = crop(var1, ext(unlist(crs(grid_type, describe=TRUE)$extent)))                                            
                          # If it does not match we need to reproject it
                          var1 = project(var1, grid_type, method="near", align = FALSE) ; gc()
                          if (std_present) { 
                              # Ensure that the extent of the input object is consistent 
                              # with the possible extent of the selected epsg
                              var2 = crop(var2, ext(unlist(crs(grid_type, describe=TRUE)$extent)))                                                
                              var2 = project(var2, grid_type, method="near", align = FALSE) ; gc() 
                          }
                      }

                      # Use extend and crop to ensure we closely match the extent of the analysis
                      var1 = extend(var1,cardamom_ext) ; var1 = crop(var1,cardamom_ext) 
                      if (std_present) { var2 = extend(var2,cardamom_ext) ; var2 = crop(var2,cardamom_ext) }

                      # Adjust spatial resolution of the datasets, this occurs in all cases
                      if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                          # Resample to correct grid.
                          # Probably should be done via aggregate function to allow for correct error propogation
                          var1 = resample(var1, cardamom_ext, method="average") ; gc() 
                          if (std_present) { var2 = resample(var2, cardamom_ext, method="average") ; gc() }
                      } # Aggrgeate to resolution

                      # Combine estimate and uncertainty variables into a stacked raster
                      if (done_first_time == FALSE) {
                          # create output raster
                          est_out_tif = var1 ; rm(var1)
                          if (std_present) { std_out_tif = var2 ; rm(var2) }
                          done_first_time = TRUE                   
                      } else {                  
                          # add to existing output raster
                          add(est_out_tif) <- var1 ; rm(var1)
                          if (std_present) { add(std_out_tif) <- var2 ; rm(var2) }
                      }                     

                 } # Looping with within year

             } # Is dataset within the analysis time period?

        } # looping available years

        # Extract spatial information just the once
        if (lat_done == FALSE) {
            # Set flag to true
            lat_done = TRUE
            # extract dimension information for the grid, 
            # note the axis switching between raster and actual array
            xdim = dim(est_out_tif)[2] ; ydim = dim(est_out_tif)[1]
            # extract the lat / long information needed
            long = crds(est_out_tif, df=TRUE, na.rm=FALSE)
            lat  = long$y ; long = long$x
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
        }

        # Sanity check for dataset
        if (lat_done == FALSE) {stop(paste('No ',est_var_name_in,' information could be found...',sep=""))}

        # Extract from the raster structure into arrays
        est_out = array(NA, dim=c(xdim,ydim,length(doy_out)))
        for (d in seq(1, length(doy_out))) {
             est_out[,,d] = array(values(subset(est_out_tif, d)), dim=c(xdim,ydim))
        }
        if (std_present) { 
            # Extract from the raster structure into arrays
            std_out = array(NA, dim=c(xdim,ydim,length(doy_out)))      
            for (d in seq(1, length(doy_out))) {
                 std_out[,,d] = array(values(subset(std_out_tif, d)), dim=c(xdim,ydim))
            }
        }

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
        if (length(missing_years) == 0) {missing_years = -9999}

        # Create output object
        output_all = list(est_out, std_out, doy_obs = doy_out, years = years_with_obs, lat = lat, long = long, missing_years = missing_years, data_available = TRUE) 
        # Update with the correct variable names
        names(output_all)[1:2]<-c(est_var_name_out,unc_var_name_out)

        # clean up variables
        rm(est_out,std_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
        return(output_all)

    } else if (data_source == " " | data_source == "site_specific") {
    
        # Do nothing as this should be not needed or will be read later
    
    } else {

        # We have a problem as something unexpected has been specified
        stop(paste("The data_source = ",data_source," is not a recognised value (Gridded_nc, Gridded_tif, site_specific or blank)",sep=""))

    } # Gridded dataset or not

} # function end load_observation_dataset_for_extraction

## Use byte compile
load_observation_dataset_for_extraction<-cmpfun(load_observation_dataset_for_extraction)
