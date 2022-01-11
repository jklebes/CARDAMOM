
###
## Function to load burned area data from various gridded datasets
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_burnt_area_fields_for_extraction<-function(latlon_in,burnt_area_source,path_to_burnt_area,start_year,end_year,cardamom_ext,spatial_type) {

    # Determine the years for which the analysis will occur
    years_to_do = as.character(seq(start_year,end_year))

    if (burnt_area_source == "GFED4") {

        # let the user know this might take some time
        print("Loading GFED4 processed burnt area estimates for subsequent sub-setting ...")

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

             if (yr == 1) {
                 # list the available files
                 available_files = list.files(path_to_burnt_area,full.names=TRUE)
                 # first check how many files we have
                 for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(years_to_do[yrr],available_files))) > 0) {
                          keepers = keepers+1
                      } else {
                          missing_years = append(missing_years,years_to_do[yrr])
                      } # missing years
                 } # loop years to check which we have
	            } # first year
              # open processed modis files
  	          input_file_1 = paste(path_to_burnt_area,"/GFED4_",years_to_do[yr],".nc",sep="")

              # check to see if file exists if it does then we read it in,
              # if not then we assume its a year we don't have data for and move on
	            if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  lat_in = ncvar_get(data1, "latitude") ; long_in = ncvar_get(data1, "longitude")
                  # read the burnt fraction estimate (units are 0-1)
                  var1 = ncvar_get(data1, "BurnedFraction")
                  # get time information (month in this case)
                  var2 = ncvar_get(data1, "time") ; time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate fire to that point
                  if (lat_done == FALSE) {
                      doy_obs = floor((var2*(365.25/12))-(365.25/24))
                  } else {
                      doy_obs = append(doy_obs,floor((var2*(365.25/12))-(365.25/24)))
                  }

                  # close files after use
                  nc_close(data1) ; rm(var2)

                  # Turn lat_in / long_in from vectors to arrays
                  lat_in = t(array(lat_in, dim=c(dim(var1)[2],dim(var1)[1])))
                  long_in = array(long_in, dim=c(dim(var1)[1],dim(var1)[2]))

                  # Convert to a raster, assuming standad WGS84 grid
                  var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                  var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
                  # Remove the input lat / long information
                  rm(lat_in,long_in)

                  # Trim the extent of the overall grid to the analysis domain
                  var1 = crop(var1,cardamom_ext)
                  # set actual missing data to 0 as missing data is actually no fire
                  var1[which(is.na(as.vector(var1)))] = 0
                  # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
                  # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
                  if (spatial_type == "grid") {
                      if (res(var1)[1] < res(cardamom_ext)[1] | res(var1)[2] < res(cardamom_ext)[2]) {
                          # Create raster with the target resolution
                          target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                          # Resample to correct grid.
                          # Probably should be done via aggregate function to allow for correct error propogation
                          var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()
                      } # Aggrgeate to resolution
                  } # spatial_type == "grid"

                  # Extract spatial information just the once
                  if (lat_done == FALSE) {
                      # extract dimension information for the grid, note the axis switching between raster and actual array
                      xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                      # extract the lat / long information needed
                      long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
                      # restructure into correct orientation
                      long = array(long, dim=c(xdim,ydim))
                      lat = array(lat, dim=c(xdim,ydim))
                  }
                  # break out from the rasters into arrays which we can manipulate
                  var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

                  # vectorise at this time
                  if (lat_done == FALSE) {
                      burnt_area = as.vector(var1)
                  } else {
                      burnt_area = append(burnt_area,as.vector(var1))
                  }

                  # update flag for lat / long load
                  if (lat_done == FALSE) {lat_done = TRUE}
                  # keep track of years actually ran
                  yrs = yrs+1

	            } # end of does file exist

  	     } # year loop

         # remove initial value
         missing_years=missing_years[-1]

         # clean up variables
         rm(var1,var2) ; gc(reset=TRUE,verbose=FALSE)

         # restructure
         burnt_area=array(burnt_area, dim=c(xdim,ydim,length(doy_obs)))

         # output variables
         return(list(burnt_area=burnt_area,doy_obs=doy_obs,lat=lat,long=long,missing_years=missing_years))

    } else if (burnt_area_source == "MCD64A1") {

         # let the user know this might take some time
         print("Loading MCD64A1 processed burnt area estimates for subsequent sub-setting ...")

         lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
         # loop for year here
         for (yr in seq(1, length(years_to_do))) {
              print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

              if (yr == 1) {
                  # list the available files
                  available_files = list.files(path_to_burnt_area,full.names=TRUE)
                  # first check how many files we have
                  for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(years_to_do[yrr],available_files))) > 0) {keepers = keepers+1} else {missing_years = append(missing_years,years_to_do[yrr])}
                  }
              }
              # open processed modis files
              input_file_1 = paste(path_to_burnt_area,"/MCD64A1_",years_to_do[yr],".nc",sep="")

              # check to see if file exists if it does then we read it in, if not then we assume its a year we don't have data for and move on
              if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  lat_in = ncvar_get(data1, "latitude") ; long_in = ncvar_get(data1, "longitude")
                  # read the burnt fraction estimate (units are 0-1)
                  var1 = ncvar_get(data1, "BurnedFraction")

                  # get time information (day of year)
                  var2 = ncvar_get(data1, "time") ; time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate fire to that point
                  if (lat_done == FALSE) {
                      doy_obs = var2
                  } else {
                      doy_obs = append(doy_obs,var2)
                  }

                  # close files after use
                  nc_close(data1) ; rm(var2)

                  # Convert to a raster, assuming standad WGS84 grid
                  var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                  var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
                  # Remove the input lat / long information
                  rm(lat_in,long_in)

                  # Trim the extent of the overall grid to the analysis domain
                  var1 = crop(var1,cardamom_ext)
                  # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
                  # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
                  if (spatial_type == "grid") {
                      if (res(var1)[1] < res(cardamom_ext)[1] | res(var1)[2] < res(cardamom_ext)[2]) {
                          # Create raster with the target resolution
                          target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                          # Resample to correct grid.
                          # Probably should be done via aggregate function to allow for correct error propogation
                          var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()
                      } # Aggrgeate to resolution
                  } # spatial_type == "grid"

                  # Extract spatial information just the once
                  if (lat_done == FALSE) {
                      # extract dimension information for the grid, note the axis switching between raster and actual array
                      xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                      # extract the lat / long information needed
                      long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
                      # restructure into correct orientation
                      long = array(long, dim=c(xdim,ydim))
                      lat = array(lat, dim=c(xdim,ydim))
                  }
                  # break out from the rasters into arrays which we can manipulate
                  var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

                  # vectorise at this time
                  if (lat_done == FALSE) {
                      burnt_area = as.vector(var1)
                  } else {
                      burnt_area = append(burnt_area,as.vector(var1))
                  }

                  # update flag for lat / long load
                  if (lat_done == FALSE) {lat_done = TRUE}
                  # keep track of years actually ran
                  yrs = yrs+1

              } # end of does file exist

          } # year loop

          # remove initial value
          missing_years = missing_years[-1]

          # clean up variables
          rm(var1,var2) ; gc(reset=TRUE,verbose=FALSE)

          # restructure
          burnt_area=array(burnt_area, dim=c(xdim,ydim,length(doy_obs)))

          # output variables
          return(list(burnt_area=burnt_area,doy_obs=doy_obs,lat=lat,long=long,missing_years=missing_years))

    } else if (burnt_area_source == " " | burnt_area_source == "site_specific"){

	        # Do nothing as this should be read directly from files or not needed

    } else {

	        stop(paste("Burnt area option (",burnt_area_source,") not valid"))

    } # if MPI biomass

} # function end

## Use byte compile
load_burnt_area_fields_for_extraction<-cmpfun(load_burnt_area_fields_for_extraction)
