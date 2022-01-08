
###
## Function to load leaf area index from global databases
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_lai_fields_for_extraction<-function(latlon_in,lai_source,years_to_load,cardamom_ext,spatial_type) {

  if (lai_source == "MODIS") {

      # let the user know this might take some time
      print("Loading processed lai fields for subsequent sub-setting ...")

    lat_done = FALSE ; missing_years=0 ; keepers=0 ; yrs=1
    # loop for year here
    for (yr in seq(1, length(years_to_load))) {
      print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

      if (yr == 1) {
        # first check how many files we have
        for (yrr in seq(1, length(years_to_load))) {
          # open processed modis files
          input_file_1=paste(path_to_lai,"modis_lai_with_lat_long_",years_to_load[yrr],".nc",sep="")
          if (file.exists(input_file_1) == TRUE) {keepers=keepers+1} else {missing_years=append(missing_years,years_to_load[yrr])}
        }
      }
      # open processed modis files
      input_file_1=paste(path_to_lai,"modis_lai_with_lat_long_",years_to_load[yr],".nc",sep="")

      # check to see if file exists if it does then we read it in, if not then we assume its a year we don't have data for and move on
      if (file.exists(input_file_1) == TRUE) {
        # open the file
        data1=nc_open(input_file_1)

        # get timing variable
        doy_in=ncvar_get(data1,"doy")

        # extract location variables
        if (lat_done == FALSE) {
          lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
          # restrict the spatial extent based on latlong ranges provided
          remove_lat=intersect(which(lat < (max(latlon_in[,1])+2)),which(lat > (min(latlon_in[,1])-2)))
          remove_long=intersect(which(long < (max(latlon_in[,2])+2)),which(long > (min(latlon_in[,2])-2)))
          # now find common where out in both contexts
          remove_lat=intersect(remove_lat,remove_long)
          # update both variables because of common matrix
          remove_long=remove_lat
          # adjust for matrix rather than vector arrangement
          remove_lat=remove_lat/dim(lat)[1]
          remove_long=(remove_long-(floor(remove_lat)*dim(lat)[1]))+1
          remove_lat=ceiling(remove_lat)
          # update new dimensions
          lat_dim=length(min(remove_lat):max(remove_lat)) ; long_dim=length(min(remove_long):max(remove_long))
          lat=lat[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)] ; long=long[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
        }

        # read the modis lai drivers
        var1=ncvar_get(data1, "modis_lai")
        # reduce spatial cover to the desired area only
        var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
        # set actual missing data to -9999
        var1[which(is.na(as.vector(var1)))] = -9999

        # close files after use
        nc_close(data1)

        # remove additional spatial information
        if (lat_done == FALSE) {
          lai_hold=array(NA, dim=c(long_dim*lat_dim*48,keepers))
          lai_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
          doy_out=doy_in
        } else {
          lai_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
          doy_out=append(doy_out,doy_in)
        }

        # update flag for lat / long load
        if (lat_done == FALSE) {lat_done = TRUE}
        # keep track of years actually ran
        yrs=yrs+1
        # clean up allocated memeory
        rm(var1) ; gc()
      } # end of does file exist

    } # year loop

    # Sanity check for LAI
    if (lat_done == FALSE) {stop('No LAI information could be found...')}

    # remove initial value
    missing_years=missing_years[-1]

    # check which ones are NA because I made them up
    not_na=is.na(as.vector(lai_hold))
    not_na=which(not_na == FALSE)

    # now remove the ones that are actual missing data
    lai_hold[which(as.vector(lai_hold) == -9999)] = NA
    # return spatial structure to data
    lai_out=array(as.vector(lai_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))
    lai_unc_out = array(-9999, dim=dim(lai_out))
    # output variables
    lai_all = list(lai_all = lai_out, lai_unc_all = lai_unc_out,
                  doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
    # clean up variables
    rm(doy_in,lai_hold,not_na,lai_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
    return(lai_all)

  } else if (lai_source == "COPERNICUS") {

    # let the user know this might take some time
    print("Loading processed lai fields for subsequent sub-setting ...")

    # check which file prefix we are using today
    # list all available files which we will then search
    avail_files = list.files(path_to_lai,full.names=TRUE)
    prefix = "c_gls(.)*_" # (.)* wildcard characters for unix standard c_gls*_

    # timing information on the number of day in a month
    month_days = rep(31,length.out=12)
    month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30

    lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1 ; doy_out = 0
    # loop for year here
    for (yr in seq(1, length(years_to_load))) {

         # Update the user as to our progress
         print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

         # If we are just starting check how many files we have
         if (yr == 1) {
             nsteps = 0
             # Loop through all the analyses years and check whether files exist for it
             for (yrr in seq(1, length(years_to_load))) {
                  # create the prefix to the files we will want for a given year
                  input_file_1=paste(prefix,years_to_load[yrr],sep="")
                  # then check whether this pattern is found in the available files
                  this_year = grepl(input_file_1, avail_files) ; this_year = which(this_year == TRUE)
                  # if we have at least one timestep for this year then we have some information otherwise it is missing!
                  if (length(this_year) > 0) {
                      keepers = keepers+1 ; nsteps = max(nsteps,length(this_year))
                  } else {
                      missing_years = append(missing_years,years_to_load[yrr])
                  }
             } # loop through possible years
             rm(yrr)
         } # first year?

         # Begin reading the files in now for real

         # Determine the unique file name pattern
         input_file_1=paste(prefix,years_to_load[yr],sep="")

         # then check whether this pattern is found in the available files
         this_year = avail_files[grepl(input_file_1, avail_files)]
         if (length(this_year) > 0) {

             # Strip out the time component entirely
             tmp = unlist(strsplit(this_year,input_file_1))
             tmp = tmp[seq(2,length(tmp),2)]
             # The re-order the files correctly
             this_year = this_year[order(tmp)]

             # now loop through the available files for the current year
             for (t in seq(1, length(this_year))) {

                  # open the file
                  data1 = nc_open(this_year[t])

                  # get timing variable
                  tmp = strsplit(this_year[t],input_file_1)[[1]][2]
                  month = as.numeric(substring(tmp,1,2)) ; doy_in = as.numeric(substring(tmp,3,4))
                  # January is correct already, so only adjust if month is >= February
                  if (month > 1) {
                      doy_in = doy_in + sum(month_days[1:(month-1)])
                  }

                  # Extract spatial information
                  lat_in = ncvar_get(data1, "lat") ; long_in = ncvar_get(data1, "lon")

                  # read the LAI observations
                  var1 = ncvar_get(data1, "LAI") # leaf area index (m2/m2)
                  # check for error variable
                  if (length(which(grepl("LAI_ERR",names(data1$var)) == TRUE)) > 0) {
                      var2 = ncvar_get(data1, "LAI_ERR") # standard error (m2/m2)
                  } else if (length(which(grepl("RMSE",names(data1$var)) == TRUE)) > 0) {
                      var2 = ncvar_get(data1, "RMSE") # standard error (m2/m2)
                  } else {
                      stop("LAI error variable cannot be found for copernicus...")
                  }

                  # Close the current file
                  nc_close(data1)

                  # Turn lat_in / long_in from vectors to arrays
                  lat_in = t(array(lat_in, dim=c(dim(var1)[2],dim(var1)[1])))
                  long_in = array(long_in, dim=c(dim(var1)[1],dim(var1)[2]))

                  # Convert to a raster, assuming standad WGS84 grid
                  var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                  var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
                  var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var2))
                  var2 = rasterFromXYZ(var2, crs = ("+init=epsg:4326"))
                  # Remove the input lat / long information
                  rm(lat_in,long_in)

                  # Trim the extent of the overall grid to the analysis domain
                  var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)
                  # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
                  # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
                  if (spatial_type == "grid") {
                      if (res(var1)[1] < res(cardamom_ext)[1] | res(var1)[2] < res(cardamom_ext)[2]) {
                          # Create raster with the target resolution
                          target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                          # Resample to correct grid.
                          # Probably should be done via aggregate function to allow for correct error propogation
                          var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()
                          var2 = resample(var2, target, method="bilinear") ; gc() ; removeTmpFiles()
                      } # Aggrgeate to resolution
                  } # spatial_type == "grid"

                  # Extract spatial information just the once
                  if (lat_done == FALSE) {
                      # Set flag to true
                      lat_done = TRUE
                      # extract dimension information for the grid, note the axis switching between raster and actual array
                      xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                      # extract the lat / long information needed
                      long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
                      # restructure into correct orientation
                      long = array(long, dim=c(xdim,ydim))
                      lat = array(lat, dim=c(xdim,ydim))
                      # create holding arrays for the lai information...
                      lai_hold = array(NA, dim=c(xdim*ydim,keepers*nsteps))
                      lai_unc_hold = array(NA, dim=c(xdim*ydim,keepers*nsteps))
                  }
                  # break out from the rasters into arrays which we can manipulate
                  var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))
                  var2 = array(as.vector(unlist(var2)), dim=c(xdim,ydim))
                  # set actual missing data to -9999
                  var1[which(is.na(as.vector(var1)))] = -9999
                  var2[which(is.na(as.vector(var2)))] = -9999

                  # begin populating the various outputs
                  lai_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
                  lai_unc_hold[1:length(as.vector(var2)),(t+((yrs-1)*nsteps))] = as.vector(var2)
                  doy_out = append(doy_out,doy_in)

             } # loop through available time steps in the current year

             # keep track of years actually ran
             yrs = yrs + 1
             # clean up allocated memeory
             rm(var1,var2) ; gc()

         } # is there information for the current year?

    } # year loop

    # Correct for initialisation
    doy_out = doy_out[-1]

    # Sanity check for LAI
    if (lat_done == FALSE) {stop('No LAI information could be found...')}

    # remove initial value
    missing_years = missing_years[-1]

    # check which ones are NA because I made them up
    not_na = is.na(as.vector(lai_hold))
    not_na = which(not_na == FALSE)

    filter = as.vector(lai_hold) == -9999 | as.vector(lai_unc_hold) == -9999
    # now remove the ones that are actual missing data
    lai_hold[filter] = NA ; lai_unc_hold[filter] = NA
    # return spatial structure to data
    lai_out = array(as.vector(lai_hold)[not_na], dim=c(xdim,ydim,length(doy_out)))
    lai_unc_out = array(as.vector(lai_unc_hold)[not_na], dim=c(xdim,ydim,length(doy_out)))

    # output variables
    lai_all = list(lai_all = lai_out, lai_unc_all = lai_unc_out,
                   doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
    # clean up variables
    rm(doy_in,lai_hold,lai_unc_hold,not_na,lai_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
    return(lai_all)

  } # if MODIS or COPERNICUS

} # function end
## Use byte compile
load_lai_fields_for_extraction<-cmpfun(load_lai_fields_for_extraction)
