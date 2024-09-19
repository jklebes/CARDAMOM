
###
## Function to load fraction of absorbed photosynthetically active radition from global databases
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_fapar_fields_for_extraction<-function(latlon_in,fapar_source,years_to_load,cardamom_ext,spatial_type) {

  if (fapar_source == "Gridded"){

      # let the user know this might take some time
      print("Loading processed fAPAR fields for subsequent sub-setting ...")

      # check which file prefix we are using today
      # list all available files which we will then search
      avail_files = list.files(path_to_fapar,full.names=TRUE)
      #prefix = "MCD15A2H_LAI_(.)*" # (.)* wildcard characters for unix standard MCD15A2H_LAI_*
      prefix = "fraction_absorbed_par_"

      # timing information on the number of day in a month
      month_days = rep(31,length.out=12)
      month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30 

      lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1 ; doy_out = 0
      # Loop through all the analyses years and check whether files exist for it
      for (yr in seq(1, length(years_to_load))) {
           # create the prefix for the observation file files we will want for a given year
           input_file_1 = paste(prefix,years_to_load[yr],sep="")
           # then check whether this pattern is found in the available files
           this_year = grepl(input_file_1, avail_files) ; this_year = which(this_year == TRUE)
           # We should have a single file per year containing both the estiate and the standard deviation
           if (length(this_year) == 1) {
               # Track the number of years we have information for
               keepers = keepers+1
           } else {
               missing_years = append(missing_years,years_to_load[yr])
           }
      } # loop through possible years
      missing_years = missing_years[-1]

      # Warn the user if there are no data found
      if (length(missing_years) == length(years_to_do)) {
          print("WARNINGS: fAPAR observations have been requested but none found for the analysis time period")
          # output variables
          fapar_all = list(fapar_all = -9999, fapar_unc_all = -9999,
                           doy_obs = -9999, lat = -9999, long = -9999, missing_years=missing_years)
          return(fapar_all)
      }

      # Flag to ensure we create output variables once
      done_first_time = FALSE
      # Loop for year here
      for (yr in seq(1, length(years_to_load))) {

           # Update the user as to our progress
           print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

           ## Begin reading the files in now for real

           # Determine the unique file name pattern
           input_file_1=paste(prefix,years_to_load[yr],sep="")

           # Then check whether this pattern is found in the available files
           this_year = avail_files[grepl(input_file_1, avail_files)]
           if (length(this_year) > 0) {

               # open the file
               data1 = nc_open(this_year)

               # Get timing variable ; and accumulate for the overall vector
               doy_in = ncvar_get(data1, "doy") ; doy_out = append(doy_out,doy_in)
               # Extract spatial information
               lat_in = ncvar_get(data1, "lat") ; long_in = ncvar_get(data1, "lon")
               # read the observations
               fapar_est_in = ncvar_get(data1, "fAPAR") # fraction absorbed PAR (0-1)
               # read error variable
               fapar_std_in = ncvar_get(data1, "fAPAR_SD") # standard deviation (0-1)
               # Close the current file
               nc_close(data1) 

               # Loop through the file to aggregate each time step in turn
               for (t in seq(1, length(doy_in))) {

                     # Convert to a raster, assuming standad WGS84 grid
                     var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(fapar_est_in[,,t]))
                     var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")
                     var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(fapar_std_in[,,t]))
                     var2 = rast(var2, crs = ("+init=epsg:4326"), type="xyz")

                     # Extend the extent of the overall grid to the analysis domain
                     var1 = extend(var1,cardamom_ext) ; var2 = extend(var2,cardamom_ext)
                     # Trim the extent of the overall grid to the analysis domain
                     var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)

                     # Adjust spatial resolution of the datasets, this occurs in all cases
                     if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                         # Create raster with the target resolution
                         target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
                         # Resample to correct grid.
                         # Probably should be done via aggregate function to allow for correct error propogation
                         var1 = resample(var1, target, method="bilinear") ; gc() 
                         var2 = resample(var2, target, method="bilinear") ; gc()
                     } # Aggrgeate to resolution

                     # Combine estimate and uncertainty variables into a stacked raster
                     if (done_first_time == FALSE) {
                         # create output raster
                         fapar_est_out = var1 ; rm(var1)
                         fapar_std_out = var2 ; rm(var2)     
                         done_first_time = TRUE                   
                     } else {                  
                         # add to existing output raster
                         add(fapar_est_out) <- var1 ; rm(var1)
                         add(fapar_std_out) <- var2 ; rm(var2)                        
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
          # extract dimension information for the grid, note the axis switching between raster and actual array
          xdim = dim(fapar_est_out)[2] ; ydim = dim(fapar_est_out)[1]
          # extract the lat / long information needed
          long = crds(fapar_est_out, df=TRUE, na.rm=FALSE)
          lat  = long$y ; long = long$x
          # restructure into correct orientation
          long = array(long, dim=c(xdim,ydim))
          lat = array(lat, dim=c(xdim,ydim))
      }

      # Correct for initialisation
      doy_out = doy_out[-1]

      # Sanity check for fAPAR
      if (lat_done == FALSE) {stop('No fAPAR information could be found...')}

      # Extract from the raster structure into arrays
      fapar_out = array(NA, dim=c(xdim,ydim,length(doy_out)))
      fapar_unc_out = array(NA, dim=c(xdim,ydim,length(doy_out)))      
      for (d in seq(1, length(doy_out))) {
           fapar_out[,,d] = array(values(subset(fapar_est_out, d)), dim=c(xdim,ydim))
           fapar_unc_out[,,d] = array(values(subset(fapar_std_out, d)), dim=c(xdim,ydim))
      }

      # Check which are NA in either the estimate or the uncertainty
      filter = which(is.na(fapar_out) | is.na(fapar_unc_out))
      # and remove those which differ
      fapar_out[filter] = NA ; fapar_unc_out[filter] = NA

      # output variables
      fapar_all = list(fapar_all = fapar_out, fapar_unc_all = fapar_unc_out,
                     doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
      # clean up variables
      rm(doy_in,fapar_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
      return(fapar_all)

  } # if gridded dataset or not

} # function end
## Use byte compile
load_fapar_fields_for_extraction<-cmpfun(load_fapar_fields_for_extraction)
