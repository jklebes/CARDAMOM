
###
## Function to load Net Biome Exchange for subsquent subsetting
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_nbe_fields_for_extraction<-function(latlon_in,nbe_source,years_to_load) {

    if (nbe_source == "Global_Combined") {

        # let the user know this might take some time
        print("Loading Global_Combined NBE estimates for subsequent sub-setting ...")

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

          if (yr == 1) {
              # list the available files
              available_files = list.files(path_to_nbe,full.names=TRUE)
              # first check how many files we have
              for (yrr in seq(1, length(years_to_do))) {
                   if (length(which(grepl(paste("Combined_NBE_OBS_",years_to_do[yrr],sep=""),available_files))) > 0) {
                       keepers = keepers+1
                   } else {
                       missing_years = append(missing_years,years_to_do[yrr])
                   } # missing years
              } # loop years to check which we have
           } # first year

           # Determine the name of the first file first file files
           input_file_1 = paste(path_to_nbe,"/Combined_NBE_OBS_",years_to_do[yr],".nc",sep="")

           # check to see if file exists if it does then we read it in,
           # if not then we assume its a year we don't have data for and move on
           if (file.exists(input_file_1) == TRUE) {

               # open the file
               data1 = nc_open(input_file_1)

               # extract location variables
               if (lat_done == FALSE) {
                   lat = ncvar_get(data1, "lat_axis") ; long = ncvar_get(data1, "long_axis")
                   lat = array(rev(lat), dim=c(length(lat),length(long))) ; lat = t(lat)
                   long = array(long, dim=dim(lat))
               }

               # read the NBE estimate (units are gC/m2/day)
               var1 = ncvar_get(data1, "NBE")
               # get some uncertainty information - in this case the max / min which will be the basis of our uncertainty
               var2 = ncvar_get(data1, "NBE_min") ; var3 = ncvar_get(data1, "NBE_max")
               var2 = (var3 - var2) * 0.5 ; rm(var3)
               # Months in a year
               time_steps_per_year = 12
               # approximate doy of the mid-month and allocate nbe to that point
               if (lat_done == FALSE) {
                   doy_obs = floor((c(1:time_steps_per_year)*(365.25/12))-(365.25/24))
               } else {
                   doy_obs = append(doy_obs,floor((c(1:time_steps_per_year)*(365.25/12))-(365.25/24)))
               }

               # close files after use
               nc_close(data1)

               # vectorise at this time
               if (lat_done == FALSE) {
                   nbe_gCm2day = as.vector(var1)
                   nbe_unc_gCm2day = as.vector(var2)
               } else {
                   nbe_gCm2day = append(nbe_gCm2day,as.vector(var1))
                   nbe_unc_gCm2day = append(nbe_unc_gCm2day,as.vector(var2))
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
      nbe_gCm2day = array(nbe_gCm2day, dim=c(dim(long)[1],dim(lat)[2],length(doy_obs)))
      nbe_unc_gCm2day = array(nbe_unc_gCm2day, dim=c(dim(long)[1],dim(lat)[2],length(doy_obs)))

      # output variables
      return(list(nbe_gCm2day = nbe_gCm2day, nbe_unc_gCm2day = nbe_unc_gCm2day,
                  doy_obs = doy_obs, lat = lat, long = long, missing_years = missing_years))

  } else if (nbe_source == "GEOSCHEM") {

      # let the user know this might take some time
      print("Loading processed NBE fields for subsequent sub-setting ...")

      # check which file prefix we are using today
      # list all available files which we will then search
      avail_files = list.files(path_to_nbe,full.names=TRUE)
      prefix = "geoschem_nbe_" # (.)* wildcard characters for unix standard c_gls*_

      # timing information on the number of day in a month
      month_days = rep(31,length.out=12)
      month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30

      lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
      # loop for year here
      for (yr in seq(1, length(years_to_load))) {
           print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

           # first check how many files we have
           if (yr == 1) {
               nsteps = 0
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

           # open processed files
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
                    month = as.numeric(substring(tmp,1,2)) #; doy_in = as.numeric(substring(tmp,3,4))
                    doy_in = ncvar_get(data1,"day")
                    # January is correct already, so only adjust if month is >= February
                    if (month > 1) {
                        doy_in = doy_in + sum(month_days[1:(month-1)])
                    }

                    # extract location variables
                    if (lat_done == FALSE) {
                        lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
                        lat = array(rev(lat), dim=c(length(lat),length(long))) ; lat = t(lat)
                        long = array(long, dim=dim(lat))
                        # restrict the spatial extent based on latlong ranges provided
                        remove_lat = intersect(which(lat < (max(latlon_in[,1])+2)),which(lat > (min(latlon_in[,1])-2)))
                        remove_long = intersect(which(long < (max(latlon_in[,2])+2)),which(long > (min(latlon_in[,2])-2)))
                        # now find common where out in both contexts
                        remove_lat = intersect(remove_lat,remove_long)
                        # update both variables because of common matrix
                        remove_long = remove_lat
                        # adjust for matrix rather than vector arrangement
                        remove_lat = remove_lat/dim(lat)[1]
                        remove_long = (remove_long-(floor(remove_lat)*dim(lat)[1]))+1
                        remove_lat = ceiling(remove_lat)
                        # update new dimensions
                        lat_dim = length(min(remove_lat):max(remove_lat)) ; long_dim = length(min(remove_long):max(remove_long))
                        lat = lat[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
                        long = long[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
                    }

                    # read the NBE observations
                    var1 = ncvar_get(data1, "NBE") # net biome exchange of CO2 (gC/m2/day)
                    # check for error variable
                    var2 = ncvar_get(data1, "NBE_unc") # NBE error  estimate(gC/m2/day)

                    # reduce spatial cover to the desired area only
                    var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
                    var2 = var2[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
                    # set actual missing data to -9999
                    var1[which(is.na(as.vector(var1)))] = -9999
                    var2[which(is.na(as.vector(var2)))] = -9999

                    # close files after use
                    nc_close(data1)

                    # remove additional spatial information
                    if (lat_done == FALSE) {
                        # create holding arrays for the nbe information...
                        nbe_hold = array(NA, dim=c(long_dim*lat_dim,keepers*nsteps))
                        nbe_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
                        # ...and its uncertainty information...
                        nbe_unc_hold = array(NA, dim=c(long_dim*lat_dim,keepers*nsteps))
                        nbe_unc_hold[1:length(as.vector(var2)),(t+((yrs-1)*nsteps))] = as.vector(var2)
                        # ...and timing
                        doy_obs = doy_in
                    } else {
                        # begin populating the various outputs
                        nbe_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
                        nbe_unc_hold[1:length(as.vector(var2)),(t+((yrs-1)*nsteps))] = as.vector(var2)
                        doy_obs = append(doy_obs,doy_in)
                    }

                    # update flag for lat / long load
                    if (lat_done == FALSE) {lat_done = TRUE}

               } # loop through available time steps in the current year

               # keep track of years actually ran
               yrs = yrs + 1
               # clean up allocated memeory
               rm(var1,var2) ; gc()

           } # is there information for the current year?

      } # year loop

      # Sanity check for NBE
      if (lat_done == FALSE) {stop('No NBE information could be found...')}

      # remove initial value
      missing_years = missing_years[-1]

      # check which ones are NA because I made them up
      not_na = is.na(as.vector(nbe_hold))
      not_na = which(not_na == FALSE)
      filter = as.vector(nbe_hold) == -9999 | as.vector(nbe_unc_hold) == -9999
      # now remove the ones that are actual missing data
      nbe_hold[filter] = NA ; nbe_unc_hold[filter] = NA

      # enforce minimum uncertainty value
      nbe_unc_hold[abs(nbe_unc_hold) < 0.01] = 0.01

      # return spatial structure to data
      nbe_out = array(as.vector(nbe_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_obs)))
      nbe_unc_out = array(as.vector(nbe_unc_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_obs)))

      # output variables
      nbe_all = list(nbe_gCm2day = nbe_out, nbe_unc_gCm2day = nbe_unc_out,
                     doy_obs = doy_obs, lat = lat, long = long, missing_years=missing_years)
      # clean up variables
      rm(doy_in,nbe_hold,nbe_unc_hold,not_na,nbe_out,doy_obs,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
      return(nbe_all)

  } # if GEOSCHEM

} # function end load_nbe_fields_for_extraction

## Use byte compile
load_nbe_fields_for_extraction<-cmpfun(load_nbe_fields_for_extraction)
