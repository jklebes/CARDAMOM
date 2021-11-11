
###
## Function to load Net Biome Exchange for subsquent subsetting
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_nbe_fields_for_extraction<-function(latlon_in,nbe_source,years_to_load) {

  if (nbe_source == "GEOSCHEM") {

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

                    # re-structure to matching orientation with the lat / long information
                    var1 = var1[,dim(var1)[2]:1] ; var2 = var2[,dim(var2)[2]:1]
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
                        doy_out = doy_in
                    } else {
                        # begin populating the various outputs
                        nbe_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
                        nbe_unc_hold[1:length(as.vector(var2)),(t+((yrs-1)*nsteps))] = as.vector(var2)
                        doy_out = append(doy_out,doy_in)
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
      nbe_out = array(as.vector(nbe_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))
      nbe_unc_out = array(as.vector(nbe_unc_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))

      # output variables
      nbe_all = list(nbe_gCm2day = nbe_out, nbe_unc_gCm2day = nbe_unc_out,
                     doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
      # clean up variables
      rm(doy_in,nbe_hold,nbe_unc_hold,not_na,nbe_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
      return(nbe_all)

  } # if GEOSCHEM

} # function end load_nbe_fields_for_extraction

## Use byte compile
load_nbe_fields_for_extraction<-cmpfun(load_nbe_fields_for_extraction)
