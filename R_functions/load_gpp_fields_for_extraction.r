
###
## Function to load GPP from various gridded datasets
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_gpp_fields_for_extraction<-function(latlon_in,gpp_source,start_year,end_year) {

    # Determine the years for which the analysis will occur
    years_to_do = as.character(seq(start_year,end_year))

    if (gpp_source == "Global_Combined") {

        # let the user know this might take some time
        print("Loading Global_Combined GPP estimates for subsequent sub-setting ...")

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

             if (yr == 1) {
                 # list the available files
                 available_files = list.files(path_to_gpp,full.names=TRUE)
                 # first check how many files we have
                 for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(paste("Combined_GPP_OBS_",years_to_do[yrr],sep=""),available_files))) > 0) {
                          keepers = keepers+1
                      } else {
                          missing_years = append(missing_years,years_to_do[yrr])
                      } # missing years
                 } # loop years to check which we have
	            } # first year

              # Determine the name of the first file first file files
  	          input_file_1 = paste(path_to_gpp,"/Combined_GPP_OBS_",years_to_do[yr],".nc",sep="")

              # check to see if file exists if it does then we read it in,
              # if not then we assume its a year we don't have data for and move on
	            if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  if (lat_done == FALSE) {
                      lat = ncvar_get(data1, "lat_axis") ; long = ncvar_get(data1, "long_axis")
                  }

                  # read the burnt fraction estimate (units are 0-1)
                  var1 = ncvar_get(data1, "GPP")
                  # get some uncertainty information - in this case the max / min which will be the basis of our uncertainty
                  var2 = ncvar_get(data1, "GPP_min") ; var3 = ncvar_get(data1, "GPP_max")
                  var2 = (var3 - var2) * 0.5 ; rm(var3)
                  # Months in a year
                  time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate gpp to that point
                  if (lat_done == FALSE) {
                      doy_obs = floor((c(1:time_steps_per_year)*(365.25/12))-(365.25/24))
                  } else {
                      doy_obs = append(doy_obs,floor((c(1:time_steps_per_year)*(365.25/12))-(365.25/24)))
                  }

                  # close files after use
                  nc_close(data1)

                  # vectorise at this time
                  if (lat_done == FALSE) {
                      gpp_gCm2day = as.vector(var1)
                      gpp_unc_gCm2day = as.vector(var2)
                  } else {
                      gpp_gCm2day = append(gpp_gCm2day,as.vector(var1))
                      gpp_unc_gCm2day = append(gpp_unc_gCm2day,as.vector(var2))
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
         gpp_gCm2day = array(gpp_gCm2day, dim=c(length(long),length(lat),length(doy_obs)))
         gpp_unc_gCm2day = array(gpp_unc_gCm2day, dim=c(length(long),length(lat),length(doy_obs)))

         # output variables
         return(list(gpp_gCm2day = gpp_gCm2day, gpp_unc_gCm2day = gpp_unc_gCm2day,
                     doy_obs = doy_obs, lat = lat, long = long, missing_years = missing_years))

    } else if (gpp_source == " " | gpp_source == "site_specific"){

	        # Do nothing as this should be read directly from files or not needed

    } else {

	        stop(paste("GPP option (",gpp_source,") not valid"))

    } # if GPP type

} # function end

## Use byte compile
load_gpp_fields_for_extraction<-cmpfun(load_gpp_fields_for_extraction)
