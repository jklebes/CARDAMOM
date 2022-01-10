
###
## Function extracts leaf area index (LAI) information from a pre-loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_lai_timeseries<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,lai_all,years_to_load) {

   # Update the user
   print(paste("LAI data extracted for current location ",Sys.time(),sep=""))

   # convert input data long to conform to what we need
   check1 = which(lai_all$long > 180) ; if (length(check1) > 0) { lai_all$long[check1]=lai_all$long[check1]-360 }
   # find the nearest location
   output = closest2d(1,lai_all$lat,lai_all$long,latlon_in[1],latlon_in[2],2)
   i1 = unlist(output, use.names=FALSE)[1] ; j1=unlist(output, use.names=FALSE)[2]

   # return long to 0-360
   if (length(check1) > 0) { lai_all$long[check1]=lai_all$long[check1]+360 }

   # Extract current location to local variable
   lai = lai_all$lai_all[i1,j1,]
   lai_unc = lai_all$lai_unc_all[i1,j1,]

   # convert missing data to -9999
   lai[which(is.na(lai))] = -9999.0 ; lai_unc[which(is.na(lai_unc))] = -9999.0

   # Determine how many days are in each year
   doy_out = 0
   for (i in seq(1, length(years_to_load))) {
        nos_days = nos_days_in_year(years_to_load[i])
        # count up days needed
        doy_out = append(doy_out,1:nos_days)
   }
   doy_out = doy_out[-1]

   # Just incase there is no missing data we best make sure there is a value which can be assessed
   if (length(lai_all$missing_years) == 0) { lai_all$missing_years=1066 }

   # declare output variable
   lai_out = array(-9999, dim=length(doy_out))
   lai_unc_out = array(-9999, dim=length(doy_out))
   # now line up the obs days with all days
   b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
   while (b <= length(lai_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != lai_all$missing_years[a]) {
          if (doy_out[i] == lai_all$doy_obs[b]) {
              lai_out[i] = lai[b] ; lai_unc_out[i] = lai_unc[b] ; b = b + 1
          } # end if doy matches
      } # end if missing year

      # but we do keep counting through the total vector length which we expect
      i = i + 1

      # each time we come back to doy_out[i]==1 we need to count on the year
      if (doy_out[i] == 1 & b <= length(lai_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == lai_all$missing_years[a]) { a = min(length(lai_all$missing_years),a+1) }
          start_year = start_year + 1
      } # end if doy_out[i] == 1

   } # end while condition

   if (length(timestep_days) == 1 & timestep_days[1] == 1) {

       # well actually we do nothing

   } else {
       # generally this now deals with time steps which are not daily.
       # However if not monthly special case
       if (length(timestep_days) == 1) {
           run_day_selector=seq(1,length(lai_out),timestep_days)
           timestep_days=rep(timestep_days, length.out=length(lai_out))
       }
       print("...calculating monthly averages for lai")
       # determine the actual daily positions
       run_day_selector=cumsum(timestep_days)
       # create needed variables
       lai_agg = array(NA,dim=length(run_day_selector))
       lai_unc_agg = array(NA,dim=length(run_day_selector))
       # Convert -9999 into NaN
       lai_out[which(lai_out == -9999)] = NA
       lai_unc_out[which(lai_unc_out == -9999)] = NA
       for (y in seq(1,length(run_day_selector))) {
            pick = (run_day_selector[y]-timestep_days[y]):run_day_selector[y]
            lai_agg[y] = mean(lai_out[pick], na.rm=TRUE)
            lai_unc_agg[y] = mean(lai_unc_out[pick], na.rm=TRUE)
       }
       # convert missing values to -9999
       lai_agg[which(is.na(lai_agg))] = -9999
       lai_unc_agg[which(is.na(lai_unc_agg))] = -9999
       # update with new output information
       lai_out = lai_agg ; lai_unc_out = lai_unc_agg
       # clean up
       rm(lai_agg,lai_unc_agg,y) ; gc()

   } # monthly aggregation etc

   # clean up
   rm(i1,j1,check1,lai,i,nos_days,doy_out,a) ; gc(reset=TRUE,verbose=FALSE)

   # CARDAMOM works best if the uncertainties are the same across each LAI observation as the framework tends towards lower LAI values
   # Therefore, to make use of the uncertainty information we take the mean for this site and apply it across each value.
   # NOTE: we put a book end the upper uncertainty linked to half the mean LAI estimate to ensure that there is some constraint
   lai_unc_out[lai_out >= 0] = max(0.25,min(mean(lai_unc_out[lai_out >= 0]), 0.5*mean(lai_out[lai_out >= 0])))

   # pass the information back
   output = list(lai = lai_out, lai_unc = lai_unc_out)
   return(output)

} # end function extract_lai_timeseries

## Use byte compile
extract_lai_timeseries<-cmpfun(extract_lai_timeseries)
