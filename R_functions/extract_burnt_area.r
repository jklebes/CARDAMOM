
###
## Function extracts locatio specific information on Burned area
## from already loaded gridded datasets
###

# This function is based by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_burnt_area_information<- function(latlon_in,timestep_days,spatial_type,grid_type,resolution,start_year,end_year,burnt_all) {

  # Update the user
  print(paste("Beginning burned fraction data extraction for current location ",Sys.time(),sep=""))

  # find the nearest location
  if (length(dim(burnt_all$long)) == 2) {
      output = closest2d(1,burnt_all$lat,burnt_all$long,latlon_in[1],latlon_in[2],2)
  } else {
      output = closest2d(1,burnt_all$lat,burnt_all$long,latlon_in[1],latlon_in[2],3)
  }
  #i1=unlist(output)[2] ; j1=unlist(output)[1]
  i1 = unlist(output)[1] ; j1=unlist(output)[2]

  # Extract location specific information to a local variable
  burnt_area = burnt_all$burnt_area[i1,j1,]
  # convert missing data to -9999
  burnt_area[which(is.na(burnt_area))] = -9999

  # next work out how many days we should have in the year
  doy_out=0
  for (i in seq(1, length(years_to_do))) {
       # is current year a leap or not
       nos_days = 365
       mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/4))*4
       if (mod == 0) {
           nos_days = 366
           mod = as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/100))*100
           if (mod == 0) {
               nos_days  = 365
               mod = as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/400))*400
               if (mod == 0) {
                   nos_days  = 366
               }
            }
       }
       # count up days needed
       doy_out = append(doy_out,1:nos_days)
  }
  # remove initial value
  doy_out = doy_out[-1]

  # just incase there is no missing data we best make sure there is a value which can be assessed
  if (length(burnt_all$missing_years) == 0) { burnt_all$missing_years=1066 }

  # declare output variable
  burnt_area_out = array(0, dim=length(doy_out))
  # now line up the obs days with all days
  b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_do[1])
  while (b <= length(burnt_all$doy_obs)) {
         # if we are in a year which is missing then we do not allow consideration of DOY
         if (start_year != burnt_all$missing_years[a]) {
             if (doy_out[i] == burnt_all$doy_obs[b]) {
                 burnt_area_out[i] = burnt_area[b] ; b = b+1
             } # end if doy matches
         } # end if missing year
         # but we do keep counting through the total vector length which we expect
         i = i + 1
         # have we just looped round the year?
         if (i != 1 & doy_out[i-1] > doy_out[i]) {
             # and if we have just been in a missing year we need to count on the missing years vector to
             if (start_year == burnt_all$missing_years[a]) {a = min(length(burnt_all$missing_years),a+1)}
             start_year = start_year+1
         } # end if doy_out[i] == 1
  } # end while condition

  if (length(timestep_days) == 1 & timestep_days[1] == 1) {
      # well actually we do nothing
  } else {
      # generally this now deals with time steps which are not daily.
      # However if not monthly special case
      if (length(timestep_days) == 1) {
          run_day_selector = seq(1,length(burnt_area_out),timestep_days)
          timestep_days = rep(timestep_days, length.out=length(burnt_area_out))
      }
      print("...calculating monthly averages for burnt area")
      # determine the actual daily positions
      run_day_selector = cumsum(timestep_days)
      # create needed variables
      burnt_area_agg = array(-9999,dim=length(run_day_selector))
      # As far is an absolute event we sum all events within the time period not average!
      for (y in seq(1,length(run_day_selector))) {
           pick = burnt_area_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
           burnt_area_agg[y] = sum(pick[which(pick != -9999)],na.rm=TRUE)
      }
      # now convert the missing values back to -9999
      burnt_area_agg[which(is.na(burnt_area_agg))] = -9999
      # update with new output information
      burnt_area_out = burnt_area_agg
      # clean up
      rm(burnt_area_agg) ; gc()
  } # monthly aggregation etc

  # pass the information back
  return(burnt_area=burnt_area_out)

} # end function extract_burnt_area_information

## Use byte compile
extract_burnt_area_information<-cmpfun(extract_burnt_area_information)
