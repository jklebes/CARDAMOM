# declare function
bintodec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), ""), use.names=FALSE) == 1))-1))
}

extract_copernicus_lai<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,lai_all,years_to_load) {

  # convert input data long to conform to what we need
  check1 = which(lai_all$long > 180) ; if (length(check1) > 0) { lai_all$long[check1]=lai_all$long[check1]-360 }
  # find the nearest location
  output = closest2d(1,lai_all$lat,lai_all$long,latlon_in[1],latlon_in[2],2)
  i1 = unlist(output, use.names=FALSE)[1] ; j1=unlist(output, use.names=FALSE)[2]
  print(paste("LAI data extracted for current location ",Sys.time(),sep=""))

  # return long to 0-360
  if (length(check1) > 0) { lai_all$long[check1]=lai_all$long[check1]+360 }

  # work out number of pixels to average over
  if (spatial_type == "grid") {
      if (grid_type == "wgs84") {
          # resolution of the product
          product_res = abs(lai_all$lat[1,2]-lai_all$lat[1,1])+abs(lai_all$long[2,1]-lai_all$long[1,1])
          product_res = (product_res * 0.5)
          # radius is ceiling of the ratio of the product vs analysis ratio
          radius = ceiling(resolution / product_res)
          max_radius = radius+4
      } else if (grid_type == "UK") {
          radius = max(0,floor(1*resolution*1e-3*0.5))
          max_radius = radius+4
      } else {
          stop("have not specified the grid used in this analysis")
      }
  } else {
      radius = 0
      max_radius = 4
  }

  answer = NA
  while (is.na(answer) == TRUE) {
    # work out average areas
    average_i = max(1,(i1-radius)):min(dim(lai_all$lai_all)[1],(i1+radius))
    average_j = max(1,(j1-radius)):min(dim(lai_all$lai_all)[2],(j1+radius))
    # carry out averaging
    lai = array(NA, dim=c(dim(lai_all$lai_all)[3])) ; lai_unc = array(NA, dim=c(dim(lai_all$lai_all)[3]))
    for (n in seq(1, dim(lai_all$lai_all)[3])) {
         lai[n] = mean(as.vector(lai_all$lai_all[average_i,average_j,n]), na.rm=TRUE)
         lai_unc[n] = mean(as.vector(lai_all$lai_unc_all[average_i,average_j,n]), na.rm=TRUE)
    }
    # are any of the my data points now filled
    answer = max(as.vector(lai),na.rm=TRUE)
    # what proportion of my data points are within a realistic range
    npoints = length(which(lai > 0 & lai < 12))/length(lai)
    # do I have at least 20 % of data points filled
    if (is.na(answer) | answer == -Inf | npoints < 0.2) {radius = radius+1 ; answer = NA}
    # restrict how far to look before giving up
    if (radius >= max_radius) {answer = 1}
  }
  # warning to the used
  print(paste("NOTE LAI averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))
  # convert missing data back to -9999
  lai[which(is.na(lai))]=-9999.0 ; lai_unc[which(is.na(lai_unc))]=-9999.0
  # next work out how many days we should have in the year
  doy_out=0
  for (i in seq(1, length(years_to_load))) {
       # is current year a leap or not
       nos_days = 365
       mod=as.numeric(years_to_load[i])-round((as.numeric(years_to_load[i])/4))*4
       if (mod == 0) {
           nos_days = 366
           mod=as.numeric(years_to_load[i])-round((as.numeric(years_to_load[i])/100))*100
           if (mod == 0) {
               nos_days  = 365
               mod=as.numeric(years_to_load[i])-round((as.numeric(years_to_load[i])/400))*400
               if (mod == 0) {
                   nos_days  = 366
               }
           }
       }
       # count up days needed
       doy_out = append(doy_out,1:nos_days)
  }
  doy_out = doy_out[-1]

  # just incase there is no missing data we best make sure there is a value which can be assessed
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
    if (doy_out[i] == 1) {
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
    for (y in seq(1,length(run_day_selector))) {
      pick = lai_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
      pick_unc = lai_unc_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
      lai_agg[y] = mean(pick[which(pick != -9999)],na.rm=TRUE)
      lai_unc_agg[y] = mean(pick_unc[which(pick_unc != -9999)],na.rm=TRUE)
    }
    # convert missing values back to -9999
    lai_agg[which(is.na(lai_agg))]=-9999
    lai_unc_agg[which(is.na(lai_unc_agg))]=-9999
    # update with new output information
    lai_out=lai_agg ; lai_unc_out=lai_unc_agg
    # clean up
    rm(lai_agg,lai_unc_agg,y) ; gc()

  } # monthly aggregation etc

  # clean up
  rm(i1,j1,check1,lai,i,nos_days,doy_out,radius,n,a) ; gc(reset=TRUE,verbose=FALSE)

  # CARDAMOM works best if the uncertainties are the same across each LAI observation as the framework tends towards lower LAI values
  # Therefore, to make use of the uncertainty information we take the mean for this site and apply it across each value
  lai_unc_out[which(lai_out > 0)] = mean(lai_unc_out[which(lai_out > 0)])

  # pass the information back
  output = list(lai = lai_out, lai_unc = lai_unc_out)
  return(output)

} # end of function
