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
# Function used to extract location specific information on gross primary productivity
# from already loaded gridded datasets
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_gpp<- function(i1,j1,timestep_days,spatial_type,resolution,grid_type,
                       latlon_in,gpp_all,years_to_load,doy_obs) {

  # Update the user
  if (use_parallel == FALSE) {print(paste("GPP data extracted for current location ",Sys.time(),sep=""))}

#  # find the nearest location
#  output = closest2d_2(1,gpp_all$lat,gpp_all$long,latlon_in[1],latlon_in[2])
#  i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

  # Extract to local variable
  gpp = gpp_all$gpp_gCm2day[i1,j1,]
  gpp_unc = gpp_all$gpp_unc_gCm2day[i1,j1,]

  # just incase there is no missing data we best make sure there is a value which can be assessed
  if (length(gpp_all$missing_years) == 0) { gpp_all$missing_years=1066 }

  # declare output variable
  gpp_out = array(NA, dim=length(doy_obs))
  gpp_unc_out = array(NA, dim=length(doy_obs))
  # now line up the obs days with all days
  b = 1 ; i = 1 ; a = 1 ; start_year=as.numeric(years_to_load[1])
  while (b <= length(gpp_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != gpp_all$missing_years[a]) {
          if (doy_obs[i] == gpp_all$doy_obs[b]) {
              gpp_out[i] = gpp[b] ; gpp_unc_out[i] = gpp_unc[b] ; b = b+1
          } # end if doy matches
      } # end if missing year
      # but we do keep counting through the total vector length which we expect
      i = i+1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      if (doy_obs[i-1] > doy_obs[i] & b <= length(gpp_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == gpp_all$missing_years[a]) {a = min(length(gpp_all$missing_years),a+1)}
          start_year=start_year+1
      } # end if doy_obs[i] == 1

  } # end while condition

  if (length(timestep_days) == 1 & timestep_days[1] == 1) {

      # well actually we do nothing

  } else {

      # Generally this now deals with time steps which are not daily.
      # However if not monthly special case
      if (length(timestep_days) == 1) {
          run_day_selector = seq(1,length(gpp_out),timestep_days)
          timestep_days = rep(timestep_days, length.out=length(gpp_out))
      }
      #print("...calculating monthly averages for GPP")
      # determine the actual daily positions
      run_day_selector = cumsum(timestep_days)
      # create needed variables
      gpp_agg = array(NA,dim=length(run_day_selector))
      gpp_unc_agg = array(NA, dim=length(run_day_selector))
      for (y in seq(1,length(run_day_selector))) {
           pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
           gpp_agg[y] = mean(gpp_out[pick],na.rm=TRUE)
           gpp_unc_agg[y] = mean(gpp_unc_out[pick],na.rm=TRUE)
      }
      # update with new output information
      gpp_out = gpp_agg ; gpp_unc_out = gpp_unc_agg
      # clean up
      rm(gpp_agg,gpp_unc_agg,y) ; gc()

  } # monthly aggregation etc

  # convert missing data back to -9999
  gpp_out[which(is.na(gpp_out))] = -9999
  gpp_unc_out[which(is.na(gpp_unc_out))] = -9999

  # pass the information back
  output = list(GPP = gpp_out, GPP_unc = gpp_unc_out)
  return(output)

} # end function extract_gpp

## Use byte compile
extract_gpp<-cmpfun(extract_gpp)
