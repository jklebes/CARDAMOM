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
# Function used to extract location specific information on leaf area index
# from already loaded gridded datasets
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_lai_timeseries<- function(i1,j1,timestep_days,spatial_type,resolution,
                                  grid_type,latlon_in,lai_all,years_to_load,doy_obs) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("LAI data extracted for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   lai = lai_all$lai_m2m2[i1,j1,]
   lai_unc = lai_all$lai_unc_m2m2[i1,j1,]

   # Just incase there is no missing data we best make sure there is a value which can be assessed
   if (length(lai_all$missing_years) == 0) { lai_all$missing_years=1066 }

   # declare output variable
   lai_out = array(NA, dim=length(doy_obs))
   lai_unc_out = array(NA, dim=length(doy_obs))
   # now line up the obs days with all days
   b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
   #print("...begin inserting LAI observations into model time steps")
   while (b <= length(lai_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != lai_all$missing_years[a]) {
          if (doy_obs[i] == lai_all$doy_obs[b]) {
              lai_out[i] = lai[b] ; lai_unc_out[i] = lai_unc[b] ; b = b + 1
          } # end if doy matches
      } # end if missing year

      # but we do keep counting through the total vector length which we expect
      i = i + 1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      if (doy_obs[i] == 1 & b <= length(lai_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == lai_all$missing_years[a]) { a = min(length(lai_all$missing_years),a+1) }
          start_year = start_year + 1
      } # end if doy_obs[i] == 1

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
       #print("...calculating monthly averages for lai")
       # determine the actual daily positions
       run_day_selector=cumsum(timestep_days)
       # create needed variables
       lai_agg = array(NA,dim=length(run_day_selector))
       lai_unc_agg = array(NA,dim=length(run_day_selector))
       # Loop through
       for (y in seq(1,length(run_day_selector))) {
            pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
            lai_agg[y] = mean(lai_out[pick], na.rm=TRUE)
            lai_unc_agg[y] = mean(lai_unc_out[pick], na.rm=TRUE)
       }
       # update with new output information
       lai_out = lai_agg ; lai_unc_out = lai_unc_agg
       # clean up
       rm(lai_agg,lai_unc_agg,y) ; gc()

   } # monthly aggregation etc

   # convert missing data to -9999
   lai_out[which(is.na(lai_out))] = -9999 ; lai_unc_out[which(is.na(lai_unc_out))] = -9999

   # clean up
   rm(i1,j1,lai,i,a) ; gc(reset=TRUE,verbose=FALSE)

   # CARDAMOM works best if the uncertainties are the same across each LAI observation as the framework tends towards lower LAI values
   # Therefore, to make use of the uncertainty information we take the mean for this site and apply it across each value.
   # NOTE: we put a book end the upper uncertainty linked to the max LAI estimate to ensure that there is some constraint
   #lai_unc_out[lai_out >= 0] = max(0.25,min(mean(lai_unc_out[lai_out >= 0]), max(lai_out[lai_out >= 0])))
   #lai_unc_out[lai_out >= 0] = pmax(0.25,lai_unc_out[lai_out >= 0])

   # pass the information back
   output = list(lai = lai_out, lai_unc = lai_unc_out)
   return(output)

} # end function extract_lai_timeseries

## Use byte compile
extract_lai_timeseries<-cmpfun(extract_lai_timeseries)
