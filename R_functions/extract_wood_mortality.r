
###
## Function extracts location specific information on the timeseries of wood stock information
## drawn from an already loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Created: 28/12/2021
# Last modified: 28/12/2021 (T. L. Smallman)

extract_wood_mortality<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,Cwood_mortality_all) {

   # Update the user
   print(paste("Cwood mortality extracted for current location ",Sys.time(),sep=""))

   # find the nearest location
   output = closest2d(1,Cwood_mortality_all$lat,Cwood_mortality_all$long,latlon_in[1],latlon_in[2],2)
   i1 = unlist(output)[1] ; j1 = unlist(output)[2]

   # Create time series output variables
   Cwood_mortality = rep(-9999, length(timestep_days))
   Cwood_mortality_unc = rep(-9999, length(timestep_days))
   Cwood_mortality_lag = rep(-9999, length(timestep_days))

   # Loop through each time step of the Cwood increment / production timeseries,
   # its associated uncertainty and period of effect
   for (t in seq(1, length(Cwood_mortality_all$place_obs_in_step))) {
        # Prodictivity estimate
        Cwood_mortality[Cwood_mortality_all$place_obs_in_step[t]] = Cwood_mortality_all$Cwood_mortality_gCm2[i1,j1,t]
        # Its uncertainty
        tmp = min(Cwood_mortality[Cwood_mortality_all$place_obs_in_step[t]], Cwood_mortality_all$Cwood_mortality_uncertainty_gCm2[i1,j1,t])
        Cwood_mortality_unc[Cwood_mortality_all$place_obs_in_step[t]] = tmp
        # Its period of effect
        Cwood_mortality_lag[Cwood_mortality_all$place_obs_in_step[t]] = Cwood_mortality_all$Cwood_mortality_lag_step[i1,j1,t]
   }

   # Set any time series values with NaN to missing data flag (-9999)
   Cwood_mortality[which(is.na(Cwood_mortality))] = -9999
   Cwood_mortality_unc[which(is.na(Cwood_mortality_unc))] = -9999
   Cwood_mortality_lag[which(is.na(Cwood_mortality_lag))] = -9999

   # pass the information back
   return(list(Cwood_mortality = Cwood_mortality, Cwood_mortality_unc = Cwood_mortality_unc, Cwood_mortality_lag = Cwood_mortality_lag))

} # end function extract_wood_productivity

## Use byte compile
extract_wood_mortality<-cmpfun(extract_wood_mortality)