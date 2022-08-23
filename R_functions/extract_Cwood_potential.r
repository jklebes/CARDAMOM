

###
## Function extracts information on the potential steady state value of the wood stock pool
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Cwood_potential<- function(i1,j1,timestep_days,spatial_type,resolution,
                                   grid_type,latlon_in,Cwood_potential_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("Cwood potential extracted for current location ",Sys.time(),sep=""))}

#   # find the nearest location
#   output = closest2d_2(1,Cwood_potential_all$lat,Cwood_potential_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract to local variable
   Cwood = Cwood_potential_all$biomass_gCm2[i1,j1]
   Cwood_unc = Cwood_potential_all$biomass_uncertainty_gCm2[i1,j1]

   # Convert any NaN to missing data flag -9999
   Cwood[which(is.na(Cwood))] = -9999 ; Cwood_unc[which(is.na(Cwood_unc))] = -9999

   # pass the information back
   return(list(Cwood_stock = Cwood, Cwood_stock_unc = Cwood_unc))

} # end function extract_Cwood_potential

## Use byte compile
extract_Cwood_potential<-cmpfun(extract_Cwood_potential)
