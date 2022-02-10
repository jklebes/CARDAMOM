
###
## Extracts location specific information on the initial wood biomass stock from a loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Cwood_initial<- function(spatial_type,resolution,grid_type,latlon_in,Cwood_initial_all) {

   # Update the user
   print(paste("Cwood initial extracted for current location ",Sys.time(),sep=""))

   # find the nearest location
   output = closest2d_2(1,Cwood_initial_all$lat,Cwood_initial_all$long,latlon_in[1],latlon_in[2])
   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Carry out averaging
   Cwood = Cwood_initial_all$biomass_gCm2[i1,j1]
   Cwood_unc = Cwood_initial_all$biomass_uncertainty_gCm2[i1,j1]

   # Convert any NaN to missing data flag -9999
   Cwood[which(is.na(Cwood))] = -9999 ; Cwood_unc[which(is.na(Cwood_unc))] = -9999

   # pass the information back
   return(list(Cwood_stock = Cwood, Cwood_stock_unc = Cwood_unc))

} # end function extract_Cwood_initial

## Use byte compile
extract_Cwood_initial<-cmpfun(extract_Cwood_initial)
