
###
## Extracts location specific information on the leaf carbon per leaf area (gCm2) from a loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_lca_prior<- function(i1,j1,spatial_type,resolution,grid_type,latlon_in,lca_all) {

   # Update the user
   print(paste("LCA prior extracted for current location ",Sys.time(),sep=""))

#   # find the nearest location
#   output = closest2d_2(1,lca_all$lat,lca_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract target location
   lca_gCm2 = lca_all$lca_gCm2[i1,j1]
   lca_unc_gCm2 = lca_all$lca_uncertainty_gCm2[i1,j1]

   # Convert any NaN to missing data flag -9999
   lca_gCm2[is.na(lca_gCm2)] = -9999 ; lca_unc_gCm2[is.na(lca_unc_gCm2)] = -9999

   # pass the information back
   return(list(lca_gCm2 = lca_gCm2, lca_unc_gCm2 = lca_unc_gCm2))

} # end function extract_lca_prior

## Use byte compile
extract_lca_prior<-cmpfun(extract_lca_prior)
