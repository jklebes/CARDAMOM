
###
## Function extracts a location specific estimate of the initial Csom value from the gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Csom_prior<-function(i1,j1,spatial_type,resolution,grid_type,
                             latlon_wanted,Csom_all) {

  # extract information on soil C content from soil grids database

  # Update the user
  print(paste("Csom prior data extracted for current location ",Sys.time(),sep=""))

#  # find desired lat / long location within the soilgrid database
#  output = closest2d_2(1,Csom_all$lat,Csom_all$long,latlon_wanted[1],latlon_wanted[2])
#  i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

  # Extract target location
  Csom = Csom_all$Csom[i1,j1]
  Csom_unc = Csom_all$Csom_unc[i1,j1]

  # Convert any NaN to missing data flag -9999
  Csom[is.na(Csom) == TRUE] = -9999 ; Csom_unc[is.na(Csom_unc) == TRUE] = -9999

  # retun back to the user
  return(list(Csom_initial = Csom, Csom_initial_unc = Csom_unc))

} # end function extract_Csom_prior

## Use byte compile
extract_Csom_prior<-cmpfun(extract_Csom_prior)
