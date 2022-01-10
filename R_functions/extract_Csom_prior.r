
###
## Function extracts a location specific estimate of the initial Csom value from the gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Csom_prior<-function(spatial_type,resolution,grid_type,latlon_wanted,Csom_all) {

  # extract information on soil C content from soil grids database

  # Update the user
  print(paste("Csom prior data extracted for current location ",Sys.time(),sep=""))

	# convert input data long to conform to what we need
	check1=which(Csom_all$long > 180) ; if (length(check1) > 0) { Csom_all$long[check1] = Csom_all$long[check1]-360 }

  # find desired lat / long location within the soilgrid database
  output = closest2d(1,Csom_all$lat,Csom_all$long,latlon_wanted[1],latlon_wanted[2],2)
  i1 = unlist(output)[1] ; j1 = unlist(output)[2]

	# return long to 0-360
	if (length(check1) > 0) { Csom_all$long[check1] = Csom_all$long[check1]+360 }

  # Extract target location
  Csom = Csom_all$Csom[i1,j1]
  Csom_unc = Csom_all$Csom_unc[i1,j1]

  # Convert any NaN to missing data flag -9999
  Csom[which(is.na(Csom))] = -9999 ; Csom_unc[which(is.na(Csom_unc))] = -9999

  # retun back to the user
  return(list(Csom_initial = Csom, Csom_initial_unc = Csom_unc))

} # end function extract_Csom_prior

## Use byte compile
extract_Csom_prior<-cmpfun(extract_Csom_prior)
