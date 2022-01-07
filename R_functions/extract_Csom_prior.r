
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
  # If resolution has been provides as single value then adjust this here
  if (length(resolution) == 1 & spatial_type == "grid") {tmp_res = resolution * c(1,1)} else {tmp_res = resolution}

  # Extract the correct value, but allow for expanding to a larger area if we pick a no data area
  radius = c(0,0) # assume precise location is known
  answer = NA
  while (is.na(answer) == TRUE) {
    # work out average areas
    average_i = (i1-radius[1]):(i1+radius[1]) ; average_j = (j1-radius[2]):(j1+radius[2])
    average_i = max(1,(i1-radius[1])):min(dim(Csom_all$Csom)[1],(i1+radius[1]))
    average_j = max(1,(j1-radius[2])):min(dim(Csom_all$Csom)[2],(j1+radius[2]))
    # carry out averaging
    Csom = mean(as.vector(Csom_all$Csom[average_i,average_j]), na.rm=TRUE)
    Csom_unc = mean(as.vector(Csom_all$Csom_unc[average_i,average_j]), na.rm=TRUE)
    # check for errors
    if (is.na(Csom) | Csom < 0) {radius = radius+1 ; answer = NA} else {answer = 1}
  }

  # retun back to the user
  return(list(Csom_initial = Csom, Csom_initial_unc = Csom_unc))

} # end function extract_Csom_prior

## Use byte compile
extract_Csom_prior<-cmpfun(extract_Csom_prior)
