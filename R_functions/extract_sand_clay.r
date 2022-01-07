
###
## Function to extract location specific information on soil texture from the gridded SoilGrids database
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_sand_clay<- function(spatial_type,resolution,grid_type,latlon_in,sand_clay_all) {

  # Update the user
	print(paste("Sand/clay data extracted for current location ",Sys.time(),sep=""))

	# convert input data long to conform to what we need
	check1=which(sand_clay_all$long > 180) ; if (length(check1) > 0) { sand_clay_all$long[check1]=sand_clay_all$long[check1]-360 }

	# find the nearest location
	output=closest2d(1,sand_clay_all$lat,sand_clay_all$long,latlon_in[1],latlon_in[2],2)
	j1=unlist(output)[2];i1=unlist(output)[1]

	# return long to 0-360
	if (length(check1) > 0) { sand_clay_all$long[check1]=sand_clay_all$long[check1]+360 }


  # Extract the correct value, but allow for expanding to a larger area if we pick a no data area
  radius = c(0,0) # assume precise location is known
  answer = NA
	while (is.na(answer) == TRUE) {
	    # work out average areas
	    average_i = (i1-radius[1]):(i1+radius[1]) ; average_j = (j1-radius[2]):(j1+radius[2])
	    average_i = max(1,(i1-radius[1])):min(dim(sand_clay_all$top_sand)[1],(i1+radius[1]))
      average_j = max(1,(j1-radius[2])):min(dim(sand_clay_all$top_sand)[2],(j1+radius[2]))
	    # carry out averaging
	    tmp1 = sand_clay_all$top_sand[average_i,average_j] ; tmp1[which(tmp1 == -9999)] = NA
	    tmp2 = sand_clay_all$bot_sand[average_i,average_j] ; tmp2[which(tmp2 == -9999)] = NA
	    tmp3 = sand_clay_all$top_clay[average_i,average_j] ; tmp3[which(tmp3 == -9999)] = NA
	    tmp4 = sand_clay_all$bot_clay[average_i,average_j] ; tmp4[which(tmp4 == -9999)] = NA
	    top_sand = mean(tmp1, na.rm=TRUE) ; bot_sand=mean(tmp2, na.rm=TRUE)
	    top_clay = mean(tmp3, na.rm=TRUE) ; bot_clay=mean(tmp4, na.rm=TRUE)
	    # error checking
	    if (is.na(top_sand) | top_sand <= 0) {radius = radius+1 ; answer = NA} else {answer = 0}
	}
  # Inform the user
	print(paste("NOTE sand/clay averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))

	# just to check because when averaging sometimes the sand / clay combinations can be > 100 %
	# 94 % chosesn as this is the highest total % found in the HWSD dataset
	if ((top_sand+top_clay) > 94) {
	    tmp1 = top_sand / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    tmp2 = top_clay / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    top_sand = tmp1*100 ; top_clay = tmp2*100
	}
	if ((bot_sand+bot_clay) > 94) {
	    tmp1 = bot_sand / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    tmp2 = bot_clay / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    top_sand = tmp1*100 ; top_clay = tmp2*100
	}

	# pass the information back
	return(list(top_sand=top_sand,bot_sand=bot_sand,top_clay=top_clay,bot_clay=bot_clay))

} # end function extract_sand_clay

## Use byte compile
extract_sand_clay<-cmpfun(extract_sand_clay)
