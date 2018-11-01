####
## Function to extract the HWSD soil carbon estimate for initial conditions
####


extract_soilwater_initial<- function(spatial_type,resolution,grid_type,latlon_in,soilwater_all) {

	# find the nearest location
	output=closest2d(1,soilwater_all$lat,soilwater_all$long,latlon_in[1],latlon_in[2],2)
	j1=unlist(output)[2];i1=unlist(output)[1]
	print(paste("Initial soil water fraction extracted for current location ",Sys.time(),sep=""))

	# work out number of pixels to average over
	if (spatial_type == "grid") {
		if (grid_type == "wgs84") {
			# resolution of the product
			product_res = abs(soilwater_all$lat[1,2]-soilwater_all$lat[1,1])+abs(soilwater_all$long[2,1]-soilwater_all$long[1,1])
			product_res = abs(product_res * 0.5)
			# radius is ceiling of the ratio of the product vs analysis ratio
			radius = ceiling(resolution / product_res)
		} else if (grid_type == "UK") {
			radius = max(0,floor(1*resolution*1e-3*0.5))
		} else {
			stop("have not specified the grid used in this analysis")
		}
	} else {
		radius = 0
	}

	answer=NA
	while (is.na(answer) == TRUE) {
		# work out average areas
		average_i=(i1-radius):(i1+radius) ; average_j=(j1-radius):(j1+radius)
		average_i=max(1,(i1-radius)):min(dim(soilwater_all$soil_water)[1],(i1+radius)) ; average_j=max(1,(j1-radius)):min(dim(soilwater_all$soil_water)[2],(j1+radius))
		# carry out averaging
		tmp=soilwater_all$soil_water[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
		soilwater=mean(tmp, na.rm=TRUE)
		tmp=soilwater_all$soil_water_unc[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
		soilwater_unc=mean(tmp, na.rm=TRUE)
		# error checking
		if (is.na(soilwater) | soilwater == 0) {radius=radius+1 ; answer=NA} else {answer=0}
	}
	print(paste("NOTE: Initial soil water averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))
	# pass the information back
	return(list(soil_water = soilwater,soil_water_unc = soilwater_unc))

} # end of function
