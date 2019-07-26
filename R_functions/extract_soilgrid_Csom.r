
extract_soilgrid_Csom<-function(spatial_type,resolution,grid_type,latlon_wanted,Csom_all) {

  # extract information on soil C content from soil grids database

  print(paste("SoilGrids Csom data extracted for current location ",Sys.time(),sep=""))
 
  # find desired lat / long location within the soilgrid database
  output = closest2d(1,Csom_all$lat,Csom_all$long,latlon_wanted[1],latlon_wanted[2],2)
  i1 = unlist(output)[1] ; j1 = unlist(output)[2]

  # work out number of pixels to average over
  if (spatial_type == "grid") {
      if (grid_type == "wgs84") {
          # resolution of the product
          product_res = abs(Csom_all$lat[1,2]-Csom_all$lat[1,1])+abs(Csom_all$long[2,1]-Csom_all$long[1,1])
          product_res = product_res * 0.5
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

  answer = NA
  while (is.na(answer) == TRUE) {
    # work out average areas
    average_i = max(1,(i1-radius)):min(dim(Csom_all$Csom)[1],(i1+radius))
    average_j = max(1,(j1-radius)):min(dim(Csom_all$Csom)[2],(j1+radius))
    # carry out averaging
    Csom = mean(as.vector(Csom_all$Csom[average_i,average_j]), na.rm=TRUE)
    Csom_unc = mean(as.vector(Csom_all$Csom_unc[average_i,average_j]), na.rm=TRUE)
    # check for errors
    if (is.na(Csom) | Csom < 0) {radius = radius+1 ; answer = NA} else {answer = 1}
  }

  # retun back to the user 
  return(list(Csom_initial = Csom, Csom_initial_unc = Csom_unc))

} # extract_soilgrid_Csom

