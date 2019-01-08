
###
## Function to extract state variables and direct the production of uncertainty plots for the key states and fluxes
###

generate_simplified_stock_and_flux_maps<-function(PROJECT) {

  print("...beginning stocks and fluxes...")

  # read in the grid simmary information
  infile = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
  if (file.exists(infile) == FALSE) {stop("grid_outputs for 'generate_simplified_stock_and_flux_maps' missing")}
  load(paste(infile))
  infile=paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep="")
  if (file.exists(infile) == FALSE) {stop("parameter_maps for 'generate_simplified_stock_and_flux_maps' missing")}
  load(paste(infile))

  # work out area matrix for the pixels in meters
  # include adjustment for g-> Tg (*1e-12)
  if (PROJECT$grid_type == "UK") {
    area_with_g_Tg = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))*1e-12
    area = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
  } else if (PROJECT$grid_type == "wgs84") {
    # generate the lat / long grid again
    output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
    # then generate the area estimates for each pixel
    area_with_g_Tg = calc_pixel_area(output$lat,output$long,PROJECT$resolution)*1e-12
    area = calc_pixel_area(output$lat,output$long,PROJECT$resolution)
    # this output is in vector form and we need matching array shapes so...
    area = array(area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
    area_with_g_Tg = array(area_with_g_Tg, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
  } else {
    stop("valid spatial grid option not selected (UK, or wgs84)")
  }

  # extract a list of all the variables stored in the output object
  par_names = names(grid_output)
  # filter for those related to the 'mean' status
  par_names = par_names[grepl("mean", par_names)]
  # determine the array value for the median,
  num_quantiles = dim(grid_output$mean_labile_gCm2)[3]
  if (num_quantiles == 5) {
    # then we assume we are dealing with 0.025, 0.25, 0.5, 0.75, 0.975 quantiles
    median_loc = 3 ; lower_loc = 1 ; upper_loc = 5
  } else {
    # otherwise we need to approximate it...
    median_loc = round(num_quantiles / 2,digits=0)
    lower_loc = ceiling(num_quantile * 0.025)
    upper_loc = floor(num_quantile * 0.975)
  }

  # calculate land mask
  landmask=array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
  # load colour palette
  colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
  # array sizes are always the same so
  colour_choices = colour_choices_upper(length(area))

  # determine correct height and widths
  fig_height=4000 ; fig_width=7200
  if (PROJECT$grid_type == "UK") { fig_height=8000 ; fig_width=7200 }

  mean_rooting_depth = NA
  if (PROJECT$model$name == "DALEC_BUCKET" | PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
    jpeg(file=paste(PROJECT$figpath,"median_root_depth_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
    mean_rooting_depth = par_array_median[,,40] * (grid_output$mean_roots_gCm2[,,median_loc]*2) / (par_array_median[,,39] + (grid_output$mean_roots_gCm2[,,median_loc]*2))
    z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
    image.plot(mean_rooting_depth,col=colour_choices, main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
    dev.off()
  } else if (PROJECT$model$name == "DALECN_BUCKET") {
    jpeg(file=paste(PROJECT$figpath,"median_root_depth_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
    par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
    mean_rooting_depth = par_array_median[,,35] * (grid_output$mean_roots_gCm2[,,median_loc]*2) / (par_array_median[,,34] + (grid_output$mean_roots_gCm2[,,median_loc]*2))
    z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
    image.plot(mean_rooting_depth,col=colour_choices, main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
    dev.off()
  }

  # loop through these mean variables, output median and CI range for each of these variables
  for (p in seq(1,length(par_names))) {

    # determine position in grid_output list which contains the variable of interest
    pp = which(names(grid_output) == par_names[p])

    jpeg(file=paste(PROJECT$figpath,"grid_mean_map_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=7200, height=4000, res=300, quality=100)
    info = " " # assume default is no header, but sometimes we add something extra...
    var1 = mean(grid_output[[pp]][,,median_loc], na.rm=TRUE)
    var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
    var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
    var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
    info = paste("Mean estimate: ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3, sep="")
    image.plot(grid_output[[pp]][,,median_loc], main=info, col = colour_choices, axes=FALSE, cex.main=2.4, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
    contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
    dev.off()

  } # loop through all "grid_output" objects

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

} # end function generate_simplified_stock_and_fluxe_maps
