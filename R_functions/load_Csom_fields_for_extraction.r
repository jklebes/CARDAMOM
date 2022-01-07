
###
## Function to load gridded dataset of soil prior information from HWSD
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_Csom_fields_for_extraction<-function(latlon_in,Csom_source,cardamom_ext,spatial_type) {

    if (Csom_source == "SoilGrids") {

        # This is a very bespoke modification so leave it here to avoid getting lost
        Csom = raster(paste(path_to_Csom,"Csom_gCm2_mean_0to1m.tif", sep=""))
        Csom_unc = raster(paste(path_to_Csom,"Csom_gCm2_sd_0to1m.tif", sep=""))

        # Create raster with the target crs
        target = raster(crs = ("+init=epsg:4326"), ext = extent(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() ; removeTmpFiles()
            Csom_unc = resample(Csom_unc, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext) ; Csom_unc = crop(Csom,cardamom_ext)
        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
        if (spatial_type == "grid") {
            if (res(Csom)[1] < res(cardamom_ext)[1] | res(Csom)[2] < res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                # Resample to correct grid
                Csom = resample(Csom, target, method="bilinear") ; gc() ; removeTmpFiles()
                Csom_unc = resample(Csom_unc, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution
        } # spatial_type == "grid"

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = coordinates(Csom)[,1] ; lat = coordinates(Csom)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))
        Csom_unc = array(as.vector(unlist(Csom_unc)), dim=c(xdim,ydim))

        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat,long = long))

    } else if (Csom_source == "HWSD") {

        # let the user know this might take some time
        print("Loading processed HWSD Csom fields for subsequent sub-setting ...")

        # open processed modis files
        input_file_1=paste(path_to_Csom,"/HWSD_Csom_with_lat_long.nc",sep="")
        data1=nc_open(input_file_1)

        # extract location variables
        lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
        # read the modis lai drivers
        hwsd_Csom=ncvar_get(data1, "HWSD_Csom")

        if (length(dim(latlon_in)) > 1) {
            max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
            min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
        } else {
            max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
            min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
        }
        keep_lat=which(lat[1,] > min_lat & lat[1,] < max_lat)
        keep_long=which(long[,1] > min_long & long[,1] < max_long)
        hwsd_Csom=hwsd_Csom[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        lat=lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        long=long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	      # close files after use
        nc_close(data1)

        # clean
        rm(keep_lat,keep_long,max_lat,max_long,min_lat,min_long) ; gc(reset=TRUE,verbose=FALSE)

        # output variables
        return(list(Csom=hwsd_Csom,Csom_unc = -9999,lat=lat,long=long))

    } else {
        # output variables
	      return(list(Csom=-9999, Csom_unc = -9999, lat=-9999,long=-9999))
    }

} # function end load_Csom_fields_for_extraction

## Use byte compile
load_Csom_fields_for_extraction<-cmpfun(load_Csom_fields_for_extraction)
