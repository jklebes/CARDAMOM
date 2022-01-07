
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
        Csom = crop(Csom,cardamom_ext) ; Csom_unc = crop(Csom_unc,cardamom_ext)
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
        # read the HWSD soil C prior
        Csom=ncvar_get(data1, "HWSD_Csom")

        # Convert to a raster, assuming standad WGS84 grid
        Csom = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(Csom))
        Csom = rasterFromXYZ(xyz, crs = ("+init=epsg:4326"))

        # Create raster with the target crs (technically this bit is not required)
        target = raster(crs = ("+init=epsg:4326"), ext = extent(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext) 
        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
        if (spatial_type == "grid") {
            if (res(Csom)[1] < res(cardamom_ext)[1] | res(Csom)[2] < res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                # Resample to correct grid
                Csom = resample(Csom, target, method="bilinear") ; gc() ; removeTmpFiles()

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

        # see papers assessing uncertainty of HWSD, ~47 %
        Csom_unc = array(Csom * 0.47, dim=c(xdim,ydim))
        # With a minium bound assumption
        Csom_unc[Csom_unc < 100] = 100

        # Return the loaded dataset
        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat,long = long))

    } else {
        # output variables
	      return(list(Csom=-9999, Csom_unc = -9999, lat=-9999,long=-9999))
    }

} # function end load_Csom_fields_for_extraction

## Use byte compile
load_Csom_fields_for_extraction<-cmpfun(load_Csom_fields_for_extraction)
