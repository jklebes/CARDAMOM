
###
## Function to load leaf carbon per unit leaf area (gC/m2) maps
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_lca_maps_for_extraction<-function(latlon_in,lca_source) {

    ###
    ## Select the correct LCA source for specific time points

    if (lca_source == "butler") {

        # let the user know this might take some time
        print("Loading Butler LCA map...")

        # Create the full file paths estimates and their uncertainty (gC/m2)
        input_file = list.files(path_to_lca)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        unc_input_file = input_file[grepl("LCA_SD",input_file) == TRUE]
        input_file = input_file[grepl("LCA_SD",input_file) == FALSE]
        # Check that we have the same number of files for both biomass and uncertainty
        if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        if (length(input_file) > 1 | length(unc_input_file) > 1) {stop("More than one file has been found for the estimate and its uncertainty, there should only be one")}

        # Read in the estimate and uncertainty rasters
        lca_gCm2 = raster(paste(path_to_lca,input_file,sep=""))
        lca_uncertainty_gCm2 = raster(paste(path_to_lca,unc_input_file,sep=""))

        # Store dimension information
        dims = dim(lca_gCm2)[1:2]
        # Extract latitude / longitude information
        lat = coordinates(lca_gCm2)
        # Split between long and lat
        long = lat[,1] ; lat = lat[,2]
        # Reconstruct the full lat / long grid and flip dimensions as needed
        long = array(long, dim=c(dims[2],dims[1]))
        lat = array(lat, dim=c(dims[2],dims[1]))
        long = long[,dim(long)[2]:1]
        lat = lat[,dim(lat)[2]:1]

        # filter around target area
        max_lat = max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
        min_lat = min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
        keep_lat_min = min(which(lat[1,] > min_lat))
        keep_lat_max = max(which(lat[1,] < max_lat))
        keep_long_min = min(which(long[,1] > min_long))
        keep_long_max = max(which(long[,1] < max_long))
        lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
        long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

        # Similarly break apart the raster and re-construct into the correct orientation
        lca_gCm2 = array(as.vector(lca_gCm2), dim=c(dims[2],dims[1]))
        lca_gCm2 = lca_gCm2[,dim(lca_gCm2)[2]:1]
        lca_uncertainty_gCm2 = array(as.vector(lca_uncertainty_gCm2), dim=c(dims[2],dims[1]))
        lca_uncertainty_gCm2 = lca_uncertainty_gCm2[,dim(lca_uncertainty_gCm2)[2]:1]

        # now remove the ones that are actual missing data
        lca_gCm2[which(as.vector(lca_gCm2) < 0)] = NA
        lca_uncertainty_gCm2[which(as.vector(lca_uncertainty_gCm2) < 0)] = NA

        # remove data outside of target area
        lca_gCm2 = lca_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
        lca_uncertainty_gCm2 = lca_uncertainty_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2]
        lca_gCm2 = array(lca_gCm2, dim=c(idim,jdim))
        lca_uncertainty_gCm2 = array(lca_uncertainty_gCm2, dim=c(idim,jdim))

        # Output variables
        return(list(lat = lat, long = long, lca_gCm2 = lca_gCm2, lca_uncertainty_gCm2 = lca_uncertainty_gCm2))

    } else {

        # Output dummy variables
        return(list(lat = -9999, long = -9999, lca_gCm2 = -9999, lca_uncertainty_gCm2 = -9999))

    } # which LCA source?

} # function end load_lca_maps_for_extraction

## Use byte compile
load_lca_maps_for_extraction<-cmpfun(load_lca_maps_for_extraction)
