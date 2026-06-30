#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Functions to aggregate spatial information to coarser resolutions
# 
# Author: T. Luke Smallman (12/11/2024)
# Exceptions states below in specific functions
#
#########################################################################################

# Griddify function (marmap library) updated to work with terra libraries
griddify <-function (xyz, nlon, nlat) {
    if (ncol(xyz) != 3) {stop("xyz must be a 3-column matrix or data.frame of longitudes, latitudes and depths/altitudes")}
    colnames(xyz) <- c("lon", "lat", "depth")
    r <- terra::rast(xmin = min(xyz$lon, na.rm = TRUE), xmax = max(xyz$lon, na.rm = TRUE),
                     ymin = min(xyz$lat, na.rm = TRUE), ymax = max(xyz$lat, na.rm = TRUE),
                     ncol = nlon, nrow = nlat)
    xy = terra::vect(xyz, geom = c("lon","lat"))
    x <- terra::rasterize(xy, r, field = "depth", fun = mean)
    s <- terra::rast(nrow = nlat, ncol = nlon)
    terra::ext(s) <- terra::ext(x)
    s <- terra::resample(x, s, method = "bilinear")
    return(s)
}

# Define local function
regrid_func<-function(var1_in, epsg_in, lat_in, long_in, cardamom_ext) {

   # Set flags
   lat_done = FALSE

   # Check dimensions of the input variable
   if (length(dim(var1_in)) == 2) {
       var1_in = array(var1_in, dim=c(dim(var1_in),1))
   }

   # Loop through each timestep in the year
   for (t in seq(1, dim(var1_in)[3])) {
        # Convert to a raster, assuming standard WGS84 grid
        var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1_in[,,t]))
        # Determine whether the grid is regular using a tolerance (all.equal,
        # ~1.5e-8 relative) rather than exact equality via unique(), so floating-point
        # round-off in the coordinates does not spuriously route a regular grid to
        # griddify() (which bilinearly interpolates and silently alters the values).
        lon_regular = isTRUE(all.equal(min(diff(long_in[,1])), max(diff(long_in[,1]))))
        lat_regular = isTRUE(all.equal(min(diff(lat_in[1,])),  max(diff(lat_in[1,]))))
        if (!lon_regular | !lat_regular) {
            var1 = griddify(var1, dim(var1_in)[1], dim(var1_in)[2])
        } else {
            var1 = rast(var1, crs = (epsg_in), type="xyz")
        }

        # Create raster with the target crs (technically this bit is not required)
        target = rast(crs = epsg_in, ext = ext(var1), resolution = res(var1))

        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(var1,target) == FALSE) {
            # Resample to correct grid
            #var1 = resample(var1, target, method="near") ; gc()
            var1 = project(var1, target, mask = FALSE, method="near")
        }

        if (ext(cardamom_ext) != ext(var1)) {
            # Extend the extent of the overall grid to the analysis domain
            #var1 = extend(var1,cardamom_ext, snap="near", fill = NA)
            var1 = extend(var1,cardamom_ext)
            # Trim the extent of the overall grid to the analysis domain
            var1 = crop(var1,cardamom_ext)
            # With different resolutions that are non-divisible with the target
            # resolution a mismatch may still be present. However, this should be
            # closer than it was and allow for the final adjustments made in the resampling
        }

        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here.
        # Despite creation of a cardamom_ext for a site run do not allow aggregation here as this will damage the fine resolution datasets
        if (res(var1)[1] != res(cardamom_ext)[1] |
            res(var1)[2] != res(cardamom_ext)[2] |
            ext(cardamom_ext) != ext(var1)) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            var1 = resample(var1, target, method="average") ; gc()

        } # Aggregate to resolution

        # If a land mask is present then also restrict to this target domain
        #if (missing("landmask") == FALSE) {var1 = crop(var1,landmask)}

        if (lat_done == FALSE) {
            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(var1)[2] ; ydim = dim(var1)[1]
            # extract the lat / long information needed
            long = crds(var1,df=TRUE, na.rm=FALSE)
            lat  = long$y ; long = long$x
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
        }
        # break out from the raster into arrays which we can manipulate
        var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

        # vectorise at this time
        if (lat_done == FALSE) {
            var_out = as.vector(var1)
        } else {
            var_out = append(var_out,as.vector(var1))
        }

        # update flag for lat / long load
        if (lat_done == FALSE) {lat_done = TRUE}

   } # Within variable time loop

   # restructure
   var_out = array(var_out, dim=c(xdim,ydim,dim(var1_in)[3]))

   # Return to user
   return(list(var = var_out, lat = lat, long = long))

} # end function regrid_func
