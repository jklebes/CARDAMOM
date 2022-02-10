
###
## Function determines all lat long coordinates needed for grid run mode
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

determine_lat_long_needed<- function(lat,long,resolution,grid_type,remove) {

    # check input data
    if (length(which(long > 180)) > 0) {stop("Long should be -180 to +180")}

    # generate UK or WGS-84 lat long grid
    if (grid_type == "UK") {
        output = generate_uk_grid(lat,long,resolution)
    } else if (grid_type=="wgs84") {
        output = generate_wgs84_grid(lat,long,resolution)
    } else {
        stop('have selected invalid grid type, the valid options are "UK" and "wgs84"')
    }
    # extract the latitude / longitude and extent/resolution information
    lat = output$lat ; long = output$long ; long_dim = output$long_dim ; lat_dim = output$lat_dim
    cardamom_ext = output$cardamom_ext

    # Create a grid specifically to be used for extracting the correct location from the gridded datasets
    obs_long_grid = array(long, dim=c(long_dim,lat_dim))
    obs_lat_grid = array(rev(lat), dim=c(long_dim,lat_dim)) # rev() accounts for the flipping of orientation conducted in the generate_*_grid()

    # remove the values we don't want
    if (length(remove) > 0) {lat = lat[-remove] ; long = long[-remove]}

    # output the result
    output = list(lat = lat, long = long, cardamom_ext = cardamom_ext,
                  obs_long_grid = obs_long_grid, obs_lat_grid = obs_lat_grid)
    rm(lat,long) ; gc(reset=TRUE, verbose=FALSE)
    return(output)

}
## Use byte compile
determine_lat_long_needed<-cmpfun(determine_lat_long_needed)
