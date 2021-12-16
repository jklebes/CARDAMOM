
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
    # extract the desired information,
    lat=output$lat ; long=output$long ; cardamom_ext = output$cardamom_ext ; rm(output)

    # remove the values we don't want
    if (length(remove) > 0) {lat = lat[-remove] ; long = long[-remove]}

    # output the result
    output=list(lat = lat, long = long, cardamom_ext = cardamom_ext)
    rm(lat,long) ; gc(reset=TRUE, verbose=FALSE)
    return(output)

}
## Use byte compile
determine_lat_long_needed<-cmpfun(determine_lat_long_needed)
