

###
## Function to determine the area of a pixel at any location in meter squared
###

## function to convert degrees to radians
deg2rad<-function(degree_in) {return(degree_in*(pi/180))}
rad2deg<-function(radian_in) {return(radian_in*(180/pi))}

calc_pixel_area<-function(latitude,longitude,resolution) {

    # resolution in degrees
    # latitude (-90,90) degrees
    # longitude (-180,180) degrees

    # mean earth radius (m)
    R = 6371e3

    pixel_area = R**2 * ( deg2rad(longitude+resolution*0.5)-deg2rad(longitude-resolution*0.5) ) * (sin(deg2rad(latitude+resolution*0.5))-sin(deg2rad(latitude-resolution*0.5)))

    # return to user in meters
    return(pixel_area)

}
