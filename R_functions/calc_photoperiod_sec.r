
###
## Function to estimate the number of seconds per day
###

# This function was coded by T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
# See function for relevant source references

calc_photoperiod_sec<-function(lat,days){

   # Function calculates the day length in hours based on day of year and latitude (degrees).
   # the output is daylength converted to seconds
   # Ref: NEEDS TO BE ADDED

   declin    = - asin ( sin ( 23.45 * ( pi / 180 ) ) * cos ( 2. * pi * ( days + 10. ) / 365. ) )
   sinld     = sin ( lat*(pi/180.) ) * sin ( declin )
   cosld     = cos ( lat*(pi/180.) ) * cos ( declin )
   aob       = sinld / cosld
   aob       = pmax(-1.0,pmin(1.0,sinld / cosld))
   daylength = 12.0 * ( 1. + 2. * asin ( aob ) / pi )
   # convert hours to seconds
   daylength=daylength*3600
   # clean up
   rm(declin,sinld,cosld,aob) ; gc()
   # now return
   return(daylength)

} # end function calc_photoperiod_sec

## Use byte compile
calc_photoperiod_sec<-cmpfun(calc_photoperiod_sec)
