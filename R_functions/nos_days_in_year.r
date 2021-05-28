
###
## Function to estimate the number of days in a given year
###

# This code was created by T. L Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)

nos_days_in_year<-function(year) {

    # is current year a leap or not
    nos_days = 365
    mod=as.numeric(year)-round((as.numeric(year)/4))*4
    if (mod == 0) {
        nos_days = 366
        mod=as.numeric(year)-round((as.numeric(year)/100))*100
        if (mod == 0) {
            nos_days  = 365
            mod=as.numeric(year)-round((as.numeric(year)/400))*400
            if (mod == 0) {
                nos_days  = 366
            }
        }
    }

    # clean up
    rm(mod) ; gc()

    # return to user
    return(nos_days)

} # end function nos_days_in_year

## Use byte compile
nos_days_in_year<-cmpfun(nos_days_in_year)
