
###
## Some low level functions for aggregating with a specified number of missing observations.
## NOTE: these are not that efficient so use roll apply() when missing data is not an issue.
###

# This function was coded by T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)

daily_mean <-function(var, interval, missing_allowed) {

   # work out how many intervals fit
   # i.e. number of days possible
   nos_days=ceiling(length(var)/interval)
   output=array(NaN, dim=c(nos_days))
   b=0
   for (i in seq(1, nos_days)) {
        if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
          output[i] = mean(var[b:(b+interval)], na.rm=T)
        } else {
          output[i] = NaN
        }
        b = b + interval
   }
   # clean up
   rm(nos_days,b,i) ; gc()
   return(output)

} # end function daily_mean

## Use byte compile
daily_mean<-cmpfun(daily_mean)

daily_sum <-function(var, interval, missing_allowed) {

   # work out how many intervals fit
   # i.e. number of days possible
   nos_days=ceiling(length(var)/interval)
   output=array(NaN, dim=c(nos_days))
   b=0
   for (i in seq(1, nos_days)) {
       if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
           output[i]=sum(var[b:(b+interval)], na.rm=T)
       } else {
         output[i]=NaN
       }
       b=b+interval
   }
   # clean up
   rm(nos_days,b,i) ; gc()
   return(output)

} # end function daily_sum

## Use byte compile
daily_sum<-cmpfun(daily_sum)
