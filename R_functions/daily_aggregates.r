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
# Some low level functions for aggregating with a specified number of missing observations.
# NOTE: these are not that efficient so use roll apply() when missing data is not an issue.
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

daily_mean <-function(var, interval, missing_allowed) {

   # work out how many intervals fit
   # i.e. number of days possible
   nos_days=ceiling(length(var)/interval)
   output=array(NaN, dim=c(nos_days))
   b=1
   for (i in seq(1, nos_days)) {
        if (length(which(is.na(var[b:(b+interval-1)]))) < missing_allowed) {
          output[i] = mean(var[b:(b+interval-1)], na.rm=T)
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
   b=1
   for (i in seq(1, nos_days)) {
       if (length(which(is.na(var[b:(b+interval-1)]))) < missing_allowed) {
           output[i]=sum(var[b:(b+interval-1)], na.rm=T)
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
