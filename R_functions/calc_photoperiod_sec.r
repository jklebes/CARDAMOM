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
# This function is to estimate the number of seconds per day
# This function was coded by T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
# See function for relevant source references
#
#########################################################################################

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
   daylength = daylength*3600
   # clean up
   rm(declin,sinld,cosld,aob) ; gc()
   # now return
   return(daylength)

} # end function calc_photoperiod_sec

## Use byte compile
calc_photoperiod_sec<-cmpfun(calc_photoperiod_sec)
