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
# Function to estimate the number of days in a given year
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

nos_days_in_year<-function(year) {

    # This function uses scalar control flow; a vector year would error under
    # R >= 4.2 ("condition has length > 1"). Enforce the scalar contract explicitly.
    if (length(year) > 1) {stop("nos_days_in_year expects a single year")}

    # is current year a leap or not
    nos_days = 365
    mod = as.numeric(year)-round((as.numeric(year)/4))*4
    if (mod == 0) {
        nos_days = 366
        mod = as.numeric(year)-round((as.numeric(year)/100))*100
        if (mod == 0) {
            nos_days  = 365
            mod = as.numeric(year)-round((as.numeric(year)/400))*400
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
