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
# Function to estimate the potential incoming shortwave radiation
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

  est_potential_radiation <- function(doy,latitude,diurnal_range,mean_diurnal_range,surf_pressure_Pa,water_vapour_Pa,hardcode_Tt) {

    # Function to calculate an estimate of potential short wave radiation (MJ.m-2.day-1).
    # Day of year and latitude (degrees) are used to first estimate potential clear skys radiation.
    # Cloudiness fraction is then estimated from diurnal range and used to scale down the estimate.
    # Ref: Thornton & Running (1999) Agri Forest Met 93, 211-228

    # The approach used here assumes vectoriation to generate 24 hour estimate

    # parameters
    So = 1360 ; deg_to_rad = pi/180.0 ; hours_of_day = seq(1,24,1) ; hours_in_day = 24 ; seconds_per_hour = 60*60
    sea_surface_Pa = 101325
    alpha = -6.1e-5 # Pa-1 impact on transmittance per Pa of water vapour
    Tnadir = 0.87 # maximum clear sky transmittance on a dry day

    # convert latitude in radians
    latitude_radians = latitude * deg_to_rad

    # calculate declination
    declination = - asin ( sin ( 23.45 * deg_to_rad ) * cos ( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    #declination = max(0.0,declination)
    # calculate angle of sun each hour
    hourangle = deg_to_rad * 15. * ( hours_of_day - 0.5 * hours_in_day ) * 24. / hours_in_day
    # estimate cosine of solar zenith angle
    cos_solar_zenith_angle = sin(latitude_radians) * sin(declination) + cos(latitude_radians) * cos(declination) * cos(hourangle)
    # calculate potential radiation (MJ.m-2.s-1)
    est_potential_radiation = So * cos_solar_zenith_angle
    est_potential_radiation = pmax(0.0,est_potential_radiation)
    if (hardcode_Tt == 0) {
        b0 = 0.031 ; b1 = 0.201 ; b2 = 0.185 ; B = b0 + b1 * exp(-b2 * mean_diurnal_range) ; C = 2.4
        # estimate optical air mass for given solar angle
        m0 = cos_solar_zenith_angle**-1
        # estimate Maximum total transmittance fraction
        Ttmax = Tnadir**m0
        if (surf_pressure_Pa != 0) { Ttmax = Ttmax ** (surf_pressure_Pa/sea_surface_Pa) }
        if (water_vapour_Pa != 0) { Ttmax = Ttmax + alpha*water_vapour_Pa }
        # estimate Total transmittance fraction
        Tt = Ttmax*(1-exp(-B * diurnal_range**C))
    } else {
        if (hardcode_Tt > Tnadir | hardcode_Tt < 0.1) {return(print("Inputted Transmittance fraction is outside of allowable bounds"))}
        Tt = hardcode_Tt
    }
    # upscale to MJ.m-2.day-1
    est_potential_radiation = sum(est_potential_radiation*Tt*seconds_per_hour)*1e-6

    # explicit return
    return(est_potential_radiation)

  } # function est_potential_radiation

## Use byte compile
est_potential_radiation<-cmpfun(est_potential_radiation)
