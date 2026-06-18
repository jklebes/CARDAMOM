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
# Function for converting between specific humidity, relative humidity and 
# vapour pressure deficit. 
# These functions were coded by T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
# See function for relevant source references
#
#########################################################################################

sp_humidity_to_vpd<-function(sp_moist,atmos_press,air_temperature) {

   # Converts specific humidity (kg/kg) into VPD (Pa)
   # Determine vapour pressure (Pa); based on specific humidity, air temperature (oC),
   # air pressure (Pa->(hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
   # a physical outline)

   # calculate vapour pressure of the air
   vpair=((sp_moist*(atmos_press*1.0e-2))/0.62197)*1.0e2
   # Saturation vapour pressure (Pa) calculation from Jones p110; uses
   # absolute air temperature (oC)
   vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
   # Difference between the vapour pressure of saturation and air, i.e. the
   # VPD (Pa)
   vpd_pa = vpsat - vpair
   # clean up
   rm(vpair,vpsat) ; gc()
   # return to user
   return(vpd_pa)

} # end function sp_humidity_to_vpd
## Use byte compile
sp_humidity_to_vpd<-cmpfun(sp_humidity_to_vpd)

vpd_to_rh<-function(vpd_in,air_temperature) {

   # Converts VPD (Pa) to rel humidity (frac)
   # Determine vapour pressure (Pa)"," based on specific humidity, air
   # pressure (Pa input)
   # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
   # a physical outline)

   # Saturation vapour pressure (Pa) calculation from Jones p110"," uses
   # absolute air temperature (oC)
   vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
   # RH is the ratio of vapour presure in the air and vapour pressure at saturation.
   # Below pressure is estimated from the saturation vapour pressure and vapour pressure deficit.
   # Units (Pa)
   rh = (vpsat - vpd_in) / vpsat
   rh[rh > 1] <- 1
   rh[rh < 0] <- 0
   # clean up
   rm(vpsat) ; gc()
   # return to user
   return(rh)

} # end function vpd_to_rh
## Use byte compile
vpd_to_rh<-cmpfun(vpd_to_rh)

dew_temp_to_sp_humidity<-function(dew_airt,airt,pressure) {

  # Ref: p95 McIlveen 1986, Basic Meteorology -
  # a physical outline

  # dew_airt (oC)
  # airt (oC)
  # pressure (Pa->hPa)
  # vapour pressues (hPa->Pa)

  ## Specific humidity (kg/kg)
  vapour_pressure = 6.11*10**((7.5*dew_airt)/(237.3+dew_airt))
  dew_temp_to_sp_humidity = 0.622*vapour_pressure/(pressure*1e-2)

} # end function dew_temp_to_sp_humidity
## Use byte compile
dew_temp_to_sp_humidity<-cmpfun(dew_temp_to_sp_humidity)

rh_to_vpd<-function(rh_in,air_temperature) {

   # Converts rel humidity (frac) into VPD (Pa)
   # Determine vapour pressure (Pa)"," based on specific humidity, air
   # pressure (Pa input)
   # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
   # a physical outline)

   # Saturation vapour pressure (Pa) calculation from Jones p110"," uses
   # absolute air temperature (oC)
   vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
   # Difference between the vapour pressure of saturation and air, i.e. the
   # VPD (Pa)
   vpair = vpsat*rh_in
   vpd_pa = vpsat-vpair

   # clean up
   rm(vpsat) ; gc()
   # return to user
   return(vpd_pa)

} # end function rh_to_vpd
## Use byte compile
rh_to_vpd<-cmpfun(rh_to_vpd)

sp_humidity_to_rh <- function(qair, temp, press = 1013.25){

   # Ref: p95 McIlveen 1986, Basic Meteorology -
   # a physical outline

   # temperature in oC, specific humidity pressure in (hPa)
   es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
   e <- qair * press / (0.378 * qair + 0.622)
   rh <- e / es
   rh[rh > 1] <- 1
   rh[rh < 0] <- 0
   # clean up
   rm(es,e) ; gc()
   return(rh)
}
## Use byte compile
sp_humidity_to_rh<-cmpfun(sp_humidity_to_rh)
