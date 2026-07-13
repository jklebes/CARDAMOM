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
# Function to extract met data for ACM recalibration
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

extract_acm_met_drivers<-function(PROJECT,latlon_wanted,site_name) {

	# construct assumed name
	infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
#  infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")
#  infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")
#  infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_avgNlessthan4.csv",sep="")
  infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_avgNlessthan4_subsample.csv.csv",sep="")
	# extract the drivers currently in use
	maxt_out=read_site_specific_obs("sat_max",infile)
	mint_out=read_site_specific_obs("sat_min",infile)
	swrad_out=read_site_specific_obs("swrad_avg",infile) ; swrad_out=swrad_out*1e-6*86400 # W.m-2 -> MJ.m-2.day-1
	co2_out=read_site_specific_obs("co2_avg",infile)
	doy_out=read_site_specific_obs("doy",infile)
	avgN_out=read_site_specific_obs("avgN",infile)
	lai_out=read_site_specific_obs("lai",infile)
	lat_out=read_site_specific_obs("lat",infile)
	vpd_out=read_site_specific_obs("vpd_avg",infile)*1e3 # convert to Pa
	wind_out=read_site_specific_obs("wind_avg",infile)
	precip_out=read_site_specific_obs("ppt_avg",infile) ; precip_out=precip_out/(60*60) # average mm.step -> kg.m-2.s-1
	Rtot_out=read_site_specific_obs("Rtot",infile)
	top_sand = read_site_specific_obs("sand_top",infile) # % sand in top soil
	bot_sand = read_site_specific_obs("sand_bot",infile) # % sand in deep soil
	top_clay = read_site_specific_obs("clay_top",infile) # % clay in top soil
	bot_clay = read_site_specific_obs("clay_bot",infile) # % clay in deep soil

  # update start / end year information
	year_out=read_site_specific_obs("year",infile)
	PROJECT$start_year=min(as.numeric(year_out))
	PROJECT$end_year=max(as.numeric(year_out))

	# not actually used but maybe one day
	lagged_precip_out=read_site_specific_obs("lagged_precip",infile)
	avgTmax_out=read_site_specific_obs("avgTmax",infile)
	photoperiod_out=read_site_specific_obs("photoperiod",infile)
	vpd_lagged_out=read_site_specific_obs("vpd_lagged",infile)

	# create day of run variable
	run_day=seq(1:length(swrad_out))

	# output variables
	return(list(run_day=run_day,mint=mint_out,maxt=maxt_out,swrad=swrad_out,co2=co2_out,doy=doy_out,lagged_precip=lagged_precip_out
             ,avgTmax=avgTmax_out,photoperiod=photoperiod_out,vpd_lagged=vpd_lagged_out,vpd=vpd_out,avgN=avgN_out
             ,lai=lai_out,lat=lat_out,wind_spd=wind_out,precip=precip_out,Rtot=Rtot_out,top_sand = top_sand
             ,bot_sand = bot_sand, top_clay  = top_clay, bot_clay = bot_clay))

} # function end extract_acm_met_drivers

## Use byte compile
extract_acm_met_drivers<-cmpfun(extract_acm_met_drivers)
