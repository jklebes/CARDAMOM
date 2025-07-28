
###
## Process CARDAMOM-DALEC output files into NetCDF files 
## consistent with the TRENDYv14 / GCP model intercomparison structure
## In constrast to the sibling script which groups variables together into different files,
## this script places everything into a single file per document consistent with the latest guidance.
### 

###
## Job specific information

print("Begin creation of Trendy v14 compatible single variable netcdf files...")

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.004_MHMCMC/global_0.5deg_dalec4_trendyv14_LCA_TWB_GPP_fAPAR"

# Specify any extra information for the filename
output_prefix = "CARDAMOM_S3_" # follow with "_"
#output_prefix = "CARDAMOM_S2_" # follow with "_"
output_suffix = "" # begin with "_"

###
## Load libraries, functions and data needed

# load needed libraries
library(ncdf4)
library(terra)
library(compiler)
library(zoo)

# Load needed functions 
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/load_all_cardamom_functions.r")

# load the CARDAMOM files
load(paste(input_dir,"/infofile.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

# Create output directory to aid storage management
out_dir = paste(PROJECT$results_processedpath,"/trendy_output/",sep="")
if (dir.exists(out_dir) == FALSE) {
    dir.create(out_dir)
}

###
## Load required functions

extract_from_met_drivers<-function(n, met_index, bias_adj, scale_adj) {

   ## Dependancies
   # 1) 'grid_output' and 'PROJECT' objects being loaded into global memory
   # 2) read_binary_file_format() loaded into global memory
   
   ## Arguments
   # n = site counter
   # met_index = which array member
   # bias_adj = a constant value which will be added to the extracted variable
   # scale_adj = a constant value which will be multiplied to the extracted variable

   # Check whether this is a viable location
   if (is.na(grid_output$i_location[n]) == FALSE) {   

       # Read in site specific drivers
       drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
       # Missing data flag
       drivers$met[drivers$met[,met_index] == -9999,met_index] = NA       
       # Extract the desired variable
       forcings[grid_output$i_location[n],grid_output$j_location[n],] <<- (drivers$met[,met_index] + bias_adj) * scale_adj
   
       return(0)

   } else {
   
       return(-1)
   
   } # valid location

} # end function extract_from_met_drivers
                          
extract_from_obs_drivers<-function(n, obs_index, bias_adj, scale_adj) {

   ## Assumptions
   # 1) 'grid_output' and 'PROJECT' objects being loaded into global memory
   # 2) read_binary_file_format() loaded into global memory
   # 3) That est, unc, lag exist in global memory
   
   ## Arguments
   # n = site counter
   # met_index = which array member
   # bias_adj = a constant value which will be added to the extracted variable
   # scale_adj = a constant value which will be multiplied to the extracted variable
   #             NOTE neither adj value is assumed to be applied to the 'lag' variable

   # Check whether this is a viable location
   if (is.na(grid_output$i_location[n]) == FALSE) {   

       # Read in site specific drivers
       drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))

       # Missing data flag
       drivers$obs[drivers$obs[,obs_index] == -9999,obs_index] = NA
       drivers$obs[drivers$obs[,obs_index+1] == -9999,obs_index] = NA
       drivers$obs[drivers$obs[,obs_index+2] == -9999,obs_index] = NA
                     
       # Extract the desired variable
       est[grid_output$i_location[n],grid_output$j_location[n],] <<- (drivers$obs[,obs_index] + bias_adj) * scale_adj
       unc[grid_output$i_location[n],grid_output$j_location[n],] <<- (drivers$obs[,obs_index+1] + bias_adj) * scale_adj
       lag[grid_output$i_location[n],grid_output$j_location[n],] <<- (drivers$obs[,obs_index+2])
   
       return(0)

   } else {
   
       return(-1)
   
   } # valid location

} # end function extract_from_obs_drivers

# Function to extract the variables from grid_output
extract_from_grid_output<-function(n, var_name, bias_adj, scale_adj) {

   ## Dependancies
   # 1) 'grid_output' and 'PROJECT' objects being loaded into global memory
   # 2) read_binary_file_format() loaded into global memory
   
   ## Arguments
   # n = site counter
   # met_index = which array member
   # bias_adj = a constant value which will be added to the extracted variable
   # scale_adj = a constant value which will be multiplied to the extracted variable

   # Check whether this is a viable location
   if (is.na(grid_output$i_location[n]) == FALSE) {   

       # Extract the desired variable
       est[grid_output$i_location[n],grid_output$j_location[n],,] <<- (grid_output[[var_name]][n,,] + bias_adj) * scale_adj
   
       return(0)

   } else {
   
       return(-1)
   
   } # valid location

} # end function extract_from_grid_output
                      
# Write forcing to netcdf file
write_to_nc_forcing<-function(var_est,var_name,var_unit,var_long) {

   ## Assumptions
   # 1) That PROJECT object has been loaded into global memory
   # 2) That area_m2 array has been loaded into global memory
   # 3) That land_fraction array has been loaded into global memory
   # 4) That analysis timestep (unit of cumulative days) has been loaded into global memory

   # Define the output file name
   output_name = paste(out_dir,output_prefix,var_name,output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   # Median
   var_new  = ncvar_def(var_name, unit=var_unit, longname = var_long, dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   
   # Create the empty file space
   new_file = nc_create(filename = output_name, vars = list(var0,var1,var2,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, analysis_timestep)
   # Grid area 
   ncvar_put(new_file, var1, area_m2)
   # Land fraction
   ncvar_put(new_file, var2, land_fraction)
   # VARIABLE
   ncvar_put(new_file, var_new, var_est)
   
   # Close the existing file to ensure its written to file
   nc_close(new_file)

   # Tidy
   rm(var_new)

} # write_to_nc_forcing

# Write forcing to netcdf file
write_to_nc_observation<-function(var_est,var_unc,var_lag,
                                  var_est_name,var_unc_name,var_lag_name,
                                  var_est_unit,var_lag_unit,
                                  var_est_long,var_unc_long,var_lag_long) {

   ## Assumptions
   # 1) That PROJECT object has been loaded into global memory
   # 2) That area_m2 array has been loaded into global memory
   # 3) That land_fraction array has been loaded into global memory
   # 4) That analysis timestep (unit of cumulative days) has been loaded into global memory
   # 5) That the units are common for the estimate and the uncertainty

   # Define the output file name
   output_name = paste(out_dir,output_prefix,var_est_name,output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   # Median
   var_est_new  = ncvar_def(var_est_name, unit=var_est_unit, longname = var_est_long, dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   var_unc_new  = ncvar_def(var_unc_name, unit=var_est_unit, longname = var_unc_long, dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   var_lag_new  = ncvar_def(var_lag_name, unit=var_lag_unit, longname = var_lag_long, dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)      
   
   # Create the empty file space
   new_file = nc_create(filename = output_name, vars = list(var0,var1,var2,var_est_new,var_unc_new,var_lag_new), force_v4 = TRUE)

   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, analysis_timestep)
   # Grid area 
   ncvar_put(new_file, var1, area_m2)
   # Land fraction
   ncvar_put(new_file, var2, land_fraction)
   # VARIABLE
   ncvar_put(new_file, var_est_new, var_est)
   ncvar_put(new_file, var_unc_new, var_unc)
   ncvar_put(new_file, var_lag_new, var_lag)
   
   # Close the existing file to ensure its written to file
   nc_close(new_file)
   
   # Tidy
   rm(var_est_new,var_unc_new,var_lag_new)

} # write_to_nc_observation

# write analysis output to netcdf file at the analysis timestep
write_to_nc_analysis_timestep<-function(var_est,var_name,var_unit,var_long) {

   # Define the output file name
   output_name = paste(out_dir,output_prefix,var_name,output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   # Median
   var_mid = ncvar_def(var_name, unit = var_unit, longname = paste(var_long," - Median estimate",sep=""), 
                      dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Quantiles for uncertainty
   var_q1 = ncvar_def(paste(var_name,"_",q1_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q1_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   var_q2 = ncvar_def(paste(var_name,"_",q2_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q2_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   var_q3 = ncvar_def(paste(var_name,"_",q3_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q3_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   var_q4 = ncvar_def(paste(var_name,"_",q4_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q4_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   var_q5 = ncvar_def(paste(var_name,"_",q5_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q5_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   var_q6 = ncvar_def(paste(var_name,"_",q6_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q6_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var1,var2,var_mid,var_q1,var_q2,var_q3,var_q4,var_q5,var_q6), force_v4 = TRUE)

   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, analysis_timestep)
   # Grid area 
   ncvar_put(new_file, var1, area_m2)
   # Land fraction
   ncvar_put(new_file, var2, land_fraction)
   # VARIABLE
   ncvar_put(new_file, var_mid, var_est[,,mid_quant,])
   ncvar_put(new_file, var_q1, var_est[,,q1_quant,])
   ncvar_put(new_file, var_q2, var_est[,,q2_quant,])
   ncvar_put(new_file, var_q3, var_est[,,q3_quant,])
   ncvar_put(new_file, var_q4, var_est[,,q4_quant,])   
   ncvar_put(new_file, var_q5, var_est[,,q5_quant,])
   ncvar_put(new_file, var_q6, var_est[,,q6_quant,])   

   # Close the existing file to ensure its written to file
   nc_close(new_file)
   
   # Tidy
   rm(var_q1,var_q2,var_q3,var_q4,var_q5,var_q6,var_mid)

} # end function write_to_nc_analysis_timestep

# write analysis output to netcdf file at the analysis timestep
write_to_nc_analysis_annual<-function(var_est,var_name,var_unit,var_long) {

   # Define the output file name
   output_name = paste(out_dir,output_prefix,"mean_annual_",var_name,output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   # Median
   var_mid = ncvar_def(var_name, unit = var_unit, longname = paste(var_long," - Median estimate",sep=""), 
                      dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   # Quantiles for uncertainty
   var_q1 = ncvar_def(paste(var_name,"_",q1_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q1_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   var_q2 = ncvar_def(paste(var_name,"_",q2_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q2_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   var_q3 = ncvar_def(paste(var_name,"_",q3_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q3_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   var_q4 = ncvar_def(paste(var_name,"_",q4_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q4_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   var_q5 = ncvar_def(paste(var_name,"_",q5_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q5_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   var_q6 = ncvar_def(paste(var_name,"_",q6_quant_lab,sep=""), unit = var_unit, longname = paste(var_long," - ",q6_quant_longlab,sep=""), 
                      dim = list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   
   # Create the empty file space
   new_file = nc_create(filename = output_name, vars = list(var1,var2,var_mid,var_q1,var_q2,var_q3,var_q4,var_q5,var_q6), force_v4 = TRUE)

   # Load first variable into the file
   ## Year information will come with the dimension being used
   # Grid area 
   ncvar_put(new_file, var1, area_m2)
   # Land fraction
   ncvar_put(new_file, var2, land_fraction)
   # VARIABLE
   ncvar_put(new_file, var_mid, var_est[,,mid_quant,])
   ncvar_put(new_file, var_q1, var_est[,,q1_quant,])
   ncvar_put(new_file, var_q2, var_est[,,q2_quant,])
   ncvar_put(new_file, var_q3, var_est[,,q3_quant,])
   ncvar_put(new_file, var_q4, var_est[,,q4_quant,])   
   ncvar_put(new_file, var_q5, var_est[,,q5_quant,])
   ncvar_put(new_file, var_q6, var_est[,,q6_quant,])   

   # Close the existing file to ensure its written to file
   nc_close(new_file)

   # Tidy
   rm(var_q1,var_q2,var_q3,var_q4,var_q5,var_q6,var_mid)

} # end function write_to_nc_analysis_annual

###
## Begin creating information for processing and subsequent saving to files

# Time information
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
steps_per_year = dim(grid_output$lai_m2m2)[3] / nos_years

# create lat / long axes, assumes regular WGS-84 grid
output = determine_lat_long_needed(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution,PROJECT$grid_type,PROJECT$waterpixels)
# NOTE: rev due to CARDAMOM grid being inverse of what comes out of raster function. Should consider changing this at some point.
longitude = output$obs_long_grid[,1] ; latitude = rev(output$obs_lat_grid[1,]) 
# Tidy up
rm(output) ; gc(reset=TRUE,verbose=FALSE)

# Extract the available quantiles 
quantiles_wanted = grid_output$num_quantiles
nos_quantiles = length(quantiles_wanted)
# Check that the quantiles we want to use are available
# Minimum quantiles
if (length(which(quantiles_wanted == 0.025)) == 1) {
    q1_quant = which(quantiles_wanted == 0.025)
    q1_quant_lab = "2.5pc"
    q1_quant_longlab = "2.5 % quantile"
} else {
    stop("Desired min quantile cannot be found")
}
# A lower quantile
if (length(which(quantiles_wanted == 0.05)) == 1) {
    q2_quant = which(quantiles_wanted == 0.05)
    q2_quant_lab = "5pc"
    q2_quant_longlab = "5 % quantile"
} else {
    stop("Desired low quantile cannot be found")
}
# A lower quartile
if (length(which(quantiles_wanted == 0.160)) == 1) {
    q3_quant = which(quantiles_wanted == 0.160)
    q3_quant_lab = "16pc"
    q3_quant_longlab = "16 % quantile"
} else {
    stop("Desired lower quartile cannot be found")
}
# The median estimate
if (length(which(quantiles_wanted == 0.5)) == 1) {
    mid_quant = which(quantiles_wanted == 0.5)
} else {
    stop("Median quantile cannot be found")
}
# A upper quartile
if (length(which(quantiles_wanted == 0.84)) == 1) {
    q4_quant = which(quantiles_wanted == 0.84)
    q4_quant_lab = "84pc"
    q4_quant_longlab = "84 % quantile"
} else {
    stop("Desired upper quartile cannot be found")
}
# A upper quantile
if (length(which(quantiles_wanted == 0.95)) == 1) {
    q5_quant = which(quantiles_wanted == 0.95)
    q5_quant_lab = "95pc"
    q5_quant_longlab = "95 % quantile"
} else {
    stop("Desired high quantile cannot be found")
}
# Maximum quantile
if (length(which(quantiles_wanted == 0.975)) == 1) {
    q6_quant = which(quantiles_wanted == 0.975)
    q6_quant_lab = "97.5pc"
    q6_quant_longlab = "97.5 % quantile"
} else {
    stop("Desired max quantile cannot be found")
}

###
## Create netcdf commons
###

# Extract to final common information
area_m2 = grid_output$area_m2
land_fraction = grid_output$land_fraction
drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[1],".bin",sep=""))
analysis_timestep = drivers$met[,1] ; rm(drivers)

## define dimension
lat_dimen <- ncdim_def( "latitude", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "longitude", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
year_dimen <- ncdim_def( "year", units="", as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
npar_dimen <- ncdim_def( "nos_parameters", units="", 1:(max(PROJECT$model$nopars)+1)) # NOTE: +1 is to account for the log-likelihood 

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), 
                 dim=list(time_dimen), missval = -99999, prec="single", compression = 9)
var1 = ncvar_def("grid_area", units = "m2", longname = paste("Pixel area",sep=""), 
                 dim=list(long_dimen,lat_dimen), missval = -99999, prec="single", compression = 9)
var2 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                 dim=list(long_dimen,lat_dimen), missval = -99999, prec="single", compression = 9)

###
## Process each variable in turn and write to netcdf file
###

###
## Forcings

## Minimum temperature with Celcius -> K adjustment
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 2, bias_adj = 273.15, scale_adj = 1)
write_to_nc_forcing(forcings,var_name = "tas_min",var_unit = "K",var_long = "Monthly average of the daily minimum air temperature") 
## Maximum temperature with Celcius -> K adjustment
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 3, bias_adj = 273.15, scale_adj = 1)
write_to_nc_forcing(forcings,var_name = "tas_max",var_unit = "K",var_long = "Monthly average of the daily maximum air temperature") 
## Average temperature with Celcius -> K adjustment
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 14, bias_adj = 273.15, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "tas",var_unit = "K",var_long = "Monthly average of the daily mean air temperature") 
## Average daily downwelling shortwave radiation with MJ/m2/day -> W/m2
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 4, bias_adj = 0, scale_adj = 1e6*(1/86400)) 
write_to_nc_forcing(forcings,var_name = "rsds",var_unit = "W.m-2",var_long = "Monthly average of the daily downwelling shortwave radiation") 
## Average atmospheric CO2 concentration
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 5, bias_adj = 0, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "co2",var_unit = "ppm",var_long = "Monthly average atmospheric CO2 concentration") 
## Average precipitation rate
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 7, bias_adj = 0, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "pr",var_unit = "kg.m-2.s",var_long = "Monthly average precipitation rate") 
## Harvest fraction
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 8, bias_adj = 0, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "FLOSS_FRAC",var_unit = "fraction",var_long = "Monthly fractional loss of biomass") 
## Burned fraction
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 9, bias_adj = 0, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "burntArea",var_unit = "fraction",var_long = "Monthly fraction burned area") 
## Wind speed
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 15, bias_adj = 0, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "wind_speed",var_unit = "m.s-1",var_long = "Monthly mean wind speed") 
## Vapour pressure deficit
# (Re-)create the object to be updated
forcings = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_met_drivers, met_index = 16, bias_adj = 0, scale_adj = 1) 
write_to_nc_forcing(forcings,var_name = "vpd",var_unit = "Pa",var_long = "Monthly mean vapour pressure deficit") 

#Tidy
rm(forcings,tmp)

###
## Assimilated observations - state variables

## LAI
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 4, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "LAI_OBS", var_unc_name = "LAI_OBS_UNC", var_lag_name = "LAI_OBS_LAG",
                        var_est_unit = "m2.m-2", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of leaf area index", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 
## Wood stocks
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 19, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "WOOD_OBS", var_unc_name = "WOOD_OBS_UNC", var_lag_name = "WOOD_OBS_LAG",
                        var_est_unit = "gC.m-2", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of wood stocks", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 
## Soil stocks
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 28, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "SOIL_OBS", var_unc_name = "SOIL_OBS_UNC", var_lag_name = "SOIL_OBS_LAG",
                        var_est_unit = "gC.m-2", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of soil organic matter stocks", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 
        
###
## Assimilated observations - fluxes variables
                                                                    
## Gross primary production
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 1, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "GPP_OBS", var_unc_name = "GPP_OBS_UNC", var_lag_name = "GPP_OBS_LAG",
                        var_est_unit = "gC.m-2.d-1", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of gross primary production", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 

## Net Ecosystem Exchange of CO2
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 7, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "NEE_OBS", var_unc_name = "NEE_OBS_UNC", var_lag_name = "NEE_OBS_LAG",
                        var_est_unit = "gC.m-2.d-1", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of net ecosystem exchange of CO2", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 

## Net Biome Exchange of CO2
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 46, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "NBE_OBS", var_unc_name = "NBE_OBS_UNC", var_lag_name = "NBE_OBS_LAG",
                        var_est_unit = "gC.m-2.d-1", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of net biomse exchange of CO2", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 
## Fire C emissions
# (Re-)create the object to be updated
est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
lag = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_obs_drivers, obs_index = 10, bias_adj = 0, scale_adj = 1) 
write_to_nc_observation(est,unc,lag,
                        var_est_name = "FIRE_OBS", var_unc_name = "FIRE_OBS_UNC", var_lag_name = "FIRE_OBS_LAG",
                        var_est_unit = "gC.m-2.d-1", var_lag_unit = "time steps",
                        var_est_long = "Assimilated estimate of carbon emissions from fire", 
                        var_unc_long = "Uncertainty of estimate", 
                        var_lag_long = "Number of model time steps over which the estmate is assumed to represent the average") 

# Tidy
rm(est,unc,lag,tmp)

###
## Analysis outputs at model timestep

## LAI
if (exists(x = "lai_m2m2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "lai_m2m2", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_timestep(est, var_name = "lai", var_unit = "m2.m-2", var_long = "Leaf Area Index") 
}
## Total C
if (exists(x = "Ctotal_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "Ctotal_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cTotal", var_unit = "kg.m-2", var_long = "Carbon in live and dead organic matter") 
}
## Labile
if (exists(x = "labile_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "labile_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cLabile", var_unit = "kg.m-2", var_long = "Carbon in labile")  
}
## Foliage
if (exists(x = "foliage_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "foliage_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cLeaf", var_unit = "kg.m-2", var_long = "Carbon in leaves")  
}
## Fine roots
if (exists(x = "roots_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "roots_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cRoot", var_unit = "kg.m-2", var_long = "Carbon in fine root")  
}
## Wood
if (exists(x = "wood_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "wood_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cWoodTotal", var_unit = "kg.m-2", var_long = "Carbon in (AGB + BGB) wood")  
}
## Fine litter (foliage + fine root)
if (exists(x = "litter_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "litter_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cLitter", var_unit = "kg.m-2", var_long = "Carbon in (Foliar + fine root)")  
}
## Soil Organic Matter
if (exists(x = "som_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "som_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cSoil", var_unit = "kg.m-2", var_long = "Carbon in soil organic matter (0-1m)")  
}
## Wood litter
if (exists(x = "woodlitter_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "woodlitter_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cCwod", var_unit = "kg.m-2", var_long = "Carbon in (wood) litter")   
}
## Dead Organic Matter
if (exists(x = "dom_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "dom_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cDOM", var_unit = "kg.m-2", var_long = "Carbon in leaf, fine root, wood litter, and soil organic matter")  
}
## Biomass
if (exists(x = "biomass_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "biomass_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "cVeg", var_unit = "kg.m-2", var_long = "Carbon in live biomass")  
} 
## Biomass change since time step 1
if (exists(x = "dCbiomass_gCm2", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "dCbiomass_gCm2", bias_adj = 0, scale_adj = 1e-3) 
    write_to_nc_analysis_timestep(est, var_name = "dcVeg", var_unit = "kg.m-2", var_long = "Change in Carbon in live biomass since t=1")  
}
## GPP gC/m2/day -> kgC/m2/s
if (exists(x = "gpp_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "gpp_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "gpp", var_unit = "kg.m-2.s-1", var_long = "Gross Primary Productivity")  
}
## Rauto gC/m2/day -> kgC/m2/s
if (exists(x = "rauto_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "rauto_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "ra", var_unit = "kg.m-2.s-1", var_long = "Autotrophic (Plant) Respiration")  
}
## Rhet gC/m2/day -> kgC/m2/s
if (exists(x = "rhet_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "rhet_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "rh", var_unit = "kg.m-2.s-1", var_long = "Heterotrophic Respiration")  
}
## NPP gC/m2/day -> kgC/m2/s
if (exists(x = "npp_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "npp_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "npp", var_unit = "kg.m-2.s-1", var_long = "Net Primary Productivity")  
}
## Fire gC/m2/day -> kgC/m2/s
if (exists(x = "fire_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "fire_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fFire", var_unit = "kg.m-2.s-1", var_long = "Fire C emission")  
}
## Harvest gC/m2/day -> kgC/m2/s
if (exists(x = "harvest_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "harvest_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fLuc", var_unit = "kg.m-2.s-1", var_long = "C extracted due to forest harvest")  
}
## Reco gC/m2/day -> kgC/m2/s
if (exists(x = "reco_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "reco_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "reco", var_unit = "kg.m-2.s-1", var_long = "Ecosystem (Ra + Rh) Respiration")  
}
## NEE gC/m2/day -> kgC/m2/s
if (exists(x = "nee_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "nee_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "nee", var_unit = "kg.m-2.s-1", var_long = "Net Ecosystem Exchange")
}
## NBE gC/m2/day -> kgC/m2/s
if (exists(x = "nbe_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "nbe_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "nbe", var_unit = "kg.m-2.s-1", var_long = "Net Biome Exchange (NEE + Fire)")
}
## NBP gC/m2/day -> kgC/m2/s
if (exists(x = "nbp_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "nbp_gCm2day", bias_adj = 0, scale_adj = 1e-3*(1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "nbp", var_unit = "kg.m-2.s-1", var_long = "Net Biome Productivity (-NEE - Fire - fLuc)")
}
## Evapotranspiration kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "ET_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "ET_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "evapotrans", var_unit = "kg.m-2.s-1", var_long = "Evapotranspiration")
}
## Transpiration kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "Etrans_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "Etrans_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "tran", var_unit = "kg.m-2.s-1", var_long = "Transpiration")
}
## Soil evaporation kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "Esoil_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "Esoil_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "evspsblsoi", var_unit = "kg.m-2.s-1", var_long = "Soil evaporation")
}
## Wet canopy evaporation kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "Ewetcanopy_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "Ewetcanopy_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "evspsblveg", var_unit = "kg.m-2.s-1", var_long = "Canopy intercepted rainfall evaporation")
}
## Total drainage (surface runoff + underflow) kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "total_drainage_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "total_drainage_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "mrro", var_unit = "kg.m-2.s-1", var_long = "Total drainage from soil surface and bottom of soil column")
}
## Surface runoff kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "runoff_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "runoff_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "runoff", var_unit = "kg.m-2.s-1", var_long = "Soil surface water runoff")
}
## Underflow kgH2O/m2/day -> kgH2O/m2/s
if (exists(x = "underflow_kgH2Om2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "underflow_kgH2Om2day", bias_adj = 0, scale_adj = (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "underflow", var_unit = "kg.m-2.s-1", var_long = "Water drainage from the bottom of the soil column")
}
## Biomass to litter gC/m2/day -> kgC/m2/s
if (exists(x = "combined_biomass_to_litter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_biomass_to_litter_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fVegLitter", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven litter creation from biomass")
}
## Labile to litter gC/m2/day -> kgC/m2/s
if (exists(x = "combined_labile_to_litter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_labile_to_litter_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fLabileLitter", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven litter creation from labile")
}
## Foliage to litter gC/m2/day -> kgC/m2/s
if (exists(x = "combined_foliage_to_litter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_foliage_to_litter_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fLeafLitter", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven litter creation from foliage")
}
## Fine roots to litter gC/m2/day -> kgC/m2/s
if (exists(x = "combined_roots_to_litter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_roots_to_litter_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fRootLitter", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven litter creation from fine root")
}
## Wood to litter gC/m2/day -> kgC/m2/s
if (exists(x = "combined_wood_to_litter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_wood_to_litter_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fVegSoil", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven litter creation from wood, which is allocated to som")
}
## Litter to soil gC/m2/day -> kgC/m2/s
if (exists(x = "combined_litter_to_som_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_litter_to_som_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fLitterSoil", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven allocation of litter to soil")
}
## Wood litter to soil gC/m2/day -> kgC/m2/s
if (exists(x = "combined_woodlitter_to_som_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_woodlitter_to_som_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fCwdSoil", var_unit = "kg.m-2.s-1", var_long = "Combined natural, fire and harvest driven allocation of wood litter to soil")
}
## C allocation to fine roots gC/m2/day -> kgC/m2/s
if (exists(x = "alloc_roots_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "alloc_roots_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fAllocRoot", var_unit = "kg.m-2.s-1", var_long = "Net Primary Productivity to fine root")
}
## C allocation to wood gC/m2/day -> kgC/m2/s
if (exists(x = "alloc_wood_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "alloc_wood_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fAllocWood", var_unit = "kg.m-2.s-1", var_long = "Net Primary Productivity to wood")
}
## C allocation to foliage gC/m2/day -> kgC/m2/s
if (exists(x = "combined_alloc_foliage_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "combined_alloc_foliage_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fAllocLeaf", var_unit = "kg.m-2.s-1", var_long = "Both direct and via labile Net Primary Productivity to foliage")
}
## Fire combusted CO2 output flux from foliar and fine root litter gC/m2/day -> kgCO2/m2/s
if (exists(x = "FIREemiss_litter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "FIREemiss_litter_gCm2day", bias_adj = 0, scale_adj = (44/12) * 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fFireLitter", var_unit = "kg.m-2.s-1", var_long = "Fire combusted CO2 output flux from foliar and fine root litter")
}
## Fire combusted CO2 output flux from wood litter gC/m2/day -> kgCO2/m2/s
if (exists(x = "FIREemiss_woodlitter_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "FIREemiss_woodlitter_gCm2day", bias_adj = 0, scale_adj = (44/12) * 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fFireCcwd", var_unit = "kg.m-2.s-1", var_long = "Fire combusted CO2 output flux from wood litter")
}
## Fire combusted CO2 output flux from soil gC/m2/day -> kgCO2/m2/s
if (exists(x = "FIREemiss_som_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "FIREemiss_som_gCm2day", bias_adj = 0, scale_adj = (44/12) * 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fFireCSoil", var_unit = "kg.m-2.s-1", var_long = "Fire combusted CO2 output flux from soil organic matter")
}
## Fire combusted CO2 output flux from biomass gC/m2/day -> kgCO2/m2/s
if (exists(x = "FIREemiss_biomass_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "FIREemiss_biomass_gCm2day", bias_adj = 0, scale_adj = (44/12) * 1e-3 * (1/86400)) 
    write_to_nc_analysis_timestep(est, var_name = "fFireCveg", var_unit = "kg.m-2.s-1", var_long = "Fire combusted CO2 output flux from vegetation")
}

###
## Analysis outputs at annual timestep

## Annual gross primary productivity gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_gpp_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_gpp_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "gpp", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Gross Primary Productivity")
}
## Annual autotrophic respiration gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_rauto_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_rauto_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "ra", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Autotrophic (plant) respiration")
}
## Annual heterotrophic respiration gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_rhet_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_rhet_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "rh", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Heterotrophic respiration")
}
## Annual NPP gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_npp_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_npp_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "npp", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Net Primary Productivity")
}
## Annual fire C emissions gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_fire_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_fire_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "fFire", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Fire C emission")
} 
## Annual harvest gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_harvest_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_harvest_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "fLuc", var_unit = "kg.m-2.s-1", var_long = "Mean Annual C extracted due to forest harvest")
}  
## Annual ecosystem respiration gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_reco_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_reco_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "reco", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Ecosystem (Ra + Rh) Respiration")
}   
## Annual net ecosystem exchange gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_nee_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_nee_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "nee", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Net Ecosystem Exchange")
}   
## Annual net biome exchange gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_nbe_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_nbe_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "nbe", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Net Biome Exchange (NEE + Fire)")
}    
## Annual net biome productivity gC/m2/day -> kgC/m2/s
if (exists(x = "mean_annual_nbp_gCm2day", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_nbp_gCm2day", bias_adj = 0, scale_adj = 1e-3 * (1/86400)) 
    write_to_nc_analysis_annual(est, var_name = "nbp", var_unit = "kg.m-2.s-1", var_long = "Mean Annual Net Biome Productivity (-NEE - Fire - fLuc)")
}    
## Annual carbon use efficiency (gC/gC)
if (exists(x = "mean_annual_cue", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "mean_annual_cue", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "cue", var_unit = "0-1", var_long = "Mean Annual Carbon Use Efficiency")
}    
## Annual mean transit times for biomass (years)
if (exists(x = "MTT_annual_biomass_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_biomass_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_biomass", var_unit = "years", var_long = "Annual mean transit (residence) time of the biomass C pool")
}    
## Annual mean transit times for labile (years)
if (exists(x = "MTT_annual_labile_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_labile_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_labile", var_unit = "years", var_long = "Annual mean transit (residence) time of the labile C pool")
}       
## Annual mean transit times for foliage (years)
if (exists(x = "MTT_annual_foliage_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_foliage_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_foliage", var_unit = "years", var_long = "Annual mean transit (residence) time of the foliar C pool")
}    
## Annual mean transit times for fine roots (years)
if (exists(x = "MTT_annual_roots_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_roots_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_roots", var_unit = "years", var_long = "Annual mean transit (residence) time of the root C pool")
}    
## Annual mean transit times for wood (years)
if (exists(x = "MTT_annual_wood_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_wood_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_wood", var_unit = "years", var_long = "Annual mean transit (residence) time of the wood C pool")
}    
## Annual mean transit times for fine litter (years)
if (exists(x = "MTT_annual_litter_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_litter_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_litter", var_unit = "years", var_long = "Annual mean transit (residence) time of the fine (foliage + fine root) litter C pool")
}    
## Annual mean transit times for wood litter (years)
if (exists(x = "MTT_annual_woodlitter_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_woodlitter_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_woodlitter", var_unit = "years", var_long = "Annual mean transit (residence) time of the coarse (wood) litter C pool")
}    
## Annual mean transit times for soil organic matter (years)
if (exists(x = "MTT_annual_som_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_som_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_som", var_unit = "years", var_long = "Annual mean transit (residence) time of the soil organic matter C pool")
}    
## Annual mean transit times for dead organic matter (years)
if (exists(x = "MTT_annual_dom_years", where = grid_output)) {
    # (Re-)create the object to be updated
    est = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
    tmp = sapply(c(1:PROJECT$nosites), FUN = extract_from_grid_output, var_name = "MTT_annual_dom_years", bias_adj = 0, scale_adj = 1) 
    write_to_nc_analysis_annual(est, var_name = "MTT_dom", var_unit = "years", var_long = "Annual mean transit (residence) time of the dead organic matter C pool")
}    

# Tidy away the main objects to save memory
rm(est,tmp,grid_output) ; gc()


