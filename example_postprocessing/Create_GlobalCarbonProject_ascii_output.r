
###
## Process CARDAMOM into Global Carbon Project C-budget domains.
## This file will create an ascii file containing the C-budget terms
## for the globe, tropics, north (>N30 degree) and south (>S30 degree).
## Terms are net biome productivity, including all simulated terms.
## Units PgC/yr, annual estimates. 
## The median and propogated uncertainty information are provided.
## Created 29/07/2022 by T. Luke Smallman
###

###
## Job specific information

print("Begin creation of GCP compatible ascii file for annual NBP / GPP / RECO / RHET / RAUTO / FIRE")

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_CsomPriorNCSCD/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_NBP_CsomPriorNCSDC3m/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Trendyv9_historical/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire_NBP/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Mexico_1deg_C7_agb_lca_gpp_fire_NBP/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Miombo_0.5deg_allWood/"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/global_1deg_C7_GCP_LCA_AGB/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/global_1deg_C7_GCP_LCA_AGB_GPP_FIRE/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/global_1deg_C7_GCP_LCA_AGB_etol_EQF_harsh/"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/global_1deg_C7_GCP_LCA_AGB_ACM2_LUE/"

# Specify any extra information for the filename
output_prefix = "CARDAMOM_S3_" # follow with "_"
#output_prefix = "CARDAMOM_S2_" # follow with "_"
output_suffix = "" # begin with "_"

###
## Load libraries, functions and data needed

# load needed libraries
library(ncdf4)
library(raster)
library(compiler)

# Load needed functions 
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/load_all_cardamom_functions.r")

# load the CARDAMOM files
load(paste(input_dir,"/infofile.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

# Create output directory to aid storage management
out_dir = paste(PROJECT$results_processedpath,"/trendy_output",sep="")
if (dir.exists(out_dir) == FALSE) {
    dir.create(out_dir)
}

###
## Begin creating information for processing and subsequent saving to files

# Time information
years = c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
nos_years = length(years)
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
if (length(which(quantiles_wanted == 0.025)) == 1) {
    min_quant = which(quantiles_wanted == 0.025)
} else {
    stop("Desired min quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.25)) == 1) {
    low_quant = which(quantiles_wanted == 0.25)
} else {
    stop("Desired low quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.5)) == 1) {
    mid_quant = which(quantiles_wanted == 0.5)
} else {
    stop("Median quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.75)) == 1) {
    high_quant = which(quantiles_wanted == 0.75)
} else {
    stop("Desired high quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.975)) == 1) {
    max_quant = which(quantiles_wanted == 0.975)
} else {
    stop("Desired max quantile cannot be found")
}

###
## Loop through grid to aggregate the net biome productivity (NBP = GPP - Ra - Rh - Fire - Harvest)

# Median estimate
NBP_global_PgCyr = rep(0, nos_years)
NBP_north_PgCyr = rep(0, nos_years)
NBP_tropics_PgCyr = rep(0, nos_years)
NBP_south_PgCyr = rep(0, nos_years)
GPP_global_PgCyr = rep(0, nos_years)
GPP_north_PgCyr = rep(0, nos_years)
GPP_tropics_PgCyr = rep(0, nos_years)
GPP_south_PgCyr = rep(0, nos_years)
RECO_global_PgCyr = rep(0, nos_years)
RECO_north_PgCyr = rep(0, nos_years)
RECO_tropics_PgCyr = rep(0, nos_years)
RECO_south_PgCyr = rep(0, nos_years)
RHET_global_PgCyr = rep(0, nos_years)
RHET_north_PgCyr = rep(0, nos_years)
RHET_tropics_PgCyr = rep(0, nos_years)
RHET_south_PgCyr = rep(0, nos_years)
RAUTO_global_PgCyr = rep(0, nos_years)
RAUTO_north_PgCyr = rep(0, nos_years)
RAUTO_tropics_PgCyr = rep(0, nos_years)
RAUTO_south_PgCyr = rep(0, nos_years)
FIRE_global_PgCyr = rep(0, nos_years)
FIRE_north_PgCyr = rep(0, nos_years)
FIRE_tropics_PgCyr = rep(0, nos_years)
FIRE_south_PgCyr = rep(0, nos_years)

# Minimum CI
NBP_global_PgCyr_minCI = rep(0, nos_years)
NBP_north_PgCyr_minCI = rep(0, nos_years)
NBP_tropics_PgCyr_minCI = rep(0, nos_years)
NBP_south_PgCyr_minCI = rep(0, nos_years)
GPP_global_PgCyr_minCI = rep(0, nos_years)
GPP_north_PgCyr_minCI = rep(0, nos_years)
GPP_tropics_PgCyr_minCI = rep(0, nos_years)
GPP_south_PgCyr_minCI = rep(0, nos_years)
RECO_global_PgCyr_minCI = rep(0, nos_years)
RECO_north_PgCyr_minCI = rep(0, nos_years)
RECO_tropics_PgCyr_minCI = rep(0, nos_years)
RECO_south_PgCyr_minCI = rep(0, nos_years)
RHET_global_PgCyr_minCI = rep(0, nos_years)
RHET_north_PgCyr_minCI = rep(0, nos_years)
RHET_tropics_PgCyr_minCI = rep(0, nos_years)
RHET_south_PgCyr_minCI = rep(0, nos_years)
RAUTO_global_PgCyr_minCI = rep(0, nos_years)
RAUTO_north_PgCyr_minCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_minCI = rep(0, nos_years)
RAUTO_south_PgCyr_minCI = rep(0, nos_years)
FIRE_global_PgCyr_minCI = rep(0, nos_years)
FIRE_north_PgCyr_minCI = rep(0, nos_years)
FIRE_tropics_PgCyr_minCI = rep(0, nos_years)
FIRE_south_PgCyr_minCI = rep(0, nos_years)
# Lower CI
NBP_global_PgCyr_lowCI = rep(0, nos_years)
NBP_north_PgCyr_lowCI = rep(0, nos_years)
NBP_tropics_PgCyr_lowCI = rep(0, nos_years)
NBP_south_PgCyr_lowCI = rep(0, nos_years)
GPP_global_PgCyr_lowCI = rep(0, nos_years)
GPP_north_PgCyr_lowCI = rep(0, nos_years)
GPP_tropics_PgCyr_lowCI = rep(0, nos_years)
GPP_south_PgCyr_lowCI = rep(0, nos_years)
RECO_global_PgCyr_lowCI = rep(0, nos_years)
RECO_north_PgCyr_lowCI = rep(0, nos_years)
RECO_tropics_PgCyr_lowCI = rep(0, nos_years)
RECO_south_PgCyr_lowCI = rep(0, nos_years)
RHET_global_PgCyr_lowCI = rep(0, nos_years)
RHET_north_PgCyr_lowCI = rep(0, nos_years)
RHET_tropics_PgCyr_lowCI = rep(0, nos_years)
RHET_south_PgCyr_lowCI = rep(0, nos_years)
RAUTO_global_PgCyr_lowCI = rep(0, nos_years)
RAUTO_north_PgCyr_lowCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_lowCI = rep(0, nos_years)
RAUTO_south_PgCyr_lowCI = rep(0, nos_years)
FIRE_global_PgCyr_lowCI = rep(0, nos_years)
FIRE_north_PgCyr_lowCI = rep(0, nos_years)
FIRE_tropics_PgCyr_lowCI = rep(0, nos_years)
FIRE_south_PgCyr_lowCI = rep(0, nos_years)
# Upper CI
NBP_global_PgCyr_highCI = rep(0, nos_years)
NBP_north_PgCyr_highCI = rep(0, nos_years)
NBP_tropics_PgCyr_highCI = rep(0, nos_years)
NBP_south_PgCyr_highCI = rep(0, nos_years)
GPP_global_PgCyr_highCI = rep(0, nos_years)
GPP_north_PgCyr_highCI = rep(0, nos_years)
GPP_tropics_PgCyr_highCI = rep(0, nos_years)
GPP_south_PgCyr_highCI = rep(0, nos_years)
RECO_global_PgCyr_highCI = rep(0, nos_years)
RECO_north_PgCyr_highCI = rep(0, nos_years)
RECO_tropics_PgCyr_highCI = rep(0, nos_years)
RECO_south_PgCyr_highCI = rep(0, nos_years)
RHET_global_PgCyr_highCI = rep(0, nos_years)
RHET_north_PgCyr_highCI = rep(0, nos_years)
RHET_tropics_PgCyr_highCI = rep(0, nos_years)
RHET_south_PgCyr_highCI = rep(0, nos_years)
RAUTO_global_PgCyr_highCI = rep(0, nos_years)
RAUTO_north_PgCyr_highCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_highCI = rep(0, nos_years)
RAUTO_south_PgCyr_highCI = rep(0, nos_years)
FIRE_global_PgCyr_highCI = rep(0, nos_years)
FIRE_north_PgCyr_highCI = rep(0, nos_years)
FIRE_tropics_PgCyr_highCI = rep(0, nos_years)
FIRE_south_PgCyr_highCI = rep(0, nos_years)
# Maximum CI
NBP_global_PgCyr_maxCI = rep(0, nos_years)
NBP_north_PgCyr_maxCI = rep(0, nos_years)
NBP_tropics_PgCyr_maxCI = rep(0, nos_years)
NBP_south_PgCyr_maxCI = rep(0, nos_years)
GPP_global_PgCyr_maxCI = rep(0, nos_years)
GPP_north_PgCyr_maxCI = rep(0, nos_years)
GPP_tropics_PgCyr_maxCI = rep(0, nos_years)
GPP_south_PgCyr_maxCI = rep(0, nos_years)
RECO_global_PgCyr_maxCI = rep(0, nos_years)
RECO_north_PgCyr_maxCI = rep(0, nos_years)
RECO_tropics_PgCyr_maxCI = rep(0, nos_years)
RECO_south_PgCyr_maxCI = rep(0, nos_years)
RHET_global_PgCyr_maxCI = rep(0, nos_years)
RHET_north_PgCyr_maxCI = rep(0, nos_years)
RHET_tropics_PgCyr_maxCI = rep(0, nos_years)
RHET_south_PgCyr_maxCI = rep(0, nos_years)
RAUTO_global_PgCyr_maxCI = rep(0, nos_years)
RAUTO_north_PgCyr_maxCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_maxCI = rep(0, nos_years)
RAUTO_south_PgCyr_maxCI = rep(0, nos_years)
FIRE_global_PgCyr_maxCI = rep(0, nos_years)
FIRE_north_PgCyr_maxCI = rep(0, nos_years)
FIRE_tropics_PgCyr_maxCI = rep(0, nos_years)
FIRE_south_PgCyr_maxCI = rep(0, nos_years)

# Fill the output arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]
         
         # Determine nos days per time step         
         deltat = 365.25

         # Accumulate global
         NBP_global_PgCyr = NBP_global_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_minCI = NBP_global_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_lowCI = NBP_global_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_highCI = NBP_global_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_maxCI = NBP_global_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         GPP_global_PgCyr = GPP_global_PgCyr + (grid_output$mean_annual_gpp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         GPP_global_PgCyr_minCI = GPP_global_PgCyr_minCI + (grid_output$mean_annual_gpp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         GPP_global_PgCyr_lowCI = GPP_global_PgCyr_lowCI + (grid_output$mean_annual_gpp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         GPP_global_PgCyr_highCI = GPP_global_PgCyr_highCI + (grid_output$mean_annual_gpp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         GPP_global_PgCyr_maxCI = GPP_global_PgCyr_maxCI + (grid_output$mean_annual_gpp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RECO_global_PgCyr = RECO_global_PgCyr + (grid_output$mean_annual_reco_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RECO_global_PgCyr_minCI = RECO_global_PgCyr_minCI + (grid_output$mean_annual_reco_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RECO_global_PgCyr_lowCI = RECO_global_PgCyr_lowCI + (grid_output$mean_annual_reco_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RECO_global_PgCyr_highCI = RECO_global_PgCyr_highCI + (grid_output$mean_annual_reco_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RECO_global_PgCyr_maxCI = RECO_global_PgCyr_maxCI + (grid_output$mean_annual_reco_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_global_PgCyr = RHET_global_PgCyr + (grid_output$mean_annual_rhet_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_global_PgCyr_minCI = RHET_global_PgCyr_minCI + (grid_output$mean_annual_rhet_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_global_PgCyr_lowCI = RHET_global_PgCyr_lowCI + (grid_output$mean_annual_rhet_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_global_PgCyr_highCI = RHET_global_PgCyr_highCI + (grid_output$mean_annual_rhet_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_global_PgCyr_maxCI = RHET_global_PgCyr_maxCI + (grid_output$mean_annual_rhet_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RAUTO_global_PgCyr = RAUTO_global_PgCyr + (grid_output$mean_annual_rauto_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RAUTO_global_PgCyr_minCI = RAUTO_global_PgCyr_minCI + (grid_output$mean_annual_rauto_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RAUTO_global_PgCyr_lowCI = RAUTO_global_PgCyr_lowCI + (grid_output$mean_annual_rauto_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RAUTO_global_PgCyr_highCI = RAUTO_global_PgCyr_highCI + (grid_output$mean_annual_rauto_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RAUTO_global_PgCyr_maxCI = RAUTO_global_PgCyr_maxCI + (grid_output$mean_annual_rauto_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FIRE_global_PgCyr = FIRE_global_PgCyr + (grid_output$mean_annual_fire_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FIRE_global_PgCyr_minCI = FIRE_global_PgCyr_minCI + (grid_output$mean_annual_fire_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FIRE_global_PgCyr_lowCI = FIRE_global_PgCyr_lowCI + (grid_output$mean_annual_fire_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FIRE_global_PgCyr_highCI = FIRE_global_PgCyr_highCI + (grid_output$mean_annual_fire_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FIRE_global_PgCyr_maxCI = FIRE_global_PgCyr_maxCI + (grid_output$mean_annual_fire_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
                  
         # Where appropriate accumulate tropical
         if (grid_output$lat[i,j] <= 30 & grid_output$lat[i,j] >= -30) {
             NBP_tropics_PgCyr = NBP_tropics_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_minCI = NBP_tropics_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_lowCI = NBP_tropics_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_highCI = NBP_tropics_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_maxCI = NBP_tropics_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_tropics_PgCyr = GPP_tropics_PgCyr + (grid_output$mean_annual_gpp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_tropics_PgCyr_minCI = GPP_tropics_PgCyr_minCI + (grid_output$mean_annual_gpp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_tropics_PgCyr_lowCI = GPP_tropics_PgCyr_lowCI + (grid_output$mean_annual_gpp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_tropics_PgCyr_highCI = GPP_tropics_PgCyr_highCI + (grid_output$mean_annual_gpp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_tropics_PgCyr_maxCI = GPP_tropics_PgCyr_maxCI + (grid_output$mean_annual_gpp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_tropics_PgCyr = RECO_tropics_PgCyr + (grid_output$mean_annual_reco_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_tropics_PgCyr_minCI = RECO_tropics_PgCyr_minCI + (grid_output$mean_annual_reco_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_tropics_PgCyr_lowCI = RECO_tropics_PgCyr_lowCI + (grid_output$mean_annual_reco_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_tropics_PgCyr_highCI = RECO_tropics_PgCyr_highCI + (grid_output$mean_annual_reco_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_tropics_PgCyr_maxCI = RECO_tropics_PgCyr_maxCI + (grid_output$mean_annual_reco_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_tropics_PgCyr = RHET_tropics_PgCyr + (grid_output$mean_annual_rhet_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_tropics_PgCyr_minCI = RHET_tropics_PgCyr_minCI + (grid_output$mean_annual_rhet_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_tropics_PgCyr_lowCI = RHET_tropics_PgCyr_lowCI + (grid_output$mean_annual_rhet_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_tropics_PgCyr_highCI = RHET_tropics_PgCyr_highCI + (grid_output$mean_annual_rhet_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_tropics_PgCyr_maxCI = RHET_tropics_PgCyr_maxCI + (grid_output$mean_annual_rhet_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_tropics_PgCyr = RAUTO_tropics_PgCyr + (grid_output$mean_annual_rauto_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_tropics_PgCyr_minCI = RAUTO_tropics_PgCyr_minCI + (grid_output$mean_annual_rauto_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_tropics_PgCyr_lowCI = RAUTO_tropics_PgCyr_lowCI + (grid_output$mean_annual_rauto_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_tropics_PgCyr_highCI = RAUTO_tropics_PgCyr_highCI + (grid_output$mean_annual_rauto_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_tropics_PgCyr_maxCI = RAUTO_tropics_PgCyr_maxCI + (grid_output$mean_annual_rauto_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_tropics_PgCyr = FIRE_tropics_PgCyr + (grid_output$mean_annual_fire_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_tropics_PgCyr_minCI = FIRE_tropics_PgCyr_minCI + (grid_output$mean_annual_fire_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_tropics_PgCyr_lowCI = FIRE_tropics_PgCyr_lowCI + (grid_output$mean_annual_fire_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_tropics_PgCyr_highCI = FIRE_tropics_PgCyr_highCI + (grid_output$mean_annual_fire_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_tropics_PgCyr_maxCI = FIRE_tropics_PgCyr_maxCI + (grid_output$mean_annual_fire_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         }

         # Where appropriate accumulate north
         if (grid_output$lat[i,j] > 30) {
             NBP_north_PgCyr = NBP_north_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_minCI = NBP_north_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_lowCI = NBP_north_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_highCI = NBP_north_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_maxCI = NBP_north_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_north_PgCyr = GPP_north_PgCyr + (grid_output$mean_annual_gpp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_north_PgCyr_minCI = GPP_north_PgCyr_minCI + (grid_output$mean_annual_gpp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_north_PgCyr_lowCI = GPP_north_PgCyr_lowCI + (grid_output$mean_annual_gpp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_north_PgCyr_highCI = GPP_north_PgCyr_highCI + (grid_output$mean_annual_gpp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_north_PgCyr_maxCI = GPP_north_PgCyr_maxCI + (grid_output$mean_annual_gpp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_north_PgCyr = RECO_north_PgCyr + (grid_output$mean_annual_reco_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_north_PgCyr_minCI = RECO_north_PgCyr_minCI + (grid_output$mean_annual_reco_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_north_PgCyr_lowCI = RECO_north_PgCyr_lowCI + (grid_output$mean_annual_reco_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_north_PgCyr_highCI = RECO_north_PgCyr_highCI + (grid_output$mean_annual_reco_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_north_PgCyr_maxCI = RECO_north_PgCyr_maxCI + (grid_output$mean_annual_reco_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_north_PgCyr = RHET_north_PgCyr + (grid_output$mean_annual_rhet_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_north_PgCyr_minCI = RHET_north_PgCyr_minCI + (grid_output$mean_annual_rhet_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_north_PgCyr_lowCI = RHET_north_PgCyr_lowCI + (grid_output$mean_annual_rhet_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_north_PgCyr_highCI = RHET_north_PgCyr_highCI + (grid_output$mean_annual_rhet_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_north_PgCyr_maxCI = RHET_north_PgCyr_maxCI + (grid_output$mean_annual_rhet_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_north_PgCyr = RAUTO_north_PgCyr + (grid_output$mean_annual_rauto_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_north_PgCyr_minCI = RAUTO_north_PgCyr_minCI + (grid_output$mean_annual_rauto_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_north_PgCyr_lowCI = RAUTO_north_PgCyr_lowCI + (grid_output$mean_annual_rauto_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_north_PgCyr_highCI = RAUTO_north_PgCyr_highCI + (grid_output$mean_annual_rauto_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_north_PgCyr_maxCI = RAUTO_north_PgCyr_maxCI + (grid_output$mean_annual_rauto_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_north_PgCyr = FIRE_north_PgCyr + (grid_output$mean_annual_fire_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_north_PgCyr_minCI = FIRE_north_PgCyr_minCI + (grid_output$mean_annual_fire_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_north_PgCyr_lowCI = FIRE_north_PgCyr_lowCI + (grid_output$mean_annual_fire_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_north_PgCyr_highCI = FIRE_north_PgCyr_highCI + (grid_output$mean_annual_fire_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_north_PgCyr_maxCI = FIRE_north_PgCyr_maxCI + (grid_output$mean_annual_fire_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)         
         }

         # Where appropriate accumulate south
         if (grid_output$lat[i,j] < -30) {
             NBP_south_PgCyr = NBP_south_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_minCI = NBP_south_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_lowCI = NBP_south_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_highCI = NBP_south_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_maxCI = NBP_south_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_south_PgCyr = GPP_south_PgCyr + (grid_output$mean_annual_gpp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_south_PgCyr_minCI = GPP_south_PgCyr_minCI + (grid_output$mean_annual_gpp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_south_PgCyr_lowCI = GPP_south_PgCyr_lowCI + (grid_output$mean_annual_gpp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_south_PgCyr_highCI = GPP_south_PgCyr_highCI + (grid_output$mean_annual_gpp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             GPP_south_PgCyr_maxCI = GPP_south_PgCyr_maxCI + (grid_output$mean_annual_gpp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_south_PgCyr = RECO_south_PgCyr + (grid_output$mean_annual_reco_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_south_PgCyr_minCI = RECO_south_PgCyr_minCI + (grid_output$mean_annual_reco_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_south_PgCyr_lowCI = RECO_south_PgCyr_lowCI + (grid_output$mean_annual_reco_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_south_PgCyr_highCI = RECO_south_PgCyr_highCI + (grid_output$mean_annual_reco_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RECO_south_PgCyr_maxCI = RECO_south_PgCyr_maxCI + (grid_output$mean_annual_reco_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_south_PgCyr = RHET_south_PgCyr + (grid_output$mean_annual_rhet_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_south_PgCyr_minCI = RHET_south_PgCyr_minCI + (grid_output$mean_annual_rhet_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_south_PgCyr_lowCI = RHET_south_PgCyr_lowCI + (grid_output$mean_annual_rhet_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_south_PgCyr_highCI = RHET_south_PgCyr_highCI + (grid_output$mean_annual_rhet_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_south_PgCyr_maxCI = RHET_south_PgCyr_maxCI + (grid_output$mean_annual_rhet_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_south_PgCyr = RAUTO_south_PgCyr + (grid_output$mean_annual_rauto_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_south_PgCyr_minCI = RAUTO_south_PgCyr_minCI + (grid_output$mean_annual_rauto_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_south_PgCyr_lowCI = RAUTO_south_PgCyr_lowCI + (grid_output$mean_annual_rauto_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_south_PgCyr_highCI = RAUTO_south_PgCyr_highCI + (grid_output$mean_annual_rauto_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RAUTO_south_PgCyr_maxCI = RAUTO_south_PgCyr_maxCI + (grid_output$mean_annual_rauto_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_south_PgCyr = FIRE_south_PgCyr + (grid_output$mean_annual_fire_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_south_PgCyr_minCI = FIRE_south_PgCyr_minCI + (grid_output$mean_annual_fire_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_south_PgCyr_lowCI = FIRE_south_PgCyr_lowCI + (grid_output$mean_annual_fire_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_south_PgCyr_highCI = FIRE_south_PgCyr_highCI + (grid_output$mean_annual_fire_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FIRE_south_PgCyr_maxCI = FIRE_south_PgCyr_maxCI + (grid_output$mean_annual_fire_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         }

     } # Does the file exist / has it been processed

} # site loop

###
## Combine NBP together the data into an output variable

output = data.frame(Year = years, Global = NBP_global_PgCyr, North = NBP_north_PgCyr, Tropics = NBP_tropics_PgCyr, South = NBP_south_PgCyr,
                                  Global_2.5pc = NBP_global_PgCyr_minCI, North_2.5pc = NBP_north_PgCyr_minCI, Tropics_2.5pc = NBP_tropics_PgCyr_minCI, South_2.5pc = NBP_south_PgCyr_minCI,
                                  Global_25pc = NBP_global_PgCyr_lowCI, North_25pc = NBP_north_PgCyr_lowCI, Tropics_25pc = NBP_tropics_PgCyr_lowCI, South_25pc = NBP_south_PgCyr_lowCI,
                                  Global_75pc = NBP_global_PgCyr_highCI, North_75pc = NBP_north_PgCyr_highCI, Tropics_75pc = NBP_tropics_PgCyr_highCI, South_75pc = NBP_south_PgCyr_highCI,
                                  Global_97.5pc = NBP_global_PgCyr_maxCI, North_97.5pc = NBP_north_PgCyr_maxCI, Tropics_97.5pc = NBP_tropics_PgCyr_maxCI, South_97.5pc = NBP_south_PgCyr_maxCI)
                                  
###
## Begin writing out NBP to file

write.table(output,file = paste(out_dir,"/",output_prefix,"NBP",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine GPP together the data into an output variable

output = data.frame(Year = years, Global = GPP_global_PgCyr, North = GPP_north_PgCyr, Tropics = GPP_tropics_PgCyr, South = GPP_south_PgCyr,
                                  Global_2.5pc = GPP_global_PgCyr_minCI, North_2.5pc = GPP_north_PgCyr_minCI, Tropics_2.5pc = GPP_tropics_PgCyr_minCI, South_2.5pc = GPP_south_PgCyr_minCI,
                                  Global_25pc = GPP_global_PgCyr_lowCI, North_25pc = GPP_north_PgCyr_lowCI, Tropics_25pc = GPP_tropics_PgCyr_lowCI, South_25pc = GPP_south_PgCyr_lowCI,
                                  Global_75pc = GPP_global_PgCyr_highCI, North_75pc = GPP_north_PgCyr_highCI, Tropics_75pc = GPP_tropics_PgCyr_highCI, South_75pc = GPP_south_PgCyr_highCI,
                                  Global_97.5pc = GPP_global_PgCyr_maxCI, North_97.5pc = GPP_north_PgCyr_maxCI, Tropics_97.5pc = GPP_tropics_PgCyr_maxCI, South_97.5pc = GPP_south_PgCyr_maxCI)

###
## Begin writing out GPP to file

write.table(output,file = paste(out_dir,"/",output_prefix,"GPP",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine RECO together the data into an output variable

output = data.frame(Year = years, Global = RECO_global_PgCyr, North = RECO_north_PgCyr, Tropics = RECO_tropics_PgCyr, South = RECO_south_PgCyr,
                                  Global_2.5pc = RECO_global_PgCyr_minCI, North_2.5pc = RECO_north_PgCyr_minCI, Tropics_2.5pc = RECO_tropics_PgCyr_minCI, South_2.5pc = RECO_south_PgCyr_minCI,
                                  Global_25pc = RECO_global_PgCyr_lowCI, North_25pc = RECO_north_PgCyr_lowCI, Tropics_25pc = RECO_tropics_PgCyr_lowCI, South_25pc = RECO_south_PgCyr_lowCI,
                                  Global_75pc = RECO_global_PgCyr_highCI, North_75pc = RECO_north_PgCyr_highCI, Tropics_75pc = RECO_tropics_PgCyr_highCI, South_75pc = RECO_south_PgCyr_highCI,
                                  Global_97.5pc = RECO_global_PgCyr_maxCI, North_97.5pc = RECO_north_PgCyr_maxCI, Tropics_97.5pc = RECO_tropics_PgCyr_maxCI, South_97.5pc = RECO_south_PgCyr_maxCI)

###
## Begin writing out RECO to file

write.table(output,file = paste(out_dir,"/",output_prefix,"RECO",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine RHET together the data into an output variable

output = data.frame(Year = years, Global = RHET_global_PgCyr, North = RHET_north_PgCyr, Tropics = RHET_tropics_PgCyr, South = RHET_south_PgCyr,
                                  Global_2.5pc = RHET_global_PgCyr_minCI, North_2.5pc = RHET_north_PgCyr_minCI, Tropics_2.5pc = RHET_tropics_PgCyr_minCI, South_2.5pc = RHET_south_PgCyr_minCI,
                                  Global_25pc = RHET_global_PgCyr_lowCI, North_25pc = RHET_north_PgCyr_lowCI, Tropics_25pc = RHET_tropics_PgCyr_lowCI, South_25pc = RHET_south_PgCyr_lowCI,
                                  Global_75pc = RHET_global_PgCyr_highCI, North_75pc = RHET_north_PgCyr_highCI, Tropics_75pc = RHET_tropics_PgCyr_highCI, South_75pc = RHET_south_PgCyr_highCI,
                                  Global_97.5pc = RHET_global_PgCyr_maxCI, North_97.5pc = RHET_north_PgCyr_maxCI, Tropics_97.5pc = RHET_tropics_PgCyr_maxCI, South_97.5pc = RHET_south_PgCyr_maxCI)

###
## Begin writing out RHET to file

write.table(output,file = paste(out_dir,"/",output_prefix,"RHET",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine RAUTO together the data into an output variable

output = data.frame(Year = years, Global = RAUTO_global_PgCyr, North = RAUTO_north_PgCyr, Tropics = RAUTO_tropics_PgCyr, South = RAUTO_south_PgCyr,
                                  Global_2.5pc = RAUTO_global_PgCyr_minCI, North_2.5pc = RAUTO_north_PgCyr_minCI, Tropics_2.5pc = RAUTO_tropics_PgCyr_minCI, South_2.5pc = RAUTO_south_PgCyr_minCI,
                                  Global_25pc = RAUTO_global_PgCyr_lowCI, North_25pc = RAUTO_north_PgCyr_lowCI, Tropics_25pc = RAUTO_tropics_PgCyr_lowCI, South_25pc = RAUTO_south_PgCyr_lowCI,
                                  Global_75pc = RAUTO_global_PgCyr_highCI, North_75pc = RAUTO_north_PgCyr_highCI, Tropics_75pc = RAUTO_tropics_PgCyr_highCI, South_75pc = RAUTO_south_PgCyr_highCI,
                                  Global_97.5pc = RAUTO_global_PgCyr_maxCI, North_97.5pc = RAUTO_north_PgCyr_maxCI, Tropics_97.5pc = RAUTO_tropics_PgCyr_maxCI, South_97.5pc = RAUTO_south_PgCyr_maxCI)

###
## Begin writing out RAUTO to file

write.table(output,file = paste(out_dir,"/",output_prefix,"RAUTO",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine FIRE together the data into an output variable

output = data.frame(Year = years, Global = FIRE_global_PgCyr, North = FIRE_north_PgCyr, Tropics = FIRE_tropics_PgCyr, South = FIRE_south_PgCyr,
                                  Global_2.5pc = FIRE_global_PgCyr_minCI, North_2.5pc = FIRE_north_PgCyr_minCI, Tropics_2.5pc = FIRE_tropics_PgCyr_minCI, South_2.5pc = FIRE_south_PgCyr_minCI,
                                  Global_25pc = FIRE_global_PgCyr_lowCI, North_25pc = FIRE_north_PgCyr_lowCI, Tropics_25pc = FIRE_tropics_PgCyr_lowCI, South_25pc = FIRE_south_PgCyr_lowCI,
                                  Global_75pc = FIRE_global_PgCyr_highCI, North_75pc = FIRE_north_PgCyr_highCI, Tropics_75pc = FIRE_tropics_PgCyr_highCI, South_75pc = FIRE_south_PgCyr_highCI,
                                  Global_97.5pc = FIRE_global_PgCyr_maxCI, North_97.5pc = FIRE_north_PgCyr_maxCI, Tropics_97.5pc = FIRE_tropics_PgCyr_maxCI, South_97.5pc = FIRE_south_PgCyr_maxCI)

###
## Begin writing out FIRE to file

write.table(output,file = paste(out_dir,"/",output_prefix,"FIRE",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)


