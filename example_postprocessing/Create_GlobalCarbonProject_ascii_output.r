
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

print("Begin creation of GCP compatible ascii file for annual NBP / NBE / NEE / NPP / GPP / RECO / RHET / RHET_LIT / RHET_SOM / RAUTO / FIRE / HARV / LAI / ET / BIO / DOM / LAB / FOL / ROOT / WOOD / LIT / SOM")

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/reccap2_permafrost_1deg_dalec2_isimip3a_agb_lca_nbe_CsomPriorNCSDC3m/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/reccap2_permafrost_1deg_dalec2_isimip3a_agb_lca_nbe_gpp_CsomPriorNCSDC3m/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/Miombo_0.5deg_allWood"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_AGB_GPP_NBE"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_oneAGB"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_AGB"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_oneAGB"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_AGB_GPP"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS//DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_AGB_NBE"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_1deg_dalec4_trendyv13_LCA_AGB"

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
NBE_global_PgCyr = rep(0, nos_years)
NBE_north_PgCyr = rep(0, nos_years)
NBE_tropics_PgCyr = rep(0, nos_years)
NBE_south_PgCyr = rep(0, nos_years)
NEE_global_PgCyr = rep(0, nos_years)
NEE_north_PgCyr = rep(0, nos_years)
NEE_tropics_PgCyr = rep(0, nos_years)
NEE_south_PgCyr = rep(0, nos_years)
NPP_global_PgCyr = rep(0, nos_years)
NPP_north_PgCyr = rep(0, nos_years)
NPP_tropics_PgCyr = rep(0, nos_years)
NPP_south_PgCyr = rep(0, nos_years)
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
RHET_LIT_global_PgC = rep(0, nos_years)
RHET_LIT_north_PgC = rep(0, nos_years)
RHET_LIT_tropics_PgC = rep(0, nos_years)
RHET_LIT_south_PgC = rep(0, nos_years)
RHET_SOM_global_PgC = rep(0, nos_years)
RHET_SOM_north_PgC = rep(0, nos_years)
RHET_SOM_tropics_PgC = rep(0, nos_years)
RHET_SOM_south_PgC = rep(0, nos_years)
RAUTO_global_PgCyr = rep(0, nos_years)
RAUTO_north_PgCyr = rep(0, nos_years)
RAUTO_tropics_PgCyr = rep(0, nos_years)
RAUTO_south_PgCyr = rep(0, nos_years)
FIRE_global_PgCyr = rep(0, nos_years)
FIRE_north_PgCyr = rep(0, nos_years)
FIRE_tropics_PgCyr = rep(0, nos_years)
FIRE_south_PgCyr = rep(0, nos_years)
HARV_global_PgCyr = rep(0, nos_years)
HARV_north_PgCyr = rep(0, nos_years)
HARV_tropics_PgCyr = rep(0, nos_years)
HARV_south_PgCyr = rep(0, nos_years)
LAI_global_m2m2 = rep(0, nos_years)
LAI_north_m2m2 = rep(0, nos_years)
LAI_tropics_m2m2 = rep(0, nos_years)
LAI_south_m2m2 = rep(0, nos_years)
BIO_global_PgC = rep(0, nos_years)
BIO_north_PgC = rep(0, nos_years)
BIO_tropics_PgC = rep(0, nos_years)
BIO_south_PgC = rep(0, nos_years)
DOM_global_PgC = rep(0, nos_years)
DOM_north_PgC = rep(0, nos_years)
DOM_tropics_PgC = rep(0, nos_years)
DOM_south_PgC = rep(0, nos_years)
ET_global_PgH2Oyr = rep(0, nos_years)
ET_north_PgH2Oyr = rep(0, nos_years)
ET_tropics_PgH2Oyr = rep(0, nos_years)
ET_south_PgH2Oyr = rep(0, nos_years)
LAB_global_PgC = rep(0, nos_years)
LAB_north_PgC = rep(0, nos_years)
LAB_tropics_PgC = rep(0, nos_years)
LAB_south_PgC = rep(0, nos_years)
FOL_global_PgC = rep(0, nos_years)
FOL_north_PgC = rep(0, nos_years)
FOL_tropics_PgC = rep(0, nos_years)
FOL_south_PgC = rep(0, nos_years)
ROOT_global_PgC = rep(0, nos_years)
ROOT_north_PgC = rep(0, nos_years)
ROOT_tropics_PgC = rep(0, nos_years)
ROOT_south_PgC = rep(0, nos_years)
WOOD_global_PgC = rep(0, nos_years)
WOOD_north_PgC = rep(0, nos_years)
WOOD_tropics_PgC = rep(0, nos_years)
WOOD_south_PgC = rep(0, nos_years)
LIT_global_PgC = rep(0, nos_years)
LIT_north_PgC = rep(0, nos_years)
LIT_tropics_PgC = rep(0, nos_years)
LIT_south_PgC = rep(0, nos_years)
SOM_global_PgC = rep(0, nos_years)
SOM_north_PgC = rep(0, nos_years)
SOM_tropics_PgC = rep(0, nos_years)
SOM_south_PgC = rep(0, nos_years)
# Minimum CI
NBP_global_PgCyr_minCI = rep(0, nos_years)
NBP_north_PgCyr_minCI = rep(0, nos_years)
NBP_tropics_PgCyr_minCI = rep(0, nos_years)
NBP_south_PgCyr_minCI = rep(0, nos_years)
NBE_global_PgCyr_minCI = rep(0, nos_years)
NBE_north_PgCyr_minCI = rep(0, nos_years)
NBE_tropics_PgCyr_minCI = rep(0, nos_years)
NBE_south_PgCyr_minCI = rep(0, nos_years)
NEE_global_PgCyr_minCI = rep(0, nos_years)
NEE_north_PgCyr_minCI = rep(0, nos_years)
NEE_tropics_PgCyr_minCI = rep(0, nos_years)
NEE_south_PgCyr_minCI = rep(0, nos_years)
NPP_global_PgCyr_minCI = rep(0, nos_years)
NPP_north_PgCyr_minCI = rep(0, nos_years)
NPP_tropics_PgCyr_minCI = rep(0, nos_years)
NPP_south_PgCyr_minCI = rep(0, nos_years)
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
RHET_LIT_global_PgC_minCI = rep(0, nos_years)
RHET_LIT_north_PgC_minCI = rep(0, nos_years)
RHET_LIT_tropics_PgC_minCI = rep(0, nos_years)
RHET_LIT_south_PgC_minCI = rep(0, nos_years)
RHET_SOM_global_PgC_minCI = rep(0, nos_years)
RHET_SOM_north_PgC_minCI = rep(0, nos_years)
RHET_SOM_tropics_PgC_minCI = rep(0, nos_years)
RHET_SOM_south_PgC_minCI = rep(0, nos_years)
RAUTO_global_PgCyr_minCI = rep(0, nos_years)
RAUTO_north_PgCyr_minCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_minCI = rep(0, nos_years)
RAUTO_south_PgCyr_minCI = rep(0, nos_years)
FIRE_global_PgCyr_minCI = rep(0, nos_years)
FIRE_north_PgCyr_minCI = rep(0, nos_years)
FIRE_tropics_PgCyr_minCI = rep(0, nos_years)
FIRE_south_PgCyr_minCI = rep(0, nos_years)
HARV_global_PgCyr_minCI = rep(0, nos_years)
HARV_north_PgCyr_minCI = rep(0, nos_years)
HARV_tropics_PgCyr_minCI = rep(0, nos_years)
HARV_south_PgCyr_minCI = rep(0, nos_years)
LAI_global_m2m2_minCI = rep(0, nos_years)
LAI_north_m2m2_minCI = rep(0, nos_years)
LAI_tropics_m2m2_minCI = rep(0, nos_years)
LAI_south_m2m2_minCI = rep(0, nos_years)
BIO_global_PgC_minCI = rep(0, nos_years)
BIO_north_PgC_minCI = rep(0, nos_years)
BIO_tropics_PgC_minCI = rep(0, nos_years)
BIO_south_PgC_minCI = rep(0, nos_years)
DOM_global_PgC_minCI = rep(0, nos_years)
DOM_north_PgC_minCI = rep(0, nos_years)
DOM_tropics_PgC_minCI = rep(0, nos_years)
DOM_south_PgC_minCI = rep(0, nos_years)
ET_global_PgH2Oyr_minCI = rep(0, nos_years)
ET_north_PgH2Oyr_minCI = rep(0, nos_years)
ET_tropics_PgH2Oyr_minCI = rep(0, nos_years)
ET_south_PgH2Oyr_minCI = rep(0, nos_years)
LAB_global_PgC_minCI = rep(0, nos_years)
LAB_north_PgC_minCI = rep(0, nos_years)
LAB_tropics_PgC_minCI = rep(0, nos_years)
LAB_south_PgC_minCI = rep(0, nos_years)
FOL_global_PgC_minCI = rep(0, nos_years)
FOL_north_PgC_minCI = rep(0, nos_years)
FOL_tropics_PgC_minCI = rep(0, nos_years)
FOL_south_PgC_minCI = rep(0, nos_years)
ROOT_global_PgC_minCI = rep(0, nos_years)
ROOT_north_PgC_minCI = rep(0, nos_years)
ROOT_tropics_PgC_minCI = rep(0, nos_years)
ROOT_south_PgC_minCI = rep(0, nos_years)
WOOD_global_PgC_minCI = rep(0, nos_years)
WOOD_north_PgC_minCI = rep(0, nos_years)
WOOD_tropics_PgC_minCI = rep(0, nos_years)
WOOD_south_PgC_minCI = rep(0, nos_years)
LIT_global_PgC_minCI = rep(0, nos_years)
LIT_north_PgC_minCI = rep(0, nos_years)
LIT_tropics_PgC_minCI = rep(0, nos_years)
LIT_south_PgC_minCI = rep(0, nos_years)
SOM_global_PgC_minCI = rep(0, nos_years)
SOM_north_PgC_minCI = rep(0, nos_years)
SOM_tropics_PgC_minCI = rep(0, nos_years)
SOM_south_PgC_minCI = rep(0, nos_years)
# Lower CI
NBP_global_PgCyr_lowCI = rep(0, nos_years)
NBP_north_PgCyr_lowCI = rep(0, nos_years)
NBP_tropics_PgCyr_lowCI = rep(0, nos_years)
NBP_south_PgCyr_lowCI = rep(0, nos_years)
NBE_global_PgCyr_lowCI = rep(0, nos_years)
NBE_north_PgCyr_lowCI = rep(0, nos_years)
NBE_tropics_PgCyr_lowCI = rep(0, nos_years)
NBE_south_PgCyr_lowCI = rep(0, nos_years)
NEE_global_PgCyr_lowCI = rep(0, nos_years)
NEE_north_PgCyr_lowCI = rep(0, nos_years)
NEE_tropics_PgCyr_lowCI = rep(0, nos_years)
NEE_south_PgCyr_lowCI = rep(0, nos_years)
NPP_global_PgCyr_lowCI = rep(0, nos_years)
NPP_north_PgCyr_lowCI = rep(0, nos_years)
NPP_tropics_PgCyr_lowCI = rep(0, nos_years)
NPP_south_PgCyr_lowCI = rep(0, nos_years)
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
RHET_LIT_global_PgC_lowCI = rep(0, nos_years)
RHET_LIT_north_PgC_lowCI = rep(0, nos_years)
RHET_LIT_tropics_PgC_lowCI = rep(0, nos_years)
RHET_LIT_south_PgC_lowCI = rep(0, nos_years)
RHET_SOM_global_PgC_lowCI = rep(0, nos_years)
RHET_SOM_north_PgC_lowCI = rep(0, nos_years)
RHET_SOM_tropics_PgC_lowCI = rep(0, nos_years)
RHET_SOM_south_PgC_lowCI = rep(0, nos_years)
RAUTO_global_PgCyr_lowCI = rep(0, nos_years)
RAUTO_north_PgCyr_lowCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_lowCI = rep(0, nos_years)
RAUTO_south_PgCyr_lowCI = rep(0, nos_years)
FIRE_global_PgCyr_lowCI = rep(0, nos_years)
FIRE_north_PgCyr_lowCI = rep(0, nos_years)
FIRE_tropics_PgCyr_lowCI = rep(0, nos_years)
FIRE_south_PgCyr_lowCI = rep(0, nos_years)
HARV_global_PgCyr_lowCI = rep(0, nos_years)
HARV_north_PgCyr_lowCI = rep(0, nos_years)
HARV_tropics_PgCyr_lowCI = rep(0, nos_years)
HARV_south_PgCyr_lowCI = rep(0, nos_years)
LAI_global_m2m2_lowCI = rep(0, nos_years)
LAI_north_m2m2_lowCI = rep(0, nos_years)
LAI_tropics_m2m2_lowCI = rep(0, nos_years)
LAI_south_m2m2_lowCI = rep(0, nos_years)
BIO_global_PgC_lowCI = rep(0, nos_years)
BIO_north_PgC_lowCI = rep(0, nos_years)
BIO_tropics_PgC_lowCI = rep(0, nos_years)
BIO_south_PgC_lowCI = rep(0, nos_years)
DOM_global_PgC_lowCI = rep(0, nos_years)
DOM_north_PgC_lowCI = rep(0, nos_years)
DOM_tropics_PgC_lowCI = rep(0, nos_years)
DOM_south_PgC_lowCI = rep(0, nos_years)
ET_global_PgH2Oyr_lowCI = rep(0, nos_years)
ET_north_PgH2Oyr_lowCI = rep(0, nos_years)
ET_tropics_PgH2Oyr_lowCI = rep(0, nos_years)
ET_south_PgH2Oyr_lowCI = rep(0, nos_years)
LAB_global_PgC_lowCI = rep(0, nos_years)
LAB_north_PgC_lowCI = rep(0, nos_years)
LAB_tropics_PgC_lowCI = rep(0, nos_years)
LAB_south_PgC_lowCI = rep(0, nos_years)
FOL_global_PgC_lowCI = rep(0, nos_years)
FOL_north_PgC_lowCI = rep(0, nos_years)
FOL_tropics_PgC_lowCI = rep(0, nos_years)
FOL_south_PgC_lowCI = rep(0, nos_years)
ROOT_global_PgC_lowCI = rep(0, nos_years)
ROOT_north_PgC_lowCI = rep(0, nos_years)
ROOT_tropics_PgC_lowCI = rep(0, nos_years)
ROOT_south_PgC_lowCI = rep(0, nos_years)
WOOD_global_PgC_lowCI = rep(0, nos_years)
WOOD_north_PgC_lowCI = rep(0, nos_years)
WOOD_tropics_PgC_lowCI = rep(0, nos_years)
WOOD_south_PgC_lowCI = rep(0, nos_years)
LIT_global_PgC_lowCI = rep(0, nos_years)
LIT_north_PgC_lowCI = rep(0, nos_years)
LIT_tropics_PgC_lowCI = rep(0, nos_years)
LIT_south_PgC_lowCI = rep(0, nos_years)
SOM_global_PgC_lowCI = rep(0, nos_years)
SOM_north_PgC_lowCI = rep(0, nos_years)
SOM_tropics_PgC_lowCI = rep(0, nos_years)
SOM_south_PgC_lowCI = rep(0, nos_years)
# Upper CI
NBP_global_PgCyr_highCI = rep(0, nos_years)
NBP_north_PgCyr_highCI = rep(0, nos_years)
NBP_tropics_PgCyr_highCI = rep(0, nos_years)
NBP_south_PgCyr_highCI = rep(0, nos_years)
NBE_global_PgCyr_highCI = rep(0, nos_years)
NBE_north_PgCyr_highCI = rep(0, nos_years)
NBE_tropics_PgCyr_highCI = rep(0, nos_years)
NBE_south_PgCyr_highCI = rep(0, nos_years)
NEE_global_PgCyr_highCI = rep(0, nos_years)
NEE_north_PgCyr_highCI = rep(0, nos_years)
NEE_tropics_PgCyr_highCI = rep(0, nos_years)
NEE_south_PgCyr_highCI = rep(0, nos_years)
NPP_global_PgCyr_highCI = rep(0, nos_years)
NPP_north_PgCyr_highCI = rep(0, nos_years)
NPP_tropics_PgCyr_highCI = rep(0, nos_years)
NPP_south_PgCyr_highCI = rep(0, nos_years)
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
RHET_LIT_global_PgC_highCI = rep(0, nos_years)
RHET_LIT_north_PgC_highCI = rep(0, nos_years)
RHET_LIT_tropics_PgC_highCI = rep(0, nos_years)
RHET_LIT_south_PgC_highCI = rep(0, nos_years)
RHET_SOM_global_PgC_highCI = rep(0, nos_years)
RHET_SOM_north_PgC_highCI = rep(0, nos_years)
RHET_SOM_tropics_PgC_highCI = rep(0, nos_years)
RHET_SOM_south_PgC_highCI = rep(0, nos_years)
RAUTO_global_PgCyr_highCI = rep(0, nos_years)
RAUTO_north_PgCyr_highCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_highCI = rep(0, nos_years)
RAUTO_south_PgCyr_highCI = rep(0, nos_years)
FIRE_global_PgCyr_highCI = rep(0, nos_years)
FIRE_north_PgCyr_highCI = rep(0, nos_years)
FIRE_tropics_PgCyr_highCI = rep(0, nos_years)
FIRE_south_PgCyr_highCI = rep(0, nos_years)
HARV_global_PgCyr_highCI = rep(0, nos_years)
HARV_north_PgCyr_highCI = rep(0, nos_years)
HARV_tropics_PgCyr_highCI = rep(0, nos_years)
HARV_south_PgCyr_highCI = rep(0, nos_years)
LAI_global_m2m2_highCI = rep(0, nos_years)
LAI_north_m2m2_highCI = rep(0, nos_years)
LAI_tropics_m2m2_highCI = rep(0, nos_years)
LAI_south_m2m2_highCI = rep(0, nos_years)
BIO_global_PgC_highCI = rep(0, nos_years)
BIO_north_PgC_highCI = rep(0, nos_years)
BIO_tropics_PgC_highCI = rep(0, nos_years)
BIO_south_PgC_highCI = rep(0, nos_years)
DOM_global_PgC_highCI = rep(0, nos_years)
DOM_north_PgC_highCI = rep(0, nos_years)
DOM_tropics_PgC_highCI = rep(0, nos_years)
DOM_south_PgC_highCI = rep(0, nos_years)
ET_global_PgH2Oyr_highCI = rep(0, nos_years)
ET_north_PgH2Oyr_highCI = rep(0, nos_years)
ET_tropics_PgH2Oyr_highCI = rep(0, nos_years)
ET_south_PgH2Oyr_highCI = rep(0, nos_years)
LAB_global_PgC_highCI = rep(0, nos_years)
LAB_north_PgC_highCI = rep(0, nos_years)
LAB_tropics_PgC_highCI = rep(0, nos_years)
LAB_south_PgC_highCI = rep(0, nos_years)
FOL_global_PgC_highCI = rep(0, nos_years)
FOL_north_PgC_highCI = rep(0, nos_years)
FOL_tropics_PgC_highCI = rep(0, nos_years)
FOL_south_PgC_highCI = rep(0, nos_years)
ROOT_global_PgC_highCI = rep(0, nos_years)
ROOT_north_PgC_highCI = rep(0, nos_years)
ROOT_tropics_PgC_highCI = rep(0, nos_years)
ROOT_south_PgC_highCI = rep(0, nos_years)
WOOD_global_PgC_highCI = rep(0, nos_years)
WOOD_north_PgC_highCI = rep(0, nos_years)
WOOD_tropics_PgC_highCI = rep(0, nos_years)
WOOD_south_PgC_highCI = rep(0, nos_years)
LIT_global_PgC_highCI = rep(0, nos_years)
LIT_north_PgC_highCI = rep(0, nos_years)
LIT_tropics_PgC_highCI = rep(0, nos_years)
LIT_south_PgC_highCI = rep(0, nos_years)
SOM_global_PgC_highCI = rep(0, nos_years)
SOM_north_PgC_highCI = rep(0, nos_years)
SOM_tropics_PgC_highCI = rep(0, nos_years)
SOM_south_PgC_highCI = rep(0, nos_years)
# Maximum CI
NBP_global_PgCyr_maxCI = rep(0, nos_years)
NBP_north_PgCyr_maxCI = rep(0, nos_years)
NBP_tropics_PgCyr_maxCI = rep(0, nos_years)
NBP_south_PgCyr_maxCI = rep(0, nos_years)
NBE_global_PgCyr_maxCI = rep(0, nos_years)
NBE_north_PgCyr_maxCI = rep(0, nos_years)
NBE_tropics_PgCyr_maxCI = rep(0, nos_years)
NBE_south_PgCyr_maxCI = rep(0, nos_years)
NEE_global_PgCyr_maxCI = rep(0, nos_years)
NEE_north_PgCyr_maxCI = rep(0, nos_years)
NEE_tropics_PgCyr_maxCI = rep(0, nos_years)
NEE_south_PgCyr_maxCI = rep(0, nos_years)
NPP_global_PgCyr_maxCI = rep(0, nos_years)
NPP_north_PgCyr_maxCI = rep(0, nos_years)
NPP_tropics_PgCyr_maxCI = rep(0, nos_years)
NPP_south_PgCyr_maxCI = rep(0, nos_years)
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
RHET_LIT_global_PgC_maxCI = rep(0, nos_years)
RHET_LIT_north_PgC_maxCI = rep(0, nos_years)
RHET_LIT_tropics_PgC_maxCI = rep(0, nos_years)
RHET_LIT_south_PgC_maxCI = rep(0, nos_years)
RHET_SOM_global_PgC_maxCI = rep(0, nos_years)
RHET_SOM_north_PgC_maxCI = rep(0, nos_years)
RHET_SOM_tropics_PgC_maxCI = rep(0, nos_years)
RHET_SOM_south_PgC_maxCI = rep(0, nos_years)
RAUTO_global_PgCyr_maxCI = rep(0, nos_years)
RAUTO_north_PgCyr_maxCI = rep(0, nos_years)
RAUTO_tropics_PgCyr_maxCI = rep(0, nos_years)
RAUTO_south_PgCyr_maxCI = rep(0, nos_years)
FIRE_global_PgCyr_maxCI = rep(0, nos_years)
FIRE_north_PgCyr_maxCI = rep(0, nos_years)
FIRE_tropics_PgCyr_maxCI = rep(0, nos_years)
FIRE_south_PgCyr_maxCI = rep(0, nos_years)
HARV_global_PgCyr_maxCI = rep(0, nos_years)
HARV_north_PgCyr_maxCI = rep(0, nos_years)
HARV_tropics_PgCyr_maxCI = rep(0, nos_years)
HARV_south_PgCyr_maxCI = rep(0, nos_years)
LAI_global_m2m2_maxCI = rep(0, nos_years)
LAI_north_m2m2_maxCI = rep(0, nos_years)
LAI_tropics_m2m2_maxCI = rep(0, nos_years)
LAI_south_m2m2_maxCI = rep(0, nos_years)
BIO_global_PgC_maxCI = rep(0, nos_years)
BIO_north_PgC_maxCI = rep(0, nos_years)
BIO_tropics_PgC_maxCI = rep(0, nos_years)
BIO_south_PgC_maxCI = rep(0, nos_years)
DOM_global_PgC_maxCI = rep(0, nos_years)
DOM_north_PgC_maxCI = rep(0, nos_years)
DOM_tropics_PgC_maxCI = rep(0, nos_years)
DOM_south_PgC_maxCI = rep(0, nos_years)
ET_global_PgH2Oyr_maxCI = rep(0, nos_years)
ET_north_PgH2Oyr_maxCI = rep(0, nos_years)
ET_tropics_PgH2Oyr_maxCI = rep(0, nos_years)
ET_south_PgH2Oyr_maxCI = rep(0, nos_years)
LAB_global_PgC_maxCI = rep(0, nos_years)
LAB_north_PgC_maxCI = rep(0, nos_years)
LAB_tropics_PgC_maxCI = rep(0, nos_years)
LAB_south_PgC_maxCI = rep(0, nos_years)
FOL_global_PgC_maxCI = rep(0, nos_years)
FOL_north_PgC_maxCI = rep(0, nos_years)
FOL_tropics_PgC_maxCI = rep(0, nos_years)
FOL_south_PgC_maxCI = rep(0, nos_years)
ROOT_global_PgC_maxCI = rep(0, nos_years)
ROOT_north_PgC_maxCI = rep(0, nos_years)
ROOT_tropics_PgC_maxCI = rep(0, nos_years)
ROOT_south_PgC_maxCI = rep(0, nos_years)
WOOD_global_PgC_maxCI = rep(0, nos_years)
WOOD_north_PgC_maxCI = rep(0, nos_years)
WOOD_tropics_PgC_maxCI = rep(0, nos_years)
WOOD_south_PgC_maxCI = rep(0, nos_years)
LIT_global_PgC_maxCI = rep(0, nos_years)
LIT_north_PgC_maxCI = rep(0, nos_years)
LIT_tropics_PgC_maxCI = rep(0, nos_years)
LIT_south_PgC_maxCI = rep(0, nos_years)
SOM_global_PgC_maxCI = rep(0, nos_years)
SOM_north_PgC_maxCI = rep(0, nos_years)
SOM_tropics_PgC_maxCI = rep(0, nos_years)
SOM_south_PgC_maxCI = rep(0, nos_years)

# Counters
nos_global = 0
nos_tropics = 0
nos_north = 0
nos_south = 0

# Fill the output arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]
         
         # Determine nos days per time step         
         deltat = 365.25

         # Update counter
         nos_global = nos_global + 1

         # Accumulate global
         NBP_global_PgCyr = NBP_global_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_minCI = NBP_global_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_lowCI = NBP_global_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_highCI = NBP_global_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_maxCI = NBP_global_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBE_global_PgCyr = NBE_global_PgCyr + (grid_output$mean_annual_nbe_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBE_global_PgCyr_minCI = NBE_global_PgCyr_minCI + (grid_output$mean_annual_nbe_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBE_global_PgCyr_lowCI = NBE_global_PgCyr_lowCI + (grid_output$mean_annual_nbe_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBE_global_PgCyr_highCI = NBE_global_PgCyr_highCI + (grid_output$mean_annual_nbe_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBE_global_PgCyr_maxCI = NBE_global_PgCyr_maxCI + (grid_output$mean_annual_nbe_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NEE_global_PgCyr = NEE_global_PgCyr + (grid_output$mean_annual_nee_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NEE_global_PgCyr_minCI = NEE_global_PgCyr_minCI + (grid_output$mean_annual_nee_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NEE_global_PgCyr_lowCI = NEE_global_PgCyr_lowCI + (grid_output$mean_annual_nee_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NEE_global_PgCyr_highCI = NEE_global_PgCyr_highCI + (grid_output$mean_annual_nee_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NEE_global_PgCyr_maxCI = NEE_global_PgCyr_maxCI + (grid_output$mean_annual_nee_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NPP_global_PgCyr = NPP_global_PgCyr + (grid_output$mean_annual_npp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NPP_global_PgCyr_minCI = NPP_global_PgCyr_minCI + (grid_output$mean_annual_npp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NPP_global_PgCyr_lowCI = NPP_global_PgCyr_lowCI + (grid_output$mean_annual_npp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NPP_global_PgCyr_highCI = NPP_global_PgCyr_highCI + (grid_output$mean_annual_npp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NPP_global_PgCyr_maxCI = NPP_global_PgCyr_maxCI + (grid_output$mean_annual_npp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
         RHET_LIT_global_PgC = RHET_LIT_global_PgC + (grid_output$mean_annual_rhet_litter_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_LIT_global_PgC_minCI = RHET_LIT_global_PgC_minCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_LIT_global_PgC_lowCI = RHET_LIT_global_PgC_lowCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_LIT_global_PgC_highCI = RHET_LIT_global_PgC_highCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_LIT_global_PgC_maxCI = RHET_LIT_global_PgC_maxCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_SOM_global_PgC = RHET_SOM_global_PgC + (grid_output$mean_annual_rhet_som_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_SOM_global_PgC_minCI = RHET_SOM_global_PgC_minCI + (grid_output$mean_annual_rhet_som_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_SOM_global_PgC_lowCI = RHET_SOM_global_PgC_lowCI + (grid_output$mean_annual_rhet_som_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_SOM_global_PgC_highCI = RHET_SOM_global_PgC_highCI + (grid_output$mean_annual_rhet_som_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         RHET_SOM_global_PgC_maxCI = RHET_SOM_global_PgC_maxCI + (grid_output$mean_annual_rhet_som_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
         HARV_global_PgCyr = HARV_global_PgCyr + (grid_output$mean_annual_harvest_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         HARV_global_PgCyr_minCI = HARV_global_PgCyr_minCI   + (grid_output$mean_annual_harvest_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         HARV_global_PgCyr_lowCI = HARV_global_PgCyr_lowCI   + (grid_output$mean_annual_harvest_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         HARV_global_PgCyr_highCI = HARV_global_PgCyr_highCI + (grid_output$mean_annual_harvest_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         HARV_global_PgCyr_maxCI = HARV_global_PgCyr_maxCI   + (grid_output$mean_annual_harvest_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LAI_global_m2m2 = LAI_global_m2m2 + (grid_output$mean_annual_lai_m2m2[n,mid_quant,])
         LAI_global_m2m2_minCI = LAI_global_m2m2_minCI + (grid_output$mean_annual_lai_m2m2[n,min_quant,])
         LAI_global_m2m2_lowCI = LAI_global_m2m2_lowCI + (grid_output$mean_annual_lai_m2m2[n,low_quant,])
         LAI_global_m2m2_highCI = LAI_global_m2m2_highCI + (grid_output$mean_annual_lai_m2m2[n,high_quant,])
         LAI_global_m2m2_maxCI = LAI_global_m2m2_maxCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])
         BIO_global_PgC = BIO_global_PgC + (grid_output$mean_annual_biomass_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         BIO_global_PgC_minCI = BIO_global_PgC_minCI + (grid_output$mean_annual_biomass_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         BIO_global_PgC_lowCI = BIO_global_PgC_lowCI + (grid_output$mean_annual_biomass_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         BIO_global_PgC_highCI = BIO_global_PgC_highCI + (grid_output$mean_annual_biomass_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         BIO_global_PgC_maxCI = BIO_global_PgC_maxCI + (grid_output$mean_annual_biomass_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         DOM_global_PgC = DOM_global_PgC + (grid_output$mean_annual_dom_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         DOM_global_PgC_minCI = DOM_global_PgC_minCI + (grid_output$mean_annual_dom_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         DOM_global_PgC_lowCI = DOM_global_PgC_lowCI + (grid_output$mean_annual_dom_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         DOM_global_PgC_highCI = DOM_global_PgC_highCI + (grid_output$mean_annual_dom_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         DOM_global_PgC_maxCI = DOM_global_PgC_maxCI + (grid_output$mean_annual_dom_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         ET_global_PgH2Oyr = ET_global_PgH2Oyr + (grid_output$mean_annual_ET_kgH2Om2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
         ET_global_PgH2Oyr_minCI = ET_global_PgH2Oyr_minCI + (grid_output$mean_annual_ET_kgH2Om2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
         ET_global_PgH2Oyr_lowCI = ET_global_PgH2Oyr_lowCI + (grid_output$mean_annual_ET_kgH2Om2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
         ET_global_PgH2Oyr_highCI = ET_global_PgH2Oyr_highCI + (grid_output$mean_annual_ET_kgH2Om2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
         ET_global_PgH2Oyr_maxCI = ET_global_PgH2Oyr_maxCI + (grid_output$mean_annual_ET_kgH2Om2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
         LAB_global_PgC = LAB_global_PgC + (grid_output$mean_annual_labile_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LAB_global_PgC_minCI = LAB_global_PgC_minCI + (grid_output$mean_annual_labile_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LAB_global_PgC_lowCI = LAB_global_PgC_lowCI + (grid_output$mean_annual_labile_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LAB_global_PgC_highCI = LAB_global_PgC_highCI + (grid_output$mean_annual_labile_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LAB_global_PgC_maxCI = LAB_global_PgC_maxCI + (grid_output$mean_annual_labile_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FOL_global_PgC = FOL_global_PgC + (grid_output$mean_annual_foliage_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FOL_global_PgC_minCI = FOL_global_PgC_minCI + (grid_output$mean_annual_foliage_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FOL_global_PgC_lowCI = FOL_global_PgC_lowCI + (grid_output$mean_annual_foliage_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FOL_global_PgC_highCI = FOL_global_PgC_highCI + (grid_output$mean_annual_foliage_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         FOL_global_PgC_maxCI = FOL_global_PgC_maxCI + (grid_output$mean_annual_foliage_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         ROOT_global_PgC = ROOT_global_PgC + (grid_output$mean_annual_roots_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         ROOT_global_PgC_minCI = ROOT_global_PgC_minCI + (grid_output$mean_annual_roots_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         ROOT_global_PgC_lowCI = ROOT_global_PgC_lowCI + (grid_output$mean_annual_roots_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         ROOT_global_PgC_highCI = ROOT_global_PgC_highCI + (grid_output$mean_annual_roots_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         ROOT_global_PgC_maxCI = ROOT_global_PgC_maxCI + (grid_output$mean_annual_roots_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         WOOD_global_PgC = WOOD_global_PgC + (grid_output$mean_annual_wood_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         WOOD_global_PgC_minCI = WOOD_global_PgC_minCI + (grid_output$mean_annual_wood_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         WOOD_global_PgC_lowCI = WOOD_global_PgC_lowCI + (grid_output$mean_annual_wood_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         WOOD_global_PgC_highCI = WOOD_global_PgC_highCI + (grid_output$mean_annual_wood_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         WOOD_global_PgC_maxCI = WOOD_global_PgC_maxCI + (grid_output$mean_annual_wood_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LIT_global_PgC = LIT_global_PgC + (grid_output$mean_annual_litter_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LIT_global_PgC_minCI = LIT_global_PgC_minCI + (grid_output$mean_annual_litter_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LIT_global_PgC_lowCI = LIT_global_PgC_lowCI + (grid_output$mean_annual_litter_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LIT_global_PgC_highCI = LIT_global_PgC_highCI + (grid_output$mean_annual_litter_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         LIT_global_PgC_maxCI = LIT_global_PgC_maxCI + (grid_output$mean_annual_litter_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         SOM_global_PgC = SOM_global_PgC + (grid_output$mean_annual_som_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         SOM_global_PgC_minCI = SOM_global_PgC_minCI + (grid_output$mean_annual_som_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         SOM_global_PgC_lowCI = SOM_global_PgC_lowCI + (grid_output$mean_annual_som_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         SOM_global_PgC_highCI = SOM_global_PgC_highCI + (grid_output$mean_annual_som_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         SOM_global_PgC_maxCI = SOM_global_PgC_maxCI + (grid_output$mean_annual_som_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
                           
         # Where appropriate accumulate tropical
         if (grid_output$lat[i,j] <= 30 & grid_output$lat[i,j] >= -30) {

             # Update counter
             nos_tropics = nos_tropics + 1

             NBP_tropics_PgCyr = NBP_tropics_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_minCI = NBP_tropics_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_lowCI = NBP_tropics_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_highCI = NBP_tropics_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_maxCI = NBP_tropics_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_tropics_PgCyr = NBE_tropics_PgCyr + (grid_output$mean_annual_nbe_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_tropics_PgCyr_minCI = NBE_tropics_PgCyr_minCI + (grid_output$mean_annual_nbe_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_tropics_PgCyr_lowCI = NBE_tropics_PgCyr_lowCI + (grid_output$mean_annual_nbe_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_tropics_PgCyr_highCI = NBE_tropics_PgCyr_highCI + (grid_output$mean_annual_nbe_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_tropics_PgCyr_maxCI = NBE_tropics_PgCyr_maxCI + (grid_output$mean_annual_nbe_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_tropics_PgCyr = NEE_tropics_PgCyr + (grid_output$mean_annual_nee_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_tropics_PgCyr_minCI = NEE_tropics_PgCyr_minCI + (grid_output$mean_annual_nee_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_tropics_PgCyr_lowCI = NEE_tropics_PgCyr_lowCI + (grid_output$mean_annual_nee_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_tropics_PgCyr_highCI = NEE_tropics_PgCyr_highCI + (grid_output$mean_annual_nee_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_tropics_PgCyr_maxCI = NEE_tropics_PgCyr_maxCI + (grid_output$mean_annual_nee_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_tropics_PgCyr = NPP_tropics_PgCyr + (grid_output$mean_annual_npp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_tropics_PgCyr_minCI = NPP_tropics_PgCyr_minCI + (grid_output$mean_annual_npp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_tropics_PgCyr_lowCI = NPP_tropics_PgCyr_lowCI + (grid_output$mean_annual_npp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_tropics_PgCyr_highCI = NPP_tropics_PgCyr_highCI + (grid_output$mean_annual_npp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_tropics_PgCyr_maxCI = NPP_tropics_PgCyr_maxCI + (grid_output$mean_annual_npp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
             RHET_LIT_tropics_PgC = RHET_LIT_tropics_PgC + (grid_output$mean_annual_rhet_litter_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_tropics_PgC_minCI = RHET_LIT_tropics_PgC_minCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_tropics_PgC_lowCI = RHET_LIT_tropics_PgC_lowCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_tropics_PgC_highCI = RHET_LIT_tropics_PgC_highCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_tropics_PgC_maxCI = RHET_LIT_tropics_PgC_maxCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_tropics_PgC = RHET_SOM_tropics_PgC + (grid_output$mean_annual_rhet_som_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_tropics_PgC_minCI = RHET_SOM_tropics_PgC_minCI + (grid_output$mean_annual_rhet_som_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_tropics_PgC_lowCI = RHET_SOM_tropics_PgC_lowCI + (grid_output$mean_annual_rhet_som_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_tropics_PgC_highCI = RHET_SOM_tropics_PgC_highCI + (grid_output$mean_annual_rhet_som_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_tropics_PgC_maxCI = RHET_SOM_tropics_PgC_maxCI + (grid_output$mean_annual_rhet_som_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
             HARV_tropics_PgCyr = HARV_tropics_PgCyr + (grid_output$mean_annual_harvest_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_tropics_PgCyr_minCI = HARV_tropics_PgCyr_minCI + (grid_output$mean_annual_harvest_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_tropics_PgCyr_lowCI = HARV_tropics_PgCyr_lowCI + (grid_output$mean_annual_harvest_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_tropics_PgCyr_highCI = HARV_tropics_PgCyr_highCI + (grid_output$mean_annual_harvest_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_tropics_PgCyr_maxCI = HARV_tropics_PgCyr_maxCI + (grid_output$mean_annual_harvest_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAI_tropics_m2m2 = LAI_tropics_m2m2 + (grid_output$mean_annual_lai_m2m2[n,mid_quant,])
             LAI_tropics_m2m2_minCI = LAI_tropics_m2m2_minCI + (grid_output$mean_annual_lai_m2m2[n,min_quant,])
             LAI_tropics_m2m2_lowCI = LAI_tropics_m2m2_lowCI + (grid_output$mean_annual_lai_m2m2[n,low_quant,])
             LAI_tropics_m2m2_highCI = LAI_tropics_m2m2_highCI + (grid_output$mean_annual_lai_m2m2[n,high_quant,])
             LAI_tropics_m2m2_maxCI = LAI_tropics_m2m2_maxCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])
             BIO_tropics_PgC = BIO_tropics_PgC + (grid_output$mean_annual_biomass_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_tropics_PgC_minCI = BIO_tropics_PgC_minCI + (grid_output$mean_annual_biomass_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_tropics_PgC_lowCI = BIO_tropics_PgC_lowCI + (grid_output$mean_annual_biomass_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_tropics_PgC_highCI = BIO_tropics_PgC_highCI + (grid_output$mean_annual_biomass_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_tropics_PgC_maxCI = BIO_tropics_PgC_maxCI + (grid_output$mean_annual_biomass_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_tropics_PgC = DOM_tropics_PgC + (grid_output$mean_annual_dom_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_tropics_PgC_minCI = DOM_tropics_PgC_minCI + (grid_output$mean_annual_dom_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_tropics_PgC_lowCI = DOM_tropics_PgC_lowCI + (grid_output$mean_annual_dom_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_tropics_PgC_highCI = DOM_tropics_PgC_highCI + (grid_output$mean_annual_dom_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_tropics_PgC_maxCI = DOM_tropics_PgC_maxCI + (grid_output$mean_annual_dom_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ET_tropics_PgH2Oyr = ET_tropics_PgH2Oyr + (grid_output$mean_annual_ET_kgH2Om2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_tropics_PgH2Oyr_minCI = ET_tropics_PgH2Oyr_minCI + (grid_output$mean_annual_ET_kgH2Om2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_tropics_PgH2Oyr_lowCI = ET_tropics_PgH2Oyr_lowCI + (grid_output$mean_annual_ET_kgH2Om2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_tropics_PgH2Oyr_highCI = ET_tropics_PgH2Oyr_highCI + (grid_output$mean_annual_ET_kgH2Om2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_tropics_PgH2Oyr_maxCI = ET_tropics_PgH2Oyr_maxCI + (grid_output$mean_annual_ET_kgH2Om2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             LAB_tropics_PgC = LAB_tropics_PgC + (grid_output$mean_annual_labile_gCm2[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_minCI = LAB_tropics_PgC_minCI + (grid_output$mean_annual_labile_gCm2[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_lowCI = LAB_tropics_PgC_lowCI + (grid_output$mean_annual_labile_gCm2[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_highCI = LAB_tropics_PgC_highCI + (grid_output$mean_annual_labile_gCm2[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_maxCI = LAB_tropics_PgC_maxCI + (grid_output$mean_annual_labile_gCm2[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC = FOL_tropics_PgC + (grid_output$mean_annual_foliage_gCm2[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_minCI = FOL_tropics_PgC_minCI + (grid_output$mean_annual_foliage_gCm2[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_lowCI = FOL_tropics_PgC_lowCI + (grid_output$mean_annual_foliage_gCm2[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_highCI = FOL_tropics_PgC_highCI + (grid_output$mean_annual_foliage_gCm2[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_maxCI = FOL_tropics_PgC_maxCI + (grid_output$mean_annual_foliage_gCm2[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC = ROOT_tropics_PgC + (grid_output$mean_annual_roots_gCm2[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_minCI = ROOT_tropics_PgC_minCI + (grid_output$mean_annual_roots_gCm2[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_lowCI = ROOT_tropics_PgC_lowCI + (grid_output$mean_annual_roots_gCm2[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_highCI = ROOT_tropics_PgC_highCI + (grid_output$mean_annual_roots_gCm2[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_maxCI = ROOT_tropics_PgC_maxCI + (grid_output$mean_annual_roots_gCm2[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC = WOOD_tropics_PgC + (grid_output$mean_annual_wood_gCm2[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_minCI = WOOD_tropics_PgC_minCI + (grid_output$mean_annual_wood_gCm2[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_lowCI = WOOD_tropics_PgC_lowCI + (grid_output$mean_annual_wood_gCm2[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_highCI = WOOD_tropics_PgC_highCI + (grid_output$mean_annual_wood_gCm2[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_maxCI = WOOD_tropics_PgC_maxCI + (grid_output$mean_annual_wood_gCm2[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC = LIT_tropics_PgC + (grid_output$mean_annual_litter_gCm2[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_minCI = LIT_tropics_PgC_minCI + (grid_output$mean_annual_litter_gCm2[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_lowCI = LIT_tropics_PgC_lowCI + (grid_output$mean_annual_litter_gCm2[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_highCI = LIT_tropics_PgC_highCI + (grid_output$mean_annual_litter_gCm2[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_maxCI = LIT_tropics_PgC_maxCI + (grid_output$mean_annual_litter_gCm2[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC = SOM_tropics_PgC + (grid_output$mean_annual_som_gCm2[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_minCI = SOM_tropics_PgC_minCI + (grid_output$mean_annual_som_gCm2[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_lowCI = SOM_tropics_PgC_lowCI + (grid_output$mean_annual_som_gCm2[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_highCI = SOM_tropics_PgC_highCI + (grid_output$mean_annual_som_gCm2[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_maxCI = SOM_tropics_PgC_maxCI + (grid_output$mean_annual_som_gCm2[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC = LAB_tropics_PgC + (grid_output$mean_annual_labile_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_minCI = LAB_tropics_PgC_minCI + (grid_output$mean_annual_labile_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_lowCI = LAB_tropics_PgC_lowCI + (grid_output$mean_annual_labile_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_highCI = LAB_tropics_PgC_highCI + (grid_output$mean_annual_labile_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_tropics_PgC_maxCI = LAB_tropics_PgC_maxCI + (grid_output$mean_annual_labile_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC = FOL_tropics_PgC + (grid_output$mean_annual_foliage_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_minCI = FOL_tropics_PgC_minCI + (grid_output$mean_annual_foliage_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_lowCI = FOL_tropics_PgC_lowCI + (grid_output$mean_annual_foliage_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_highCI = FOL_tropics_PgC_highCI + (grid_output$mean_annual_foliage_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_tropics_PgC_maxCI = FOL_tropics_PgC_maxCI + (grid_output$mean_annual_foliage_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC = ROOT_tropics_PgC + (grid_output$mean_annual_roots_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_minCI = ROOT_tropics_PgC_minCI + (grid_output$mean_annual_roots_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_lowCI = ROOT_tropics_PgC_lowCI + (grid_output$mean_annual_roots_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_highCI = ROOT_tropics_PgC_highCI + (grid_output$mean_annual_roots_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_tropics_PgC_maxCI = ROOT_tropics_PgC_maxCI + (grid_output$mean_annual_roots_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC = WOOD_tropics_PgC + (grid_output$mean_annual_wood_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_minCI = WOOD_tropics_PgC_minCI + (grid_output$mean_annual_wood_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_lowCI = WOOD_tropics_PgC_lowCI + (grid_output$mean_annual_wood_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_highCI = WOOD_tropics_PgC_highCI + (grid_output$mean_annual_wood_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_tropics_PgC_maxCI = WOOD_tropics_PgC_maxCI + (grid_output$mean_annual_wood_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC = LIT_tropics_PgC + (grid_output$mean_annual_litter_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_minCI = LIT_tropics_PgC_minCI + (grid_output$mean_annual_litter_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_lowCI = LIT_tropics_PgC_lowCI + (grid_output$mean_annual_litter_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_highCI = LIT_tropics_PgC_highCI + (grid_output$mean_annual_litter_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_tropics_PgC_maxCI = LIT_tropics_PgC_maxCI + (grid_output$mean_annual_litter_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC = SOM_tropics_PgC + (grid_output$mean_annual_som_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_minCI = SOM_tropics_PgC_minCI + (grid_output$mean_annual_som_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_lowCI = SOM_tropics_PgC_lowCI + (grid_output$mean_annual_som_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_highCI = SOM_tropics_PgC_highCI + (grid_output$mean_annual_som_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_tropics_PgC_maxCI = SOM_tropics_PgC_maxCI + (grid_output$mean_annual_som_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         }

         # Where appropriate accumulate north
         if (grid_output$lat[i,j] > 30) {
         
             # Update counter
             nos_north = nos_north + 1
             
             NBP_north_PgCyr = NBP_north_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_minCI = NBP_north_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_lowCI = NBP_north_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_highCI = NBP_north_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_maxCI = NBP_north_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_north_PgCyr = NBE_north_PgCyr + (grid_output$mean_annual_nbe_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_north_PgCyr_minCI = NBE_north_PgCyr_minCI + (grid_output$mean_annual_nbe_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_north_PgCyr_lowCI = NBE_north_PgCyr_lowCI + (grid_output$mean_annual_nbe_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_north_PgCyr_highCI = NBE_north_PgCyr_highCI + (grid_output$mean_annual_nbe_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_north_PgCyr_maxCI = NBE_north_PgCyr_maxCI + (grid_output$mean_annual_nbe_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_north_PgCyr = NEE_north_PgCyr + (grid_output$mean_annual_nee_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_north_PgCyr_minCI = NEE_north_PgCyr_minCI + (grid_output$mean_annual_nee_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_north_PgCyr_lowCI = NEE_north_PgCyr_lowCI + (grid_output$mean_annual_nee_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_north_PgCyr_highCI = NEE_north_PgCyr_highCI + (grid_output$mean_annual_nee_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_north_PgCyr_maxCI = NEE_north_PgCyr_maxCI + (grid_output$mean_annual_nee_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_north_PgCyr = NPP_north_PgCyr + (grid_output$mean_annual_npp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_north_PgCyr_minCI = NPP_north_PgCyr_minCI + (grid_output$mean_annual_npp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_north_PgCyr_lowCI = NPP_north_PgCyr_lowCI + (grid_output$mean_annual_npp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_north_PgCyr_highCI = NPP_north_PgCyr_highCI + (grid_output$mean_annual_npp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_north_PgCyr_maxCI = NPP_north_PgCyr_maxCI + (grid_output$mean_annual_npp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
             RHET_LIT_north_PgC = RHET_LIT_north_PgC + (grid_output$mean_annual_rhet_litter_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_north_PgC_minCI = RHET_LIT_north_PgC_minCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_north_PgC_lowCI = RHET_LIT_north_PgC_lowCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_north_PgC_highCI = RHET_LIT_north_PgC_highCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_north_PgC_maxCI = RHET_LIT_north_PgC_maxCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_north_PgC = RHET_SOM_north_PgC + (grid_output$mean_annual_rhet_som_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_north_PgC_minCI = RHET_SOM_north_PgC_minCI + (grid_output$mean_annual_rhet_som_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_north_PgC_lowCI = RHET_SOM_north_PgC_lowCI + (grid_output$mean_annual_rhet_som_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_north_PgC_highCI = RHET_SOM_north_PgC_highCI + (grid_output$mean_annual_rhet_som_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_north_PgC_maxCI = RHET_SOM_north_PgC_maxCI + (grid_output$mean_annual_rhet_som_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
             HARV_north_PgCyr = HARV_north_PgCyr + (grid_output$mean_annual_harvest_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_north_PgCyr_minCI = HARV_north_PgCyr_minCI + (grid_output$mean_annual_harvest_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_north_PgCyr_lowCI = HARV_north_PgCyr_lowCI + (grid_output$mean_annual_harvest_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_north_PgCyr_highCI = HARV_north_PgCyr_highCI + (grid_output$mean_annual_harvest_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_north_PgCyr_maxCI = HARV_north_PgCyr_maxCI + (grid_output$mean_annual_harvest_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAI_north_m2m2 = LAI_north_m2m2 + (grid_output$mean_annual_lai_m2m2[n,mid_quant,])
             LAI_north_m2m2_minCI = LAI_north_m2m2_minCI + (grid_output$mean_annual_lai_m2m2[n,min_quant,])
             LAI_north_m2m2_lowCI = LAI_north_m2m2_lowCI + (grid_output$mean_annual_lai_m2m2[n,low_quant,])
             LAI_north_m2m2_highCI = LAI_north_m2m2_highCI + (grid_output$mean_annual_lai_m2m2[n,high_quant,])
             LAI_north_m2m2_maxCI = LAI_north_m2m2_maxCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])        
             BIO_north_PgC = BIO_north_PgC + (grid_output$mean_annual_biomass_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_north_PgC_minCI = BIO_north_PgC_minCI + (grid_output$mean_annual_biomass_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_north_PgC_lowCI = BIO_north_PgC_lowCI + (grid_output$mean_annual_biomass_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_north_PgC_highCI = BIO_north_PgC_highCI + (grid_output$mean_annual_biomass_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_north_PgC_maxCI = BIO_north_PgC_maxCI + (grid_output$mean_annual_biomass_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_north_PgC = DOM_north_PgC + (grid_output$mean_annual_dom_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_north_PgC_minCI = DOM_north_PgC_minCI + (grid_output$mean_annual_dom_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_north_PgC_lowCI = DOM_north_PgC_lowCI + (grid_output$mean_annual_dom_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_north_PgC_highCI = DOM_north_PgC_highCI + (grid_output$mean_annual_dom_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_north_PgC_maxCI = DOM_north_PgC_maxCI + (grid_output$mean_annual_dom_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)              
             ET_north_PgH2Oyr = ET_north_PgH2Oyr + (grid_output$mean_annual_ET_kgH2Om2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_north_PgH2Oyr_minCI = ET_north_PgH2Oyr_minCI + (grid_output$mean_annual_ET_kgH2Om2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_north_PgH2Oyr_lowCI = ET_north_PgH2Oyr_lowCI + (grid_output$mean_annual_ET_kgH2Om2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_north_PgH2Oyr_highCI = ET_north_PgH2Oyr_highCI + (grid_output$mean_annual_ET_kgH2Om2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_north_PgH2Oyr_maxCI = ET_north_PgH2Oyr_maxCI + (grid_output$mean_annual_ET_kgH2Om2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             LAB_north_PgC = LAB_north_PgC + (grid_output$mean_annual_labile_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_north_PgC_minCI = LAB_north_PgC_minCI + (grid_output$mean_annual_labile_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_north_PgC_lowCI = LAB_north_PgC_lowCI + (grid_output$mean_annual_labile_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_north_PgC_highCI = LAB_north_PgC_highCI + (grid_output$mean_annual_labile_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_north_PgC_maxCI = LAB_north_PgC_maxCI + (grid_output$mean_annual_labile_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_north_PgC = FOL_north_PgC + (grid_output$mean_annual_foliage_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_north_PgC_minCI = FOL_north_PgC_minCI + (grid_output$mean_annual_foliage_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_north_PgC_lowCI = FOL_north_PgC_lowCI + (grid_output$mean_annual_foliage_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_north_PgC_highCI = FOL_north_PgC_highCI + (grid_output$mean_annual_foliage_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_north_PgC_maxCI = FOL_north_PgC_maxCI + (grid_output$mean_annual_foliage_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_north_PgC = ROOT_north_PgC + (grid_output$mean_annual_roots_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_north_PgC_minCI = ROOT_north_PgC_minCI + (grid_output$mean_annual_roots_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_north_PgC_lowCI = ROOT_north_PgC_lowCI + (grid_output$mean_annual_roots_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_north_PgC_highCI = ROOT_north_PgC_highCI + (grid_output$mean_annual_roots_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_north_PgC_maxCI = ROOT_north_PgC_maxCI + (grid_output$mean_annual_roots_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_north_PgC = WOOD_north_PgC + (grid_output$mean_annual_wood_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_north_PgC_minCI = WOOD_north_PgC_minCI + (grid_output$mean_annual_wood_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_north_PgC_lowCI = WOOD_north_PgC_lowCI + (grid_output$mean_annual_wood_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_north_PgC_highCI = WOOD_north_PgC_highCI + (grid_output$mean_annual_wood_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_north_PgC_maxCI = WOOD_north_PgC_maxCI + (grid_output$mean_annual_wood_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_north_PgC = LIT_north_PgC + (grid_output$mean_annual_litter_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_north_PgC_minCI = LIT_north_PgC_minCI + (grid_output$mean_annual_litter_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_north_PgC_lowCI = LIT_north_PgC_lowCI + (grid_output$mean_annual_litter_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_north_PgC_highCI = LIT_north_PgC_highCI + (grid_output$mean_annual_litter_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_north_PgC_maxCI = LIT_north_PgC_maxCI + (grid_output$mean_annual_litter_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_north_PgC = SOM_north_PgC + (grid_output$mean_annual_som_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_north_PgC_minCI = SOM_north_PgC_minCI + (grid_output$mean_annual_som_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_north_PgC_lowCI = SOM_north_PgC_lowCI + (grid_output$mean_annual_som_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_north_PgC_highCI = SOM_north_PgC_highCI + (grid_output$mean_annual_som_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_north_PgC_maxCI = SOM_north_PgC_maxCI + (grid_output$mean_annual_som_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         }

         # Where appropriate accumulate south
         if (grid_output$lat[i,j] < -30) {

             # Update counter
             nos_south = nos_south + 1         
         
             NBP_south_PgCyr = NBP_south_PgCyr + (grid_output$mean_annual_nbp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_minCI = NBP_south_PgCyr_minCI + (grid_output$mean_annual_nbp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_lowCI = NBP_south_PgCyr_lowCI + (grid_output$mean_annual_nbp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_highCI = NBP_south_PgCyr_highCI + (grid_output$mean_annual_nbp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_maxCI = NBP_south_PgCyr_maxCI + (grid_output$mean_annual_nbp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_south_PgCyr = NBE_south_PgCyr + (grid_output$mean_annual_nbe_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_south_PgCyr_minCI = NBE_south_PgCyr_minCI + (grid_output$mean_annual_nbe_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_south_PgCyr_lowCI = NBE_south_PgCyr_lowCI + (grid_output$mean_annual_nbe_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_south_PgCyr_highCI = NBE_south_PgCyr_highCI + (grid_output$mean_annual_nbe_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBE_south_PgCyr_maxCI = NBE_south_PgCyr_maxCI + (grid_output$mean_annual_nbe_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_south_PgCyr = NEE_south_PgCyr + (grid_output$mean_annual_nee_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_south_PgCyr_minCI = NEE_south_PgCyr_minCI + (grid_output$mean_annual_nee_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_south_PgCyr_lowCI = NEE_south_PgCyr_lowCI + (grid_output$mean_annual_nee_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_south_PgCyr_highCI = NEE_south_PgCyr_highCI + (grid_output$mean_annual_nee_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NEE_south_PgCyr_maxCI = NEE_south_PgCyr_maxCI + (grid_output$mean_annual_nee_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_south_PgCyr = NPP_south_PgCyr + (grid_output$mean_annual_npp_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_south_PgCyr_minCI = NPP_south_PgCyr_minCI + (grid_output$mean_annual_npp_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_south_PgCyr_lowCI = NPP_south_PgCyr_lowCI + (grid_output$mean_annual_npp_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_south_PgCyr_highCI = NPP_south_PgCyr_highCI + (grid_output$mean_annual_npp_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NPP_south_PgCyr_maxCI = NPP_south_PgCyr_maxCI + (grid_output$mean_annual_npp_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
             RHET_LIT_south_PgC = RHET_LIT_south_PgC + (grid_output$mean_annual_rhet_litter_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_south_PgC_minCI = RHET_LIT_south_PgC_minCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_south_PgC_lowCI = RHET_LIT_south_PgC_lowCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_south_PgC_highCI = RHET_LIT_south_PgC_highCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_LIT_south_PgC_maxCI = RHET_LIT_south_PgC_maxCI + (grid_output$mean_annual_rhet_litter_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_south_PgC = RHET_SOM_south_PgC + (grid_output$mean_annual_rhet_som_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_south_PgC_minCI = RHET_SOM_south_PgC_minCI + (grid_output$mean_annual_rhet_som_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_south_PgC_lowCI = RHET_SOM_south_PgC_lowCI + (grid_output$mean_annual_rhet_som_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_south_PgC_highCI = RHET_SOM_south_PgC_highCI + (grid_output$mean_annual_rhet_som_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             RHET_SOM_south_PgC_maxCI = RHET_SOM_south_PgC_maxCI + (grid_output$mean_annual_rhet_som_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
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
             HARV_south_PgCyr = HARV_south_PgCyr + (grid_output$mean_annual_harvest_gCm2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_south_PgCyr_minCI = HARV_south_PgCyr_minCI + (grid_output$mean_annual_harvest_gCm2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_south_PgCyr_lowCI = HARV_south_PgCyr_lowCI + (grid_output$mean_annual_harvest_gCm2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_south_PgCyr_highCI = HARV_south_PgCyr_highCI + (grid_output$mean_annual_harvest_gCm2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             HARV_south_PgCyr_maxCI = HARV_south_PgCyr_maxCI + (grid_output$mean_annual_harvest_gCm2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAI_south_m2m2 = LAI_south_m2m2 + (grid_output$mean_annual_lai_m2m2[n,max_quant,])         
             LAI_south_m2m2_minCI = LAI_south_m2m2_minCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])         
             LAI_south_m2m2_lowCI = LAI_south_m2m2_lowCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])         
             LAI_south_m2m2_highCI = LAI_south_m2m2_highCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])         
             LAI_south_m2m2_maxCI = LAI_south_m2m2_maxCI + (grid_output$mean_annual_lai_m2m2[n,max_quant,])         
             BIO_south_PgC = BIO_south_PgC + (grid_output$mean_annual_biomass_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_south_PgC_minCI = BIO_south_PgC_minCI + (grid_output$mean_annual_biomass_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_south_PgC_lowCI = BIO_south_PgC_lowCI + (grid_output$mean_annual_biomass_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_south_PgC_highCI = BIO_south_PgC_highCI + (grid_output$mean_annual_biomass_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             BIO_south_PgC_maxCI = BIO_south_PgC_maxCI + (grid_output$mean_annual_biomass_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_south_PgC = DOM_south_PgC + (grid_output$mean_annual_dom_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_south_PgC_minCI = DOM_south_PgC_minCI + (grid_output$mean_annual_dom_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_south_PgC_lowCI = DOM_south_PgC_lowCI + (grid_output$mean_annual_dom_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_south_PgC_highCI = DOM_south_PgC_highCI + (grid_output$mean_annual_dom_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             DOM_south_PgC_maxCI = DOM_south_PgC_maxCI + (grid_output$mean_annual_dom_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ET_south_PgH2Oyr = ET_south_PgH2Oyr + (grid_output$mean_annual_ET_kgH2Om2day[n,mid_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_south_PgH2Oyr_minCI = ET_south_PgH2Oyr_minCI + (grid_output$mean_annual_ET_kgH2Om2day[n,min_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_south_PgH2Oyr_lowCI = ET_south_PgH2Oyr_lowCI + (grid_output$mean_annual_ET_kgH2Om2day[n,low_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_south_PgH2Oyr_highCI = ET_south_PgH2Oyr_highCI + (grid_output$mean_annual_ET_kgH2Om2day[n,high_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             ET_south_PgH2Oyr_maxCI = ET_south_PgH2Oyr_maxCI + (grid_output$mean_annual_ET_kgH2Om2day[n,max_quant,]*deltat*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-12)
             LAB_south_PgC = LAB_south_PgC + (grid_output$mean_annual_labile_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_south_PgC_minCI = LAB_south_PgC_minCI + (grid_output$mean_annual_labile_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_south_PgC_lowCI = LAB_south_PgC_lowCI + (grid_output$mean_annual_labile_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_south_PgC_highCI = LAB_south_PgC_highCI + (grid_output$mean_annual_labile_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LAB_south_PgC_maxCI = LAB_south_PgC_maxCI + (grid_output$mean_annual_labile_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_south_PgC = FOL_south_PgC + (grid_output$mean_annual_foliage_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_south_PgC_minCI = FOL_south_PgC_minCI + (grid_output$mean_annual_foliage_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_south_PgC_lowCI = FOL_south_PgC_lowCI + (grid_output$mean_annual_foliage_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_south_PgC_highCI = FOL_south_PgC_highCI + (grid_output$mean_annual_foliage_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             FOL_south_PgC_maxCI = FOL_south_PgC_maxCI + (grid_output$mean_annual_foliage_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_south_PgC = ROOT_south_PgC + (grid_output$mean_annual_roots_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_south_PgC_minCI = ROOT_south_PgC_minCI + (grid_output$mean_annual_roots_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_south_PgC_lowCI = ROOT_south_PgC_lowCI + (grid_output$mean_annual_roots_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_south_PgC_highCI = ROOT_south_PgC_highCI + (grid_output$mean_annual_roots_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             ROOT_south_PgC_maxCI = ROOT_south_PgC_maxCI + (grid_output$mean_annual_roots_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_south_PgC = WOOD_south_PgC + (grid_output$mean_annual_wood_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_south_PgC_minCI = WOOD_south_PgC_minCI + (grid_output$mean_annual_wood_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_south_PgC_lowCI = WOOD_south_PgC_lowCI + (grid_output$mean_annual_wood_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_south_PgC_highCI = WOOD_south_PgC_highCI + (grid_output$mean_annual_wood_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             WOOD_south_PgC_maxCI = WOOD_south_PgC_maxCI + (grid_output$mean_annual_wood_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_south_PgC = LIT_south_PgC + (grid_output$mean_annual_litter_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_south_PgC_minCI = LIT_south_PgC_minCI + (grid_output$mean_annual_litter_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_south_PgC_lowCI = LIT_south_PgC_lowCI + (grid_output$mean_annual_litter_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_south_PgC_highCI = LIT_south_PgC_highCI + (grid_output$mean_annual_litter_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             LIT_south_PgC_maxCI = LIT_south_PgC_maxCI + (grid_output$mean_annual_litter_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_south_PgC = SOM_south_PgC + (grid_output$mean_annual_som_gCm2[n,mid_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_south_PgC_minCI = SOM_south_PgC_minCI + (grid_output$mean_annual_som_gCm2[n,min_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_south_PgC_lowCI = SOM_south_PgC_lowCI + (grid_output$mean_annual_som_gCm2[n,low_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_south_PgC_highCI = SOM_south_PgC_highCI + (grid_output$mean_annual_som_gCm2[n,high_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             SOM_south_PgC_maxCI = SOM_south_PgC_maxCI + (grid_output$mean_annual_som_gCm2[n,max_quant,]*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         }

     } # Does the file exist / has it been processed

} # site loop

# Median
LAI_global_m2m2 = LAI_global_m2m2 / nos_global
LAI_north_m2m2 = LAI_north_m2m2 / nos_north
LAI_tropics_m2m2 = LAI_tropics_m2m2 / nos_tropics
LAI_south_m2m2 = LAI_south_m2m2 / nos_south
# Minimum CI
LAI_global_m2m2_minCI = LAI_global_m2m2_minCI / nos_global
LAI_north_m2m2_minCI = LAI_north_m2m2_minCI / nos_north
LAI_tropics_m2m2_minCI = LAI_tropics_m2m2_minCI / nos_tropics
LAI_south_m2m2_minCI = LAI_south_m2m2_minCI / nos_south
# Lower CI
LAI_global_m2m2_lowCI = LAI_global_m2m2_lowCI / nos_global
LAI_north_m2m2_lowCI = LAI_north_m2m2_lowCI / nos_north
LAI_tropics_m2m2_lowCI = LAI_tropics_m2m2_lowCI / nos_tropics
LAI_south_m2m2_lowCI = LAI_south_m2m2_lowCI / nos_south
# Upper CI
LAI_global_m2m2_highCI = LAI_global_m2m2_highCI / nos_global
LAI_north_m2m2_highCI = LAI_north_m2m2_highCI / nos_north
LAI_tropics_m2m2_highCI = LAI_tropics_m2m2_highCI / nos_tropics
LAI_south_m2m2_highCI = LAI_south_m2m2_highCI / nos_south
# Maximum CI
LAI_global_m2m2_maxCI = LAI_global_m2m2_maxCI / nos_global
LAI_north_m2m2_maxCI = LAI_north_m2m2_maxCI / nos_north
LAI_tropics_m2m2_maxCI = LAI_tropics_m2m2_maxCI / nos_tropics
LAI_south_m2m2_maxCI = LAI_south_m2m2_maxCI / nos_south

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
## Combine NBE together the data into an output variable

output = data.frame(Year = years, Global = NBE_global_PgCyr, North = NBE_north_PgCyr, Tropics = NBE_tropics_PgCyr, South = NBE_south_PgCyr,
                                  Global_2.5pc = NBE_global_PgCyr_minCI, North_2.5pc = NBE_north_PgCyr_minCI, Tropics_2.5pc = NBE_tropics_PgCyr_minCI, South_2.5pc = NBE_south_PgCyr_minCI,
                                  Global_25pc = NBE_global_PgCyr_lowCI, North_25pc = NBE_north_PgCyr_lowCI, Tropics_25pc = NBE_tropics_PgCyr_lowCI, South_25pc = NBE_south_PgCyr_lowCI,
                                  Global_75pc = NBE_global_PgCyr_highCI, North_75pc = NBE_north_PgCyr_highCI, Tropics_75pc = NBE_tropics_PgCyr_highCI, South_75pc = NBE_south_PgCyr_highCI,
                                  Global_97.5pc = NBE_global_PgCyr_maxCI, North_97.5pc = NBE_north_PgCyr_maxCI, Tropics_97.5pc = NBE_tropics_PgCyr_maxCI, South_97.5pc = NBE_south_PgCyr_maxCI)
                                  
###
## Begin writing out NBE to file

write.table(output,file = paste(out_dir,"/",output_prefix,"NBE",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine NEE together the data into an output variable

output = data.frame(Year = years, Global = NEE_global_PgCyr, North = NEE_north_PgCyr, Tropics = NEE_tropics_PgCyr, South = NEE_south_PgCyr,
                                  Global_2.5pc = NEE_global_PgCyr_minCI, North_2.5pc = NEE_north_PgCyr_minCI, Tropics_2.5pc = NEE_tropics_PgCyr_minCI, South_2.5pc = NEE_south_PgCyr_minCI,
                                  Global_25pc = NEE_global_PgCyr_lowCI, North_25pc = NEE_north_PgCyr_lowCI, Tropics_25pc = NEE_tropics_PgCyr_lowCI, South_25pc = NEE_south_PgCyr_lowCI,
                                  Global_75pc = NEE_global_PgCyr_highCI, North_75pc = NEE_north_PgCyr_highCI, Tropics_75pc = NEE_tropics_PgCyr_highCI, South_75pc = NEE_south_PgCyr_highCI,
                                  Global_97.5pc = NEE_global_PgCyr_maxCI, North_97.5pc = NEE_north_PgCyr_maxCI, Tropics_97.5pc = NEE_tropics_PgCyr_maxCI, South_97.5pc = NEE_south_PgCyr_maxCI)
                                  
###
## Begin writing out NEE to file

write.table(output,file = paste(out_dir,"/",output_prefix,"NEE",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine NPP together the data into an output variable

output = data.frame(Year = years, Global = NPP_global_PgCyr, North = NPP_north_PgCyr, Tropics = NPP_tropics_PgCyr, South = NPP_south_PgCyr,
                                  Global_2.5pc = NPP_global_PgCyr_minCI, North_2.5pc = NPP_north_PgCyr_minCI, Tropics_2.5pc = NPP_tropics_PgCyr_minCI, South_2.5pc = NPP_south_PgCyr_minCI,
                                  Global_25pc = NPP_global_PgCyr_lowCI, North_25pc = NPP_north_PgCyr_lowCI, Tropics_25pc = NPP_tropics_PgCyr_lowCI, South_25pc = NPP_south_PgCyr_lowCI,
                                  Global_75pc = NPP_global_PgCyr_highCI, North_75pc = NPP_north_PgCyr_highCI, Tropics_75pc = NPP_tropics_PgCyr_highCI, South_75pc = NPP_south_PgCyr_highCI,
                                  Global_97.5pc = NPP_global_PgCyr_maxCI, North_97.5pc = NPP_north_PgCyr_maxCI, Tropics_97.5pc = NPP_tropics_PgCyr_maxCI, South_97.5pc = NPP_south_PgCyr_maxCI)
                                  
###
## Begin writing out NPP to file

write.table(output,file = paste(out_dir,"/",output_prefix,"NPP",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

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
## Combine RHET_LIT together the data into an output variable

output = data.frame(Year = years, Global = RHET_LIT_global_PgC, North = RHET_LIT_north_PgC, Tropics = RHET_LIT_tropics_PgC, South = RHET_LIT_south_PgC,
                                  Global_2.5pc = RHET_LIT_global_PgC_minCI, North_2.5pc = RHET_LIT_north_PgC_minCI, Tropics_2.5pc = RHET_LIT_tropics_PgC_minCI, South_2.5pc = RHET_LIT_south_PgC_minCI,
                                  Global_25pc = RHET_LIT_global_PgC_lowCI, North_25pc = RHET_LIT_north_PgC_lowCI, Tropics_25pc = RHET_LIT_tropics_PgC_lowCI, South_25pc = RHET_LIT_south_PgC_lowCI,
                                  Global_75pc = RHET_LIT_global_PgC_highCI, North_75pc = RHET_LIT_north_PgC_highCI, Tropics_75pc = RHET_LIT_tropics_PgC_highCI, South_75pc = RHET_LIT_south_PgC_highCI,
                                  Global_97.5pc = RHET_LIT_global_PgC_maxCI, North_97.5pc = RHET_LIT_north_PgC_maxCI, Tropics_97.5pc = RHET_LIT_tropics_PgC_maxCI, South_97.5pc = RHET_LIT_south_PgC_maxCI)
                                  
###
## Begin writing out RHET_LIT to file

write.table(output,file = paste(out_dir,"/",output_prefix,"RHET_LIT",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine RHET_SOM together the data into an output variable

output = data.frame(Year = years, Global = RHET_SOM_global_PgC, North = RHET_SOM_north_PgC, Tropics = RHET_SOM_tropics_PgC, South = RHET_SOM_south_PgC,
                                  Global_2.5pc = RHET_SOM_global_PgC_minCI, North_2.5pc = RHET_SOM_north_PgC_minCI, Tropics_2.5pc = RHET_SOM_tropics_PgC_minCI, South_2.5pc = RHET_SOM_south_PgC_minCI,
                                  Global_25pc = RHET_SOM_global_PgC_lowCI, North_25pc = RHET_SOM_north_PgC_lowCI, Tropics_25pc = RHET_SOM_tropics_PgC_lowCI, South_25pc = RHET_SOM_south_PgC_lowCI,
                                  Global_75pc = RHET_SOM_global_PgC_highCI, North_75pc = RHET_SOM_north_PgC_highCI, Tropics_75pc = RHET_SOM_tropics_PgC_highCI, South_75pc = RHET_SOM_south_PgC_highCI,
                                  Global_97.5pc = RHET_SOM_global_PgC_maxCI, North_97.5pc = RHET_SOM_north_PgC_maxCI, Tropics_97.5pc = RHET_SOM_tropics_PgC_maxCI, South_97.5pc = RHET_SOM_south_PgC_maxCI)
                                  
###
## Begin writing out RHET_SOM to file

write.table(output,file = paste(out_dir,"/",output_prefix,"RHET_SOM",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

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

###
## Combine HARV together the data into an output variable

output = data.frame(Year = years, Global = HARV_global_PgCyr, North = HARV_north_PgCyr, Tropics = HARV_tropics_PgCyr, South = HARV_south_PgCyr,
                                  Global_2.5pc = HARV_global_PgCyr_minCI, North_2.5pc = HARV_north_PgCyr_minCI, Tropics_2.5pc = HARV_tropics_PgCyr_minCI, South_2.5pc = HARV_south_PgCyr_minCI,
                                  Global_25pc = HARV_global_PgCyr_lowCI, North_25pc = HARV_north_PgCyr_lowCI, Tropics_25pc = HARV_tropics_PgCyr_lowCI, South_25pc = HARV_south_PgCyr_lowCI,
                                  Global_75pc = HARV_global_PgCyr_highCI, North_75pc = HARV_north_PgCyr_highCI, Tropics_75pc = HARV_tropics_PgCyr_highCI, South_75pc = HARV_south_PgCyr_highCI,
                                  Global_97.5pc = HARV_global_PgCyr_maxCI, North_97.5pc = HARV_north_PgCyr_maxCI, Tropics_97.5pc = HARV_tropics_PgCyr_maxCI, South_97.5pc = HARV_south_PgCyr_maxCI)

###
## Begin writing out HARV to file

write.table(output,file = paste(out_dir,"/",output_prefix,"HARV",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine LAI together the data into an output variable

output = data.frame(Year = years, Global = LAI_global_m2m2, North = LAI_north_m2m2, Tropics = LAI_tropics_m2m2, South = LAI_south_m2m2,
                                  Global_2.5pc = LAI_global_m2m2_minCI, North_2.5pc = LAI_north_m2m2_minCI, Tropics_2.5pc = LAI_tropics_m2m2_minCI, South_2.5pc = LAI_south_m2m2_minCI,
                                  Global_25pc = LAI_global_m2m2_lowCI, North_25pc = LAI_north_m2m2_lowCI, Tropics_25pc = LAI_tropics_m2m2_lowCI, South_25pc = LAI_south_m2m2_lowCI,
                                  Global_75pc = LAI_global_m2m2_highCI, North_75pc = LAI_north_m2m2_highCI, Tropics_75pc = LAI_tropics_m2m2_highCI, South_75pc = LAI_south_m2m2_highCI,
                                  Global_97.5pc = LAI_global_m2m2_maxCI, North_97.5pc = LAI_north_m2m2_maxCI, Tropics_97.5pc = LAI_tropics_m2m2_maxCI, South_97.5pc = LAI_south_m2m2_maxCI)
                                  
###
## Begin writing out LAI to file

write.table(output,file = paste(out_dir,"/",output_prefix,"LAI",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine BIO together the data into an output variable

output = data.frame(Year = years, Global = BIO_global_PgC, North = BIO_north_PgC, Tropics = BIO_tropics_PgC, South = BIO_south_PgC,
                                  Global_2.5pc = BIO_global_PgC_minCI, North_2.5pc = BIO_north_PgC_minCI, Tropics_2.5pc = BIO_tropics_PgC_minCI, South_2.5pc = BIO_south_PgC_minCI,
                                  Global_25pc = BIO_global_PgC_lowCI, North_25pc = BIO_north_PgC_lowCI, Tropics_25pc = BIO_tropics_PgC_lowCI, South_25pc = BIO_south_PgC_lowCI,
                                  Global_75pc = BIO_global_PgC_highCI, North_75pc = BIO_north_PgC_highCI, Tropics_75pc = BIO_tropics_PgC_highCI, South_75pc = BIO_south_PgC_highCI,
                                  Global_97.5pc = BIO_global_PgC_maxCI, North_97.5pc = BIO_north_PgC_maxCI, Tropics_97.5pc = BIO_tropics_PgC_maxCI, South_97.5pc = BIO_south_PgC_maxCI)

###
## Begin writing out BIO to file

write.table(output,file = paste(out_dir,"/",output_prefix,"BIO",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine DOM together the data into an output variable

output = data.frame(Year = years, Global = DOM_global_PgC, North = DOM_north_PgC, Tropics = DOM_tropics_PgC, South = DOM_south_PgC,
                                  Global_2.5pc = DOM_global_PgC_minCI, North_2.5pc = DOM_north_PgC_minCI, Tropics_2.5pc = DOM_tropics_PgC_minCI, South_2.5pc = DOM_south_PgC_minCI,
                                  Global_25pc = DOM_global_PgC_lowCI, North_25pc = DOM_north_PgC_lowCI, Tropics_25pc = DOM_tropics_PgC_lowCI, South_25pc = DOM_south_PgC_lowCI,
                                  Global_75pc = DOM_global_PgC_highCI, North_75pc = DOM_north_PgC_highCI, Tropics_75pc = DOM_tropics_PgC_highCI, South_75pc = DOM_south_PgC_highCI,
                                  Global_97.5pc = DOM_global_PgC_maxCI, North_97.5pc = DOM_north_PgC_maxCI, Tropics_97.5pc = DOM_tropics_PgC_maxCI, South_97.5pc = DOM_south_PgC_maxCI)

###
## Begin writing out DOM to file

write.table(output,file = paste(out_dir,"/",output_prefix,"DOM",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine ET together the data into an output variable

output = data.frame(Year = years, Global = ET_global_PgH2Oyr, North = ET_north_PgH2Oyr, Tropics = ET_tropics_PgH2Oyr, South = ET_south_PgH2Oyr,
                                  Global_2.5pc = ET_global_PgH2Oyr_minCI, North_2.5pc = ET_north_PgH2Oyr_minCI, Tropics_2.5pc = ET_tropics_PgH2Oyr_minCI, South_2.5pc = ET_south_PgH2Oyr_minCI,
                                  Global_25pc = ET_global_PgH2Oyr_lowCI, North_25pc = ET_north_PgH2Oyr_lowCI, Tropics_25pc = ET_tropics_PgH2Oyr_lowCI, South_25pc = ET_south_PgH2Oyr_lowCI,
                                  Global_75pc = ET_global_PgH2Oyr_highCI, North_75pc = ET_north_PgH2Oyr_highCI, Tropics_75pc = ET_tropics_PgH2Oyr_highCI, South_75pc = ET_south_PgH2Oyr_highCI,
                                  Global_97.5pc = ET_global_PgH2Oyr_maxCI, North_97.5pc = ET_north_PgH2Oyr_maxCI, Tropics_97.5pc = ET_tropics_PgH2Oyr_maxCI, South_97.5pc = ET_south_PgH2Oyr_maxCI)

###
## Begin writing out ET to file

write.table(output,file = paste(out_dir,"/",output_prefix,"ET",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine LAB together the data into an output variable

output = data.frame(Year = years, Global = LAB_global_PgC, North = LAB_north_PgC, Tropics = LAB_tropics_PgC, South = LAB_south_PgC,
                                  Global_2.5pc = LAB_global_PgC_minCI, North_2.5pc = LAB_north_PgC_minCI, Tropics_2.5pc = LAB_tropics_PgC_minCI, South_2.5pc = LAB_south_PgC_minCI,
                                  Global_25pc = LAB_global_PgC_lowCI, North_25pc = LAB_north_PgC_lowCI, Tropics_25pc = LAB_tropics_PgC_lowCI, South_25pc = LAB_south_PgC_lowCI,
                                  Global_75pc = LAB_global_PgC_highCI, North_75pc = LAB_north_PgC_highCI, Tropics_75pc = LAB_tropics_PgC_highCI, South_75pc = LAB_south_PgC_highCI,
                                  Global_97.5pc = LAB_global_PgC_maxCI, North_97.5pc = LAB_north_PgC_maxCI, Tropics_97.5pc = LAB_tropics_PgC_maxCI, South_97.5pc = LAB_south_PgC_maxCI)
                                  
###
## Begin writing out LAB to file

write.table(output,file = paste(out_dir,"/",output_prefix,"LAB",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine FOL together the data into an output variable

output = data.frame(Year = years, Global = FOL_global_PgC, North = FOL_north_PgC, Tropics = FOL_tropics_PgC, South = FOL_south_PgC,
                                  Global_2.5pc = FOL_global_PgC_minCI, North_2.5pc = FOL_north_PgC_minCI, Tropics_2.5pc = FOL_tropics_PgC_minCI, South_2.5pc = FOL_south_PgC_minCI,
                                  Global_25pc = FOL_global_PgC_lowCI, North_25pc = FOL_north_PgC_lowCI, Tropics_25pc = FOL_tropics_PgC_lowCI, South_25pc = FOL_south_PgC_lowCI,
                                  Global_75pc = FOL_global_PgC_highCI, North_75pc = FOL_north_PgC_highCI, Tropics_75pc = FOL_tropics_PgC_highCI, South_75pc = FOL_south_PgC_highCI,
                                  Global_97.5pc = FOL_global_PgC_maxCI, North_97.5pc = FOL_north_PgC_maxCI, Tropics_97.5pc = FOL_tropics_PgC_maxCI, South_97.5pc = FOL_south_PgC_maxCI)
                                  
###
## Begin writing out FOL to file

write.table(output,file = paste(out_dir,"/",output_prefix,"FOL",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine ROOT together the data into an output variable

output = data.frame(Year = years, Global = ROOT_global_PgC, North = ROOT_north_PgC, Tropics = ROOT_tropics_PgC, South = ROOT_south_PgC,
                                  Global_2.5pc = ROOT_global_PgC_minCI, North_2.5pc = ROOT_north_PgC_minCI, Tropics_2.5pc = ROOT_tropics_PgC_minCI, South_2.5pc = ROOT_south_PgC_minCI,
                                  Global_25pc = ROOT_global_PgC_lowCI, North_25pc = ROOT_north_PgC_lowCI, Tropics_25pc = ROOT_tropics_PgC_lowCI, South_25pc = ROOT_south_PgC_lowCI,
                                  Global_75pc = ROOT_global_PgC_highCI, North_75pc = ROOT_north_PgC_highCI, Tropics_75pc = ROOT_tropics_PgC_highCI, South_75pc = ROOT_south_PgC_highCI,
                                  Global_97.5pc = ROOT_global_PgC_maxCI, North_97.5pc = ROOT_north_PgC_maxCI, Tropics_97.5pc = ROOT_tropics_PgC_maxCI, South_97.5pc = ROOT_south_PgC_maxCI)
                                  
###
## Begin writing out ROOT to file

write.table(output,file = paste(out_dir,"/",output_prefix,"ROOT",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine WOOD together the data into an output variable

output = data.frame(Year = years, Global = WOOD_global_PgC, North = WOOD_north_PgC, Tropics = WOOD_tropics_PgC, South = WOOD_south_PgC,
                                  Global_2.5pc = WOOD_global_PgC_minCI, North_2.5pc = WOOD_north_PgC_minCI, Tropics_2.5pc = WOOD_tropics_PgC_minCI, South_2.5pc = WOOD_south_PgC_minCI,
                                  Global_25pc = WOOD_global_PgC_lowCI, North_25pc = WOOD_north_PgC_lowCI, Tropics_25pc = WOOD_tropics_PgC_lowCI, South_25pc = WOOD_south_PgC_lowCI,
                                  Global_75pc = WOOD_global_PgC_highCI, North_75pc = WOOD_north_PgC_highCI, Tropics_75pc = WOOD_tropics_PgC_highCI, South_75pc = WOOD_south_PgC_highCI,
                                  Global_97.5pc = WOOD_global_PgC_maxCI, North_97.5pc = WOOD_north_PgC_maxCI, Tropics_97.5pc = WOOD_tropics_PgC_maxCI, South_97.5pc = WOOD_south_PgC_maxCI)
                                  
###
## Begin writing out WOOD to file

write.table(output,file = paste(out_dir,"/",output_prefix,"WOOD",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine LIT together the data into an output variable

output = data.frame(Year = years, Global = LIT_global_PgC, North = LIT_north_PgC, Tropics = LIT_tropics_PgC, South = LIT_south_PgC,
                                  Global_2.5pc = LIT_global_PgC_minCI, North_2.5pc = LIT_north_PgC_minCI, Tropics_2.5pc = LIT_tropics_PgC_minCI, South_2.5pc = LIT_south_PgC_minCI,
                                  Global_25pc = LIT_global_PgC_lowCI, North_25pc = LIT_north_PgC_lowCI, Tropics_25pc = LIT_tropics_PgC_lowCI, South_25pc = LIT_south_PgC_lowCI,
                                  Global_75pc = LIT_global_PgC_highCI, North_75pc = LIT_north_PgC_highCI, Tropics_75pc = LIT_tropics_PgC_highCI, South_75pc = LIT_south_PgC_highCI,
                                  Global_97.5pc = LIT_global_PgC_maxCI, North_97.5pc = LIT_north_PgC_maxCI, Tropics_97.5pc = LIT_tropics_PgC_maxCI, South_97.5pc = LIT_south_PgC_maxCI)
                                  
###
## Begin writing out LIT to file

write.table(output,file = paste(out_dir,"/",output_prefix,"LIT",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

###
## Combine SOM together the data into an output variable

output = data.frame(Year = years, Global = SOM_global_PgC, North = SOM_north_PgC, Tropics = SOM_tropics_PgC, South = SOM_south_PgC,
                                  Global_2.5pc = SOM_global_PgC_minCI, North_2.5pc = SOM_north_PgC_minCI, Tropics_2.5pc = SOM_tropics_PgC_minCI, South_2.5pc = SOM_south_PgC_minCI,
                                  Global_25pc = SOM_global_PgC_lowCI, North_25pc = SOM_north_PgC_lowCI, Tropics_25pc = SOM_tropics_PgC_lowCI, South_25pc = SOM_south_PgC_lowCI,
                                  Global_75pc = SOM_global_PgC_highCI, North_75pc = SOM_north_PgC_highCI, Tropics_75pc = SOM_tropics_PgC_highCI, South_75pc = SOM_south_PgC_highCI,
                                  Global_97.5pc = SOM_global_PgC_maxCI, North_97.5pc = SOM_north_PgC_maxCI, Tropics_97.5pc = SOM_tropics_PgC_maxCI, South_97.5pc = SOM_south_PgC_maxCI)
                                  
###
## Begin writing out SOM to file

write.table(output,file = paste(out_dir,"/",output_prefix,"SOM",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)

