
###
## Generate some basic plots resulting from the gridded future simulation analyses
## Author: T. Luke Smallman (t.l.smallman@ed.ac.uk)
## Created: 09/10/23
## Last modified: 09/10/23
###

# Set the location of cardamom outputs
cardamom_output_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/"
# Specify the infofile for the project to be used
infofile = paste(cardamom_output_dir,"/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_AGB/infofile.RData",sep="")
#infofile = paste(cardamom_output_dir,"/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_oneAGB/infofile.RData",sep="")
# Load project for climate
load(infofile) 
orig_nos_years = length(c(as.vector(PROJECT$start_year):as.vector(PROJECT$end_year)))
  
# Suffix for the output files
output_suffix = "_steadystate" # must include "_" at the beginning
# Location to place the outputs of this script
outdir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_AGB/RESULTS_PROCESSED/"
#outdir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_oneAGB/RESULTS_PROCESSED/"

# Which Shared Socioeconomic Pathways (SSPs) do you want to use?
# Currently available in GCEL are: "ssp119","ssp126","ssp434","ssp245","ssp370","ssp585"
#ssp_scenarios = c("ssp119","ssp126","ssp434","ssp245","ssp370","ssp585")
#ssp_scenarios = c("ssp126","ssp245","ssp585")
ssp_scenarios = c("ssp245")

# Which ESM to extract climate change from?
# Currently available are: "MOHC" 
ESM = "MOHC"

# Load libraries
library(colorspace)
library(RColorBrewer)
library(zoo)
library(fields)
library(compiler)
library(terra)
source("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/read_binary_file_format.r")
source("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/generate_wgs_grid.r")
source("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/calc_pixel_area.r")

###
## Load spatial information

# generate UK or WGS-84 lat long grid
output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
grid_lat = array(output$lat, dim=c(output$long_dim,output$lat_dim))
grid_long = array(output$long,dim=c(output$long_dim,output$lat_dim))
# then generate the area estimates for each pixel
area = calc_pixel_area(grid_long,grid_lat)
# Extract useful information (output re-used later)
lat = output$lat ; long = output$long
lat_dim = output$lat_dim ; long_dim = output$long_dim
cardamom_ext = output$cardamom_ext

# Load the landmask
landmask = vect("./R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
# just to be sure enforce the projection to WGS-84
landmask = project(landmask,"EPSG:4326")
# Clip to the extent of the CARDAMOM analysis
landmask = crop(landmask, cardamom_ext)

###
## Begin loading into global time series grids
###

# Load the specific ssp file
load(paste(outdir,PROJECT$name,"_",ssp_scenarios[1],"_",ESM,output_suffix,".RData",sep=""))
# Extract number of years
nos_years = grid_output$nos_years
# Extract how many time steps per year
steps_per_year = grid_output$steps_per_year
# Determine the years being simulated
run_years = as.numeric(PROJECT$start_year) : (as.numeric(PROJECT$start_year)+nos_years-1)
land_fraction = grid_output$land_fraction

# Extract gridded information on the observations
dims = dim(grid_output$mean_lai_m2m2)
# Soil prior
SoilCPrior = array(NA, dim=c(dims[1], dims[2]))
# Mean annual LAI obs
LAIobs_m2m2 = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
LAIobs_unc_m2m2 = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
LAIcount = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
# Observed mean annual NBE information
NBECobs_gCm2day = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
NBECobs_unc_gCm2day = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
# Observed wood trends information
WoodCobs_gCm2 = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
WoodCobs_unc_gCm2 = array(NA, dim=c(dims[1],dims[2],orig_nos_years))
WoodCobs_trend_gCm2yr = array(NA, dim=c(dims[1],dims[2]))
WoodCmod_trend_gCm2yr = array(NA, dim=c(length(ssp_scenarios),dims[1],dims[2]))
# Disturbance
HarvestFraction = array(NA, dim=c(dims[1], dims[2]))
BurnedFraction = array(NA, dim=c(dims[1], dims[2]))
FireFreq = array(NA, dim=c(dims[1],dims[2]))
# Meteorology trends
mint_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
maxt_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
swrad_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
co2_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
precip_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
vpd_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
# Initialise variables for aggregate time series
lai_m2m2 = array(0, dim=c(length(ssp_scenarios),nos_years)) ; lai_lower_m2m2 = array(0, dim=c(length(ssp_scenarios),nos_years)) ; lai_upper_m2m2 = array(0, dim=c(length(ssp_scenarios),nos_years))
cica_ratio = array(0, dim=c(length(ssp_scenarios),nos_years)) ; cica_lower_ratio = array(0, dim=c(length(ssp_scenarios),nos_years)) ; cica_upper_ratio = array(0, dim=c(length(ssp_scenarios),nos_years))  
SurfWater_mm = array(0, dim=c(length(ssp_scenarios),nos_years)) ; SurfWater_lower_mm = array(0, dim=c(length(ssp_scenarios),nos_years)) ; SurfWater_upper_mm = array(0, dim=c(length(ssp_scenarios),nos_years))  
wSWP_MPa = array(0, dim=c(length(ssp_scenarios),nos_years)) ; wSWP_lower_MPa = array(0, dim=c(length(ssp_scenarios),nos_years)) ; wSWP_upper_MPa = array(0, dim=c(length(ssp_scenarios),nos_years))  
et_PgH2Oyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; et_lower_PgH2Oyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; et_upper_PgH2Oyr = array(0, dim=c(length(ssp_scenarios),nos_years))
gpp_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; gpp_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; gpp_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
npp_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; npp_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; npp_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
rauto_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; rauto_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; rauto_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
rhet_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; rhet_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; rhet_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
nee_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; nee_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; nee_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
nbe_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; nbe_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; nbe_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
fire_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; fire_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; fire_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
harvest_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; harvest_lower_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years)) ; harvest_upper_TgCyr = array(0, dim=c(length(ssp_scenarios),nos_years))
# Pool totals
wood_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))       ; wood_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))       ; wood_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
litter_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))     ; litter_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))     ; litter_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
woodlitter_TgC = array(0, dim=c(length(ssp_scenarios),nos_years)) ; woodlitter_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years)) ; woodlitter_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
soil_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))       ; soil_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))       ; soil_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
bio_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))        ; bio_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))        ; bio_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dom_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))        ; dom_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))        ; dom_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
Ctotal_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))     ; Ctotal_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))     ; Ctotal_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCtotal_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))    ; dCtotal_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))    ; dCtotal_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCdom_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))      ; dCdom_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))      ; dCdom_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCbio_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))      ; dCbio_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))      ; dCbio_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dClabile_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))   ; dClabile_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))   ; dClabile_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCfoliage_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))  ; dCfoliage_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))  ; dCfoliage_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCroots_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))    ; dCroots_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))    ; dCroots_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCwood_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))     ; dCwood_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))     ; dCwood_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dClitter_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))   ; dClitter_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))   ; dClitter_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
dCsom_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))      ; dCsom_lower_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))      ; dCsom_upper_TgC = array(0, dim=c(length(ssp_scenarios),nos_years))
# Flux trends across whole simulated period
gpp_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
rauto_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
rhet_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
lai_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
et_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
wood_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
som_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
gs_DemandSupply_trend = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
# Ensemble obs+CI intersection statistics
gpp_assim_data_overlap_fraction  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
lai_assim_data_overlap_fraction  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
nee_assim_data_overlap_fraction  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
wood_assim_data_overlap_fraction = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
soil_assim_data_overlap_fraction = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
et_assim_data_overlap_fraction   = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
nbe_assim_data_overlap_fraction  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
fire_assim_data_overlap_fraction = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
# Gridded MRTs
MTT_labile_years        = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_labile_lower_years  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_labile_upper_years  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_foliage_years       = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_foliage_lower_years = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_foliage_upper_years = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_roots_years         = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_roots_lower_years   = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_roots_upper_years   = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_wood_years          = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_wood_lower_years    = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_wood_upper_years    = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_litter_years        = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_litter_lower_years  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_litter_upper_years  = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_som_years           = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_som_lower_years     = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
MTT_som_upper_years     = array(NA, dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))

## Gridded mean annual totals
# State variables
lai_m2m2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
wood_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
litter_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
woodlitter_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
soil_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
bio_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
dom_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
Ctotal_gCm2_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
cica_ratio_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
SurfWater_mm_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
wSWP_MPa_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
# Fluxes
et_kgH2Om2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
gpp_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
npp_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
rauto_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
rhet_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
nee_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
nbe_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
fire_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
harvest_gCm2day_grid = array(NA,dim=c(length(ssp_scenarios),dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))

# Flux normalised by first year trends
#gpp_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
#rauto_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
#rhet_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
#lai_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
# Timing variable needed
time_vector = seq(0,nos_years, length.out = dim(grid_output$nee_gCm2day)[3])


###
## Create plotting colour schemes needed

# set up colour scheme
smoothScatter_colours=colorRampPalette(c("white",rep(rev(brewer.pal(11,"Spectral")),each=3)))
colour_choices_default = colorRampPalette(brewer.pal(11,"Spectral")) 
colour_choices_sign = colorRampPalette(brewer.pal(11,"PRGn"))
colour_choices_gain = colorRampPalette(brewer.pal(9,"YlGnBu"))
colour_choices_loss = colorRampPalette(brewer.pal(9,"YlOrRd"))
colour_choices_CI = colorRampPalette(brewer.pal(9,"Purples"))
# Model specific colour choices
scenario_colours = colorRampPalette(brewer.pal(8,"Dark2"))
scenario_colours = scenario_colours(4)
model_colours = colorRampPalette(brewer.pal(12,"Paired"))
model_colours = model_colours(5)
obs_colours = colorRampPalette(brewer.pal(8,"Dark2"))
obs_colours = obs_colours(4)

# array sizes are always the same so
colour_choices_default = colour_choices_default(100)
colour_choices_sign = colour_choices_sign(100)
colour_choices_gain = colour_choices_gain(100)
colour_choices_loss = colour_choices_loss(100)
colour_choices_CI = colour_choices_CI(100)

# Define counting function
counting<-function(var) {return(length(which(is.na(var) == FALSE)))}

do_loads = TRUE
if (do_loads) {
# Loop through each ssp scenario
for (ssp in seq(1, length(ssp_scenarios))) {
     # Loop through all sites
     nos_sites_inc = 0 ; cumarea = 0
     # Load the specific ssp file
     load(paste(outdir,PROJECT$name,"_",ssp_scenarios[ssp],"_",ESM,output_suffix,".RData",sep=""))

     # Specify the position within the stored ensemble for the median estimate and the desired uncertainty bands
     mid_quant = which(grid_output$num_quantiles == 0.5) 
     low_quant = which(grid_output$num_quantiles == 0.025) 
     high_quant = which(grid_output$num_quantiles == 0.975)
     wanted_quant = c(low_quant,3,mid_quant,5,high_quant)

     # Whole grid extractions
     wood_trend[ssp,,] = grid_output$final_dCwood_gCm2[,,mid_quant] / nos_years
     som_trend[ssp,,]  = grid_output$final_dCsom_gCm2[,,mid_quant] / nos_years

     gpp_assim_data_overlap_fraction[ssp,,]  = grid_output$gpp_assim_data_overlap_fraction
     lai_assim_data_overlap_fraction[ssp,,]  = grid_output$lai_assim_data_overlap_fraction
     nee_assim_data_overlap_fraction[ssp,,]  = grid_output$nee_assim_data_overlap_fraction
     wood_assim_data_overlap_fraction[ssp,,] = grid_output$wood_assim_data_overlap_fraction
     soil_assim_data_overlap_fraction[ssp,,] = grid_output$soil_assim_data_overlap_fraction
     et_assim_data_overlap_fraction[ssp,,]   = grid_output$et_assim_data_overlap_fraction
     nbe_assim_data_overlap_fraction[ssp,,]  = grid_output$nbe_assim_data_overlap_fraction
     fire_assim_data_overlap_fraction[ssp,,] = grid_output$fire_assim_data_overlap_fraction

     MTT_labile_years[ssp,,] = grid_output$MTT_labile_years[,,mid_quant]
     MTT_labile_lower_years[ssp,,] = grid_output$MTT_labile_years[,,low_quant]
     MTT_labile_upper_years[ssp,,] = grid_output$MTT_labile_years[,,high_quant]
     MTT_foliage_years[ssp,,] = grid_output$MTT_foliage_years[,,mid_quant]
     MTT_foliage_lower_years[ssp,,] = grid_output$MTT_foliage_years[,,low_quant]
     MTT_foliage_upper_years[ssp,,] = grid_output$MTT_foliage_years[,,high_quant]
     MTT_roots_years[ssp,,] = grid_output$MTT_roots_years[,,mid_quant]
     MTT_roots_lower_years[ssp,,] = grid_output$MTT_roots_years[,,low_quant]
     MTT_roots_upper_years[ssp,,] = grid_output$MTT_roots_years[,,high_quant]
     MTT_wood_years[ssp,,] = grid_output$MTT_wood_years[,,mid_quant]
     MTT_wood_lower_years[ssp,,] = grid_output$MTT_wood_years[,,low_quant]
     MTT_wood_upper_years[ssp,,] = grid_output$MTT_wood_years[,,high_quant]
     MTT_litter_years[ssp,,] = grid_output$MTT_litter_years[,,mid_quant]
     MTT_litter_lower_years[ssp,,] = grid_output$MTT_litter_years[,,low_quant]
     MTT_litter_upper_years[ssp,,] = grid_output$MTT_litter_years[,,high_quant]
     MTT_som_years[ssp,,] = grid_output$MTT_som_years[,,mid_quant]
     MTT_som_lower_years[ssp,,] = grid_output$MTT_som_years[,,low_quant]
     MTT_som_upper_years[ssp,,] = grid_output$MTT_som_years[,,high_quant]
      
     # Loop through every site to generate TIME VARYING estimates
     for (n in seq(1, PROJECT$nosites)) {

          # Extract each sites location within the grid
          i_loc = grid_output$i_location[n] ; j_loc = grid_output$j_location[n]
     
          # Check that location has run
          if (is.na(i_loc) == FALSE & is.na(j_loc) == FALSE) {

              ###
              ## Aggregating data
              
              nos_sites_inc = nos_sites_inc + 1
              # Estimate pixel level trends
              tmp1 = rollapply(grid_output$gpp_gCm2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
              tmp2 = rollapply(grid_output$rauto_gCm2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
              tmp3 = rollapply(grid_output$rhet_gCm2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
              tmp4 = rollapply(grid_output$lai_m2m2[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)     
              tmp5 = rollapply(grid_output$ET_kgH2Om2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
              tmp6 = rollapply(grid_output$gs_demand_supply_ratio[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              gpp_trend[ssp,i_loc,j_loc]   = coef(lm(tmp1 ~ c(1:nos_years)))[2] 
              rauto_trend[ssp,i_loc,j_loc] = coef(lm(tmp2 ~ c(1:nos_years)))[2] 
              rhet_trend[ssp,i_loc,j_loc]  = coef(lm(tmp3 ~ c(1:nos_years)))[2] 
              lai_trend[ssp,i_loc,j_loc]   = coef(lm(tmp4 ~ c(1:nos_years)))[2] 
              et_trend[ssp,i_loc,j_loc]    = coef(lm(tmp5 ~ c(1:nos_years)))[2] 
              gs_DemandSupply_trend[ssp,i_loc,j_loc] = coef(lm(tmp6 ~ c(1:nos_years)))[2] 

              # Grid mean annual states and fluxes
              lai_m2m2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$lai_m2m2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              wood_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$wood_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              litter_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$litter_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              if (length(which(names(grid_output) == "woodlitter_gCm2")) > 0) {
                  woodlitter_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$woodlitter_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              }
              soil_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$som_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              bio_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$biomass_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              dom_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$dom_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              Ctotal_gCm2_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$Ctotal_gCm2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              cica_ratio_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$CiCa[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              SurfWater_mm_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$SurfWater_kgH2Om2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              wSWP_MPa_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$wSWP_MPa[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              et_kgH2Om2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$ET_kgH2Om2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              gpp_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$gpp_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              npp_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$npp_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              rauto_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$rauto_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              rhet_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$rhet_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              nee_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$nee_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              nbe_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$nbe_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              fire_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$fire_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              harvest_gCm2day_grid[ssp,i_loc,j_loc,] = rollapply(grid_output$harvest_gCm2day[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
                     
              # Cumulate the total grid_output$area actually used in the analysis
              cumarea = cumarea + grid_output$area[i_loc,j_loc]
              lai_m2m2[ssp,]               = lai_m2m2[ssp,]       + lai_m2m2_grid[ssp,i_loc,j_loc,]
              lai_lower_m2m2[ssp,]         = lai_lower_m2m2[ssp,] + rollapply(grid_output$lai_m2m2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              lai_upper_m2m2[ssp,]         = lai_upper_m2m2[ssp,] + rollapply(grid_output$lai_m2m2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              # grid_output$area averaged states
              cica_ratio[ssp,]             = cica_ratio[ssp,]         + rollapply(grid_output$CiCa[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              cica_lower_ratio[ssp,]       = cica_lower_ratio[ssp,]   + rollapply(grid_output$CiCa[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              cica_upper_ratio[ssp,]       = cica_upper_ratio[ssp,]   + rollapply(grid_output$CiCa[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              SurfWater_mm[ssp,]           = SurfWater_mm[ssp,]       + rollapply(grid_output$SurfWater_kgH2Om2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              SurfWater_lower_mm[ssp,]     = SurfWater_lower_mm[ssp,] + rollapply(grid_output$SurfWater_kgH2Om2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              SurfWater_upper_mm[ssp,]     = SurfWater_upper_mm[ssp,] + rollapply(grid_output$SurfWater_kgH2Om2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              wSWP_MPa[ssp,]               = wSWP_MPa[ssp,]           + rollapply(grid_output$wSWP_MPa[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              wSWP_lower_MPa[ssp,]         = wSWP_lower_MPa[ssp,]     + rollapply(grid_output$wSWP_MPa[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              wSWP_upper_MPa[ssp,]         = wSWP_upper_MPa[ssp,]     + rollapply(grid_output$wSWP_MPa[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
              # Stocks
              wood_TgC[ssp,]               = wood_TgC[ssp,]          + rollapply(grid_output$wood_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
              wood_lower_TgC[ssp,]         = wood_lower_TgC[ssp,]    + rollapply(grid_output$wood_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              wood_upper_TgC[ssp,]         = wood_upper_TgC[ssp,]    + rollapply(grid_output$wood_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              litter_TgC[ssp,]             = litter_TgC[ssp,]        + rollapply(grid_output$litter_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
              litter_lower_TgC[ssp,]       = litter_lower_TgC[ssp,]  + rollapply(grid_output$litter_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              litter_upper_TgC[ssp,]       = litter_upper_TgC[ssp,]  + rollapply(grid_output$litter_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              if (length(which(names(grid_output) == "woodlitter_gCm2")) > 0) {
                  woodlitter_TgC[ssp,]        = woodlitter_TgC[ssp,]       + rollapply(grid_output$woodlitter_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
                  woodlitter_lower_TgC[ssp,]  = woodlitter_lower_TgC[ssp,] + rollapply(grid_output$woodlitter_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
                  woodlitter_upper_TgC[ssp,]  = woodlitter_upper_TgC[ssp,] + rollapply(grid_output$woodlitter_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              }
              soil_TgC[ssp,]               = soil_TgC[ssp,]          + rollapply(grid_output$som_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
              soil_lower_TgC[ssp,]         = soil_lower_TgC[ssp,]    + rollapply(grid_output$som_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              soil_upper_TgC[ssp,]         = soil_upper_TgC[ssp,]    + rollapply(grid_output$som_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              bio_TgC[ssp,]                = bio_TgC[ssp,]           + rollapply(grid_output$biomass_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
              bio_lower_TgC[ssp,]          = bio_lower_TgC[ssp,]     + rollapply(grid_output$biomass_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              bio_upper_TgC[ssp,]          = bio_upper_TgC[ssp,]     + rollapply(grid_output$biomass_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              dom_TgC[ssp,]                = dom_TgC[ssp,]           + rollapply(grid_output$dom_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
              dom_lower_TgC[ssp,]          = dom_lower_TgC[ssp,]     + rollapply(grid_output$dom_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              dom_upper_TgC[ssp,]          = dom_upper_TgC[ssp,]     + rollapply(grid_output$dom_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              Ctotal_TgC[ssp,]             = Ctotal_TgC[ssp,]        + rollapply(grid_output$Ctotal_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
              Ctotal_lower_TgC[ssp,]       = Ctotal_lower_TgC[ssp,]  + rollapply(grid_output$Ctotal_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              Ctotal_upper_TgC[ssp,]       = Ctotal_upper_TgC[ssp,]  + rollapply(grid_output$Ctotal_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
              # For the C stock change variables, determine the correct sub-sampling we want to achieve
              tmp = cumsum(rep(steps_per_year,grid_output$nos_year))
              dCtotal_TgC[ssp,]            = dCtotal_TgC[ssp,]        + (grid_output$dCtotal_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCtotal_lower_TgC[ssp,]      = dCtotal_lower_TgC[ssp,]  + (grid_output$dCtotal_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCtotal_upper_TgC[ssp,]      = dCtotal_upper_TgC[ssp,]  + (grid_output$dCtotal_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCdom_TgC[ssp,]              = dCdom_TgC[ssp,]          + (grid_output$dCdom_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCdom_lower_TgC[ssp,]        = dCdom_lower_TgC[ssp,]    + (grid_output$dCdom_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCdom_upper_TgC[ssp,]        = dCdom_upper_TgC[ssp,]    + (grid_output$dCdom_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCbio_TgC[ssp,]              = dCbio_TgC[ssp,]          + (grid_output$dCbiomass_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCbio_lower_TgC[ssp,]        = dCbio_lower_TgC[ssp,]    + (grid_output$dCbiomass_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCbio_upper_TgC[ssp,]        = dCbio_upper_TgC[ssp,]    + (grid_output$dCbiomass_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dClabile_TgC[ssp,]           = dClabile_TgC[ssp,]       + (grid_output$dClabile_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dClabile_lower_TgC[ssp,]     = dClabile_lower_TgC[ssp,] + (grid_output$dClabile_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dClabile_upper_TgC[ssp,]     = dClabile_upper_TgC[ssp,] + (grid_output$dClabile_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCfoliage_TgC[ssp,]          = dCfoliage_TgC[ssp,]      + (grid_output$dCfoliage_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCfoliage_lower_TgC[ssp,]    = dCfoliage_lower_TgC[ssp,]+ (grid_output$dCfoliage_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCfoliage_upper_TgC[ssp,]    = dCfoliage_upper_TgC[ssp,]+ (grid_output$dCfoliage_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCroots_TgC[ssp,]            = dCroots_TgC[ssp,]        + (grid_output$dCroots_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCroots_lower_TgC[ssp,]      = dCroots_lower_TgC[ssp,]  + (grid_output$dCroots_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCroots_upper_TgC[ssp,]      = dCroots_upper_TgC[ssp,]  + (grid_output$dCroots_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCwood_TgC[ssp,]             = dCwood_TgC[ssp,]         + (grid_output$dCwood_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCwood_lower_TgC[ssp,]       = dCwood_lower_TgC[ssp,]   + (grid_output$dCwood_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCwood_upper_TgC[ssp,]       = dCwood_upper_TgC[ssp,]   + (grid_output$dCwood_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dClitter_TgC[ssp,]           = dClitter_TgC[ssp,]       + (grid_output$dClitter_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dClitter_lower_TgC[ssp,]     = dClitter_lower_TgC[ssp,] + (grid_output$dClitter_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dClitter_upper_TgC[ssp,]     = dClitter_upper_TgC[ssp,] + (grid_output$dClitter_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCsom_TgC[ssp,]              = dCsom_TgC[ssp,]          + (grid_output$dCsom_gCm2[n,mid_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCsom_lower_TgC[ssp,]        = dCsom_lower_TgC[ssp,]    + (grid_output$dCsom_gCm2[n,low_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              dCsom_upper_TgC[ssp,]        = dCsom_upper_TgC[ssp,]    + (grid_output$dCsom_gCm2[n,high_quant,tmp]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc])
              # Fluxes
              et_PgH2Oyr[ssp,]             = et_PgH2Oyr[ssp,]     + (rollapply(grid_output$ET_kgH2Om2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              et_lower_PgH2Oyr[ssp,]       = et_lower_PgH2Oyr[ssp,]    + (rollapply(grid_output$ET_kgH2Om2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              et_upper_PgH2Oyr[ssp,]       = et_upper_PgH2Oyr[ssp,]    + (rollapply(grid_output$ET_kgH2Om2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              gpp_TgCyr[ssp,]              = gpp_TgCyr[ssp,]      + (rollapply(grid_output$gpp_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              gpp_lower_TgCyr[ssp,]        = gpp_lower_TgCyr[ssp,]     + (rollapply(grid_output$gpp_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              gpp_upper_TgCyr[ssp,]        = gpp_upper_TgCyr[ssp,]     + (rollapply(grid_output$gpp_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              npp_TgCyr[ssp,]              = npp_TgCyr[ssp,]      + (rollapply(grid_output$npp_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              npp_lower_TgCyr[ssp,]        = npp_lower_TgCyr[ssp,]     + (rollapply(grid_output$npp_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              npp_upper_TgCyr[ssp,]        = npp_upper_TgCyr[ssp,]     + (rollapply(grid_output$npp_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              rauto_TgCyr[ssp,]            = rauto_TgCyr[ssp,]    + (rollapply(grid_output$rauto_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              rauto_lower_TgCyr[ssp,]      = rauto_lower_TgCyr[ssp,]   + (rollapply(grid_output$rauto_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              rauto_upper_TgCyr[ssp,]      = rauto_upper_TgCyr[ssp,]   + (rollapply(grid_output$rauto_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)         
              rhet_TgCyr[ssp,]             = rhet_TgCyr[ssp,]     + (rollapply(grid_output$rhet_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              rhet_lower_TgCyr[ssp,]       = rhet_lower_TgCyr[ssp,]    + (rollapply(grid_output$rhet_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              rhet_upper_TgCyr[ssp,]       = rhet_upper_TgCyr[ssp,]    + (rollapply(grid_output$rhet_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              nee_TgCyr[ssp,]              = nee_TgCyr[ssp,]      + (rollapply(grid_output$nee_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              nee_lower_TgCyr[ssp,]        = nee_lower_TgCyr[ssp,]     + (rollapply(grid_output$nee_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              nee_upper_TgCyr[ssp,]        = nee_upper_TgCyr[ssp,]     + (rollapply(grid_output$nee_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              nbe_TgCyr[ssp,]              = nbe_TgCyr[ssp,]      + (rollapply(grid_output$nbe_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              nbe_lower_TgCyr[ssp,]        = nbe_lower_TgCyr[ssp,]     + (rollapply(grid_output$nbe_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              nbe_upper_TgCyr[ssp,]        = nbe_upper_TgCyr[ssp,]     + (rollapply(grid_output$nbe_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              fire_TgCyr[ssp,]             = fire_TgCyr[ssp,]     + (rollapply(grid_output$fire_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              fire_lower_TgCyr[ssp,]       = fire_lower_TgCyr[ssp,]    + (rollapply(grid_output$fire_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              fire_upper_TgCyr[ssp,]       = fire_upper_TgCyr[ssp,]    + (rollapply(grid_output$fire_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              harvest_TgCyr[ssp,]          = harvest_TgCyr[ssp,]  + (rollapply(grid_output$harvest_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              harvest_lower_TgCyr[ssp,]    = harvest_lower_TgCyr[ssp,] + (rollapply(grid_output$harvest_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
              harvest_upper_TgCyr[ssp,]    = harvest_upper_TgCyr[ssp,] + (rollapply(grid_output$harvest_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*grid_output$area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)

              ###
              ## Determining Drivers and trends

              # Read in pixel driving data
              drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
              orig_nos_years = dim(drivers$met)[1] / steps_per_year
              # Determine forest harvest intensity
              HarvestFraction[i_loc,j_loc] = sum(drivers$met[,8]) / orig_nos_years
              # Determine mean annual fire intensity and frequency
              BurnedFraction[i_loc,j_loc] = sum(drivers$met[,9]) / orig_nos_years
              FireFreq[i_loc,j_loc] = length(which(drivers$met[,9] > 0)) / orig_nos_years
              # Meteorological trends for min temperature, max temperature, SW radiation, co2 precipitation, VPD
              tmp1 = rollapply(drivers$met[,2], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              mint_trend[ssp,i_loc,j_loc]   = coef(lm(tmp1 ~ c(1:orig_nos_years)))[2] 
              tmp1 = rollapply(drivers$met[,3], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              maxt_trend[ssp,i_loc,j_loc]   = coef(lm(tmp1 ~ c(1:orig_nos_years)))[2] 
              tmp1 = rollapply(drivers$met[,4], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              swrad_trend[ssp,i_loc,j_loc]  = coef(lm(tmp1 ~ c(1:orig_nos_years)))[2] 
              tmp1 = rollapply(drivers$met[,5], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              co2_trend[ssp,i_loc,j_loc]    = coef(lm(tmp1 ~ c(1:orig_nos_years)))[2] 
              tmp1 = rollapply(drivers$met[,7], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              precip_trend[ssp,i_loc,j_loc] = coef(lm(tmp1 ~ c(1:orig_nos_years)))[2] 
              tmp1 = rollapply(drivers$met[,16], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
              vpd_trend[ssp,i_loc,j_loc]    = coef(lm(tmp1 ~ c(1:orig_nos_years)))[2]              
              # Load any priors
              SoilCPrior[i_loc,j_loc] = drivers$parpriors[23] ; if (SoilCPrior[i_loc,j_loc] == -9999) {SoilCPrior[i_loc,j_loc] = NA}
              # Clear missing data from and extract observed LAI
              drivers$obs[which(drivers$obs[,3] == -9999),3] = NA ; drivers$obs[which(drivers$obs[,4] == -9999),4] = NA
              # Extract mean annual information from LAI observations
              LAIcount[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,3], width = steps_per_year, by = steps_per_year, counting)
              LAIobs_m2m2[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,3], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
              LAIobs_unc_m2m2[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,4], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
              # Clear missing information from NBE
              drivers$obs[which(drivers$obs[,35] == -9999),35] = NA ; drivers$obs[which(drivers$obs[,36] == -9999),36] = NA
              NBECobs_gCm2day[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,35], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
              NBECobs_unc_gCm2day[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,36], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
              # Clear missing data from and extract observed LAI
              drivers$obs[which(drivers$obs[,13] == -9999),13] = NA ; drivers$obs[which(drivers$obs[,14] == -9999),14] = NA
              WoodCobs_gCm2[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,13], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
              WoodCobs_unc_gCm2[i_loc,j_loc,1:orig_nos_years] = rollapply(drivers$obs[,14], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)

              # If wood stock estimate available get that too
              tmp = which(WoodCobs_gCm2[i_loc,j_loc,] > 0)
              if (length(tmp) > 0) {
                  # Wood stock trends
                  obs_period_start = tmp[1] ; obs_period_end = tmp[length(tmp)] ; obs_period_years = length(c(obs_period_start:obs_period_end))
                  WoodCobs_trend_gCm2yr[i_loc,j_loc] = coef(lm(WoodCobs_gCm2[i_loc,j_loc,tmp] ~ c(1:obs_period_years)))[2] * 12 # *12 is month to yr adjustment
                  tmp = rollapply(grid_output$wood_gCm2[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
                  WoodCmod_trend_gCm2yr[ssp,i_loc,j_loc] = coef(lm(tmp[obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12 # *12 is month to yr adjustment                                
              } # we have more than zero obs              
          } # Did this location run
     } # Site loop

     # Tidy away unwanted variables
     rm(grid_output)

     # LAI averaging
     lai_m2m2[ssp,] = lai_m2m2[ssp,] / nos_sites_inc
     lai_lower_m2m2[ssp,] = lai_lower_m2m2[ssp,] / nos_sites_inc
     lai_upper_m2m2[ssp,] = lai_upper_m2m2[ssp,] / nos_sites_inc
     # grid_output$area averaging
     cica_ratio[ssp,] = cica_ratio[ssp,] / nos_sites_inc
     cica_lower_ratio[ssp,] = cica_lower_ratio[ssp,] / nos_sites_inc
     cica_upper_ratio[ssp,] = cica_upper_ratio[ssp,] / nos_sites_inc
     SurfWater_mm[ssp,] = SurfWater_mm[ssp,] / nos_sites_inc
     SurfWater_lower_mm[ssp,] = SurfWater_lower_mm[ssp,] / nos_sites_inc
     SurfWater_upper_mm[ssp,] = SurfWater_upper_mm[ssp,] / nos_sites_inc
     wSWP_MPa[ssp,] = wSWP_MPa[ssp,] / nos_sites_inc
     wSWP_lower_MPa[ssp,] = wSWP_lower_MPa[ssp,]  / nos_sites_inc
     wSWP_upper_MPa[ssp,] = wSWP_upper_MPa[ssp,] / nos_sites_inc
     # Now adjust units gC/yr -> TgC/yr
     gpp_TgCyr[ssp,]            = gpp_TgCyr[ssp,] * 1e-12
     npp_TgCyr[ssp,]            = npp_TgCyr[ssp,] * 1e-12
     rauto_TgCyr[ssp,]          = rauto_TgCyr[ssp,] * 1e-12
     rhet_TgCyr[ssp,]           = rhet_TgCyr[ssp,] * 1e-12
     nee_TgCyr[ssp,]            = nee_TgCyr[ssp,] * 1e-12
     nbe_TgCyr[ssp,]            = nbe_TgCyr[ssp,] * 1e-12
     fire_TgCyr[ssp,]           = fire_TgCyr[ssp,] * 1e-12
     harvest_TgCyr[ssp,]        = harvest_TgCyr[ssp,] * 1e-12
     litter_TgC[ssp,]           = litter_TgC[ssp,] * 1e-12
     woodlitter_TgC[ssp,]       = woodlitter_TgC[ssp,] * 1e-12
     wood_TgC[ssp,]             = wood_TgC[ssp,] * 1e-12
     soil_TgC[ssp,]             = soil_TgC[ssp,] * 1e-12
     bio_TgC[ssp,]              = bio_TgC[ssp,] * 1e-12
     dom_TgC[ssp,]              = dom_TgC[ssp,] * 1e-12
     Ctotal_TgC[ssp,]           = Ctotal_TgC[ssp,] * 1e-12
     dCtotal_TgC[ssp,]          = dCtotal_TgC[ssp,] * 1e-12
     dCbio_TgC[ssp,]            = dCbio_TgC[ssp,] * 1e-12
     dCdom_TgC[ssp,]            = dCdom_TgC[ssp,] * 1e-12
     dClabile_TgC[ssp,]         = dClabile_TgC[ssp,] * 1e-12
     dCfoliage_TgC[ssp,]        = dCfoliage_TgC[ssp,] * 1e-12
     dCroots_TgC[ssp,]          = dCroots_TgC[ssp,] * 1e-12
     dCwood_TgC[ssp,]           = dCwood_TgC[ssp,] * 1e-12
     dClitter_TgC[ssp,]         = dClitter_TgC[ssp,] * 1e-12
     dCsom_TgC[ssp,]            = dCsom_TgC[ssp,] * 1e-12          
     # lower
     gpp_lower_TgCyr[ssp,]      = gpp_lower_TgCyr[ssp,] * 1e-12
     npp_lower_TgCyr[ssp,]      = npp_lower_TgCyr[ssp,] * 1e-12
     rauto_lower_TgCyr[ssp,]    = rauto_lower_TgCyr[ssp,] * 1e-12
     rhet_lower_TgCyr[ssp,]     = rhet_lower_TgCyr[ssp,] * 1e-12
     nee_lower_TgCyr[ssp,]      = nee_lower_TgCyr[ssp,] * 1e-12
     nbe_lower_TgCyr[ssp,]      = nbe_lower_TgCyr[ssp,] * 1e-12
     fire_lower_TgCyr[ssp,]     = fire_lower_TgCyr[ssp,] * 1e-12
     harvest_lower_TgCyr[ssp,]  = harvest_lower_TgCyr[ssp,] * 1e-12
     litter_lower_TgC[ssp,]     = litter_lower_TgC[ssp,] * 1e-12
     woodlitter_lower_TgC[ssp,] = woodlitter_lower_TgC[ssp,] * 1e-12
     wood_lower_TgC[ssp,]       = wood_lower_TgC[ssp,] * 1e-12
     soil_lower_TgC[ssp,]       = soil_lower_TgC[ssp,] * 1e-12
     bio_lower_TgC[ssp,]        = bio_lower_TgC[ssp,] * 1e-12
     dom_lower_TgC[ssp,]        = dom_lower_TgC[ssp,] * 1e-12
     Ctotal_lower_TgC[ssp,]     = Ctotal_lower_TgC[ssp,] * 1e-12
     dCtotal_lower_TgC[ssp,]    = dCtotal_lower_TgC[ssp,] * 1e-12
     dCbio_lower_TgC[ssp,]      = dCbio_lower_TgC[ssp,] * 1e-12
     dCdom_lower_TgC[ssp,]      = dCdom_lower_TgC[ssp,] * 1e-12
     dClabile_lower_TgC[ssp,]   = dClabile_lower_TgC[ssp,] * 1e-12
     dCfoliage_lower_TgC[ssp,]  = dCfoliage_lower_TgC[ssp,] * 1e-12
     dCroots_lower_TgC[ssp,]    = dCroots_lower_TgC[ssp,] * 1e-12
     dCwood_lower_TgC[ssp,]     = dCwood_lower_TgC[ssp,] * 1e-12
     dClitter_lower_TgC[ssp,]   = dClitter_lower_TgC[ssp,] * 1e-12
     dCsom_lower_TgC[ssp,]      = dCsom_lower_TgC[ssp,] * 1e-12          
     # upper
     gpp_upper_TgCyr[ssp,]      = gpp_upper_TgCyr[ssp,] * 1e-12
     npp_upper_TgCyr[ssp,]      = npp_upper_TgCyr[ssp,] * 1e-12
     rauto_upper_TgCyr[ssp,]    = rauto_upper_TgCyr[ssp,] * 1e-12
     rhet_upper_TgCyr[ssp,]     = rhet_upper_TgCyr[ssp,] * 1e-12
     nee_upper_TgCyr[ssp,]      = nee_upper_TgCyr[ssp,] * 1e-12
     nbe_upper_TgCyr[ssp,]      = nbe_upper_TgCyr[ssp,] * 1e-12
     fire_upper_TgCyr[ssp,]     = fire_upper_TgCyr[ssp,] * 1e-12
     harvest_upper_TgCyr[ssp,]  = harvest_upper_TgCyr[ssp,] * 1e-12
     litter_upper_TgC[ssp,]     = litter_upper_TgC[ssp,] * 1e-12
     woodlitter_upper_TgC[ssp,] = woodlitter_upper_TgC[ssp,] * 1e-12
     wood_upper_TgC[ssp,]       = wood_upper_TgC[ssp,] * 1e-12
     soil_upper_TgC[ssp,]       = soil_upper_TgC[ssp,] * 1e-12
     bio_upper_TgC[ssp,]        = bio_upper_TgC[ssp,] * 1e-12
     dom_upper_TgC[ssp,]        = dom_upper_TgC[ssp,] * 1e-12
     Ctotal_upper_TgC[ssp,]     = Ctotal_upper_TgC[ssp,] * 1e-12
     dCtotal_upper_TgC[ssp,]    = dCtotal_upper_TgC[ssp,] * 1e-12
     dCbio_upper_TgC[ssp,]      = dCbio_upper_TgC[ssp,] * 1e-12
     dCdom_upper_TgC[ssp,]      = dCdom_upper_TgC[ssp,] * 1e-12
     dClabile_upper_TgC[ssp,]   = dClabile_upper_TgC[ssp,] * 1e-12
     dCfoliage_upper_TgC[ssp,]  = dCfoliage_upper_TgC[ssp,] * 1e-12
     dCroots_upper_TgC[ssp,]    = dCroots_upper_TgC[ssp,] * 1e-12
     dCwood_upper_TgC[ssp,]     = dCwood_upper_TgC[ssp,] * 1e-12
     dClitter_upper_TgC[ssp,]   = dClitter_upper_TgC[ssp,] * 1e-12
     dCsom_upper_TgC[ssp,]      = dCsom_upper_TgC[ssp,] * 1e-12     
     # Now adjust units kgH2O -> PgH2O
     et_PgH2Oyr[ssp,]           = et_PgH2Oyr[ssp,] * 1e-12
     et_lower_PgH2Oyr[ssp,]     = et_lower_PgH2Oyr[ssp,] * 1e-12
     et_upper_PgH2Oyr[ssp,]     = et_upper_PgH2Oyr[ssp,] * 1e-12

} # ssp_scenario loop

###
## Assign all variables into a list object for easier manipulation

# Plotting code below assumes that the first analysis will be added to a list called "orig" and the second will be "alt"
# alt = 
# orig = 
alt = list(SoilCPrior = SoilCPrior,
            LAIobs_m2m2 = LAIobs_m2m2,
            LAIobs_unc_m2m2 = LAIobs_unc_m2m2,
            LAIcount = LAIcount,
            HarvestFraction = HarvestFraction,
            BurnedFraction = BurnedFraction,
            FireFreq = FireFreq,
            mint_trend = mint_trend, 
            maxt_trend = maxt_trend, 
            swrad_trend = swrad_trend, 
            co2_trend = co2_trend,
            precip_trend = precip_trend,
            vpd_trend = vpd_trend,
            WoodCobs_gCm2 = WoodCobs_gCm2,
            WoodCobs_unc_gCm2 = WoodCobs_unc_gCm2,
            WoodCobs_trend_gCm2yr = WoodCobs_trend_gCm2yr,
            WoodCmod_trend_gCm2yr = WoodCmod_trend_gCm2yr,             
            NBECobs_gCm2day = NBECobs_gCm2day,
            NBECobs_unc_gCm2day = NBECobs_unc_gCm2day,
            gpp_assim_data_overlap_fraction = gpp_assim_data_overlap_fraction,
            lai_assim_data_overlap_fraction = lai_assim_data_overlap_fraction,
            nee_assim_data_overlap_fraction = nee_assim_data_overlap_fraction,
            wood_assim_data_overlap_fraction = wood_assim_data_overlap_fraction,
            soil_assim_data_overlap_fraction = soil_assim_data_overlap_fraction,
            et_assim_data_overlap_fraction = et_assim_data_overlap_fraction,
            nbe_assim_data_overlap_fraction = nbe_assim_data_overlap_fraction,
            fire_assim_data_overlap_fraction = fire_assim_data_overlap_fraction,
            MTT_labile_years = MTT_labile_years, MTT_labile_lower_years = MTT_labile_lower_years, MTT_labile_upper_years = MTT_labile_upper_years,
            MTT_foliage_years = MTT_foliage_years, MTT_foliage_lower_years = MTT_foliage_lower_years, MTT_foliage_upper_years = MTT_foliage_upper_years,
            MTT_roots_years = MTT_roots_years, MTT_roots_lower_years = MTT_roots_lower_years, MTT_roots_upper_years = MTT_roots_upper_years,
            MTT_wood_years = MTT_wood_years, MTT_wood_lower_years = MTT_wood_lower_years, MTT_wood_upper_years = MTT_wood_upper_years,
            MTT_litter_years = MTT_litter_years, MTT_litter_lower_years = MTT_litter_lower_years, MTT_litter_upper_years = MTT_litter_upper_years,
            MTT_som_years = MTT_som_years, MTT_som_lower_years = MTT_som_lower_years, MTT_som_upper_years = MTT_som_upper_years,
            lai_m2m2_grid = lai_m2m2_grid,
            wood_gCm2_grid = wood_gCm2_grid,
            litter_gCm2_grid = litter_gCm2_grid,
            woodlitter_gCm2_grid = woodlitter_gCm2_grid,
            soil_gCm2_grid = soil_gCm2_grid,
            bio_gCm2_grid = bio_gCm2_grid,
            dom_gCm2_grid = dom_gCm2_grid,
            Ctotal_gCm2_grid = Ctotal_gCm2_grid,
            cica_ratio_grid = cica_ratio_grid,
            SurfWater_mm_grid = SurfWater_mm_grid,
            wSWP_MPa_grid = wSWP_MPa_grid,
            et_kgH2Om2day_grid = et_kgH2Om2day_grid,
            gpp_gCm2day_grid = gpp_gCm2day_grid,
            npp_gCm2day_grid = npp_gCm2day_grid,
            rauto_gCm2day_grid = rauto_gCm2day_grid,
            rhet_gCm2day_grid = rhet_gCm2day_grid,
            nee_gCm2day_grid = nee_gCm2day_grid,
            nbe_gCm2day_grid = nbe_gCm2day_grid,
            fire_gCm2day_grid = fire_gCm2day_grid,
            harvest_gCm2day_grid = harvest_gCm2day_grid,
            lai_m2m2 = lai_m2m2,                lai_lower_m2m2 = lai_lower_m2m2,             lai_upper_m2m2 = lai_upper_m2m2,
            cica_ratio = cica_ratio,            cica_lower_ratio = cica_lower_ratio,         cica_upper_ratio = cica_upper_ratio,
            SurfWater_mm = SurfWater_mm,        SurfWater_lower_mm = SurfWater_lower_mm,     SurfWater_upper_mm = SurfWater_upper_mm,
            wSWP_MPa = wSWP_MPa,                wSWP_lower_MPa = wSWP_lower_MPa,             wSWP_upper_MPa = wSWP_upper_MPa,
            et_PgH2Oyr = et_PgH2Oyr,            et_lower_PgH2Oyr = et_lower_PgH2Oyr,         et_upper_PgH2Oyr = et_upper_PgH2Oyr,
            gpp_TgCyr = gpp_TgCyr,              gpp_lower_TgCyr = gpp_lower_TgCyr,           gpp_upper_TgCyr = gpp_upper_TgCyr,
            npp_TgCyr = npp_TgCyr,              npp_lower_TgCyr = npp_lower_TgCyr,           npp_upper_TgCyr = npp_upper_TgCyr,
            rauto_TgCyr = rauto_TgCyr,          rauto_lower_TgCyr = rauto_lower_TgCyr,       rauto_upper_TgCyr = rauto_upper_TgCyr,
            rhet_TgCyr = rhet_TgCyr,            rhet_lower_TgCyr = rhet_lower_TgCyr,         rhet_upper_TgCyr = rhet_upper_TgCyr,
            nee_TgCyr = nee_TgCyr,              nee_lower_TgCyr = nee_lower_TgCyr,           nee_upper_TgCyr = nee_upper_TgCyr,
            nbe_TgCyr = nbe_TgCyr,              nbe_lower_TgCyr = nbe_lower_TgCyr,           nbe_upper_TgCyr = nbe_upper_TgCyr,
            fire_TgCyr = fire_TgCyr,            fire_lower_TgCyr = fire_lower_TgCyr,         fire_upper_TgCyr = fire_upper_TgCyr,
            harvest_TgCyr = harvest_TgCyr,      harvest_lower_TgCyr = harvest_lower_TgCyr,   harvest_upper_TgCyr = harvest_upper_TgCyr,
            wood_TgC = wood_TgC,                wood_lower_TgC = wood_lower_TgC,             wood_upper_TgC = wood_upper_TgC,
            litter_TgC = litter_TgC,            litter_lower_TgC = litter_lower_TgC,         litter_upper_TgC = litter_upper_TgC,
            woodlitter_TgC = woodlitter_TgC,    woodlitter_lower_TgC = woodlitter_lower_TgC, woodlitter_upper_TgC = woodlitter_upper_TgC,
            soil_TgC = soil_TgC,                soil_lower_TgC = soil_lower_TgC,             soil_upper_TgC = soil_upper_TgC,
            bio_TgC = bio_TgC,                  bio_lower_TgC = bio_lower_TgC,               bio_upper_TgC = bio_upper_TgC,
            dom_TgC = dom_TgC,                  dom_lower_TgC = dom_lower_TgC,               dom_upper_TgC = soil_upper_TgC,
            Ctotal_TgC = Ctotal_TgC,            Ctotal_lower_TgC = Ctotal_lower_TgC,         Ctotal_upper_TgC = Ctotal_upper_TgC,
            dCtotal_TgC = dCtotal_TgC,          dCtotal_lower_TgC = dCtotal_lower_TgC,       dCtotal_upper_TgC = dCtotal_upper_TgC,
            dCbio_TgC = dCbio_TgC,              dCbio_lower_TgC = dCbio_lower_TgC,           dCbio_upper_TgC = dCbio_upper_TgC,
            dCdom_TgC = dCdom_TgC,              dCdom_lower_TgC = dCdom_lower_TgC,           dCdom_upper_TgC = dCdom_upper_TgC,
            dClabile_TgC = dClabile_TgC,        dClabile_lower_TgC = dClabile_lower_TgC,     dClabile_upper_TgC = dClabile_upper_TgC,
            dCfoliage_TgC = dCfoliage_TgC,      dCfoliage_lower_TgC = dCfoliage_lower_TgC,   dCfoliage_upper_TgC = dCfoliage_upper_TgC,                                    
            dCroots_TgC = dCroots_TgC,          dCroots_lower_TgC = dCroots_lower_TgC,       dCroots_upper_TgC = dCroots_upper_TgC,            
            dCwood_TgC = dCwood_TgC,            dCwood_lower_TgC = dCwood_lower_TgC,         dCwood_upper_TgC = dCwood_upper_TgC,
            dClitter_TgC = dClitter_TgC,        dClitter_lower_TgC = dClitter_lower_TgC,     dClitter_upper_TgC = dClitter_upper_TgC,
            dCsom_TgC = dCsom_TgC,              dCsom_lower_TgC = dCsom_lower_TgC,           dCsom_upper_TgC = dCsom_upper_TgC,                                    
            gpp_trend = gpp_trend,
            rauto_trend = rauto_trend,
            rhet_trend = rhet_trend,
            lai_trend = lai_trend,
            et_trend = et_trend,
            gs_DemandSupply_trend = gs_DemandSupply_trend,
            wood_trend = wood_trend,
            som_trend = som_trend
            #gpp_trend_normalised = 
            #rauto_trend_normalised = 
            #rhet_trend_normalised = 
            #lai_trend_normalised = 
            )

}
   
do_plots = TRUE
if (do_plots) {            
###
## Do some timeseries plots...

### CREATE A CUMULATIVE NBE PLOT, + WOOD/BIO AND SOIL / DOM
### What would the figures of the paper be?
### Pulling timeseries plots for specific target regions of interest?

# Specify the labels for the two different analyses being compared
# Note that if length(alt_name) == 0 these will not be plotted
orig_name = "Single_TWC"
alt_name = "Repeat_TWC"

# Specify some colours for each SSP and correspondings between different model runs
colours_orig = rev(rainbow(length(ssp_scenarios)))
if (alt_name != "") {colours_alt = desaturate(colours_orig, amount = 0.85)}

## Plot major ecosystem fluxes over time
fig_height = 4200 ; fig_width = 2500 ; fig_res = 300
if (alt_name != "") {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_and_",alt_name,"_NBE_GPP_Fire_ET_timeseries_comparison_plusCI",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res)
} else {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_NBE_GPP_Fire_ET_timeseries_comparison_plusCI",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res) 
}
par(mfrow=c(4,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
## Plot NBE
# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$nbe_TgCyr,orig$nbe_lower_TgCyr,orig$nbe_upper_TgCyr,
                     alt$nbe_TgCyr,alt$nbe_lower_TgCyr,alt$nbe_upper_TgCyr)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$nbe_TgCyr,orig$nbe_lower_TgCyr,orig$nbe_upper_TgCyr)*1e-3, na.rm=TRUE)
}
zrange[2] = zrange[2] + 0.5

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$nbe_TgCyr[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$nbe_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$nbe_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$nbe_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$nbe_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$nbe_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$nbe_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$nbe_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$nbe_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(paste(orig_name,ssp_scenarios),paste(alt_name,ssp_scenarios)), 
       col = c(colours_orig,colours_alt), lty = rep(1,2*length(ssp_scenarios)), 
       pch=rep(NA,2*length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Net Biome Exchange (PgC y",r^-1,")",sep="")), side=2, padj=-1.60,cex=1.2)
#mtext("Year", side=1, padj=2.0,cex=1.6)

## plot GPP

# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$gpp_TgCyr,orig$gpp_lower_TgCyr,orig$gpp_upper_TgCyr,
                     alt$gpp_TgCyr,alt$gpp_lower_TgCyr,alt$gpp_upper_TgCyr)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$gpp_TgCyr,orig$gpp_lower_TgCyr,orig$gpp_upper_TgCyr)*1e-3, na.rm=TRUE)
}
zrange = zrange * c(0.9,1.0)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$gpp_TgCyr[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$gpp_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$gpp_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$gpp_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$gpp_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$gpp_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$gpp_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$gpp_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$gpp_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Gross Primary Productivity (PgC y",r^-1,")",sep="")), side=2, padj=-1.60, cex=1.2)

## Plot fire

# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$fire_TgCyr,orig$fire_lower_TgCyr,orig$fire_upper_TgCyr,
                     alt$fire_TgCyr,alt$fire_lower_TgCyr,alt$fire_upper_TgCyr)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$fire_TgCyr,orig$fire_lower_TgCyr,orig$fire_upper_TgCyr)*1e-3, na.rm=TRUE)
}
zrange = zrange * c(0.9,1.1)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$fire_TgCyr[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$fire_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$fire_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$fire_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$fire_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$fire_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$fire_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$fire_lower_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$fire_upper_TgCyr[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Fire Emissions (PgC y",r^-1,")",sep="")), side=2, padj=-1.60,cex=1.2)

## Plot Evapotranspiration

# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$et_PgH2Oyr,orig$et_lower_PgH2Oyr,orig$et_upper_PgH2Oyr,
                     alt$et_PgH2Oyr,alt$et_lower_PgH2Oyr,alt$et_upper_PgH2Oyr), na.rm=TRUE) 
} else {
    zrange = range(c(orig$et_PgH2Oyr,orig$et_lower_PgH2Oyr,orig$et_upper_PgH2Oyr), na.rm=TRUE)
}
zrange = zrange * c(0.9,1.1)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$et_PgH2Oyr[ssp,]~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$et_lower_PgH2Oyr[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$et_upper_PgH2Oyr[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$et_PgH2Oyr[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$et_lower_PgH2Oyr[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$et_upper_PgH2Oyr[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$et_PgH2Oyr[ssp,]~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$et_lower_PgH2Oyr[ssp,]~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$et_upper_PgH2Oyr[ssp,]~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Evapotranspiration (PgC y",r^-1,")",sep="")), side=2, padj=-1.60,cex=1.2)
mtext("Year", side=1, padj=2.0,cex=1.6)

dev.off()

## Plot major ecosystem pools over time
fig_height = 4200 ; fig_width = 2500 ; fig_res = 300
if (alt_name != "") {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_and_",alt_name,"_lai_wood_litter_som_timeseries_comparison_plusCI",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res)
} else {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_lai_wood_litter_som_timeseries_comparison_plusCI",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res) 
}
par(mfrow=c(4,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
## Plot leaf area index

# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$lai_m2m2,orig$lai_lower_m2m2,orig$lai_upper_m2m2,
                     alt$lai_m2m2,alt$lai_lower_m2m2,alt$lai_upper_m2m2), na.rm=TRUE) 
} else {
    zrange = range(c(orig$lai_m2m2,orig$lai_lower_m2m2,orig$lai_upper_m2m2), na.rm=TRUE)
}
zrange[2] = zrange[2] + 0.5

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$lai_m2m2[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$lai_lower_m2m2[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$lai_upper_m2m2[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$lai_m2m2[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$lai_lower_m2m2[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$lai_upper_m2m2[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$lai_m2m2[ssp,]~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$lai_lower_m2m2[ssp,]~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$lai_upper_m2m2[ssp,]~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
abline(0,0,col="grey", lwd=2)
#legend("topleft", legend = c(ssp_scenarios), col = colours_orig, lty = rep(1,length(ssp_scenarios)), 
#       pch=rep(NA,length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
legend("topleft", legend = c(paste(orig_name,ssp_scenarios),paste(alt_name,ssp_scenarios)), 
       col = c(colours_orig,colours_alt), lty = rep(1,2*length(ssp_scenarios)), 
       pch=rep(NA,2*length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)

mtext(expression(paste("Leaf area index (",m^2,m^-2,")",sep="")), side=2, padj=-1.60,cex=1.2)
#mtext("Year", side=1, padj=2.0,cex=1.6)

## plot wood carbon stock

# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$wood_TgC,orig$wood_lower_TgC,orig$wood_upper_TgC,
                     alt$wood_TgC,alt$wood_lower_TgC,alt$wood_upper_TgC)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$wood_TgC,orig$wood_lower_TgC,orig$wood_upper_TgC)*1e-3, na.rm=TRUE)
}
zrange = zrange * c(0.9,1.0)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$wood_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$wood_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$wood_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$wood_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$wood_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$wood_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$wood_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$wood_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$wood_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Wood (PgC)",sep="")), side=2, padj=-1.60, cex=1.2)

## Plot fine litter carbon stock (foliage + fine root)

# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$litter_TgC,orig$litter_lower_TgC,orig$litter_upper_TgC,
                     alt$litter_TgC,alt$litter_lower_TgC,alt$litter_upper_TgC)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$litter_TgC,orig$litter_lower_TgC,orig$litter_upper_TgC)*1e-3, na.rm=TRUE)
}
zrange = zrange * c(0.9,1.1)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$litter_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$litter_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$litter_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$litter_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$litter_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$litter_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$litter_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$litter_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$litter_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Litter (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)

## Plot soil carbon stock

# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$soil_TgC,orig$soil_lower_TgC,orig$soil_upper_TgC,
                     alt$soil_TgC,alt$soil_lower_TgC,alt$soil_upper_TgC), na.rm=TRUE)*1e-3
} else {
    zrange = range(c(orig$soil_TgC,orig$soil_lower_TgC,orig$soil_upper_TgC), na.rm=TRUE)*1e-3
}
zrange = zrange * c(0.9,1.1)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$soil_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$soil_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$soil_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
    for (ssp in seq(2, length(ssp_scenarios))) {
         lines(orig$soil_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
         lines(orig$soil_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
         lines(orig$soil_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
    }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$soil_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$soil_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$soil_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Soil (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)
mtext("Year", side=1, padj=2.0,cex=1.6)

dev.off()

if (alt_name != "") {

    ## Plot major ecosystem pools over time
    fig_height = 4200 ; fig_width = 2500 ; fig_res = 300
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_and_",alt_name,"_lai_wood_litter_som_timeseries_CI_change",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res)

    par(mfrow=c(4,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))

    ## Plot the difference in uncertainty propagated overtime
    ## Leaf area index
    ci_differences = (alt$lai_upper_m2m2-alt$lai_lower_m2m2) - (orig$lai_upper_m2m2-orig$lai_lower_m2m2)
    zrange = range(ci_differences, na.rm=TRUE)
    zrange = max(abs(zrange)) * c(-1,1)
    zrange[1] = min(-0.25,zrange[1]) ; zrange[2] = max(0.25,zrange[2])

    # Start with first SSP
    ssp = 1
    # Plot initial curve
    plot(ci_differences[ssp,]~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
         col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
    # Add the other SSPs
    if (length(ssp_scenarios) > 1) {
        for (ssp in seq(2, length(ssp_scenarios))) {
             lines(ci_differences[ssp,]~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
        }
    }
    abline(0,0,col="grey", lwd=2)
    mtext(expression(paste(Delta,"LAI CI (",m^2,m^-2,")",sep="")), side=2, padj=-1.60,cex=1.2)
#    legend("topleft", legend = c(ssp_scenarios), col = colours_orig, lty = rep(1,length(ssp_scenarios)), 
#           pch=rep(NA,length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
    legend("topleft", legend = c(paste(orig_name,ssp_scenarios),paste(alt_name,ssp_scenarios)), 
           col = c(colours_orig,colours_alt), lty = rep(1,2*length(ssp_scenarios)), 
           pch=rep(NA,2*length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
           
    ## Wood stock
    ci_differences = (alt$wood_upper_TgC-alt$wood_lower_TgC) - (orig$wood_upper_TgC-orig$wood_lower_TgC)
    zrange = range(ci_differences, na.rm=TRUE)*1e-3
    zrange = max(abs(zrange)) * c(-1,1)

    # Start with first SSP
    ssp = 1
    # Plot initial curve
    plot(ci_differences[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
         col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
    # Add the other SSPs
    if (length(ssp_scenarios) > 1) {
        for (ssp in seq(2, length(ssp_scenarios))) {
             lines(ci_differences[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
        }
    }
    abline(0,0,col="grey", lwd=2)    
    mtext(expression(paste(Delta,"Wood CI (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)

    ## Litter stock
    ci_differences = (alt$litter_upper_TgC-alt$litter_lower_TgC) - (orig$litter_upper_TgC-orig$litter_lower_TgC)
    zrange = range(ci_differences, na.rm=TRUE)*1e-3
    zrange = max(abs(zrange)) * c(-1,1)

    # Start with first SSP
    ssp = 1
    # Plot initial curve
    plot(ci_differences[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
         col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
    # Add the other SSPs
    if (length(ssp_scenarios) > 1) {
        for (ssp in seq(2, length(ssp_scenarios))) {
             lines(ci_differences[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
        }
    }
    abline(0,0,col="grey", lwd=2)    
    mtext(expression(paste(Delta,"Litter CI (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)

    ## Soil stock
    ci_differences = (alt$soil_upper_TgC-alt$soil_lower_TgC) - (orig$soil_upper_TgC-orig$soil_lower_TgC)
    zrange = range(ci_differences, na.rm=TRUE)*1e-3
    zrange = max(abs(zrange)) * c(-1,1)

    # Start with first SSP
    ssp = 1
    # Plot initial curve
    plot(ci_differences[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
         col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
    # Add the other SSPs
    if (length(ssp_scenarios) > 1) {
        for (ssp in seq(2, length(ssp_scenarios))) {
             lines(ci_differences[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
        }  
    }
    abline(0,0,col="grey", lwd=2)    
    mtext(expression(paste(Delta,"Soil CI (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)
    mtext("Year", side=1, padj=2.0,cex=1.6)

    dev.off()
} # alt_name != ""

## Plot major ecosystem pools over time
fig_height = 4200 ; fig_width = 2500 ; fig_res = 300
if (alt_name != "") {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_and_",alt_name,"_dtotal_total_biomass_dom_timeseries_comparison_plusCI",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res)
} else {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_dtotal_total_biomass_dom_timeseries_comparison_plusCI",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res) 
}
par(mfrow=c(4,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
## Change in total NBE, not we used total C stock, which requires sign chnage for NBE

# Define y-axis range
if (alt_name != "") {
    zrange = range(c(-orig$dCtotal_TgC,-orig$dCtotal_lower_TgC,-orig$dCtotal_upper_TgC,
                     -alt$dCtotal_TgC,-alt$dCtotal_lower_TgC,-alt$dCtotal_upper_TgC)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(-orig$dCtotal_TgC,-orig$dCtotal_lower_TgC,-orig$dCtotal_upper_TgC)*1e-3, na.rm=TRUE)
}
zrange[2] = zrange[2] + 200

# Start with first SSP
ssp = 1
# Plot initial curve
plot(-orig$dCtotal_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(-orig$dCtotal_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(-orig$dCtotal_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
     for (ssp in seq(2, length(ssp_scenarios))) {
          lines(-orig$dCtotal_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
          lines(-orig$dCtotal_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
          lines(-orig$dCtotal_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
     }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(-alt$dCtotal_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(-alt$dCtotal_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(-alt$dCtotal_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
abline(0,0,col="grey", lwd=2)
#legend("topleft", legend = c(ssp_scenarios), col = colours_orig, lty = rep(1,length(ssp_scenarios)), 
#       pch=rep(NA,length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
legend("topleft", legend = c(paste(orig_name,ssp_scenarios),paste(alt_name,ssp_scenarios)), 
       col = c(colours_orig,colours_alt), lty = rep(1,2*length(ssp_scenarios)), 
       pch=rep(NA,2*length(ssp_scenarios)), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)

mtext(expression(paste(Sigma,"NBE (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)
#mtext("Year", side=1, padj=2.0,cex=1.6)

## plot total carbon stock

# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$Ctotal_TgC,orig$Ctotal_lower_TgC,orig$Ctotal_upper_TgC,
                     alt$Ctotal_TgC,alt$Ctotal_lower_TgC,alt$Ctotal_upper_TgC)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$Ctotal_TgC,orig$Ctotal_lower_TgC,orig$Ctotal_upper_TgC)*1e-3, na.rm=TRUE)
}
zrange = zrange * c(0.9,1.0)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$Ctotal_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$Ctotal_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$Ctotal_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
     for (ssp in seq(2, length(ssp_scenarios))) {
          lines(orig$Ctotal_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
          lines(orig$Ctotal_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
          lines(orig$Ctotal_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
     }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$Ctotal_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$Ctotal_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$Ctotal_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Total C (PgC)",sep="")), side=2, padj=-1.60, cex=1.2)

## Plot biomass (foliage + fine root + wood + labile)

# Units conversion, annual time series TgC/yr -> PgC/yr
# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$bio_TgC,orig$bio_lower_TgC,orig$bio_upper_TgC,
                     alt$bio_TgC,alt$bio_lower_TgC,alt$bio_upper_TgC)*1e-3, na.rm=TRUE) 
} else {
    zrange = range(c(orig$bio_TgC,orig$bio_lower_TgC,orig$bio_upper_TgC)*1e-3, na.rm=TRUE)
}
zrange = zrange * c(0.9,1.1)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$bio_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$bio_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$bio_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
     for (ssp in seq(2, length(ssp_scenarios))) {
          lines(orig$bio_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
          lines(orig$bio_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
          lines(orig$bio_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
     }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$bio_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$bio_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$bio_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("Biomass (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)

## Plot dead organic matter carbon stock (litter + soil)

# Define y-axis range
if (alt_name != "") {
    zrange = range(c(orig$dom_TgC,orig$dom_lower_TgC,orig$dom_upper_TgC,
                     alt$dom_TgC,alt$dom_lower_TgC,alt$dom_upper_TgC), na.rm=TRUE)*1e-3
} else {
    zrange = range(c(orig$dom_TgC,orig$dom_lower_TgC,orig$dom_upper_TgC), na.rm=TRUE)*1e-3
}
zrange = zrange * c(0.9,1.1)

# Start with first SSP
ssp = 1
# Plot initial curve
plot(orig$dom_TgC[ssp,]*1e-3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=colours_orig[ssp], type="l", lwd=4, ylab="", xlab="", lty=1)
# Overlay the uncertainty information
lines(orig$dom_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
lines(orig$dom_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
# Add the other SSPs
if (length(ssp_scenarios) > 1) {
     for (ssp in seq(2, length(ssp_scenarios))) {
          lines(orig$dom_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 1) 
          lines(orig$dom_lower_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
          lines(orig$dom_upper_TgC[ssp,]*1e-3~run_years, col=colours_orig[ssp], lwd=3, lty = 2) 
     }
}
# Add second model if available
if (alt_name != "") {
    # Add the other SSPs
    for (ssp in seq(1, length(ssp_scenarios))) {
         lines(alt$dom_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 1) 
         lines(alt$dom_lower_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
         lines(alt$dom_upper_TgC[ssp,]*1e-3~run_years, col=colours_alt[ssp], lwd=3, lty = 2) 
    }
}
mtext(expression(paste("DOM (PgC)",sep="")), side=2, padj=-1.60,cex=1.2)
mtext("Year", side=1, padj=2.0,cex=1.6)

dev.off()

## Plot trend maps for key fluxes

# Figure 3. Median maps of Mean transit times (top row) and domain average time series of wood and SOM stock changes (bottom row)
## Plot major ecosystem fluxes over time
fig_height = 2400 ; fig_width = 4200 ; fig_res = 300
ssp = 1 # select SSP
if (alt_name != "") {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_and_",alt_name,"_",ssp_scenarios[ssp],"_Wood_SOM_change",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res)
    
    par(mfrow=c(2,2),mai=c(0.3,0.35,0.3,0.7),omi=c(0.2,0.2,0.3,0.01))
    # Climate variables
    # Assign variables, note unit conversion from gCm2yr-> MgC/ha/yr
    var1 = orig$wood_trend[ssp,,]*1e-2
    var2 = orig$som_trend[ssp,,]*1e-2
    var3 = alt$wood_trend[ssp,,]*1e-2
    var4 = alt$som_trend[ssp,,]*1e-2

    # Convert to raster
    var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
    var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    # Trim to data area
    var1 = trim(var1)
    var2 = trim(var2)
    var3 = trim(var3)  
    var4 = trim(var4)
    # legend position
    ee = ext(var1) ; e = rep(NA, 4)
    e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
    e[3] = ee[3] ; e[4] = ee[4]
    # Determine ranges
    zrange1 = c(-1,1)*max(abs(range(c(values(var1),values(var3)),na.rm=TRUE)))
    zrange2 = c(-1,1)*max(abs(range(c(values(var2),values(var4)),na.rm=TRUE)))
    zrange3 = zrange1
    zrange4 = zrange2

    # Mean annual median estimates
    plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_sign))
    mtext(expression(paste("Wood change (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.9, padj = 0.9)     
    mtext(gsub("_"," ",orig_name), side=2, cex = 1.9, padj = -0.9)     
    plot(crop(landmask,var1), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_sign))
    mtext(expression(paste("SOM change (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.9, padj = 0.9)     
    plot(crop(landmask,var2), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_sign))
    #mtext(expression(paste("Wood change (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.9, padj = 0.9)     
    mtext(gsub("_"," ",alt_name), side=2, cex = 1.9, padj = -0.9)     
    plot(crop(landmask,var3), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_sign))
    #mtext(expression(paste("SOM change (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.9, padj = 0.9)     
    plot(crop(landmask,var4), add=TRUE, lwd=0.5, border="grey")

    dev.off()
}

## Plot major ecosystem fluxes over time
fig_height = 1000 ; fig_width = 5000 ; fig_res = 300
ssp = 1 # select SSP
if (alt_name != "") {
    png(file = paste(outdir,"/",gsub("%","_",PROJECT$name),"_",orig_name,"_and_",alt_name,"_",ssp_scenarios[ssp],"_Wood_LAI_NBE_GPP_Rauto_Rhet_fire_harvest",output_suffix,".png",sep=""), 
        height=fig_height, width=fig_width, res=fig_res) 
    
    par(mfrow=c(2,5),mai=c(0.15,0.20,0.15,0.55),omi=c(0.2,0.2,0.3,0.01))
    # Assign variables (gCm2day -> MgC/ha/yr)
    var1 = apply(orig$nbe_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var2 = apply(orig$gpp_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var3 = apply(orig$rauto_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var4 = apply(orig$rhet_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var5 = apply(orig$fire_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var6 = apply(alt$nbe_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var7 = apply(alt$gpp_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var8 = apply(alt$rauto_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var9 = apply(alt$rhet_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2
    var10 = apply(alt$fire_gCm2day_grid[ssp,,,],c(1,2),mean,na.rm=TRUE) * 365.25 * 1e-2

    # Convert to raster
    var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
    var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var7 = rast(vals = t((var7)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var8 = rast(vals = t((var8)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var9 = rast(vals = t((var9)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var10 = rast(vals = t((var10)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    # Trim to data area
    var1 = trim(var1)
    var2 = trim(var2)
    var3 = trim(var3)  
    var4 = trim(var4)
    var5 = trim(var5)
    var6 = trim(var6)
    var7 = trim(var7)
    var8 = trim(var8)
    var9 = trim(var9)
    var10 = trim(var10)
    # legend position
    ee = ext(var1) ; e = rep(NA, 4)
    e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
    e[3] = ee[3] ; e[4] = ee[4]
    # Determine ranges
    zrange1 = c(-1,1) * max(abs(range(c(values(var1),values(var6)),na.rm=TRUE)))
    zrange2 = c(0,max(c(values(var2),values(var3),values(var4),values(var7),values(var8),values(var9)),na.rm=TRUE))
    zrange3 = zrange2
    zrange4 = zrange2
    zrange5 = c(0,max(c(values(var5),values(var10)),na.rm=TRUE))
    zrange6 = zrange1
    zrange7 = zrange2
    zrange8 = zrange2
    zrange9 = zrange2
    zrange10 = zrange5

    ### continue updating from here
    ### also in previous figure add a third row showin the differences between analyses
    ### and finally what has happened to the Ctotal_gCm2 plots, where are the values there?
    # Mean annual median estimates
    plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=rev(colour_choices_sign))
    mtext(expression(paste("NBE (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)     
    mtext(gsub("_"," ",orig_name), side=2, cex = 1.4, padj = -0.9)     
    plot(crop(landmask,var1), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_gain))
    mtext(expression(paste("GPP (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)     
    plot(crop(landmask,var2), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_loss))
    mtext(expression(paste("Ra (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)          
    plot(crop(landmask,var3), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_loss))
    mtext(expression(paste("Rh (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)          
    plot(crop(landmask,var4), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_loss))
    mtext(expression(paste("Fire emission (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)          
    plot(crop(landmask,var5), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=rev(colour_choices_sign))
    #mtext(expression(paste("NBE (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)     
    mtext(gsub("_"," ",alt_name), side=2, cex = 1.4, padj = -0.9)     
    plot(crop(landmask,var6), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var7, range=zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_gain))
    #mtext(expression(paste("GPP (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)     
    plot(crop(landmask,var7), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var8, range=zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_loss))
    #mtext(expression(paste("Ra (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)          
    plot(crop(landmask,var8), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var9, range=zrange9, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_loss))
    #mtext(expression(paste("Rh (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)          
    plot(crop(landmask,var9), add=TRUE, lwd=0.5, border="grey")
    # Mean annual median estimates
    plot(var10, range=zrange10, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
         cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
         main = "", col=(colour_choices_loss))
    #mtext(expression(paste("Fire emission (MgC h",a^-1,y^-1,")",sep="")), side=3, cex = 1.4, padj = -0.35)          
    plot(crop(landmask,var10), add=TRUE, lwd=0.5, border="grey")

    dev.off()

}

## Summary values for the figures
# Create a local target region mask that can be used to cancel out from the data constrain
for (ssp in seq(1, length(ssp_scenarios))) {
     print(paste("...Summary information for ",ssp_scenarios[ssp],"...",sep=""))
     print(paste("......",gsub("_"," ",orig_name),sep=""))
     local_mask = land_fraction 
     local_mask = local_mask * ((area*local_mask) / sum(area*local_mask, na.rm=TRUE))
     # Mean and range of the histogram consistency for assimilated NBE and wood likelihoods
#     tmp1 = round(sum(orig$wood_assim_data_overlap_fraction[ssp,,]*local_mask,na.rm=TRUE),digits=3)
#     tmp2 = round(quantile(orig$wood_assim_data_overlap_fraction[ssp,,],prob=c(0.025),na.rm=TRUE),digits=3)
#     tmp3 = round(quantile(orig$wood_assim_data_overlap_fraction[ssp,,],prob=c(0.975),na.rm=TRUE),digits=3)
#     print(paste("Assimilated wood: Mean (2.5 % / 97.5 %) histogram consistency = ",tmp1," (",tmp2,"/",tmp3,")", sep=""))
#     tmp1 = round(sum(orig$nbe_assim_data_overlap_fraction[ssp,,]*local_mask,na.rm=TRUE),digits=3)
#     tmp2 = round(quantile(orig$nbe_assim_data_overlap_fraction[ssp,,],prob=c(0.025),na.rm=TRUE),digits=3)
#     tmp3 = round(quantile(orig$nbe_assim_data_overlap_fraction[ssp,,],prob=c(0.975),na.rm=TRUE),digits=3)
#     print(paste("Assimilated NBE: Mean (2.5 % / 97.5 %) histogram consistency = ",tmp1," (",tmp2,"/",tmp3,")", sep=""))
#     tmp1 = round(sum(orig$lai_assim_data_overlap_fraction[ssp,,]*local_mask,na.rm=TRUE),digits=3)
#     tmp2 = round(quantile(orig$lai_assim_data_overlap_fraction[ssp,,],prob=c(0.025),na.rm=TRUE),digits=3)
#     tmp3 = round(quantile(orig$lai_assim_data_overlap_fraction[ssp,,],prob=c(0.975),na.rm=TRUE),digits=3)
#     print(paste("Assimilated LAI: Mean (2.5 % / 97.5 %) histogram consistency = ",tmp1," (",tmp2,"/",tmp3,")", sep=""))
     # Mean residence times for wood, litter, som
     print("Median estimates:")
     print(paste("...Wood MTT   = ",round(sum(as.vector(orig$MTT_wood_years[ssp,,]*local_mask),na.rm=TRUE),digits = 3)," years",sep=""))
     print(paste("...Litter MTT = ",round(sum(as.vector(orig$MTT_litter_years[ssp,,]*local_mask),na.rm=TRUE),digits = 3)," years",sep=""))
     print(paste("...SOM MTT    = ",round(sum(as.vector(orig$MTT_som_years[ssp,,]*local_mask),na.rm=TRUE),digits = 3)," years",sep=""))
     # Mean accumulation rates for wood, litter, som
     print(paste("...Total accumulation rate  = ",round((orig$dCtotal_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     print(paste("...Wood accumulation rate   = ",round((orig$dCwood_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     print(paste("...Litter accumulation rate = ",round((orig$dClitter_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     print(paste("...SOM accumulation rate    = ",round((orig$dCsom_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     # Fraction accumulated in wood vs soil
     tmp1 = orig$dCtotal_TgC[ssp,nos_years]
     tmp2 = orig$dCwood_TgC[ssp,nos_years]
     tmp3 = orig$dCsom_TgC[ssp,nos_years]
     tmp = round(tmp2 / tmp1,digits = 3)
     print(paste("Ratio (0-1) of C accumulation in wood:total = ",tmp,sep=""))
     tmp = round(tmp3 / tmp1,digits = 3)
     print(paste("Ratio (0-1) of C accumulation in  som:total = ",tmp,sep=""))
     print(paste("......",gsub("_"," ",alt_name),sep=""))
     # Mean and range of the histogram consistency for assimilated NBE and wood likelihoods
#     tmp1 = round(sum(alt$wood_assim_data_overlap_fraction[ssp,,]*local_mask,na.rm=TRUE),digits=3)
#     tmp2 = round(quantile(alt$wood_assim_data_overlap_fraction[ssp,,],prob=c(0.025),na.rm=TRUE),digits=3)
#     tmp3 = round(quantile(alt$wood_assim_data_overlap_fraction[ssp,,],prob=c(0.975),na.rm=TRUE),digits=3)
#     print(paste("Assimilated wood: Mean (2.5 % / 97.5 %) histogram consistency = ",tmp1," (",tmp2,"/",tmp3,")", sep=""))
#     tmp1 = round(sum(alt$nbe_assim_data_overlap_fraction[ssp,,]*local_mask,na.rm=TRUE),digits=3)
#     tmp2 = round(quantile(alt$nbe_assim_data_overlap_fraction[ssp,,],prob=c(0.025),na.rm=TRUE),digits=3)
#     tmp3 = round(quantile(alt$nbe_assim_data_overlap_fraction[ssp,,],prob=c(0.975),na.rm=TRUE),digits=3)
#     print(paste("Assimilated NBE: Mean (2.5 % / 97.5 %) histogram consistency = ",tmp1," (",tmp2,"/",tmp3,")", sep=""))
#     tmp1 = round(sum(alt$lai_assim_data_overlap_fraction[ssp,,]*local_mask,na.rm=TRUE),digits=3)
#     tmp2 = round(quantile(alt$lai_assim_data_overlap_fraction[ssp,,],prob=c(0.025),na.rm=TRUE),digits=3)
#     tmp3 = round(quantile(alt$lai_assim_data_overlap_fraction[ssp,,],prob=c(0.975),na.rm=TRUE),digits=3)
#     print(paste("Assimilated LAI: Mean (2.5 % / 97.5 %) histogram consistency = ",tmp1," (",tmp2,"/",tmp3,")", sep=""))
     # Mean residence times for wood, litter, som
     print("Median estimates:")
     print(paste("...Wood MTT   = ",round(sum(as.vector(alt$MTT_wood_years[ssp,,]*local_mask),na.rm=TRUE),digits = 3)," years",sep=""))
     print(paste("...Litter MTT = ",round(sum(as.vector(alt$MTT_litter_years[ssp,,]*local_mask),na.rm=TRUE),digits = 3)," years",sep=""))
     print(paste("...SOM MTT    = ",round(sum(as.vector(alt$MTT_som_years[ssp,,]*local_mask),na.rm=TRUE),digits = 3)," years",sep=""))
     # Mean accumulation rates for wood, litter, som
     print(paste("...Total accumulation rate  = ",round((alt$dCtotal_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     print(paste("...Wood accumulation rate   = ",round((alt$dCwood_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     print(paste("...Litter accumulation rate = ",round((alt$dClitter_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     print(paste("...SOM accumulation rate    = ",round((alt$dCsom_TgC[ssp,nos_years]*1e-3)/nos_years,digits = 3)," PgC/y",sep=""))
     # Fraction accumulated in wood vs soil
     tmp1 = alt$dCtotal_TgC[ssp,nos_years]
     tmp2 = alt$dCwood_TgC[ssp,nos_years]
     tmp3 = alt$dCsom_TgC[ssp,nos_years]
     tmp = round(tmp2 / tmp1,digits = 3)
     print(paste("Ratio (0-1) of C accumulation in wood:total = ",tmp,sep=""))
     tmp = round(tmp3 / tmp1,digits = 3)
     print(paste("Ratio (0-1) of C accumulation in  som:total = ",tmp,sep=""))

     print("CI reduction (alt-orig)")
     tmp = sum((orig$MTT_wood_upper_years[ssp,,] - orig$MTT_wood_lower_years[ssp,,]) * local_mask, na.rm=TRUE)
     tmp = ((sum((alt$MTT_wood_upper_years[ssp,,] - alt$MTT_wood_lower_years[ssp,,]) * local_mask, na.rm=TRUE) - tmp) / tmp) * 100
     print(paste("...Wood MTT   = ",round(tmp,digits = 3)," %",sep=""))
     tmp = sum((orig$MTT_litter_upper_years[ssp,,] - orig$MTT_litter_lower_years[ssp,,]) * local_mask, na.rm=TRUE)
     tmp = ((sum((alt$MTT_litter_upper_years[ssp,,] - alt$MTT_litter_lower_years[ssp,,]) * local_mask, na.rm=TRUE) - tmp) / tmp) * 100
     print(paste("...Litter MTT = ",round(tmp,digits = 3)," %",sep=""))
     tmp = sum((orig$MTT_som_upper_years[ssp,,] - orig$MTT_som_lower_years[ssp,,]) * local_mask, na.rm=TRUE)
     tmp = ((sum((alt$MTT_som_upper_years[ssp,,] - alt$MTT_som_lower_years[ssp,,]) * local_mask, na.rm=TRUE) - tmp) / tmp) * 100
     print(paste("...SOM MTT    = ",round(tmp,digits = 3)," %",sep=""))
     tmp = orig$dCtotal_upper_TgC[ssp,nos_years] - orig$dCtotal_lower_TgC[ssp,nos_years]
     tmp = (((alt$dCtotal_upper_TgC[ssp,nos_years] - alt$dCtotal_lower_TgC[ssp,nos_years]) - tmp) / tmp) * 100
     print(paste("...Total accumulation rate  = ",round(tmp,digits = 3)," %",sep=""))
     tmp = orig$dCwood_upper_TgC[ssp,nos_years] - orig$dCwood_lower_TgC[ssp,nos_years]
     tmp = (((alt$dCwood_upper_TgC[ssp,nos_years] - alt$dCwood_lower_TgC[ssp,nos_years]) - tmp) / tmp) * 100
     print(paste("...Wood accumulation rate   = ",round(tmp,digits = 3)," %",sep=""))
     tmp = orig$dClitter_upper_TgC[ssp,nos_years] - orig$dClitter_lower_TgC[ssp,nos_years]
     tmp = (((alt$dClitter_upper_TgC[ssp,nos_years] - alt$dClitter_lower_TgC[ssp,nos_years]) - tmp) / tmp) * 100
     print(paste("...Litter accumulation rate = ",round(tmp,digits = 3)," %",sep=""))
     tmp = orig$dCsom_upper_TgC[ssp,nos_years] - orig$dCsom_lower_TgC[ssp,nos_years]
     tmp = (((alt$dCsom_upper_TgC[ssp,nos_years] - alt$dCsom_lower_TgC[ssp,nos_years]) - tmp) / tmp) * 100
     print(paste("...SOM accumulation rate    = ",round(tmp,digits = 3)," %",sep=""))
}

ssp = 1
par(mfrow=c(2,3))
# Original
image.plot(orig$gpp_trend[ssp,,], main="GPP trend") 
image.plot(orig$lai_trend[ssp,,], main="LAI trend")
image.plot(orig$et_trend[ssp,,], main="ET trend") 
# Alternate
image.plot(alt$gpp_trend[ssp,,]) 
image.plot(alt$lai_trend[ssp,,]) 
image.plot(alt$et_trend[ssp,,])

par(mfrow=c(2,3))
image.plot(orig$mint_trend[ssp,,], main="Min T trend") 
image.plot(orig$maxt_trend[ssp,,], main="Max T trend") 
image.plot(orig$swrad_trend[ssp,,], main="SW trend") 
image.plot(orig$co2_trend[ssp,,], main="CO2 trend") 
image.plot(orig$precip_trend[ssp,,], main="Precipitation trend") 
image.plot(orig$vpd_trend[ssp,,], main="VPD trend")

}

