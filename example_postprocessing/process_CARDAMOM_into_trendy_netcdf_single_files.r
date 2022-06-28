
###
## Process CARDAMOM-DALEC output files into NetCDF files 
## consistent with the TRENDYv11 / GCP model intercomparison structure
## In constrast to the sibling script which groups variables together into different files,
## this script places everything into a single file per document consistent with the latest guidance.
### 

###
## Job specific information

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_CsomPriorNCSCD/"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_nbe_CsomPriorNCSDC3m/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Trendyv9_historical/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire_nbe/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Mexico_1deg_C7_agb_lca_gpp_fire_nbe/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Miombo_0.5deg_allWood/"

# Specify any extra information for the filename
output_prefix = "" # follow with "_"
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

###
## Begin defining variables
###

###
## Variables we assume always exist

## At model time step variables
# DRIVERS
AIRT_MIN = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
AIRT_MAX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SWRAD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
CO2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
DOY = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
PRECIP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FLOSS_FRAC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days))) 
BURNT_FRAC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WINDSPD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
VPD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
# STATE OBSERVATIONS (gC/m2, except LAI = m2/m2)
LAI_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
LAI_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WOOD_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WOOD_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SOIL_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SOIL_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
# FLUX OBSERVATIONS (gC/m2/day)
GPP_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
GPP_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NEE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NEE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NBE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NBE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FIRE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FIRE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))

###
## Variables based on their presence

# C STATES
if (exists(x = "lai_m2m2", where = grid_output)) {LAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Ctotal_gCm2", where = grid_output)) {TOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "labile_gCm2", where = grid_output)) {LAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "foliage_gCm2", where = grid_output)) {FOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "roots_gCm2", where = grid_output)) {ROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "wood_gCm2", where = grid_output)) {WOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "litter_gCm2", where = grid_output)) {LIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "som_gCm2", where = grid_output)) {SOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "woodlitter_gCm2", where = grid_output)) {WLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "dom_gCm2", where = grid_output)) {DOM = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "biomass_gCm2", where = grid_output)) {BIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# C FLUXES
if (exists(x = "gpp_gCm2day", where = grid_output)) {GPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "rauto_gCm2day", where = grid_output)) {RAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "rhet_gCm2day", where = grid_output)) {RHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "npp_gCm2day", where = grid_output)) {NPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "fire_gCm2day", where = grid_output)) {FIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "harvest_gCm2day", where = grid_output)) {HARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "reco_gCm2day", where = grid_output)) {RECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "nee_gCm2day", where = grid_output)) {NEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "nbe_gCm2day", where = grid_output)) {NBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "nbp_gCm2day", where = grid_output)) {NBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# H2O FLUXES
if (exists(x = "ET_kgH2Om2day", where = grid_output)) {ET = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# H2O STATES
if (exists(x = "SurfWater_kgH2Om2", where = grid_output)) {SurfWater = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "snow_kgH2Om2", where = grid_output)) {SNOW = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# BIOPHYSICAL
if (exists(x = "CiCa", where = grid_output)) {CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "wSWP_MPa", where = grid_output)) {wSWP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "APAR_MJm2day", where = grid_output)) {APAR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "gs_demand_supply_ratio", where = grid_output)) {gs_demand_supply_ratio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "gs_mmolH2Om2day", where = grid_output)) {gs = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "gb_mmolH2Om2day", where = grid_output)) {gb = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Direct allocation of NPP (labile, foliar, fine root, wood, gC/m2/day)
if (exists(x = "alloc_labile_gCm2day", where = grid_output)) {NPP_labile_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_foliage_gCm2day", where = grid_output)) {NPP_foliage_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_wood_gCm2day", where = grid_output)) {NPP_wood_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_roots_gCm2day", where = grid_output)) {NPP_root_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Combined (direct / indirect) allocation of C to foliage (gC/m2/day)
if (exists(x = "combined_alloc_foliage_gCm2day", where = grid_output)) {NPP_combinedfoliage_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Allocation fraction of expressed NPP (foliar, fine root, wood)
if (exists(x = "NPP_foliage_fraction", where = grid_output)) {fNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "NPP_wood_fraction", where = grid_output)) {wNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "NPP_roots_fraction", where = grid_output)) {rNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# MRT (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; years)
if (exists(x = "MTT_annual_labile_years", where = grid_output)) {labMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_foliage_years", where = grid_output)) {folMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_wood_years", where = grid_output)) {wooMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_roots_years", where = grid_output)) {rooMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_litter_years", where = grid_output)) {litMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_woodlitter_years", where = grid_output)) {wlitMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_som_years", where = grid_output)) {somMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_biomass_years", where = grid_output)) {bioMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_dom_years", where = grid_output)) {domMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}

# Total outflux from these pools (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "outflux_labile_gCm2day", where = grid_output)) {outflux_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_foliage_gCm2day", where = grid_output)) {outflux_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_wood_gCm2day", where = grid_output)) {outflux_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_roots_gCm2day", where = grid_output)) {outflux_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_litter_gCm2day", where = grid_output)) {outflux_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_woodlitter_gCm2day", where = grid_output)) {outflux_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_som_gCm2day", where = grid_output)) {outflux_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_biomass_gCm2day", where = grid_output)) {outflux_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_dom_gCm2day", where = grid_output)) {outflux_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Fire mortality outflux (labile, foliar, wood, fine root, litter, biomass, dead organic matter; gC/m2/day)
if (exists(x = "FIRElitter_labile_gCm2day", where = grid_output)) {FIRElitter_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_foliage_gCm2day", where = grid_output)) {FIRElitter_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_wood_gCm2day", where = grid_output)) {FIRElitter_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_roots_gCm2day", where = grid_output)) {FIRElitter_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_litter_gCm2day", where = grid_output)) {FIRElitter_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_woodlitter_gCm2day", where = grid_output)) {FIRElitter_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_biomass_gCm2day", where = grid_output)) {FIRElitter_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Fire combustion outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "FIREemiss_labile_gCm2day", where = grid_output)) {FIREemiss_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_foliage_gCm2day", where = grid_output)) {FIREemiss_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_wood_gCm2day", where = grid_output)) {FIREemiss_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_roots_gCm2day", where = grid_output)) {FIREemiss_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_litter_gCm2day", where = grid_output)) {FIREemiss_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_woodlitter_gCm2day", where = grid_output)) {FIREemiss_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_som_gCm2day", where = grid_output)) {FIREemiss_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_biomass_gCm2day", where = grid_output)) {FIREemiss_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_dom_gCm2day", where = grid_output)) {FIREemiss_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Harvest mortality outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "HARVESTlitter_labile_gCm2day", where = grid_output)) {HARVESTlitter_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_foliage_gCm2day", where = grid_output)) {HARVESTlitter_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_wood_gCm2day", where = grid_output)) {HARVESTlitter_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_roots_gCm2day", where = grid_output)) {HARVESTlitter_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_biomass_gCm2day", where = grid_output)) {HARVESTlitter_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Harvest extracted outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "HARVESTextracted_labile_gCm2day", where = grid_output)) {HARVESTextracted_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_foliage_gCm2day", where = grid_output)) {HARVESTextracted_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_wood_gCm2day", where = grid_output)) {HARVESTextracted_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_roots_gCm2day", where = grid_output)) {HARVESTextracted_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_litter_gCm2day", where = grid_output)) {HARVESTextracted_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_woodlitter_gCm2day", where = grid_output)) {HARVESTextracted_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_som_gCm2day", where = grid_output)) {HARVESTextracted_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_biomass_gCm2day", where = grid_output)) {HARVESTextracted_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_dom_gCm2day", where = grid_output)) {HARVESTextracted_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}

# Fill the output arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]

         # Read in site specific drivers
         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
    
         # DRIVERS
         AIRT_MIN[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,2]+273.15     # mint C -> K
         AIRT_MAX[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,3]+273.15     # maxt C -> K
         SWRAD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,4]*1e6*(1/86400) # SWRAD MJ/m2/day -> W/m2
         CO2[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,5]                 # CO2 ppm
         DOY[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,6]                 # Julian day of year
         PRECIP[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,7]              # Precipitation kgH2O/m2/s
         FLOSS_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,8]          # Forest loss fraction
         BURNT_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,9]          # Burned fraction
         WINDSPD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,15]            # Wind speed m/s
         VPD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,16]                # Vapour pressure deficit Pa

         # STATE OBSERVATIONS
         LAI_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,3]     # LAI m2/m2
         LAI_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,4] # LAI UNC m2/m2
         # ...first assign priors to the 1st time step to simplify our storage
         WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriors[21] # Wood gC/m2
         WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriorunc[21] # Wood UNC gC/m2
         SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriors[23] # SOM gC/m2
         SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriorunc[23] # SOM UNC gC/m2
         # ...second assign time series information if any exists
         WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),13]
         WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),14]
         SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),19]
         SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),20]
         # FLUX OBSERVATIONS (gC/m2/day)
         GPP_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,1] 
         GPP_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,2]
         NEE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,5] 
         NEE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,6]
         NBE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,35]
         NBE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,36]
         FIRE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,7]
         FIRE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,8]

         ## At model time step
         # STATES (NOTE: unit conversions gC/m2 -> kgC/m2, except LAI)
         if (exists("LAI")) {LAI[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lai_m2m2[n,,]}
         if (exists("TOT")) {TOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Ctotal_gCm2[n,,]*1e-3}
         if (exists("LAB")) {LAB[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$labile_gCm2[n,,]*1e-3}
         if (exists("FOL")) {FOL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$foliage_gCm2[n,,]*1e-3}
         if (exists("ROOT")) {ROOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$roots_gCm2[n,,]*1e-3}
         if (exists("WOOD")) {WOOD[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wood_gCm2[n,,]*1e-3}
         if (exists("LIT")) {LIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$litter_gCm2[n,,]*1e-3}
         if (exists("SOIL")) {SOIL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$som_gCm2[n,,]*1e-3}
         if (exists("WLIT")) {WLIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$woodlitter_gCm2[n,,]*1e-3}
         if (exists("DOM")) {DOM[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dom_gCm2[n,,]*1e-3}
         if (exists("BIO")) {BIO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$biomass_gCm2[n,,]*1e-3}
         # FLUXES (NOTE; unit conversion gC/m2/day -> kgC/m2/s)
         if (exists("GPP")) {GPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gpp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("RAU")) {RAU[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rauto_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("RHE")) {RHE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rhet_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP")) {NPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$npp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("FIR")) {FIR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$fire_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("HARV")) {HARV[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$harvest_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("RECO")) {RECO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$reco_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NEE")) {NEE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nee_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NBE")) {NBE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbe_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NBP")) {NBP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbp_gCm2day[n,,]* 1e-3 * (1/86400)}
         # H2O FLUXES
         if (exists(x = "ET_kgH2Om2day", where = grid_output)) {ET[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$ET_kgH2Om2day[n,,] * (1/86400)}
         # H2O STATES
         if (exists(x = "SurfWater_kgH2Om2", where = grid_output)) {SurfWater[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$SurfWater_kgH2Om2[n,,]}
         if (exists(x = "snow_kgH2Om2", where = grid_output)) {SNOW[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$snow_kgH2Om2[n,,]}
         # BIOPHYSICAL
         if (exists("CiCa")) {CiCa[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$CiCa[n,,]}
         if (exists("wSWP")) {wSWP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wSWP_MPa[n,,]}
         if (exists("APAR")) {APAR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$APAR_MJm2day[n,,]}
         if (exists("gs_demand_supply_ratio")) {gs_demand_supply_ratio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gs_demand_supply_ratio[n,,]}
         if (exists("gs")) {gs[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gs_mmolH2Om2day[n,,]}
         if (exists("gb")) {gb[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gb_mmolH2Om2day[n,,]}
         # NPP (foliar, root, wood; gC/m2/day -> kgC/m2/s)
         if (exists("NPP_labile_FLX")) {NPP_labile_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_labile_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_root_FLX")) {NPP_root_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_roots_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_wood_FLX")) {NPP_wood_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_wood_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_combinedfoliage_FLX")) {NPP_combinedfoliage_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_alloc_foliage_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_foliage_FLX")) {NPP_foliage_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_foliage_gCm2day[n,,]* 1e-3 * (1/86400)}

         # NPP (fraction) and MRT years are requested to have same number of time steps as stocks and fluxes
         # This is awkward as no easy way to repeat specific elements without loop for variables which have no meaningful value at sub-annual timescales
         # (and are therefore calculated as annuals)
         for (q in seq(1, nos_quantiles)) {
              # MRT
              if (exists("labMRT")) {labMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_labile_years[n,q,], each = steps_per_year)}
              if (exists("folMRT")) {folMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_foliage_years[n,q,], each = steps_per_year)} 
              if (exists("rooMRT")) {rooMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_roots_years[n,q,], each = steps_per_year)} 
              if (exists("wooMRT")) {wooMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_wood_years[n,q,], each = steps_per_year)}
              if (exists("litMRT")) {litMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_litter_years[n,q,], each = steps_per_year)}
              if (exists("wlitMRT")) {wlitMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_woodlitter_years[n,q,], each = steps_per_year)}
              if (exists("somMRT")) {somMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_som_years[n,q,], each = steps_per_year)}
              if (exists("bioMRT")) {bioMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_biomass_years[n,q,], each = steps_per_year)}
              if (exists("domMRT")) {domMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_dom_years[n,q,], each = steps_per_year)}
              # NPP fractional allocation
              if (exists("fNPP")) {fNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$NPP_foliage_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))}
              if (exists("rNPP")) {rNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$NPP_roots_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))}
              if (exists("wNPP")) {wNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$NPP_wood_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))}
         } # loop quantiles
              
         # Total outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("outflux_labile")) {outflux_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_foliage")) {outflux_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_wood")) {outflux_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_root")) {outflux_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_litter")) {outflux_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_woodlitter")) {outflux_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_som")) {outflux_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_som_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_bio")) {outflux_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_dom")) {outflux_dom[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_dom_gCm2day[n,,] * 1e-3 * (1/86400)}

         # Fire mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIRElitter_labile")) {FIRElitter_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_foliage")) {FIRElitter_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_wood")) {FIRElitter_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_root")) {FIRElitter_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_litter")) {FIRElitter_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_woodlitter")) {FIRElitter_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_bio")) {FIRElitter_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}

         # Fire combustion outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIREemiss_labile")) {FIREemiss_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_foliage")) {FIREemiss_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_wood")) {FIREemiss_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_root")) {FIREemiss_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_litter")) {FIREemiss_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_woodlitter")) {FIREemiss_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_som")) {FIREemiss_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_som_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_bio")) {FIREemiss_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_dom")) {FIREemiss_dom[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_dom_gCm2day[n,,] * 1e-3 * (1/86400)}
         
         # HARVEST mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIRElitter_labile")) {HARVESTlitter_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_foliage")) {HARVESTlitter_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_wood")) {HARVESTlitter_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_root")) {HARVESTlitter_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_bio")) {HARVESTlitter_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}

         # HARVEST extraction outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("HARVESTextracted_labile")) {HARVESTextracted_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_foliage")) {HARVESTextracted_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_wood")) {HARVESTextracted_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_root")) {HARVESTextracted_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_litter")) {HARVESTextracted_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_woodlitter")) {HARVESTextracted_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_som")) {HARVESTextracted_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_som_gCm2day[n,,] * 1e-3 * (1/86400)} 
         if (exists("HARVESTextracted_bio")) {HARVESTextracted_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_dom")) {HARVESTextracted_dom[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_dom_gCm2day[n,,] * 1e-3 * (1/86400)}

     } # Does the file exist / has it been processed

} # site loop

###
## Define dimensions which will be used across all files

## define dimension
lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
year_dimen <- ncdim_def( "year", units="", 1:nos_years)
npar_dimen <- ncdim_def( "nos_parameters", units="", 1:(max(PROJECT$model$nopars)+1)) # NOTE: +1 is to account for the log-likelihood 

###
## Begin writing out individual files per variable
###

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), 
                 dim=list(time_dimen), missval = -99999, prec="double", compression = 9)

# LAI
if(exists("LAI")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"lai_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("lai_ensemble", unit="m2.m-2", longname = "Leaf Area Index - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  LAI)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}


# Labile
if(exists("LAB")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cLabile_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("cLabile_ensemble", unit="kg.m-2", longname = "Carbon in labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  LAB)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Foliar
if(exists("FOL")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cLeaf_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cLeaf_ensemble", unit="kg.m-2", longname = "Carbon in leaves - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FOL)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fine root
if(exists("ROOT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cFineRoot_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cFineRoot_ensemble", unit="kg.m-2", longname = "Carbon in fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  ROOT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Wood
if(exists("WOOD")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cWoodTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cWoodTotal_ensemble", unit="kg.m-2", longname = "Carbon in (AGB + BGB) wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  WOOD)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Foliar + fine root litter
if(exists("LIT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cLeafFineRootLitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cLeafFineRootLitter_ensemble", unit="kg.m-2", longname = "Carbon in (Foliar + fine root) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  LIT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Wood litter
if(exists("WLIT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cWoodLitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cWoodLitter_ensemble", unit="kg.m-2", longname = "Carbon in (wood) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  WLIT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Soil organic matter
if(exists("SOIL")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cSOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cSOM_ensemble", unit="kg.m-2", longname = "Carbon in soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  SOIL)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Dead organic matter
if(exists("DOM")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cDOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cDOM_ensemble", unit="kg.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  DOM)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Biomass
if(exists("BIO")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cVeg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cVeg_ensemble", unit="kg.m-2", longname = "Carbon in live biomass - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  BIO)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# TotalC
if(exists("TOT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"cTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("cTotal_ensemble", unit="kg.m-2", longname = "Carbon in live and dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  TOT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# GPP
if(exists("GPP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"gpp_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("gpp_ensemble", unit="kg.m-2.s-1", longname = "Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  GPP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Autotrophic respiration
if(exists("RAU")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"ra_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("ra_ensemble", unit="kg.m-2.s-1", longname = "Autotrophic (Plant) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  RAU)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Heterotrophic respiration
if(exists("RHE")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"rh_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("rh_ensemble", unit="kg.m-2.s-1", longname = "Heterotrophic Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  RHE)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Ecosystem respiration
if(exists("RECO")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"reco_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("reco_ensemble", unit="kg.m-2.s-1", longname = "Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  RECO)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Net Primary Productivity
if(exists("NPP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"npp_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("npp_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NPP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Net Ecosystem Exchange
if(exists("NEE")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"nee_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("nee_ensemble", unit="kg.m-2.s-1", longname = "Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NEE)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Net Biome Exchange
if(exists("NBE")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"nbe_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("nbe_ensemble", unit="kg.m-2.s-1", longname = "Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NBE)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Net Biome Productivity
if(exists("NBP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"nbp_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("nbp_ensemble", unit="kg.m-2.s-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NBP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire emissions
if(exists("FIR")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"fFire_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("fFire_ensemble", unit="kg.m-2.s-1", longname = "Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIR)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Flux from forest loss
if(exists("HARV")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"fLuc_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("fLuc_ensemble", unit="kg.m-2.s-1", longname = "Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARV)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from labile
if(exists("outflux_labile")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cLabile_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_labile)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from foliage
if(exists("outflux_foliage")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cLeaf_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_foliage)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from root
if(exists("outflux_root")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cFineRoot_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_root)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from wood
if(exists("outflux_wood")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cWoodTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_wood)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from litter
if(exists("outflux_litter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cLeafFineRootlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_litter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from wood litter
if(exists("outflux_woodlitter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cWoodlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new, outflux_woodlitter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Total outflux from som
if(exists("outflux_som")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cSOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_som)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Biomass
if(exists("outflux_bio")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cVeg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_bio)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Dead Organic Matter
if(exists("outflux_dom")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"outflux_cDOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("outflux_cDOM_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  outflux_dom)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from labile
if(exists("FIRElitter_labile")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIRElitter_cLabile_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_labile)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from foliage
if(exists("FIRElitter_foliage")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIRElitter_cLeaf_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_foliage)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from root
if(exists("FIRElitter_root")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIRElitter_cFineRoot_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_root)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from wood
if(exists("FIRElitter_wood")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIRElitter_cWoodTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_wood)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from litter
if(exists("FIRElitter_litter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIRElitter_cLeafFineRootlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_litter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from litter
if(exists("FIRElitter_woodlitter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIRElitter_cWoodlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_woodlitter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from som
if(exists("FIRElitter_bio")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cVeg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIRElitter_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from Vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIRElitter_bio)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from labile
if(exists("FIREemiss_labile")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cLabile_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new, FIREemiss_labile)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from foliage
if(exists("FIREemiss_foliage")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cLeaf_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_foliage)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from root
if(exists("FIREemiss_root")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cFineRoot_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_root)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from wood
if(exists("FIREemiss_wood")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cWoodTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_wood)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from litter
if(exists("FIREemiss_litterLAI")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cLeafFineRootlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_litter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from wood litter
if(exists("FIREemiss_woodlitter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cWoodlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_woodlitter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from som
if(exists("FIREemiss_som")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cSOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_som)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from vegetation
if(exists("FIREemiss_bio")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cVeg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_bio)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fire mortality outflux from dom
if(exists("FIREemiss_dom")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"FIREemiss_cDOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("FIREemiss_cDOM_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  FIREemiss_dom)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest litter outflux from labile
if(exists("HARVESTlitter_labile")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTlitter_cLabile_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTlitter_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new, HARVESTlitter_labile)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest litter outflux from foliage
if(exists("HARVESTlitter_foliage")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTlitter_cLeaf_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTlitter_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new, HARVESTlitter_foliage)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest litter outflux from root
if(exists("HARVESTlitter_root")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTlitter_cFineRoot_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTlitter_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new, HARVESTlitter_root)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest litter outflux from wood
if(exists("HARVESTlitter_wood")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTlitter_cWoodTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTlitter_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTlitter_wood)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest litter outflux from vegetation
if(exists("HARVESTlitter_bio")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTlitter_cVeg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTlitter_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from Vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTlitter_bio)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from labile
if(exists("HARVESTextracted_labile")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cLabile_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_labile)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from foliage
if(exists("HARVESTextracted_foliage")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cLeaf_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_foliage)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from root
if(exists("HARVESTextracted_root")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cFineRoot_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_root)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from wood
if(exists("HARVESTextracted_wood")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cWoodTotal_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_wood)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from litter
if(exists("HARVESTextracted_litter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cLeafFineRootlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_litter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from wood litter
if(exists("HARVESTextracted_woodlitter")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cWoodlitter_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_woodlitter)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from som
if(exists("HARVESTextracted_som")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cSOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_som)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from vegetation
if(exists("HARVESTextracted_bio")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cVeg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_bio)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Harvest extracted outflux from dom
if(exists("HARVESTextracted_dom")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"HARVESTextracted_cDOM_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("HARVESTextracted_cDOM_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  HARVESTextracted_dom)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

###
## Mean residence times and NPP allocation / fluxes
###
              
## Mean Residence Times 
# Labile
if(exists("labMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_lab_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_lab_ensemble", unit="year", longname = "Mean Labile Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  labMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Foliar
if(exists("folMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_fol_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_fol_ensemble", unit="year", longname = "Mean Foliar Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  folMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fine root
if(exists("rooMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_root_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_root_ensemble", unit="year", longname = "Mean fine root Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  rooMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Wood
if(exists("wooMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_wood_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_wood_ensemble", unit="year", longname = "Mean wood Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  wooMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fine litter (fol + fine root)
if(exists("litMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_lit_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_lit_ensemble", unit="year", longname = "Mean foliar + fine root litter Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  litMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Wood litter 
if(exists("wlitMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_wlit_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_wlit_ensemble", unit="year", longname = "Mean wood litter Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  wlitMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Soil
if(exists("somMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_som_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_som_ensemble", unit="year", longname = "Mean Soil Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  somMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Dead Organic Matter
if(exists("domMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_dom_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_dom_ensemble", unit="year", longname = "Mean Soil + litter Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  domMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Biomass / vegetation
if(exists("bioMRT")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"MTT_veg_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("MTT_veg_ensemble", unit="year", longname = "Mean biomass / vegetation Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  bioMRT)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

## NPP allocation fractions
# Foliar
if(exists("fNPP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_fol_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("NPP_fol_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  fNPP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fine root
if(exists("rNPP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_root_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("NPP_root_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  rNPP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Wood
if(exists("wNPP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_wood_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("NPP_wood_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  wNPP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}


## NPP allocation fluxes
# Labile
if(exists("NPP_labile_FLX")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_labile_flx_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new= ncvar_def("NPP_labile_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NPP_labile_FLX)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Direct allocation to Foliar
if(exists("NPP_foliage_FLX")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_direct_fol_flx_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("NPP_direct_fol_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity direct to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NPP_foliage_FLX)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Combined direct and via labile allocation of NPP to foliage
if(exists("NPP_combinedfoliage_FLX")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_fol_flx_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("NPP_fol_flx_ensemble", unit="kg.m-2.s-1", longname = "Both direct and via labile Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NPP_combinedfoliage_FLX)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Fine root
if(exists("NPP_root_FLX")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_root_flx_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new= ncvar_def("NPP_root_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NPP_root_FLX)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Wood
if(exists("NPP_wood_FLX")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_wood_flx_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("NPP_wood_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  NPP_wood_FLX)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

###
## H2O FLUXES
###

# Evapotranspiration
if(exists("ET")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"evapotrans_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("evapotrans_ensemble", unit="kg.m-2.s-1", longname = "Evapotranspiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  ET)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

###
## H2O STATES
###

# Soil surface water content
if(exists("SurfWater")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"SurfWater_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("SurfWater_ensemble", unit="kg.m-2", longname = "Soil water content (0-30 cm) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  SurfWater)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Soil surface snow cover
if(exists("SNOW")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"SNOW_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("SNOW_ensemble", unit="kg.m-2", longname = "Soil surface snow cover - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  SNOW)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

###
## Biophysical diagnostics
###

# CiCa
if(exists("CiCa")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"CiCa_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new = ncvar_def("CiCa_ensemble", unit="1", longname = "Internal:Ambiant CO2 ratio - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  CiCa)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Soil water potential weighted by root access
if(exists("wSWP")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"wSWP_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("wSWP_ensemble", unit="MPa", longname = "Soil water potential weighted by root access - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  wSWP)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
   }

# Canopy absorbed photosynthetically active radiation
if(exists("APAR")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"APAR_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("APAR_ensemble", unit="MJ.m2.d-1", longname = "Canopy absorbed photosynthatically active radiation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  APAR)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Ratio of actual stomatal conductance to that allowed by water supply to canopy
if(exists("gs_demand_supply_ratio")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"gs_demand_supply_ratio_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("gs_demand_supply_ratio_ensemble", unit="1", longname = "Ratio of actual stomatal conductance to potential allowed by water supply to canopy - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  gs_demand_supply_ratio)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Stomatal conductance for water
if(exists("gs")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"gs_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("gs_ensemble", unit="mmol.m-2.d-1", longname = "Stomatal conductance for water - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  gs)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}

# Boundary layer conductance for water
if(exists("gb")) {
   # Define the output file name
   output_name = paste(PROJECT$results_processedpath,output_prefix,"gb_ensemble",output_suffix,".nc",sep="")
   # Delete if the file currently exists
   if (file.exists(output_name)) {file.remove(output_name)}
   # Define the new variable
   var_new  = ncvar_def("gb_ensemble", unit="mmol.m-2.d-1", longname = "Boundary conductance for water - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   # Create the empty file space
   new_file=nc_create(filename=output_name, vars=list(var0,var_new), force_v4 = TRUE)
   # Load first variable into the file
   # TIMING
   ncvar_put(new_file, var0, drivers$met[,1])
   # VARIABLE
   ncvar_put(new_file, var_new,  gb)
   # Close the existing file to ensure its written to file
   nc_close(new_file)
}


