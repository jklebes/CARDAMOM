
###
## Process CARDAMOM-DALEC output files into NetCDF files 
## consistent with the TRENDYv11 / GCP model intercomparison structure.
## NOTE: unlike the sibling script this script does not put each variable into a single file
## and is thus not consistent with the latest guidance.
### 

##TO DO:
# 1) Add total soil water variable and output
# 2) Determine which internal dynamics variables we want to keep
# 3) Update the single file version with the corresponding variables
# 4) Update the generic to netcdf version which uses our own file / variable names
 
###
## Job specific information

print("Begin creation of Trendy v11 compatible grouped netcdf files...")

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/reccap2_permafrost_1deg_dalec2_isimip3a_agb_lca_nbe_CsomPriorNCSDC3m"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/Miombo_0.5deg_allWood"

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
if (length(which(quantiles_wanted == 0.025)) == 1) {
    low_quant = which(quantiles_wanted == 0.025)
} else {
    stop("Desired low quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.5)) == 1) {
    mid_quant = which(quantiles_wanted == 0.5)
} else {
    stop("Median quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.975)) == 1) {
    high_quant = which(quantiles_wanted == 0.975)
} else {
    stop("Desired high quantile cannot be found")
}

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
# C STATE CHANGE ESTIMATES
if (exists(x = "dCbiomass_gCm2", where = grid_output)) {dBIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
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
# C FLUXES - ANNUAL
if (exists(x = "mean_annual_gpp_gCm2day", where = grid_output)) {AGPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_rauto_gCm2day", where = grid_output)) {ARAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_rhet_gCm2day", where = grid_output)) {ARHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_npp_gCm2day", where = grid_output)) {ANPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_fire_gCm2day", where = grid_output)) {AFIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_harvest_gCm2day", where = grid_output)) {AHARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_reco_gCm2day", where = grid_output)) {ARECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_nee_gCm2day", where = grid_output)) {ANEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_nbe_gCm2day", where = grid_output)) {ANBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_nbp_gCm2day", where = grid_output)) {ANBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
# H2O FLUXES
if (exists(x = "ET_kgH2Om2day", where = grid_output)) {ET = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Etrans_kgH2Om2day", where = grid_output)) {Etrans = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Esoil_kgH2Om2day", where = grid_output)) {Esoil = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Ewetcanopy_kgH2Om2day", where = grid_output)) {Ewetcanopy = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "total_drainage_kgH2Om2day", where = grid_output)) {total_drainage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "runoff_kgH2Om2day", where = grid_output)) {runoff = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "underflow_kgH2Om2day", where = grid_output)) {underflow = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Combined natural, fire and harvest driven litter creation
if (exists(x = "combined_biomass_to_litter_gCm2day", where = grid_output)) {Combined_bio_litter_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "combined_labile_to_litter_gCm2day", where = grid_output)) {Combined_labile_litter_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "combined_foliage_to_litter_gCm2day", where = grid_output)) {Combined_foliage_litter_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "combined_roots_to_litter_gCm2day", where = grid_output)) {Combined_roots_litter_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "combined_wood_to_litter_gCm2day", where = grid_output)) {Combined_wood_litter_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Combined natural, fire and harvest driven litter to som creation
if (exists(x = "combined_litter_to_som_gCm2day", where = grid_output)) {Combined_litter_som_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "combined_woodlitter_to_som_gCm2day", where = grid_output)) {Combined_woodlitter_som_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Direct allocation of NPP (labile, foliar, fine root, wood, gC/m2/day)
if (exists(x = "alloc_wood_gCm2day", where = grid_output)) {NPP_wood_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_roots_gCm2day", where = grid_output)) {NPP_root_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Combined (direct / indirect) allocation of C to foliage (gC/m2/day)
if (exists(x = "combined_alloc_foliage_gCm2day", where = grid_output)) {NPP_combinedfoliage_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Fire combustion outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "FIREemiss_litter_gCm2day", where = grid_output)) {FIREemiss_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_woodlitter_gCm2day", where = grid_output)) {FIREemiss_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_som_gCm2day", where = grid_output)) {FIREemiss_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_biomass_gCm2day", where = grid_output)) {FIREemiss_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}

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
         # States at time step
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
         # Change in stocks
         if (exists("dBIO")) {dBIO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dCbiomass_gCm2[n,,]*1e-3}
 

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
         # FLUXES - ANNUAL (NOTE; unit conversion gC/m2/day -> kgC/m2/s)
         if (exists("AGPP")) {AGPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_gpp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ARAU")) {ARAU[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_rauto_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ARHE")) {ARHE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_rhet_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANPP")) {ANPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_npp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("AFIR")) {AFIR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_fire_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("AHARV")) {AHARV[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_harvest_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ARECO")) {ARECO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_reco_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANEE")) {ANEE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_nee_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANBE")) {ANBE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_nbe_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANBP")) {ANBP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_nbp_gCm2day[n,,]* 1e-3 * (1/86400)}
         # H2O FLUXES
         if (exists("ET")) {ET[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$ET_kgH2Om2day[n,,] * (1/86400)}
         if (exists("Etrans")) {Etrans[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Etrans_kgH2Om2day[n,,] * (1/86400)}
         if (exists("Esoil")) {Esoil[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Esoil_kgH2Om2day[n,,] * (1/86400)}
         if (exists("Ewetcanopy")) {Ewetcanopy[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Ewetcanopy_kgH2Om2day[n,,] * (1/86400)}
         if (exists("total_drainage")) {total_drainage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$total_drainage_kgH2Om2day[n,,] * (1/86400)}         
         if (exists("runoff")) {runoff[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$runoff_kgH2Om2day[n,,] * (1/86400)}
         if (exists("underflow")) {underflow[grid_output$i_location[n],grid_output$j_location[n],,]= grid_output$underflow_kgH2Om2day[n,,] * (1/86400)}
         # Combined natural, fire and harvest driven litter creation
         if (exists(x = "Combined_bio_litter_FLX")) {Combined_bio_litter_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_biomass_to_litter_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists(x = "Combined_labile_litter_FLX")) {Combined_labile_litter_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_labile_to_litter_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists(x = "Combined_foliage_litter_FLX")) {Combined_foliage_litter_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_foliage_to_litter_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists(x = "Combined_roots_litter_FLX")) {Combined_roots_litter_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_roots_to_litter_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists(x = "Combined_wood_litter_FLX")) {Combined_wood_litter_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_wood_to_litter_gCm2day[n,,]* 1e-3 * (1/86400)}
         # Combined natural, fire and harvest driven litter to som creation
         if (exists(x = "Combined_litter_som_FLX")) {Combined_litter_som_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_litter_to_som_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists(x = "Combined_woodlitter_som_FLX")) {Combined_woodlitter_som_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_woodlitter_to_som_gCm2day[n,,]* 1e-3 * (1/86400)}
         # NPP (foliar, root, wood; gC/m2/day -> kgC/m2/s)
         if (exists("NPP_root_FLX")) {NPP_root_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_roots_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_wood_FLX")) {NPP_wood_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_wood_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_combinedfoliage_FLX")) {NPP_combinedfoliage_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_alloc_foliage_gCm2day[n,,]* 1e-3 * (1/86400)}
               
         # Fire combustion outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIREemiss_litter")) {FIREemiss_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_woodlitter")) {FIREemiss_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_som")) {FIREemiss_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_som_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_bio")) {FIREemiss_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}

     } # Does the file exist / has it been processed

} # site loop

###
## Define dimensions which will be used across all files

## define dimension
lat_dimen <- ncdim_def( "latitude", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "longitude", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
year_dimen <- ncdim_def( "year", units="", as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
npar_dimen <- ncdim_def( "nos_parameters", units="", 1:(max(PROJECT$model$nopars)+1)) # NOTE: +1 is to account for the log-likelihood 

###
## Create C STATES file with timing information, landsea fraction and grid area
###

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), 
                 dim=list(time_dimen), missval = -99999, prec="single", compression = 9)
var1 = ncvar_def("grid_area", units = "m2", longname = paste("Pixel area",sep=""), 
                 dim=list(long_dimen,lat_dimen), missval = -99999, prec="single", compression = 9)
var2 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                 dim=list(long_dimen,lat_dimen), missval = -99999, prec="single", compression = 9)
# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area 
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###                 
## ADD STATE VARIABLES
###

# LAI
if(exists("LAI")) {
   # Median
   var_new  = ncvar_def("lai", unit="m2.m-2", longname = "Leaf Area Index - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LAI[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("lai_2.5pc", unit="m2.m-2", longname = "Leaf Area Index - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LAI[,,low_quant,])
   # High CI
   var_new  = ncvar_def("lai_97.5pc", unit="m2.m-2", longname = "Leaf Area Index - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LAI[,,high_quant,])
}

# Labile
if(exists("LAB")) {
   var_new  = ncvar_def("cLabile", unit="kg.m-2", longname = "Carbon in labile - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new,  LAB[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cLabile_2.5pc", unit="kg.m-2", longname = "Carbon in labile - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LAB[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cLabile_97.5pc", unit="kg.m-2", longname = "Carbon in labile - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LAB[,,high_quant,])
}

# Foliar
if(exists("FOL")) {
   var_new = ncvar_def("cLeaf", unit="kg.m-2", longname = "Carbon in leaves - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new,  FOL[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cLeaf_2.5pc", unit="kg.m-2", longname = "Carbon in leaves - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FOL[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cLeaf_97.5pc", unit="kg.m-2", longname = "Carbon in leaves - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FOL[,,high_quant,])
}

# Fine root
if(exists("ROOT")) {
   var_new = ncvar_def("cFineRoot", unit="kg.m-2", longname = "Carbon in fine root - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  ROOT[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cFineRoot_2.5pc", unit="kg.m-2", longname = "Carbon in fine root - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ROOT[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cFineRoot_97.5pc", unit="kg.m-2", longname = "Carbon in fine root - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ROOT[,,high_quant,])
}

# Wood
if(exists("WOOD")) {
   var_new = ncvar_def("cWoodTotal", unit="kg.m-2", longname = "Carbon in (AGB + BGB) wood - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  WOOD[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cWoodTotal_2.5pc", unit="kg.m-2", longname = "Carbon in (AGB + BGB) wood - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  WOOD[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cWoodTotal_97.5pc", unit="kg.m-2", longname = "Carbon in (AGB + BGB) wood - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  WOOD[,,high_quant,])
}

# Foliar + fine root litter
if(exists("LIT")) {
   var_new = ncvar_def("cLitter", unit="kg.m-2", longname = "Carbon in (Foliar + fine root) litter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new,  LIT[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cLitter_2.5pc", unit="kg.m-2", longname = "Carbon in (Foliar + fine root) litter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LIT[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cLitter_97.5pc", unit="kg.m-2", longname = "Carbon in (Foliar + fine root) litter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LIT[,,high_quant,])
}

# Wood litter
if(exists("WLIT")) {
   var_new = ncvar_def("cCwd", unit="kg.m-2", longname = "Carbon in (wood) litter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new,  WLIT[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cCwd_2.5pc", unit="kg.m-2", longname = "Carbon in (wood) litter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  WLIT[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cCwd_97.5pc", unit="kg.m-2", longname = "Carbon in (wood) litter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  WLIT[,,high_quant,])
}

# Soil organic matter
if(exists("SOIL")) {
   var_new = ncvar_def("cSOM", unit="kg.m-2", longname = "Carbon in soil organic matter (0-1m) - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new,  SOIL[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cSOM_2.5pc", unit="kg.m-2", longname = "Carbon in soil organic matter (0-1m) - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  SOIL[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cSOM_97.5pc", unit="kg.m-2", longname = "Carbon in soil organic matter (0-1m) - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  SOIL[,,high_quant,])
}

# Dead organic matter
if(exists("DOM")) {
   var_new = ncvar_def("cDOM", unit="kg.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new,  DOM[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cDOM_2.5pc", unit="kg.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  DOM[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cDOM_97.5pc", unit="kg.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  DOM[,,high_quant,])
}

# Biomass
if(exists("BIO")) {
   var_new = ncvar_def("cVeg", unit="kg.m-2", longname = "Carbon in live biomass - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, BIO[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cVeg_2.5pc", unit="kg.m-2", longname = "Carbon in live biomass - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  BIO[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cVeg_97.5pc", unit="kg.m-2", longname = "Carbon in live biomass - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  BIO[,,high_quant,])
}

# TotalC
if(exists("TOT")) {
   var_new = ncvar_def("cTotal", unit="kg.m-2", longname = "Carbon in live and dead organic matter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	
   ncvar_put(new_file, var_new, TOT[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("cTotal_2.5pc", unit="kg.m-2", longname = "Carbon in live and dead organic matter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  TOT[,,low_quant,])
   # High CI
   var_new  = ncvar_def("cTotal_97.5pc", unit="kg.m-2", longname = "Carbon in live and dead organic matter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  TOT[,,high_quant,])
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for C FLUXES with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area 
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## FLUXES
###

## Start with the per-timestep values, annual averages to come later

# GPP
if(exists("GPP")) {
   var_new  = ncvar_def("gpp", unit="kg.m-2.s-1", longname = "Gross Primary Productivity - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, GPP[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("gpp_2.5pc", unit="kg.m-2.s-1", longname = "Gross Primary Productivity - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  GPP[,,low_quant,])
   # High CI
   var_new  = ncvar_def("gpp_97.5pc", unit="kg.m-2.s-1", longname = "Gross Primary Productivity - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  GPP[,,high_quant,])
}

# Autotrophic respiration
if(exists("RAU")) {
   var_new  = ncvar_def("ra", unit="kg.m-2.s-1", longname = "Autotrophic (Plant) Respiration - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, RAU[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("ra_2.5pc", unit="kg.m-2.s-1", longname = "Autotrophic (Plant) Respiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  RAU[,,low_quant,])
   # High CI
   var_new  = ncvar_def("ra_97.5pc", unit="kg.m-2.s-1", longname = "Autotrophic (Plant) Respiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  RAU[,,high_quant,])
}

# Heterotrophic respiration
if(exists("RHE")) {
   var_new  = ncvar_def("rh", unit="kg.m-2.s-1", longname = "Heterotrophic Respiration - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, RHE[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("rh_2.5pc", unit="kg.m-2.s-1", longname = "Heterotrophic Respiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  RHE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("rh_97.5pc", unit="kg.m-2.s-1", longname = "Heterotrophic Respiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  RHE[,,high_quant,])
}

# Ecosystem respiration
if(exists("RECO")) {
   var_new  = ncvar_def("reco", unit="kg.m-2.s-1", longname = "Ecosystem (Ra + Rh) Respiration - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, RECO[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("reco_2.5pc", unit="kg.m-2.s-1", longname = "Ecosystem (Ra + Rh) Respiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  RECO[,,low_quant,])
   # High CI
   var_new  = ncvar_def("reco_97.5pc", unit="kg.m-2.s-1", longname = "Ecosystem (Ra + Rh) Respiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  RECO[,,high_quant,])
}

# Net Primary Productivity
if(exists("NPP")) {
   var_new  = ncvar_def("npp", unit="kg.m-2.s-1", longname = "Net Primary Productivity - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("npp_2.5pc", unit="kg.m-2.s-1", longname = "Net Primary Productivity - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP[,,low_quant,])
   # High CI
   var_new  = ncvar_def("npp_97.5pc", unit="kg.m-2.s-1", longname = "Net Primary Productivity - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP[,,high_quant,])
}

# Net Ecosystem Exchange
if(exists("NEE")) {
   var_new  = ncvar_def("nee", unit="kg.m-2.s-1", longname = "Net Ecosystem Exchange - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NEE[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("nee_2.5pc", unit="kg.m-2.s-1", longname = "Net Ecosystem Exchange - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NEE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("nee_97.5pc", unit="kg.m-2.s-1", longname = "Net Ecosystem Exchange - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NEE[,,high_quant,])
}

# Net Biome Exchange
if(exists("NBE")) {
   var_new  = ncvar_def("nbe", unit="kg.m-2.s-1", longname = "Net Biome Exchange (NEE + Fire) - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NBE[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("nbe_2.5pc", unit="kg.m-2.s-1", longname = "Net Biome Exchange (NEE + Fire) - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NBE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("nbe_97.5pc", unit="kg.m-2.s-1", longname = "Net Biome Exchange (NEE + Fire) - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NBE[,,high_quant,])
}

# Net Biome Productivity
if(exists("NBP")) {
   var_new  = ncvar_def("nbp", unit="kg.m-2.s-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NBP[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("nbp_2.5pc", unit="kg.m-2.s-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NBE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("nbp_97.5pc", unit="kg.m-2.s-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NBE[,,high_quant,])
}

# Fire emissions
if(exists("FIR")) {
   var_new  = ncvar_def("fFire", unit="kg.m-2.s-1", longname = "Fire C emission - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIR[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fFire_2.5pc", unit="kg.m-2.s-1", longname = "Fire C emission - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIR[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fFire_97.5pc", unit="kg.m-2.s-1", longname = "Fire C emission - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIR[,,high_quant,])
}

# Extracted C due to harvest
if(exists("HARV")) {
   var_new = ncvar_def("fLuc", unit="kg.m-2.s-1", longname = "C extracted due to forest harvest - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARV[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fLuc_2.5pc", unit="kg.m-2.s-1", longname = "C extracted due to forest harvest - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  HARV[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fLuc_97.5pc", unit="kg.m-2.s-1", longname = "C extracted due to forest harvest - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  HARV[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from biomass to litter
if(exists("Combined_bio_litter_FLX")) {
   var_new = ncvar_def("fVegLitter", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from biomass - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_bio_litter_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fVegLitter_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from biomass - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_bio_litter_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fVegLitter_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from biomass - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_bio_litter_FLX[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from foliage to litter
if(exists("Combined_labile_litter_FLX")) {
   var_new = ncvar_def("fLabileLitter", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from labile - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_labile_litter_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fLabileLitter_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from labile - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_labile_litter_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fLabileLitter_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from labile - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_labile_litter_FLX[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from labile to litter
if(exists("Combined_foliage_litter_FLX")) {
   var_new = ncvar_def("fLeafLitter", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from foliage - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_foliage_litter_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fLeafLitter_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from foliage - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_foliage_litter_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fLeafLitter_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from foliage - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_foliage_litter_FLX[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from roots to litter
if(exists("Combined_roots_litter_FLX")) {
   var_new = ncvar_def("fRootLitter", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from fine root - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_roots_litter_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fRootLitter_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from fine root - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_roots_litter_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fRootLitter_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from fine root - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_roots_litter_FLX[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from wood to som
# NOTE due to available variable names in TRENDYv11, this variable is used to provide fVegSoil
if(exists("Combined_wood_litter_FLX")) {
   var_new = ncvar_def("fVegSoil", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from wood, which is allocated to som - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_wood_litter_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fVegSoil_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from wood, which is allocated to som - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_wood_litter_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fVegSoil_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven litter creation from wood, which is allocaed to som - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_wood_litter_FLX[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from litter to som
if(exists("Combined_litter_som_FLX")) {
   var_new = ncvar_def("fLitterSoil", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven allocation of litter to soil - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_litter_som_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fLitterSoil_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest allocation of litter to soil - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_litter_som_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fLitterSoil_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven allocation of litter to soil - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_litter_som_FLX[,,high_quant,])
}

# Combined natural, fire and harvest driven flux from wood litter to som
if(exists("Combined_woodlitter_som_FLX")) {
   var_new = ncvar_def("fCwdSoil", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven allocation of wood litter to soil - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Combined_woodlitter_som_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fCwdSoil_2.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest allocation of wood litter to soil - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_woodlitter_som_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fCwdSoil_97.5pc", unit="kg.m-2.s-1", longname = "Combined natural, fire and harvest driven allocation of wood litter to soil - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Combined_woodlitter_som_FLX[,,high_quant,])
}
         
# Fire mortality outflux from litter
if(exists("FIREemiss_litter")) {
   var_new  = ncvar_def("fFireLitter", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from foliar and fine root litter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_litter[,,mid_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # low CI
   var_new  = ncvar_def("fFireLitter_2.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from foliar and fine root litter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_litter[,,low_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # High CI
   var_new  = ncvar_def("fFireLitter_97.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from foliar and fine root litter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_litter[,,high_quant,]*(44/12)) # NOTE: unit change from C -> CO2
}

# Fire mortality outflux from wood litter
if(exists("FIREemiss_woodlitter")) {
   var_new  = ncvar_def("fFireCcwd", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from wood litter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_woodlitter[,,mid_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # low CI
   var_new  = ncvar_def("fFireCcwd_2.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from wood litter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_woodlitter[,,low_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # High CI
   var_new  = ncvar_def("fFireCcwd_97.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from wood litter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_woodlitter[,,high_quant,]*(44/12)) # NOTE: unit change from C -> CO2
}

# Fire mortality outflux from som
if(exists("FIREemiss_som")) {
   var_new  = ncvar_def("fFireCsoil", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from soil organic matter - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_som[,,mid_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # low CI
   var_new  = ncvar_def("fFireCsoil_2.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from soil organic matter - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_som[,,low_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # High CI
   var_new  = ncvar_def("fFireCsoil_97.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from soil organic matter - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_som[,,high_quant,]*(44/12)) # NOTE: unit change from C -> CO2
}

# Fire mortality outflux from vegetation
if(exists("FIREemiss_bio")) {
   var_new  = ncvar_def("fFireCveg", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from vegetation - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_bio[,,mid_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # low CI
   var_new  = ncvar_def("fFireCveg_2.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from vegetation - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_bio[,,low_quant,]*(44/12)) # NOTE: unit change from C -> CO2
   # High CI
   var_new  = ncvar_def("fFireCveg_97.5pc", unit="kg.m-2.s-1", longname = "Fire combusted CO2 output flux from vegetation - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  FIREemiss_bio[,,high_quant,]*(44/12)) # NOTE: unit change from C -> CO2
}

## Annual gross fluxes 

# GPP - mean annual
if(exists("AGPP")) {
   var_new  = ncvar_def("annual_gpp", unit="kg.m-2.s-1", longname = "Mean Annual Gross Primary Productivity - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, AGPP[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_gpp_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Gross Primary Productivity - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  AGPP[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_gpp_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Gross Primary Productivity - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  AGPP[,,high_quant,])
}

# Autotrophic respiration - mean annual
if(exists("ARAU")) {
   var_new  = ncvar_def("annual_ra", unit="kg.m-2.s-1", longname = "Mean Annual Autotrophic (Plant) Respiration - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ARAU[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_ra_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Autotrophic (Plant) Respiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ARAU[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_ra_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Autotrophic (Plant) Respiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ARAU[,,high_quant,])
}

# Heterotrophic respiration - mean annual
if(exists("ARHE")) {
   var_new  = ncvar_def("annual_rh", unit="kg.m-2.s-1", longname = "Mean Annual Heterotrophic Respiration - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ARHE[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_rh_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Heterotrophic Respiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ARHE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_rh_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Heterotrophic Respiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ARHE[,,high_quant,])
}

# Ecosystem respiration - mean annual
if(exists("ARECO")) {
   var_new  = ncvar_def("annual_reco", unit="kg.m-2.s-1", longname = "Mean Annual Ecosystem (Ra + Rh) Respiration - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ARECO[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_reco_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Ecosystem (Ra + Rh) Respiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ARECO[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_reco_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Ecosystem (Ra + Rh) Respiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ARECO[,,high_quant,])
}

# Net Primary Productivity - mean annual
if(exists("ANPP")) {
   var_new  = ncvar_def("annual_npp", unit="kg.m-2.s-1", longname = "Mean Annual Net Primary Productivity - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANPP[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_npp_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Primary Productivity - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANPP[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_npp_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Primary Productivity - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANPP[,,high_quant,])
}

# Net Ecosystem Exchange - mean annual
if(exists("ANEE")) {
   var_new  = ncvar_def("annual_nee", unit="kg.m-2.s-1", longname = "Mean Annual Net Ecosystem Exchange - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANEE[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_nee_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Ecosystem Exchange - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANEE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_nee_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Ecosystem Exchange - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANEE[,,high_quant,])
}

# Net Biome Exchange - mean annual
if(exists("ANBE")) {
   var_new  = ncvar_def("annual_nbe", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Exchange (NEE + Fire) - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANBE[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_nbe_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Exchange (NEE + Fire) - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANBE[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_nbe_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Exchange (NEE + Fire) - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANBE[,,high_quant,])
}

# Net Biome Exchange - mean annual
if(exists("ANBP")) {
   var_new  = ncvar_def("annual_nbp", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Productivity (-NEE - Fire - fLuc) - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANBP[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_nbp_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Productivity (-NEE - Fire - fLuc) - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANBP[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_nbp_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Productivity (-NEE - Fire - fLuc) - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ANBP[,,high_quant,])
}

# Fire emissions - mean annual
if(exists("AFIR")) {
   var_new  = ncvar_def("annual_fFire", unit="kg.m-2.s-1", longname = "Mean Annual Fire C emission - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, AFIR[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_fFire_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Fire C emission - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  AFIR[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_fFire_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual Fire C emission - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  AFIR[,,high_quant,])
}

# Flux from forest loss - mean annual
if(exists("AHARV")) {
   var_new  = ncvar_def("annual_fLuc", unit="kg.m-2.s-1", longname = "Mean Annual C extracted due to forest harvest - Median estimate", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, AHARV[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("annual_fLuc_2.5pc", unit="kg.m-2.s-1", longname = "Mean Annual C extracted due to forest harvest - 2.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  AHARV[,,low_quant,])
   # High CI
   var_new  = ncvar_def("annual_fLuc_97.5pc", unit="kg.m-2.s-1", longname = "Mean Annual C extracted due to forest harvest - 97.5% quantile", dim=list(long_dimen,lat_dimen,year_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  AHARV[,,high_quant,])
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for NPP with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area 
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## NPP allocation / fluxes
###
              
## NPP allocation fluxes

# Combined direct and via labile allocation of NPP to foliage
if(exists("NPP_combinedfoliage_FLX")) {
   var_new = ncvar_def("fAllocLeaf", unit="kg.m-2.s-1", longname = "Both direct and via labile Net Primary Productivity to foliage - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_combinedfoliage_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fAllocLeaf_2.5pc", unit="kg.m-2.s-1", longname = "Both direct and via labile Net Primary Productivity to foliage - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP_combinedfoliage_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fAllocLeaf_97.5pc", unit="kg.m-2.s-1", longname = "Both direct and via labile Net Primary Productivity to foliage - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP_combinedfoliage_FLX[,,high_quant,])
}

# Fine root
if(exists("NPP_root_FLX")) {
   var_new= ncvar_def("fAllocRoot", unit="kg.m-2.s-1", longname = "Net Primary Productivity to fine root - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_root_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fAllocRoot_2.5pc", unit="kg.m-2.s-1", longname = "Net Primary Productivity to fine root - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP_root_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fAllocRoot_97.5pc", unit="kg.m-2.s-1", longname = "Net Primary Productivity to fine root - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP_root_FLX[,,high_quant,])
}

# Wood
if(exists("NPP_wood_FLX")) {
   var_new = ncvar_def("fAllocWood", unit="kg.m-2.s-1", longname = "Net Primary Productivity to wood - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_wood_FLX[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("fAllocWood_2.5pc", unit="kg.m-2.s-1", longname = "Net Primary Productivity to wood - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP_wood_FLX[,,low_quant,])
   # High CI
   var_new  = ncvar_def("fAllocWood_97.5pc", unit="kg.m-2.s-1", longname = "Net Primary Productivity to wood - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  NPP_wood_FLX[,,high_quant,])
}
   
###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for parameters and clusters
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"PARS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the parameter variable definition
var_pars = ncvar_def("parameters", unit="-", longname = "Parameter array + log-likelihood, see specific model version for desccription and units - Ensemble", dim=list(long_dimen,lat_dimen,npar_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
# Create the empty file space, using the already defined timing variable
new_file=nc_create(filename=output_name, vars=list(var_pars), force_v4 = TRUE)
# Load first variable into the file
# PARAMETERS
ncvar_put(new_file, var_pars, grid_output$parameters)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## Other parameter or cluster related information
###
              
## Cluster by affinity propogation process
# parameters only, i.e. initial conditions not included
if(exists(x = "pars_clusters", where = grid_output)) {
   var_new = ncvar_def("pars_clusters", unit="-", longname = "Affinity propogation clustering based on process parameters (i.e. not including initial conditions)", dim=list(long_dimen,lat_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, grid_output$pars_clusters)
}

## Cluster by affinity propogation process
# parameters only, i.e. initial conditions not included
if(exists(x = "clusters", where = grid_output)) {
   var_new = ncvar_def("clusters", unit="-", longname = "Affinity propogation clustering based on process parameters and initial conditions", dim=list(long_dimen,lat_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, grid_output$clusters)
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for DRIVERS and OBSERVATIONS with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"DRIVERS_OBS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area 
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

## DRIVERS
# Minimum air temperature
var_new = ncvar_def("tas_min", unit="K", longname = "Mean daily minimum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, AIRT_MIN)

# Maximum air temperature
var_new = ncvar_def("tas_max", unit="K", longname = "Mean daily maximum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, AIRT_MAX)

# Short Wave radiation
var_new = ncvar_def("rsds", unit="W.m-2", longname = "Mean downwelling short wave radiation", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, SWRAD)

# Atmospheric CO2 concentration
var_new = ncvar_def("co2", unit="ppm", longname = "Mean atmospheric CO2 concentration", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, CO2)

# Precipitation
var_new = ncvar_def("pr", unit="kg.m-2.s-1", longname = "Mean precipitation - combined liquid and solid phase", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, PRECIP)

# Forest loss fraction
var_new = ncvar_def("forest_loss_fraction", unit="1", longname = "Forest loss fraction", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, FLOSS_FRAC)

# Burnt fraction
var_new = ncvar_def("burntArea", unit="1", longname = "Burnt fraction", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, BURNT_FRAC)

# Wind Speed
var_new = ncvar_def("wsp", unit="m.s-1", longname = "Mean wind speed", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, WINDSPD)

# Vapour pressure deficit
var_new = ncvar_def("vpd", unit="Pa", longname = "Mean vapour pressure deficit", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, VPD)

## STATE OBSERVATIONS
# Leaf area index
var_new = ncvar_def("LAI_OBS", unit="m-2.m-2", longname = "Observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, LAI_OBS)

# Leaf area index Uncertainty
var_new = ncvar_def("LAI_UNC_OBS", unit="m-2.m-2", longname = "Uncertainty on observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, LAI_UNC_OBS)

# Wood stock
var_new = ncvar_def("WOOD_OBS", unit="g.m-2", longname = "Observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, WOOD_OBS)

# Wood stock uncertainty
var_new = ncvar_def("WOOD_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, WOOD_UNC_OBS)

# Soil stock
var_new = ncvar_def("SOIL_OBS", unit="g.m-2", longname = "Observed soil stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, SOIL_OBS)

# Soil stock uncertainty
var_new = ncvar_def("SOIL_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, SOIL_UNC_OBS)

## FLUX OBSERVATIONS
# GPP
var_new = ncvar_def("GPP_OBS", unit="g.m-2.d-1", longname = "Observed gross primary productivity of C (GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, GPP_OBS)

# GPP uncertainty
var_new = ncvar_def("GPP_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed gross primary productivity of C (GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, GPP_UNC_OBS)

# NEE 
var_new = ncvar_def("NEE_OBS", unit="g.m-2.d-1", longname = "Observed net ecosystem exchange of C (NEE = Reco - GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NEE_OBS)

# NEE uncertainty
var_new = ncvar_def("NEE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed net ecosystem exchange of C (NEE = Reco - GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NEE_UNC_OBS)

# NBE 
var_new = ncvar_def("NBE_OBS", unit="g.m-2.d-1", longname = "Observed net biome exchange of C (NEE = Reco - GPP + FIRE)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NBE_OBS)

# NBE uncertainty
var_new = ncvar_def("NBE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed net biome exchange of C (NEE = Reco - GPP + FIRE)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NBE_UNC_OBS)

# Fire 
var_new = ncvar_def("FIRE_OBS", unit="g.m-2.d-1", longname = "Observed C emission due to fire combustion", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, FIRE_OBS)

# Fire uncertainty
var_new = ncvar_def("FIRE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed C emission due to fire combustion", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, FIRE_UNC_OBS)

## close the file to write to disk
nc_close(new_file)

###
## Create file for H2O FLUXES with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"H2OFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area 
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## H2O FLUXES
###

# Evapotranspiration
if(exists("ET")) {
   var_new  = ncvar_def("evapotrans", unit="kg.m-2.s-1", longname = "Evapotranspiration - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ET[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("evapotrans_2.5pc", unit="kg.m-2.s-1", longname = "Evapotranspiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ET[,,low_quant,])
   # High CI
   var_new  = ncvar_def("evapotrans_97.5pc", unit="kg.m-2.s-1", longname = "Evapotranspiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  ET[,,high_quant,])
}

# Transpiration
if(exists("Etrans")) {
   var_new  = ncvar_def("tran", unit="kg.m-2.s-1", longname = "Transpiration - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Etrans[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("tran_2.5pc", unit="kg.m-2.s-1", longname = "Transpiration - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Etrans[,,low_quant,])
   # High CI
   var_new  = ncvar_def("tran_97.5pc", unit="kg.m-2.s-1", longname = "Transpiration - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Etrans[,,high_quant,])
}

# Soil evaporation
if(exists("Esoil")) {
   var_new  = ncvar_def("evspsblsoi", unit="kg.m-2.s-1", longname = "Soil evaporation - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Esoil[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("evspsblsoi_2.5pc", unit="kg.m-2.s-1", longname = "Soil evaporation - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Esoil[,,low_quant,])
   # High CI
   var_new  = ncvar_def("evspsblsoi_97.5pc", unit="kg.m-2.s-1", longname = "Soil evaporation - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Esoil[,,high_quant,])
}

# Canopy intercepted rainfall evaporation
if(exists("Ewetevap")) {
   var_new  = ncvar_def("evspsblveg", unit="kg.m-2.s-1", longname = "Canopy intercepted rainfall evaporation - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Ewetcanopy[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("evspsblveg_2.5pc", unit="kg.m-2.s-1", longname = "Canopy intercepted rainfall evaporation - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Ewetcanopy[,,low_quant,])
   # High CI
   var_new  = ncvar_def("evspsblveg_97.5pc", unit="kg.m-2.s-1", longname = "Canopy intercepted rainfall evaporation - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  Ewetcanopy[,,high_quant,])   
}

# Total drainage from soil surface and bottom of soil column
if(exists("total_drainage")) {
   var_new  = ncvar_def("mrro", unit="kg.m-2.s-1", longname = "Total drainage from soil surface and bottom of soil column - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, total_drainage[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("mrro_2.5pc", unit="kg.m-2.s-1", longname = "Total drainage from soil surface and bottom of soil column - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  total_drainage[,,low_quant,])
   # High CI
   var_new  = ncvar_def("mrro_97.5pc", unit="kg.m-2.s-1", longname = "Total drainage from soil surface and bottom of soil column - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  total_drainage[,,high_quant,])   
}

# Runoff from soil surface
if(exists("runoff")) {
   var_new  = ncvar_def("runoff", unit="kg.m-2.s-1", longname = "Soil surface water runoff - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, runoff[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("runoff_2.5pc", unit="kg.m-2.s-1", longname = "Soil surface water runoff - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  runoff[,,low_quant,])
   # High CI
   var_new  = ncvar_def("runoff_97.5pc", unit="kg.m-2.s-1", longname = "Soil surface water runoff - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  runoff[,,high_quant,])   
}

# Runoff from soil surface
if(exists("underflow")) {
   var_new  = ncvar_def("underflow", unit="kg.m-2.s-1", longname = "Water drainage from the bottom of the soil column - Median estimate", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, underflow[,,mid_quant,])
   # low CI
   var_new  = ncvar_def("underflow_2.5pc", unit="kg.m-2.s-1", longname = "Water drainage from the bottom of the soil column - 2.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  underflow[,,low_quant,])
   # High CI
   var_new  = ncvar_def("underflow_97.5pc", unit="kg.m-2.s-1", longname = "Water drainage from the bottom of the soil column - 97.5% quantile", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="single",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  underflow[,,high_quant,])   
}

###
## close the file to write to disk
###

nc_close(new_file)

