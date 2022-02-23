
###
## Process Default CARDAMOM output files into NetCDF
### 

# load needed libraries
library(ncdf4)
library(raster)
library(compiler)

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# Load needed functions 
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/load_all_cardamom_functions.r")

# set input and output directories
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire_nbe/"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Mexico_1deg_C7_agb_lca_gpp_fire_nbe/"
# Specify any extra information for the filename
output_prefix = "" # follow with "_"
output_suffix = "" # begin with "_"

# load the CARDAMOM files
load(paste(input_dir,"/infofile.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))

# Time information
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
steps_per_year = dim(grid_output$lai_m2m2)[3] / nos_years

# create lat / long axes, assumes regular WGS-84 grid
output = determine_lat_long_needed(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution,PROJECT$grid_type,PROJECT$waterpixels)
# NOTE: rev due to CARDAMOM grid being inverse of what comes out of raster function. Should consider changing this at some point.
longitude = output$obs_long_grid[,1] ; latitude = rev(output$obs_lat_grid[1,]) 
# Tidy up
rm(output) ; gc(reset=TRUE,verbose=FALSE)

# restricture each of the variable into a complete grid
quantiles_wanted = c(0.025,0.05,0.25,0.50,0.75,0.95,0.975)
nos_quantiles = length(quantiles_wanted)

## At model time step varaibles
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
NEE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NEE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NBE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NBE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FIRE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FIRE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
# STATES
LAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
TOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
LAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
FOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
ROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
WOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
LIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
SOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
WLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
DOM = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
BIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
# FLUXES
GPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
RAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
RHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
NPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
FIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
HARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
RECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
NEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
NBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
NBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
# BIOPHYSICAL
CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
# NPP allocation (labile, foliar, fine root, wood, gC/m2/day)
fNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
wNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
rNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
# NPP allocation (labile, foliar, fine root, wood, fraction)
fNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
wNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
rNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
# MRT (foliar, wood, fine root, litter, soil; years)
fMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
wMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
rMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
lMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
sMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))

## Annual time step variables
# Natural outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day)
NAToutflux_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
NAToutflux_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
NAToutflux_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
NAToutflux_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
NAToutflux_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
NAToutflux_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
# Fire mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day)
FIRElitter_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIRElitter_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIRElitter_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIRElitter_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIRElitter_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIRElitter_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
# Fire combustion outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day)
FIREemiss_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIREemiss_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIREemiss_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIREemiss_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIREemiss_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))
FIREemiss_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))

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
         AIRT_MIN[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,2] # mint C
         AIRT_MAX[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,3] # maxt C
         SWRAD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,4] # SWRAD MJ/m2/day
         CO2[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,5] # CO2 ppm
         DOY[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,6] # Julian day of year
         PRECIP[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,7] # precipitation kgH2O/m2/s
         FLOSS_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,8] # Forest loss fraction
         BURNT_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,9] # Burned fraction
         WINDSPD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,15] # Wind speed m/s
         VPD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,16] # Vapour pressure deficit Pa

         # STATE OBSERVATIONS
         LAI_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,3] # LAI m2/m2
         LAI_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,4] # LAI UNC m2/m2
         # ...first assign priors to the 1st time step to simplify our storage
         WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriors[21]
         WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriorunc[21]
         SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriors[23]
         SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriorunc[23]
         # ...second assign time series information if any exists
         WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),13]
         WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),14]
         SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),19]
         SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),20]
         # FLUX OBSERVATIONS (gC/m2/day)
         NEE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,5] 
         NEE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,6]
         NBE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,35]
         NBE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,36]
         FIRE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,7]
         FIRE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,8]

         ## At model time step
         # STATES
         LAI[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lai_m2m2[n,,] 
         TOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$totalC_gCm2[n,,]
         LAB[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$labile_gCm2[n,,]
         FOL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$foliage_gCm2[n,,] 
         ROOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$roots_gCm2[n,,]
         WOOD[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wood_gCm2[n,,]
         LIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lit_gCm2[n,,]
         SOIL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$som_gCm2[n,,]
         if (length(which(names(grid_output) == "litwood_gCm2")) > 0) {
             WLIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$litwood_gCm2[n,,]
         }
         DOM[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dom_gCm2[n,,]
         BIO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$biomass_gCm2[n,,]
         # FLUXES
         GPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gpp_gCm2day[n,,]
         RAU[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rauto_gCm2day[n,,]
         RHE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rhet_gCm2day[n,,]
         NPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$npp_gCm2day[n,,]
         FIR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$fire_gCm2day[n,,]
         HARV[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$harvest_gCm2day[n,,]
         RECO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$reco_gCm2day[n,,]
         NEE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nee_gCm2day[n,,]
         NBE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbe_gCm2day[n,,]
         NBP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbp_gCm2day[n,,]
         # BIOPHYSICAL
         CiCa[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$CiCa[n,,]
         # NPP (foliar, root, wood; gC/m2/day)
         fNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$fnpp_gCm2day[n,,] 
         rNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rnpp_gCm2day[n,,]
         wNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wnpp_gCm2day[n,,]
        
         # NPP (fraction) and MRT years are requested to have same number of time steps as stocks and fluxes
         # This is awkward as no easy way to repeat specific elements without loop for variables which have no meaningful value at sub-annual timescales
         # (and are therefore calculated as annuals)
         for (q in seq(1, nos_quantiles)) {
              # MRT
              fMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$aMTT_foliar_years[grid_output$i_location[n],grid_output$j_location[n],,q], each = steps_per_year)
              rMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$aMTT_root_years[grid_output$i_location[n],grid_output$j_location[n],,q], each = steps_per_year)
              wMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$aMTT_wood_years[grid_output$i_location[n],grid_output$j_location[n],,q], each = steps_per_year)
              lMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$aMTT_DeadOrg_years[grid_output$i_location[n],grid_output$j_location[n],,q], each = steps_per_year)
              sMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$aMTT_som_years[grid_output$i_location[n],grid_output$j_location[n],,q], each = steps_per_year)
              # NPP
              fNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$NPP_foliar_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))
              rNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$NPP_root_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))
              wNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$NPP_wood_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))
         } # loop quantiles

         # At annual time step
         # Natural outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day)
         NAToutflux_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$NAToutflux_labile_gCm2yr[n,,]
         NAToutflux_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$NAToutflux_foliar_gCm2yr[n,,]
         NAToutflux_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$NAToutflux_wood_gCm2yr[n,,]
         NAToutflux_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$NAToutflux_root_gCm2yr[n,,]
         NAToutflux_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$NAToutflux_litter_gCm2yr[n,,]
         NAToutflux_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$NAToutflux_som_gCm2yr[n,,]
         # Fire mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day)
         FIRElitter_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_labile_gCm2yr[n,,]
         FIRElitter_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_foliar_gCm2yr[n,,]
         FIRElitter_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_wood_gCm2yr[n,,]
         FIRElitter_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_root_gCm2yr[n,,]
         FIRElitter_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_litter_gCm2yr[n,,]
         FIRElitter_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_som_gCm2yr[n,,]
         # Fire combustion outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day)
         FIREemiss_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_labile_gCm2yr[n,,]
         FIREemiss_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_foliar_gCm2yr[n,,]
         FIREemiss_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_wood_gCm2yr[n,,]
         FIREemiss_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_root_gCm2yr[n,,]
         FIREemiss_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_litter_gCm2yr[n,,]
         FIREemiss_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_som_gCm2yr[n,,]

     } # Does the file exist / has it been processed

} # site loop

# Unit conversions
# STOCKS (gC/m2 -> kgC/m2)
TOT = TOT * 1e-3
BIO = BIO * 1e-3
DOM = DOM * 1e-3
LAB = LAB * 1e-3
FOL = FOL * 1e-3
ROOT = ROOT * 1e-3
WOOD = WOOD * 1e-3
LIT = LIT * 1e-3
WLIT = WLIT * 1e-3
SOIL = SOIL * 1e-3
# FLUXES (gC/m2/day -> kgC/m2/s)
GPP = GPP * 1e-3 * (1/86400)
RAU = RAU * 1e-3 * (1/86400)
RHE = RHE * 1e-3 * (1/86400)
NPP = NPP * 1e-3 * (1/86400)
FIR = FIR * 1e-3 * (1/86400)
HARV = HARV * 1e-3 * (1/86400)
RECO = RECO * 1e-3 * (1/86400)
NEE = NEE * 1e-3 * (1/86400)
NBE = NBE * 1e-3 * (1/86400)
NBP = NBP * 1e-3 * (1/86400)
# NPP (foliar, root, wood; gC/m2/day -> kgC/m2/s)
fNPP_FLX = fNPP_FLX * 1e-3 * (1/86400)
rNPP_FLX = rNPP_FLX * 1e-3 * (1/86400)
wNPP_FLX = wNPP_FLX * 1e-3 * (1/86400)
# Natural outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
NAToutflux_labile = NAToutflux_labile * 1e-3 * (1/86400)
NAToutflux_foliage = NAToutflux_foliage * 1e-3 * (1/86400)
NAToutflux_wood = NAToutflux_wood * 1e-3 * (1/86400)
NAToutflux_root = NAToutflux_root * 1e-3 * (1/86400)
NAToutflux_litter = NAToutflux_litter * 1e-3 * (1/86400)
NAToutflux_som = NAToutflux_som * 1e-3 * (1/86400)
# Fire mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
FIRElitter_labile = FIRElitter_labile * 1e-3 * (1/86400)
FIRElitter_foliage = FIRElitter_foliage * 1e-3 * (1/86400)
FIRElitter_wood = FIRElitter_wood * 1e-3 * (1/86400)
FIRElitter_root = FIRElitter_root * 1e-3 * (1/86400)
FIRElitter_litter = FIRElitter_litter * 1e-3 * (1/86400)
FIRElitter_som = FIRElitter_som * 1e-3 * (1/86400)
# Fire combustion outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
FIREemiss_labile = FIREemiss_labile * 1e-3 * (1/86400)
FIREemiss_foliage = FIREemiss_foliage * 1e-3 * (1/86400)
FIREemiss_wood = FIREemiss_wood * 1e-3 * (1/86400)
FIREemiss_root = FIREemiss_root * 1e-3 * (1/86400)
FIREemiss_litter = FIREemiss_litter * 1e-3 * (1/86400)
FIREemiss_som = FIREemiss_som * 1e-3 * (1/86400)

## define dimension
lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
year_dimen <- ncdim_def( "year", units="", 1:nos_years)

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), 
                 dim=list(time_dimen), missval = -99999, prec="double", compression = 9)
## STATES
# LAI
var1  = ncvar_def("lai_ensemble",       unit="m2.m-2", longname = "Leaf Area Index - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Labile
var2  = ncvar_def("cLabile_ensemble",   unit="kg.m-2", longname = "Carbon in labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Foliar
var3 = ncvar_def("cLeaf_ensemble",     unit="kg.m-2", longname = "Carbon in leaves - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var4 = ncvar_def("cFineRoot_ensemble", unit="kg.m-2", longname = "Carbon in fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var5 = ncvar_def("cWoodTotal_ensemble",unit="kg.m-2", longname = "Carbon in (AGB + BGB) wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Foliar + fine root litter
var6 = ncvar_def("cLeafFineRootLitter_ensemble",   unit="kg.m-2", longname = "Carbon in (Foliar + fine root) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
## Wood litter
var7 = ncvar_def("cWoodLitter_ensemble", unit="kg.m-2", longname = "Carbon in (wood) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Soil organic matter
var8 = ncvar_def("cSOM_ensemble",      unit="kg.m-2", longname = "Carbon in soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Dead organic matter
var9 = ncvar_def("cDOM_ensemble",      unit="kg.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Biomass
var10 = ncvar_def("cVeg_ensemble",     unit="kg.m-2", longname = "Carbon in live biomass - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# TotalC
var11 = ncvar_def("cTotal_ensemble",   unit="kg.m-2", longname = "Carbon in live and dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                   var1,var2,var3,var4,var5,var6,
                                                   var7,var8,var9,var10,var11),
                                                   force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])
## STATES
# LAI
ncvar_put(new_file, var1,  LAI)
# LAB
ncvar_put(new_file, var2,  LAB)
# FOL
ncvar_put(new_file, var3,  FOL)
# ROOT
ncvar_put(new_file, var4,  ROOT)
# WOOD
ncvar_put(new_file, var5,  WOOD)
# LIT
ncvar_put(new_file, var6,  LIT)
# Wood LIT
ncvar_put(new_file, var7,  WLIT)
# SOIL
ncvar_put(new_file, var8,  SOIL)
# DOM
ncvar_put(new_file, var9,  DOM)
# Biomass
ncvar_put(new_file, var10, BIO)
# Total C
ncvar_put(new_file, var11, TOT)

## close the file to write to disk
nc_close(new_file)

## FLUXES
# GPP
var1  = ncvar_def("gpp_ensemble",   unit="kg.m-2.s-1", longname = "Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Autotrophic respiration
var2  = ncvar_def("ra_ensemble",    unit="kg.m-2.s-1", longname = "Autotrophic (Plant) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Heterotrophic respiration
var3  = ncvar_def("rh_ensemble",    unit="kg.m-2.s-1", longname = "Heterotrophic Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Ecosystem respiration
var4  = ncvar_def("reco_ensemble",  unit="kg.m-2.s-1", longname = "Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Primary Productivity
var5  = ncvar_def("npp_ensemble",   unit="kg.m-2.s-1", longname = "Net Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Ecosystem Exchange
var6  = ncvar_def("nee_ensemble",   unit="kg.m-2.s-1", longname = "Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Biome Exchange
var7  = ncvar_def("nbe_ensemble",   unit="kg.m-2.s-1", longname = "Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Biome Exchange
var8  = ncvar_def("nbp_ensemble",   unit="kg.m-2.s-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fire emissions
var9  = ncvar_def("fFire_ensemble", unit="kg.m-2.s-1", longname = "Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Flux from forest loss
var10 = ncvar_def("fLuc_ensemble",  unit="kg.m-2.s-1", longname = "Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# CiCa
var11 = ncvar_def("CiCa_ensemble",  unit="1", longname = "Internal:Ambiant CO2 ratio - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
## Annual fluxes
# Natural outflux from labile
var12  = ncvar_def("NAToutflux_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Natural C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Natural outflux from foliage
var13  = ncvar_def("NAToutflux_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Natural C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Natural outflux from root
var14  = ncvar_def("NAToutflux_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Natural C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Natural outflux from wood
var15  = ncvar_def("NAToutflux_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Natural C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Natural outflux from litter
var16  = ncvar_def("NAToutflux_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Natural C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Natural outflux from som
var17  = ncvar_def("NAToutflux_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Natural C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from labile
var18  = ncvar_def("FIRElitter_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from foliage
var19  = ncvar_def("FIRElitter_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from root
var20  = ncvar_def("FIRElitter_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from wood
var21  = ncvar_def("FIRElitter_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from litter
var22  = ncvar_def("FIRElitter_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from som
var23  = ncvar_def("FIRElitter_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from labile
var24  = ncvar_def("FIREemiss_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from foliage
var25  = ncvar_def("FIREemiss_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from root
var26  = ncvar_def("FIREemiss_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from wood
var27  = ncvar_def("FIREemiss_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from litter
var28  = ncvar_def("FIREemiss_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
# Fire mortality outflux from som
var29  = ncvar_def("FIREemiss_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                   var1,var2,var3,var4,var5,var6,var7,        
                                                   var8,var9,var10,var11,var12,var13,var14,
                                                   var15,var16,var17,var18,var19,var20,
                                                   var21,var22,var23,var24,var25,var26,var27,
                                                   var28,var29),
                                                   force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])

## FLUXES
# GPP
ncvar_put(new_file, var1, GPP)
# RAU
ncvar_put(new_file, var2, RAU)
# RHE
ncvar_put(new_file, var3, RHE)
# RECO
ncvar_put(new_file, var4, RECO)
# NPP
ncvar_put(new_file, var5, NPP)
# NEE
ncvar_put(new_file, var6, NEE)
# NBE
ncvar_put(new_file, var7, NBE)
# NBP
ncvar_put(new_file, var8, NBP)
# FIR
ncvar_put(new_file, var9, FIR)
# Forest Harvest
ncvar_put(new_file, var10, HARV)
# CiCa
ncvar_put(new_file, var11, CiCa)
## Annual fluxes
# Natural output labile
ncvar_put(new_file, var12, NAToutflux_labile)
# Natural output foliage
ncvar_put(new_file, var13, NAToutflux_foliage)
# Natural output fine root
ncvar_put(new_file, var14, NAToutflux_root)
# Natural output wood
ncvar_put(new_file, var15, NAToutflux_wood)
# Natural output litter
ncvar_put(new_file, var16, NAToutflux_litter)
# Natural output som
ncvar_put(new_file, var17, NAToutflux_som)
# Fire mortality output labile
ncvar_put(new_file, var18, FIRElitter_labile)
# Fire mortality output foliage
ncvar_put(new_file, var19, FIRElitter_foliage)
# Fire mortality output fine root
ncvar_put(new_file, var20, FIRElitter_root)
# Fire mortality output wood
ncvar_put(new_file, var21, FIRElitter_wood)
# Fire mortality output litter
ncvar_put(new_file, var22, FIRElitter_litter)
# Fire mortality output som
ncvar_put(new_file, var23, FIRElitter_som)
# Fire combustion output labile
ncvar_put(new_file, var24, FIREemiss_labile)
# Fire combustion output foliage
ncvar_put(new_file, var25, FIREemiss_foliage)
# Fire combustion output fine root
ncvar_put(new_file, var26, FIREemiss_root)
# Fire combustion output wood
ncvar_put(new_file, var27, FIREemiss_wood)
# Fire combustion output litter
ncvar_put(new_file, var28, FIREemiss_litter)
# Fire combustion output som
ncvar_put(new_file, var29, FIREemiss_som)

## close the file to write to disk
nc_close(new_file)

## Mean Residence Times 
# Foliar
var1 = ncvar_def("MTT_fol_ensemble", unit="year", longname = "Mean Foliar Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var2 = ncvar_def("MTT_root_ensemble", unit="year", longname = "Mean fine root Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var3 = ncvar_def("MTT_wood_ensemble", unit="year", longname = "Mean wood Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine litter (fol + fine root)
var4 = ncvar_def("MTT_lit_ensemble", unit="year", longname = "Mean lit+litwood Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Soil
var5 = ncvar_def("MTT_som_ensemble", unit="year", longname = "Mean Soil Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

## NPP allocation fractions
# Foliar
var6 = ncvar_def("NPP_fol_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var7 = ncvar_def("NPP_root_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var8 = ncvar_def("NPP_wood_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

## NPP allocation fluxes
# Foliar
var9 = ncvar_def("NPP_fol_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var10 = ncvar_def("NPP_root_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var11 = ncvar_def("NPP_wood_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   
# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_MRT_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                     
                                                   var1,var2,var3,var4,var5,var6,var7,       
                                                   var8,var9,var10,var11),
                                                   force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])

## MTT - time series
# FOL
ncvar_put(new_file, var1, fMRT)
# ROOT
ncvar_put(new_file, var2, rMRT)
# WOOD
ncvar_put(new_file, var3, wMRT)
# LIT
ncvar_put(new_file, var4, lMRT)
# SOIL
ncvar_put(new_file, var5, sMRT)
## NPP fractions
# FOL
ncvar_put(new_file, var6, fNPP)
# ROOT
ncvar_put(new_file, var7, rNPP)
# WOOD
ncvar_put(new_file, var8, wNPP)
## NPP fluxes
# FOL
ncvar_put(new_file, var9, fNPP_FLX)
# ROOT
ncvar_put(new_file, var10, rNPP_FLX)
# WOOD
ncvar_put(new_file, var11, wNPP_FLX)

## close the file to write to disk
nc_close(new_file)

## DRIVERS
# Minimum air temperature
var1 = ncvar_def("tas_min", unit="K", longname = "Mean daily minimum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Maximum air temperature
var2 = ncvar_def("tas_max", unit="K", longname = "Mean daily maximum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Short Wave radiation
var3 = ncvar_def("rsds", unit="W.m-2", longname = "Mean downwelling short wave radiation", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Atmospheric CO2 concentration
var4 = ncvar_def("co2", unit="ppm", longname = "Mean atmospheric CO2 concentration", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Precipitation
var5 = ncvar_def("pr", unit="kg.m-2.s-1", longname = "Mean precipitation - combined liquid and solid phase", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Forest loss fraction
var6 = ncvar_def("forest_loss_fraction", unit="1", longname = "Forest loss fraction", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Burnt fraction
var7 = ncvar_def("burnt_fraction", unit="1", longname = "Burnt fraction", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wind Speed
var8 = ncvar_def("wsp", unit="m.s-1", longname = "Mean wind speed", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Vapour pressure deficit
var9 = ncvar_def("vpd", unit="Pa", longname = "Mean vapour pressure deficit", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

## STATE OBSERVATIONS
# Leaf area index
var10 = ncvar_def("LAI_OBS", unit="m-2.m-2", longname = "Observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Leaf area index Uncertainty
var11 = ncvar_def("LAI_UNC_OBS", unit="m-2.m-2", longname = "Uncertainty on observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock
var12 = ncvar_def("WOOD_OBS", unit="g.m-2", longname = "Observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock uncertainty
var13 = ncvar_def("WOOD_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Soil stock
var14 = ncvar_def("SOIL_OBS", unit="g.m-2", longname = "Observed soil stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Soil stock uncertainty
var15 = ncvar_def("SOIL_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
## FLUX OBSERVATIONS
# NEE 
var16 = ncvar_def("NEE_OBS", unit="g.m-2.d-1", longname = "Observed net ecosystem exchange of C (NEE = Reco - GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# NEE uncertainty
var17 = ncvar_def("NEE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed net ecosystem exchange of C (NEE = Reco - GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# NBE 
var18 = ncvar_def("NBE_OBS", unit="g.m-2.d-1", longname = "Observed net biome exchange of C (NEE = Reco - GPP + FIRE)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# NBE uncertainty
var19 = ncvar_def("NBE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed net biome exchange of C (NEE = Reco - GPP + FIRE)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fire 
var20 = ncvar_def("FIRE_OBS", unit="g.m-2.d-1", longname = "Observed C emission due to fire combustion", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fire uncertainty
var21 = ncvar_def("FIRE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed C emission due to fire combustion", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"DRIVERS_OBS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                     
                                                   var1,var2,var3,var4,var5,var6,var7,       
                                                   var8,var9,var10,var11,var12,var13,
                                                   var14,var15,var16,var17,var18,var19,
                                                   var20,var21),
                                                   force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])

## DRIVERS
# min air temperature
ncvar_put(new_file, var1, AIRT_MIN)
# max air temperature
ncvar_put(new_file, var2, AIRT_MAX)
# SW radiation
ncvar_put(new_file, var3, SWRAD)
# atmospheri CO2
ncvar_put(new_file, var4, CO2)
# precipitation
ncvar_put(new_file, var5, PRECIP)
# Forest loss fraction
ncvar_put(new_file, var6, FLOSS_FRAC)
# Burnt fraction
ncvar_put(new_file, var7, BURNT_FRAC)
# Wind speed
ncvar_put(new_file, var8, WINDSPD)
# Vapour pressure deficit
ncvar_put(new_file, var9, VPD)
## STATE OBSERVATIONS
# LAI
ncvar_put(new_file, var10, LAI_OBS)
# LAI UNC
ncvar_put(new_file, var11, LAI_UNC_OBS)
# WOOD
ncvar_put(new_file, var12, WOOD_OBS)
# WOOD UNC
ncvar_put(new_file, var13, WOOD_UNC_OBS)
# SOIL
ncvar_put(new_file, var14, SOIL_OBS)
# SOIL UNC
ncvar_put(new_file, var15, SOIL_UNC_OBS)
## FLUX OBSERVATIONS
# NEE 
ncvar_put(new_file, var16, NEE_OBS)
# NEE uncertainty
ncvar_put(new_file, var17, NEE_UNC_OBS)
# NBE 
ncvar_put(new_file, var18, NBE_OBS)
# NBE uncertainty
ncvar_put(new_file, var19, NBE_UNC_OBS)
# Fire 
ncvar_put(new_file, var20, FIRE_OBS)
# Fire uncertainty
ncvar_put(new_file, var21, FIRE_UNC_OBS)

## close the file to write to disk
nc_close(new_file)

