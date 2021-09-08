
###
## Process Default CARDAMOM output files into NetCDF
### 

# load needed libraries
library(ncdf4)

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# Load needed functions 
library(zoo) ; library(compiler) ; library(raster)
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/read_binary_file_format.r")

# load the CARDAMOM files
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/ODA_Kenya/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_MHMCMC/UK_5km_baseline/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/NCEO_ODA_Kenya_Forest/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/majadas_0.25deg_C7/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/sodankyla_0.25deg_C7/infofile.RData")
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))

# Specify any extra information for the filename
output_prefix = paste(PROJECT$name,"_",sep="") # follow with "_"
output_suffix = "" # begin with "_"

# Time information
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
steps_per_year = dim(grid_output$lai_m2m2)[3] / nos_years

# create lat / long axes, assumes WGS-84 grid
latitude = seq(PROJECT$latitude[1]+(PROJECT$resolution*0.5),PROJECT$latitude[2]-(PROJECT$resolution*0.5), length.out = PROJECT$lat_dim)
longitude = seq(PROJECT$longitude[1]+(PROJECT$resolution*0.5),PROJECT$longitude[2]-(PROJECT$resolution*0.5), length.out = PROJECT$long_dim)

# restricture each of the variable into a complete grid
quantiles_wanted = c(0.025,0.05,0.25,0.50,0.75,0.95,0.975)
#quantiles_wanted = c(0.025,0.25,0.50,0.75,0.975)
nos_quantiles = length(quantiles_wanted)
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
# OBSERVATIONS
LAI_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
LAI_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WOOD_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WOOD_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SOIL_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SOIL_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
# STATES
LAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
TOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
LAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
FOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
ROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
WOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
LIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
SOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
#WLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))
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
# PARAMETERS
PARS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,max(PROJECT$model$nopars)))

# Fill the arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]

         # Read in site specific drivers and parameters
#         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")) 
         load(paste(PROJECT$results_processedpath,"/",PROJECT$site[n],"_parameters.RData",sep=""))
         rm(MTT,aMTT,aNPP,parameter_covariance,SS,site_ctessel_pft) # tidy away unwanted information

         # Replace any -9999 missing values within the observation object with NA
         drivers$obs[which(drivers$obs == -9999)] = NA
         # Restructure the parameter array
         dims = dim(parameters) ; parameters = array(parameters, dim=c(dims[1],prod(dims[2:3])))

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
         # PARAMETERS
         PARS[grid_output$i_location[n],grid_output$j_location[n],,1:PROJECT$model$nopars[n]] = apply(parameters,1,quantile, prob=quantiles_wanted)[,1:PROJECT$model$nopars[n]]
         # OBSERVATIONS
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
         # STATES
         LAI[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lai_m2m2[n,,] 
         TOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$totalC_gCm2[n,,]
         LAB[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$labile_gCm2[n,,]
         FOL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$foliage_gCm2[n,,] 
         ROOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$roots_gCm2[n,,]
         WOOD[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wood_gCm2[n,,]
         LIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lit_gCm2[n,,]
         SOIL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$som_gCm2[n,,]
#         WLIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$litwood_gCm2[n,,]
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

     } # Does the file exist / has it been processed

} # site loop

# Adjust units were needed
AIRT_MIN = AIRT_MIN + 273.15 # C -> K
AIRT_MAX = AIRT_MAX + 273.15 # C -> K
SWRAD = SWRAD * 1e6 * (1/86400) # MJ/m2/day -> W/m2

## define dimension
lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
par_dimen <- ncdim_def( "parameter", units="-", 1:max(PROJECT$model$nopars))
year_dimen <- ncdim_def( "year", units="", 1:nos_years)

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), 
                 dim=list(time_dimen), missval = -99999, prec="double", compression = 9)
## STATES
# LAI
var1  = ncvar_def("lai_ensemble",       unit="m2.m-2", longname = "Leaf Area Index - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Labile
var2  = ncvar_def("cLabile_ensemble",   unit="gC.m-2", longname = "Carbon in labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Foliar
var3 = ncvar_def("cLeaf_ensemble",     unit="gC.m-2", longname = "Carbon in leaves - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var4 = ncvar_def("cFineRoot_ensemble", unit="gC.m-2", longname = "Carbon in fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var5 = ncvar_def("cWoodTotal_ensemble",unit="gC.m-2", longname = "Carbon in (AGB + BGB) wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Foliar + fine root litter
var6 = ncvar_def("cLeafFineRootLitter_ensemble",   unit="gC.m-2", longname = "Carbon in (Foliar + fine root) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
## Wood litter
#var7 = ncvar_def("cWoodLitter_ensemble", unit="gC.m-2", longname = "Carbon in (wood) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Soil organic matter
var8 = ncvar_def("cSOM_ensemble",      unit="gC.m-2", longname = "Carbon in soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Dead organic matter
var9 = ncvar_def("cDOM_ensemble",      unit="gC.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Biomass
var10 = ncvar_def("cVeg_ensemble",     unit="gC.m-2", longname = "Carbon in live biomass - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# TotalC
var11 = ncvar_def("cTotal_ensemble",   unit="gC.m-2", longname = "Carbon in live and dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                   var1,var2,var3,var4,var5,var6,
                                                  #var7,
                                                   var8,var9,var10,var11),
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
## Wood LIT
#ncvar_put(new_file, var7,  WLIT)
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
var1  = ncvar_def("gpp_ensemble",   unit="gC.m-2.d-1", longname = "Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Autotrophic respiration
var2  = ncvar_def("ra_ensemble",    unit="gC.m-2.d-1", longname = "Autotrophic (Plant) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Heterotrophic respiration
var3  = ncvar_def("rh_ensemble",    unit="gC.m-2.d-1", longname = "Heterotrophic Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Ecosystem respiration
var4  = ncvar_def("reco_ensemble",  unit="gC.m-2.d-1", longname = "Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Primary Productivity
var5  = ncvar_def("npp_ensemble",   unit="gC.m-2.d-1", longname = "Net Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Ecosystem Exchange
var6  = ncvar_def("nee_ensemble",   unit="gC.m-2.d-1", longname = "Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Biome Exchange
var7  = ncvar_def("nbe_ensemble",   unit="gC.m-2.d-1", longname = "Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Biome Exchange
var8  = ncvar_def("nbp_ensemble",   unit="gC.m-2.d-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fire emissions
var9  = ncvar_def("fFire_ensemble", unit="gC.m-2.d-1", longname = "Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Flux from forest loss
var10 = ncvar_def("fLuc_ensemble",  unit="gC.m-2.d-1", longname = "Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# CiCa
var11 = ncvar_def("CiCa_ensemble",  unit="1", longname = "Internal:Ambiant CO2 ratio - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                   var1,var2,var3,var4,var5,var6,var7,        
                                                   var8,var9,var10,var11),
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
var9 = ncvar_def("NPP_fol_flx_ensemble", unit="gC.m-2.d-1", longname = "Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var10 = ncvar_def("NPP_root_flx_ensemble", unit="gC.m-2.d-1", longname = "Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var11 = ncvar_def("NPP_wood_flx_ensemble", unit="gC.m-2.d-1", longname = "Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   
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
ncvar_put(new_file, var1,  fMRT)
# ROOT
ncvar_put(new_file, var2,  rMRT)
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

## OBSERVATIONS
# Leaf area index
var10 = ncvar_def("LAI_OBS", unit="m-2.m-2", longname = "Observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Leaf area index Uncertainty
var11 = ncvar_def("LAI_UNC_OBS", unit="m-2.m-2", longname = "Uncertainty on observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock
var12 = ncvar_def("WOOD_OBS", unit="g.m-2", longname = "Observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock uncertainty
var13 = ncvar_def("WOOD_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock
var14 = ncvar_def("SOIL_OBS", unit="g.m-2", longname = "Observed soil stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock uncertainty
var15 = ncvar_def("SOIL_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   
## PARAMETERS
# Parameter posteriors
var16 = ncvar_def("PARS", unit="-", longname = "Model / location specific parameter distributions", dim=list(long_dimen,lat_dimen,quantile_dimen,par_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"DRIVERS_OBS_PARS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,                                     
                                                   var1,var2,var3,var4,var5,var6,var7,       
                                                   var8,var9,var10,var11,var12,var13,var14,var15,
                                                   var16),
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
## OBSERVATIONS
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
## PARAMETERS
# Parameter posteriors
ncvar_put(new_file, var16, PARS)
## close the file to write to disk
nc_close(new_file)

