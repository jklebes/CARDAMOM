
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
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/quincy_400_disturbance/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/quincy_400_nodisturbance/infofile.RData")

# Specify any extra information for the filename
output_prefix = paste(PROJECT$name,"_",sep="") # follow with "_"
output_suffix = "" # begin with "_"

# Time information
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
steps_per_year = length(PROJECT$model$timestep_days) / nos_years

# Expected ensemble size
nos_ensemble = PROJECT$nochains * 100

# DRIVERS
AIRT_MIN = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
AIRT_MAX = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
SWRAD = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
CO2 = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
DOY = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
PRECIP = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
FLOSS_FRAC = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days))) 
BURNT_FRAC = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
WINDSPD = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
VPD = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
# OBSERVATIONS
LAI_OBS = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
LAI_UNC_OBS = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
WOOD_OBS = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
WOOD_UNC_OBS = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
SOIL_OBS = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
SOIL_UNC_OBS = array(NA, dim=c(PROJECT$nosites,length(PROJECT$model$timestep_days)))
# STATES
LAI = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
LAB = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
FOL = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
ROOT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
WOOD = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
LIT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
SOIL = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
WLIT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
BIO = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
# FLUXES
GPP = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
RAU = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
RHE = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
FIR = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
HARV = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
RECO = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
NEE = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
# BIOPHYSICAL
CiCa = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
# NPP allocation (labile, foliar, fine root, wood, fraction)
fNPP = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
wNPP = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
rNPP = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
# MRT (foliar, wood, fine root, litter, soil; years)
fMRT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
wMRT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
rMRT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
lMRT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))
sMRT = array(NA, dim=c(PROJECT$nosites,nos_ensemble,length(PROJECT$model$timestep_days)))

# Set initial flag to determine whether wood litter pool is used
WLIT_USED = FALSE

# Fill the arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Create the expected file name
     infile = paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")

     # Ensure the site has been processed
     if (file.exists(infile)) {

         # Load the file
         load(infile)

         # Determine the number of iterations stored in this file
         niter = dim(states_all$lai_m2m2)[1]

         # Replace any -9999 missing values within the observation object with NA
         drivers$obs[which(drivers$obs == -9999)] = NA

         # DRIVERS
         AIRT_MIN[n,] = drivers$met[,2] # mint C
         AIRT_MAX[n,] = drivers$met[,3] # maxt C
         SWRAD[n,] = drivers$met[,4] # SWRAD MJ/m2/day
         CO2[n,] = drivers$met[,5] # CO2 ppm
         DOY[n,] = drivers$met[,6] # Julian day of year
         PRECIP[n,] = drivers$met[,7] # precipitation kgH2O/m2/s
         FLOSS_FRAC[n,] = drivers$met[,8] # Forest loss fraction
         BURNT_FRAC[n,] = drivers$met[,9] # Burned fraction
         WINDSPD[n,] = drivers$met[,15] # Wind speed m/s
         VPD[n,] = drivers$met[,16] # Vapour pressure deficit Pa
         # OBSERVATIONS
         LAI_OBS[n,] = drivers$obs[,3] # LAI m2/m2
         LAI_UNC_OBS[n,] = drivers$obs[,4] # LAI UNC m2/m2
         # ...first assign priors to the 1st time step to simplify our storage
         WOOD_OBS[n,1] = drivers$parpriors[21]
         WOOD_UNC_OBS[n,1] = drivers$parpriorunc[21]
         SOIL_OBS[n,1] = drivers$parpriors[23]
         SOIL_UNC_OBS[n,1] = drivers$parpriorunc[23]
         # ...second assign time series information if any exists
         WOOD_OBS[n,2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),13]
         WOOD_UNC_OBS[n,2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),14]
         SOIL_OBS[n,2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),19]
         SOIL_UNC_OBS[n,2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),20]
         # STATES
         LAI[n,1:niter,] = states_all$lai_m2m2 
         LAB[n,1:niter,] = states_all$lab_gCm2
         FOL[n,1:niter,] = states_all$fol_gCm2
         ROOT[n,1:niter,] = states_all$root_gCm2
         WOOD[n,1:niter,] = states_all$wood_gCm2
         LIT[n,1:niter,] = states_all$lit_gCm2
         SOIL[n,1:niter,] = states_all$som_gCm2
         if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
             WLIT[n,1:niter,] = states_all$litwood_gCm2
             WLIT_USED = TRUE
         }
         BIO[n,1:niter,] = states_all$bio_gCm2
         # FLUXES
         GPP[n,1:niter,] = states_all$gpp_gCm2day
         RAU[n,1:niter,] = states_all$rauto_gCm2day
         RHE[n,1:niter,] = states_all$rhet_gCm2day
         FIR[n,1:niter,] = states_all$fire_gCm2day
         HARV[n,1:niter,] = states_all$harvest_C_gCm2day
         RECO[n,1:niter,] = states_all$reco_gCm2day
         NEE[n,1:niter,] = states_all$nee_gCm2day
         # BIOPHYSICAL
         CiCa[n,1:niter,] = states_all$CiCa
        
         # NPP (fraction) and MRT years are requested to have same number of time steps as stocks and fluxes
         # This is awkward as no easy way to repeat specific elements without loop for variables which have no meaningful value at sub-annual timescales
         # (and are therefore calculated as annuals)
         for (q in seq(1, niter)) {
              # MRT
              fMRT[n,q,] = rep(states_all$aMTT[q,1,], each = steps_per_year)
              rMRT[n,q,] = rep(states_all$aMTT[q,2,], each = steps_per_year)
              wMRT[n,q,] = rep(states_all$aMTT[q,3,], each = steps_per_year)
              lMRT[n,q,] = rep(states_all$aMTT[q,4,], each = steps_per_year)
              sMRT[n,q,] = rep(states_all$aMTT[q,5,], each = steps_per_year)
              # NPP
              fNPP[n,q,] = rep(states_all$aNPP[q,1], each = length(PROJECT$model$timestep_days))
              rNPP[n,q,] = rep(states_all$aNPP[q,2], each = length(PROJECT$model$timestep_days))
              wNPP[n,q,] = rep(states_all$aNPP[q,3], each = length(PROJECT$model$timestep_days))
         } # loop quantiles

     } # Does the file exist / has it been processed

} # site loop

# Adjust units were needed
AIRT_MIN = AIRT_MIN + 273.15 # C -> K
AIRT_MAX = AIRT_MAX + 273.15 # C -> K
SWRAD = SWRAD * 1e6 * (1/86400) # MJ/m2/day -> W/m2

## define dimension
site_dimen <- ncdim_def( "site", units="Site Name", 1:PROJECT$nosites )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
ensemble_dimen <- ncdim_def( "ensemble", units="-", 1:nos_ensemble)

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), 
                 dim=list(time_dimen), missval = -99999, prec="double", compression = 9)
varName = ncvar_def("SiteName", units = "-", longname = "Name of site", 
                    dim=list(site_dimen), prec="char", compression = 9)
## STATES
# LAI
var1  = ncvar_def("lai_ensemble",       unit="m2.m-2", longname = "Leaf Area Index - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Labile
var2  = ncvar_def("cLabile_ensemble",   unit="gC.m-2", longname = "Carbon in labile - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Foliar
var3 = ncvar_def("cLeaf_ensemble",     unit="gC.m-2", longname = "Carbon in leaves - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var4 = ncvar_def("cFineRoot_ensemble", unit="gC.m-2", longname = "Carbon in fine root - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var5 = ncvar_def("cWoodTotal_ensemble",unit="gC.m-2", longname = "Carbon in (AGB + BGB) wood - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Foliar + fine root litter
var6 = ncvar_def("cLeafFineRootLitter_ensemble",   unit="gC.m-2", longname = "Carbon in (Foliar + fine root) litter - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# var7 see below
# Soil organic matter
var8 = ncvar_def("cSOM_ensemble",      unit="gC.m-2", longname = "Carbon in soil organic matter - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Biomass
var9 = ncvar_def("cVeg_ensemble",     unit="gC.m-2", longname = "Carbon in live biomass - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

# Depending on whether an explicit wood litter pool exists
if (WLIT_USED) {
    # Wood litter
    var7 = ncvar_def("cWoodLitter_ensemble", unit="gC.m-2", longname = "Carbon in (wood) litter - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

    # create the empty file
    output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
    new_file=nc_create(filename=output_name, vars=list(var0,varName,                                     
                                                       var1,var2,var3,var4,var5,var6,
                                                       var7,var8,var9),
                                                       force_v4 = TRUE)
} else {
    # create the empty file
    output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
    new_file=nc_create(filename=output_name, vars=list(var0,varName,var1,var2,var3,var4,var5,var6,var8,var9),force_v4 = TRUE)
} # WLIT_USED

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])
## SiteName
ncvar_put(new_file, varName, PROJECT$sites)
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
if (WLIT_USED) {ncvar_put(new_file, var7,  WLIT)}
# SOIL
ncvar_put(new_file, var8,  SOIL)
# Biomass
ncvar_put(new_file, var9, BIO)

## close the file to write to disk
nc_close(new_file)

## FLUXES
# GPP
var1  = ncvar_def("gpp_ensemble",   unit="gC.m-2.d-1", longname = "Gross Primary Productivity - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Autotrophic respiration
var2  = ncvar_def("ra_ensemble",    unit="gC.m-2.d-1", longname = "Autotrophic (Plant) Respiration - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Heterotrophic respiration
var3  = ncvar_def("rh_ensemble",    unit="gC.m-2.d-1", longname = "Heterotrophic Respiration - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Ecosystem respiration
var4  = ncvar_def("reco_ensemble",  unit="gC.m-2.d-1", longname = "Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Net Ecosystem Exchange
var5  = ncvar_def("nee_ensemble",   unit="gC.m-2.d-1", longname = "Net Ecosystem Exchange - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fire emissions
var6  = ncvar_def("fFire_ensemble", unit="gC.m-2.d-1", longname = "Fire - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Flux from forest loss
var7 = ncvar_def("fLuc_ensemble",  unit="gC.m-2.d-1", longname = "Forest harvest - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# CiCa
var8 = ncvar_def("CiCa_ensemble",  unit="1", longname = "Internal:Ambiant CO2 ratio - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,varName,                                  
                                                   var1,var2,var3,var4,var5,var6,var7,        
                                                   var8),
                                                   force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])
## SiteName
ncvar_put(new_file, varName, PROJECT$sites)
## FLUXES
# GPP
ncvar_put(new_file, var1, GPP)
# RAU
ncvar_put(new_file, var2, RAU)
# RHE
ncvar_put(new_file, var3, RHE)
# RECO
ncvar_put(new_file, var4, RECO)
# NEE
ncvar_put(new_file, var5, NEE)
# FIR
ncvar_put(new_file, var6, FIR)
# Forest Harvest
ncvar_put(new_file, var7, HARV)
# CiCa
ncvar_put(new_file, var8, CiCa)

## close the file to write to disk
nc_close(new_file)

## Mean Residence Times 
# Foliar
var1 = ncvar_def("MTT_fol_ensemble", unit="year", longname = "Mean Foliar Transit Time - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var2 = ncvar_def("MTT_root_ensemble", unit="year", longname = "Mean fine root Transit Time - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var3 = ncvar_def("MTT_wood_ensemble", unit="year", longname = "Mean wood Transit Time - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine litter (fol + fine root)
var4 = ncvar_def("MTT_lit_ensemble", unit="year", longname = "Mean lit+litwood Transit Time - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Soil
var5 = ncvar_def("MTT_som_ensemble", unit="year", longname = "Mean Soil Transit Time - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

## NPP allocation fractions
# Foliar
var6 = ncvar_def("NPP_fol_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to foliage - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Fine root
var7 = ncvar_def("NPP_root_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to fine root - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood
var8 = ncvar_def("NPP_wood_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to wood - Ensemble", dim=list(site_dimen,ensemble_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   
# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_MRT_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,varName,                                  
                                                   var1,var2,var3,var4,var5,var6,var7,var8),
                                                   force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, drivers$met[,1])
## SiteName
ncvar_put(new_file, varName, PROJECT$sites)
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

## close the file to write to disk
nc_close(new_file)

## DRIVERS
# Minimum air temperature
var1 = ncvar_def("tas_min", unit="K", longname = "Mean daily minimum near surface air temperature", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Maximum air temperature
var2 = ncvar_def("tas_max", unit="K", longname = "Mean daily maximum near surface air temperature", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Short Wave radiation
var3 = ncvar_def("rsds", unit="W.m-2", longname = "Mean downwelling short wave radiation", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Atmospheric CO2 concentration
var4 = ncvar_def("co2", unit="ppm", longname = "Mean atmospheric CO2 concentration", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Precipitation
var5 = ncvar_def("pr", unit="kg.m-2.s-1", longname = "Mean precipitation - combined liquid and solid phase", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Forest loss fraction
var6 = ncvar_def("forest_loss_fraction", unit="1", longname = "Forest loss fraction", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Burnt fraction
var7 = ncvar_def("burnt_fraction", unit="1", longname = "Burnt fraction", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wind Speed
var8 = ncvar_def("wsp", unit="m.s-1", longname = "Mean wind speed", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Vapour pressure deficit
var9 = ncvar_def("vpd", unit="Pa", longname = "Mean vapour pressure deficit", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

## OBSERVATIONS
# Leaf area index
var10 = ncvar_def("LAI_OBS", unit="m-2.m-2", longname = "Observed Leaf area index", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Leaf area index Uncertainty
var11 = ncvar_def("LAI_UNC_OBS", unit="m-2.m-2", longname = "Uncertainty on observed Leaf area index", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock
var12 = ncvar_def("WOOD_OBS", unit="g.m-2", longname = "Observed wood stock C (above + below + coarse root)", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock uncertainty
var13 = ncvar_def("WOOD_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (above + below + coarse root)", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock
var14 = ncvar_def("SOIL_OBS", unit="g.m-2", longname = "Observed soil stock C (assumed to include wood litter)", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
# Wood stock uncertainty
var15 = ncvar_def("SOIL_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (assumed to include wood litter)", dim=list(site_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   
# create the empty file
output_name = paste(PROJECT$results_processedpath,output_prefix,"DRIVERS_OBS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
new_file=nc_create(filename=output_name, vars=list(var0,varName,                                    
                                                   var1,var2,var3,var4,var5,var6,var7,       
                                                   var8,var9,var10,var11,var12,var13,var14,var15),
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

## close the file to write to disk
nc_close(new_file)

