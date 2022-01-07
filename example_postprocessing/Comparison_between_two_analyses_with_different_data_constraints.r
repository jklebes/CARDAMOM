
###
## Compare CARDAMOM analyses using the same model version 
## but different observational constraints
## Key outcomes are how does uncertainty and estimated parameters / traits vary between.
###

#### TO DO
# Compare correlations between pixels, how have they changed to improve constraint on wood
# Ensure all relative plots are restricted +/- 1 
# Include parameter correlations check
#orig_mean_parameter_correlation = array(NA, dim=dim(grid_output$mean_lai_m2m2)[1:2])
#alt_mean_parameter_correlation = array(NA, dim=dim(grid_output$mean_lai_m2m2)[1:2])
## This time run through and find the parameter files
#for (n in seq(1, orig_PROJECT$nosites)) {
#
#     # determine the lat / long location within the grid
#     j_loc=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
#     i_loc=as.numeric(PROJECT$sites[n])-(floor(j_loc)*PROJECT$long_dim)
#     if(i_loc == 0) {i_loc = PROJECT$long_dim} ; j_loc=ceiling(j_loc)
# 
##     # Extract each sites location within the grid
##     i_loc = grid_output$i_location[n] ; j_loc = grid_output$j_location[n]
#        
#     # Check that location has run
#     if (is.na(i_loc) == FALSE & is.na(j_loc) == FALSE) {
#              
#         ###
#         ## Determining Drivers and trends
#
#         # Read in pixel paramater information
#         infile1 = paste(orig_PROJECT$results_processedpath,orig_PROJECT$sites[n],"_parameters.RData",sep="")   
#         infile2 = paste(alt_PROJECT$results_processedpath,alt_PROJECT$sites[n],"_parameters.RData",sep="")   
#         if (file.exists(infile1) & file.exists(infile2)) {
#
#             load(infile1)
#             # Rearrange the parameter matrix for the correlation analysis
#             parameters = array(parameters, dim = c(dim(parameters)[1], prod(dim(parameters)[2:3])))
#             # Calculate the mean correlation magnitude for each pixel
#             tmp = cor(t(parameters)) ; tmp = tmp[lower.tri(tmp) == 1] ; orig_mean_parameter_correlation[i_loc,j_loc] = mean(abs(tmp), na.rm=TRUE)
#
#             load(infile2)
#             # Rearrange the parameter matrix for the correlation analysis
#             parameters = array(parameters, dim = c(dim(parameters)[1], prod(dim(parameters)[2:3])))
#             # Calculate the mean correlation magnitude for each pixel
#             tmp = cor(t(parameters)) ; tmp = tmp[lower.tri(tmp) == 1] ; alt_mean_parameter_correlation[i_loc,j_loc] = mean(abs(tmp), na.rm=TRUE)
#         } # both files exist
#         
#     } # Did this location run
#} # Site loop
#
#par(mfrow=c(2,2))
#image.plot(alt_mean_parameter_correlation) ; image.plot(orig_mean_parameter_correlation) 
#image.plot(alt_mean_parameter_correlation-orig_mean_parameter_correlation) ; image.plot((orig_mean_parameter_correlation-orig_mean_parameter_correlation)/orig_mean_parameter_correlation)

#summary(as.vector(alt_mean_parameter_correlation)) ; summary(as.vector(orig_mean_parameter_correlation) )
#summary(as.vector(alt_mean_parameter_correlation-orig_mean_parameters_correlations))


# Read library
library(fields)
library(compiler)
library(RColorBrewer)
library(plotrix)
library(zoo)
library(raster)
library(abind)
#library(maptools)
library(ncdf4)

# Set CARDAMOM directory
setwd("~/WORK/GREENHOUSE/models/CARDAMOM/")
# Load any CARDAMOM functions which might be useful
source("./R_functions/generate_wgs_grid.r")
source("./R_functions/calc_pixel_area.r")
source("./R_functions/read_binary_file_format.r")
source("./R_functions/function_closest2d.r")
source("./R_functions/read_src_model_priors.r")
source("./R_functions/nos_days_in_year.r")
source("./R_functions/plotconfidence.r")

# Set quantiles we want to consider for uncertainty 
# Here assuming a total of 7 quantiles (2.5 %, 5 %, 25 %, 50 %, 75 %, 95 %, 97.5 %)
mid_quant = 4 ; low_quant = 2 ; high_quant = 6
# Set output directory
out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/LTSS_CARBON_INTEGRATION/figures_africa_one_vs_all/"
#out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/LTSS_CARBON_INTEGRATION/figures_africa_nbe_vs_plusRaGPPhighConf/"
#out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/cssp_brazil_2/figures_productivity_without_vs_with/"
outsuffix = "_original_vs_alternate"

# Assign the baseline analysis - the original
# Original AGB assimilated (2007)
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/ODA_extension_Africa_one_AGB/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/ODA_extension_Africa_nbe/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/NoRainfor_woody_productivity_mortality/infofile.RData")
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))
orig_PROJECT = PROJECT ; orig_grid_parameters = grid_parameters ; orig_grid_output = grid_output
#orig_name = "Baseline"
orig_name = "One AGB" # used in labelling figures
#orig_name = "NBE" # used in labelling figures
# Assign the alternate analysis - the new data constraint
# Repeat AGB assimilated (2007-2010)
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/ODA_extension_Africa_actualCI_agb/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/ODA_extension_Africa_nbe_RaGPPhighConf/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_productivity_mortality/infofile.RData")
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))
alt_PROJECT = PROJECT ; alt_grid_parameters = grid_parameters ; alt_grid_output = grid_output
alt_name = "Repeat" # used in labelling figures
#alt_name = "+Productivity" # used in labelling figures
#alt_name = "+RaGPPhighConf" # used in labelling figures

# Tidy
rm(PROJECT,grid_parameters,grid_output)

###
## Determine needed spatial and temporal information

# generate the lat / long grid again
output = generate_wgs84_grid(orig_PROJECT$latitude,orig_PROJECT$longitude,orig_PROJECT$resolution)
grid_lat = array(output$lat, dim=c(orig_PROJECT$long_dim,orig_PROJECT$lat_dim))
grid_long = array(output$long,dim=c(orig_PROJECT$long_dim,orig_PROJECT$lat_dim))
cardamom_ext = output$cardamom_ext
# then generate the area estimates for each pixel (m)
area = calc_pixel_area(output$lat,output$long,orig_PROJECT$resolution)
# this output is in vector form and we need matching array shapes so...
area = array(area, dim=c(orig_PROJECT$long_dim,orig_PROJECT$lat_dim))

# Useful information
run_years = as.numeric(orig_PROJECT$start_year) : as.numeric(orig_PROJECT$end_year)
nos_years = length(as.numeric(orig_PROJECT$start_year) : as.numeric(orig_PROJECT$end_year))
steps_per_year = length(orig_PROJECT$model$timestep_days) / nos_years

# Load any map information you might want to overlay, such as biome or national boundaries
# if you don't want to have any set as false

# PointsOfChange

# Load a baseline land mask
# load global shape file for land sea mask
landmask = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
# just to be sure enforce the projection to WGS-84
landmask = spTransform(landmask,crs(cardamom_ext))
# subset by continent (could also do by country)
#landmask = subset(landmask, CONTINENT == "South America") # Change continent to target area or comment out if spanning zones
#landmask = subset(landmask, CONTINENT == "Africa") # Change continent to target area or comment out if spanning zones
# Clip to the extent of the CARDAMOM analysis
landmask = crop(landmask, cardamom_ext)

#add_biomes = " "
#add_biomes = "ssa_wwf"
add_biomes = "wwf_ecoregions"
if (add_biomes == "ssa_wwf") {
    # Read in shape file for boundaries
    biomes = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/SECO/analysis/ssa_wwf_dissolve/ssa_wwf_dissolve.shp")
    # just to be sure enforce the projection to WGS-84
    biomes = spTransform(biomes, crs(cardamom_ext))
    # Clip to the extent of the CARDAMOM analysis
    biomes = crop(biomes, cardamom_ext)
    # Clip to the area covered bythe landmask in case of geographic restriction
    biomes = crop(biomes, landmask)    
    # Overwrite the existing landmask
    landmask = biomes
    # Extract the current biome names for raster code
    biome_names = levels(factor(biomes$ECO_NAME))
    # Turn biomes object into a raster now in case we want to use it later
    biomes = rasterize(biomes, cardamom_ext, factor(biomes$ECO_NAME, fun="last"))
} else if (add_biomes == "wwf_ecoregions") {
    # Read in shape file for boundaries
    biomes = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/LTSS_CARBON_INTEGRATION/wwf_biomes_map/Ecoregions2017.shp")
    # just to be sure enforce the projection to WGS-84
    biomes = spTransform(biomes, crs(cardamom_ext))
    # Clip to the extent of the CARDAMOM analysis
    biomes = crop(biomes, cardamom_ext)
    # Clip to the area covered bythe landmask in case of geographic restriction
    biomes = crop(biomes, landmask)
    # Aggregate to the biome regions rather than the more complex ecoregions
    biomes <- aggregate(biomes, by='BIOME_NAME', dissolve = TRUE)

    # Overwrite the existing landmask
    landmask = biomes

    # Extract the current biome names for raster code
    biome_names = levels(factor(biomes$BIOME_NAME))
    # Turn the Biomes object into a raster now
    biomes = rasterize(biomes,cardamom_ext, factor(biomes$BIOME_NAME), fun="last")
} else {
    # DO NOTHING
    biomes = NA
} 


# This will be used to filter the analysis to include specific locations only
use_filter = TRUE
if (use_filter) {
    #  Design a user created / loaded filter 
    landfilter = array(NA,dim=dim(orig_grid_parameters$obs_wood_gCm2))
    landfilter[which(orig_grid_parameters$obs_wood_gCm2 > 0)] = 1
} else { 
    # Use this option if you don't want to filter
    landfilter = array(1,dim=dim(orig_grid_parameters$obs_wood_gCm2)) 
}

# Update the land filter with the mask information
landfilter = raster(vals = t(landfilter[,dim(landfilter)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
landfilter = mask(landfilter, landmask, updatevalue = NA)
# Reconstruct back into an array
landfilter = (array(as.vector(landfilter), dim=c(dim(orig_grid_parameters$obs_wood_gCm2)[1],dim(orig_grid_parameters$obs_wood_gCm2)[2])))
landfilter = landfilter[,dim(landfilter)[2]:1]

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

###
## Aggregate information needed for calibration data comparison

## Extract gridded information on the observations
dims = dim(orig_grid_output$mean_lai_m2m2) ; timestep_dim = length(orig_PROJECT$model$timestep_days)
## Soil prior
SoilCPrior = array(NA, dim=c(dims[1], dims[2]))
## Mean annual LAI obs
LAIobs = array(NA, dim=c(dims[1],dims[2],nos_years))
## Disturbance
HarvestFraction = array(NA, dim=c(dims[1], dims[2]))
BurnedFraction = array(NA, dim=c(dims[1], dims[2]))
FireFreq = array(NA, dim=c(dims[1],dims[2]))
## Observed wood trends information
WoodCobs = array(NA, dim=c(dims[1], dims[2],timestep_dim))
WoodCobs_CI = array(NA, dim=c(dims[1], dims[2],timestep_dim))
WoodCobs_trend_map = array(NA, dim=c(dims[1], dims[2]))
WoodCobs_trend = rep(NA, orig_PROJECT$nosites)
mean_obs_wood = rep(NA, orig_PROJECT$nosites)
WoodCobs_mean_CI = rep(0, orig_PROJECT$nosites)
## Modelled wood trends information
# Original
orig_WoodC = array(NA, dim=c(dims[1], dims[2],timestep_dim))
orig_wood_trend_map = array(NA, dim=c(dims[1], dims[2]))
orig_wood_trend = rep(NA, orig_PROJECT$nosites)
orig_mean_wood = rep(NA, orig_PROJECT$nosites)
# Alternate
alt_WoodC = array(NA, dim=c(dims[1], dims[2],timestep_dim))
alt_wood_trend_map = array(NA, dim=c(dims[1], dims[2]))
alt_wood_trend = rep(NA, orig_PROJECT$nosites)
alt_mean_wood = rep(NA, orig_PROJECT$nosites)
## Initialise stock and flux variables for aggregate time series
cumarea = 0 
# Original
orig_lai_grid = array(NA,dim=c(dims[1],dims[2],nos_years))
orig_lai_m2m2 = rep(0,nos_years)      ; orig_lai_lower_m2m2 = rep(0,nos_years)      ; orig_lai_upper_m2m2 = rep(0,nos_years)
orig_gpp_TgCyr = rep(0,nos_years)     ; orig_gpp_lower_TgCyr = rep(0,nos_years)     ; orig_gpp_upper_TgCyr = rep(0,nos_years)
orig_rauto_TgCyr = rep(0,nos_years)   ; orig_rauto_lower_TgCyr = rep(0,nos_years)   ; orig_rauto_upper_TgCyr = rep(0,nos_years)
orig_rhet_TgCyr = rep(0,nos_years)    ; orig_rhet_lower_TgCyr = rep(0,nos_years)    ; orig_rhet_upper_TgCyr = rep(0,nos_years)
orig_nee_TgCyr = rep(0,nos_years)     ; orig_nee_lower_TgCyr = rep(0,nos_years)     ; orig_nee_upper_TgCyr = rep(0,nos_years)
orig_nbe_TgCyr = rep(0,nos_years)     ; orig_nbe_lower_TgCyr = rep(0,nos_years)     ; orig_nbe_upper_TgCyr = rep(0,nos_years)
orig_fire_TgCyr = rep(0,nos_years)    ; orig_fire_lower_TgCyr = rep(0,nos_years)    ; orig_fire_upper_TgCyr = rep(0,nos_years)
orig_harvest_TgCyr = rep(0,nos_years) ; orig_harvest_lower_TgCyr = rep(0,nos_years) ; orig_harvest_upper_TgCyr = rep(0,nos_years)
orig_wood_TgC = rep(0,nos_years)      ; orig_wood_lower_TgC = rep(0,nos_years)      ; orig_wood_upper_TgC = rep(0,nos_years)
orig_lit_TgC = rep(0,nos_years)       ; orig_lit_lower_TgC = rep(0,nos_years)       ; orig_lit_upper_TgC = rep(0,nos_years)
orig_litwood_TgC = rep(0,nos_years)   ; orig_litwood_lower_TgC = rep(0,nos_years)   ; orig_litwood_upper_TgC = rep(0,nos_years)
orig_soil_TgC = rep(0,nos_years)      ; orig_soil_lower_TgC = rep(0,nos_years)      ; orig_soil_upper_TgC = rep(0,nos_years)
# Alternate
alt_lai_grid = array(NA,dim=c(dims[1],dims[2],nos_years))
alt_lai_m2m2 = rep(0,nos_years)      ; alt_lai_lower_m2m2 = rep(0,nos_years)      ; alt_lai_upper_m2m2 = rep(0,nos_years)
alt_gpp_TgCyr = rep(0,nos_years)     ; alt_gpp_lower_TgCyr = rep(0,nos_years)     ; alt_gpp_upper_TgCyr = rep(0,nos_years)
alt_rauto_TgCyr = rep(0,nos_years)   ; alt_rauto_lower_TgCyr = rep(0,nos_years)   ; alt_rauto_upper_TgCyr = rep(0,nos_years)
alt_rhet_TgCyr = rep(0,nos_years)    ; alt_rhet_lower_TgCyr = rep(0,nos_years)    ; alt_rhet_upper_TgCyr = rep(0,nos_years)
alt_nee_TgCyr = rep(0,nos_years)     ; alt_nee_lower_TgCyr = rep(0,nos_years)     ; alt_nee_upper_TgCyr = rep(0,nos_years)
alt_nbe_TgCyr = rep(0,nos_years)     ; alt_nbe_lower_TgCyr = rep(0,nos_years)     ; alt_nbe_upper_TgCyr = rep(0,nos_years)
alt_fire_TgCyr = rep(0,nos_years)    ; alt_fire_lower_TgCyr = rep(0,nos_years)    ; alt_fire_upper_TgCyr = rep(0,nos_years)
alt_harvest_TgCyr = rep(0,nos_years) ; alt_harvest_lower_TgCyr = rep(0,nos_years) ; alt_harvest_upper_TgCyr = rep(0,nos_years)
alt_wood_TgC = rep(0,nos_years)      ; alt_wood_lower_TgC = rep(0,nos_years)      ; alt_wood_upper_TgC = rep(0,nos_years)
alt_lit_TgC = rep(0,nos_years)       ; alt_lit_lower_TgC = rep(0,nos_years)       ; alt_lit_upper_TgC = rep(0,nos_years)
alt_litwood_TgC = rep(0,nos_years)   ; alt_litwood_lower_TgC = rep(0,nos_years)   ; alt_litwood_upper_TgC = rep(0,nos_years)
alt_soil_TgC = rep(0,nos_years)      ; alt_soil_lower_TgC = rep(0,nos_years)      ; alt_soil_upper_TgC = rep(0,nos_years)
## Initialise trend variables for flux trends
# Original
orig_gpp_trend = array(NA, dim=c(dims[1],dims[2]))
orig_rauto_trend = array(NA, dim=c(dims[1],dims[2]))
orig_rhet_trend = array(NA, dim=c(dims[1],dims[2]))
orig_lai_trend = array(NA, dim=c(dims[1],dims[2]))
# Alternate
alt_gpp_trend = array(NA, dim=c(dims[1],dims[2]))
alt_rauto_trend = array(NA, dim=c(dims[1],dims[2]))
alt_rhet_trend = array(NA, dim=c(dims[1],dims[2]))
alt_lai_trend = array(NA, dim=c(dims[1],dims[2]))
# Timing variable needed
time_vector = seq(0,nos_years, length.out = dim(orig_grid_output$nee_gCm2day)[3])

# Loop through all sites
nos_sites_inc = 0
# Loop through every site
for (n in seq(1, orig_PROJECT$nosites)) {
     
     # Check that location has run
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(alt_grid_output$i_location[n]) == FALSE  & is.na(alt_grid_output$j_location[n]) == FALSE ) {
         
         # Check that the location is always witin a location we want
         if (is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {

             # Extract each sites location within the grid
             i_loc = orig_grid_output$i_location[n] ; j_loc = orig_grid_output$j_location[n]
     
             ###
             ## Aggregating data
         
             ## Keep track of the number of site we are aggregating over
             nos_sites_inc = nos_sites_inc + 1
             # Cumulate the total area actually used in the analysis
             cumarea = cumarea + area[i_loc,j_loc]
             ## Estimate pixel level trends
             # Original
             orig_gpp_trend[i_loc,j_loc]   = coef(lm(orig_grid_output$gpp_gCm2day[n,mid_quant,] ~ time_vector))[2]     # median selected
             orig_rauto_trend[i_loc,j_loc] = coef(lm(orig_grid_output$rauto_gCm2day[n,mid_quant,] ~ time_vector))[2] # median selected
             orig_rhet_trend[i_loc,j_loc]  = coef(lm(orig_grid_output$rhet_gCm2day[n,mid_quant,] ~ time_vector))[2]   # median selected
             orig_lai_trend[i_loc,j_loc]   = coef(lm(orig_grid_output$lai_m2m2[n,mid_quant,] ~ time_vector))[2]   # median selected
             # Alternate
             alt_gpp_trend[i_loc,j_loc]   = coef(lm(alt_grid_output$gpp_gCm2day[n,mid_quant,] ~ time_vector))[2]     # median selected
             alt_rauto_trend[i_loc,j_loc] = coef(lm(alt_grid_output$rauto_gCm2day[n,mid_quant,] ~ time_vector))[2] # median selected
             alt_rhet_trend[i_loc,j_loc]  = coef(lm(alt_grid_output$rhet_gCm2day[n,mid_quant,] ~ time_vector))[2]   # median selected
             alt_lai_trend[i_loc,j_loc]   = coef(lm(alt_grid_output$lai_m2m2[n,mid_quant,] ~ time_vector))[2]   # median selected
             ## Extract time series aggregate information
             # Original
             orig_lai_grid[i_loc,j_loc,] = rollapply(orig_grid_output$lai_m2m2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
             orig_lai_m2m2               = orig_lai_m2m2            + orig_lai_grid[i_loc,j_loc,]
             orig_lai_lower_m2m2         = orig_lai_lower_m2m2      + rollapply(orig_grid_output$lai_m2m2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
             orig_lai_upper_m2m2         = orig_lai_upper_m2m2      + rollapply(orig_grid_output$lai_m2m2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
             orig_wood_TgC               = orig_wood_TgC            + rollapply(orig_grid_output$wood_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             orig_wood_lower_TgC         = orig_wood_lower_TgC      + rollapply(orig_grid_output$wood_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             orig_wood_upper_TgC         = orig_wood_upper_TgC      + rollapply(orig_grid_output$wood_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             orig_lit_TgC                = orig_lit_TgC             + rollapply(orig_grid_output$lit_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             orig_lit_lower_TgC          = orig_lit_lower_TgC       + rollapply(orig_grid_output$lit_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             orig_lit_upper_TgC          = orig_lit_upper_TgC       + rollapply(orig_grid_output$lit_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             if (length(which(names(orig_grid_output) == "litwood_gCm2")) > 0) {
                 orig_litwood_TgC        = orig_litwood_TgC         + rollapply(orig_grid_output$litwood_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
                 orig_litwood_lower_TgC  = orig_litwood_lower_TgC   + rollapply(orig_grid_output$litwood_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
                 orig_litwood_upper_TgC  = orig_litwood_upper_TgC   + rollapply(orig_grid_output$litwood_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             }
             orig_soil_TgC               = orig_soil_TgC            + rollapply(orig_grid_output$som_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             orig_soil_lower_TgC         = orig_soil_lower_TgC      + rollapply(orig_grid_output$som_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             orig_soil_upper_TgC         = orig_soil_upper_TgC      + rollapply(orig_grid_output$som_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             orig_gpp_TgCyr              = orig_gpp_TgCyr           + (rollapply(orig_grid_output$gpp_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_gpp_lower_TgCyr        = orig_gpp_lower_TgCyr     + (rollapply(orig_grid_output$gpp_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_gpp_upper_TgCyr        = orig_gpp_upper_TgCyr     + (rollapply(orig_grid_output$gpp_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_rauto_TgCyr            = orig_rauto_TgCyr         + (rollapply(orig_grid_output$rauto_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_rauto_lower_TgCyr      = orig_rauto_lower_TgCyr   + (rollapply(orig_grid_output$rauto_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_rauto_upper_TgCyr      = orig_rauto_upper_TgCyr   + (rollapply(orig_grid_output$rauto_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)         
             orig_rhet_TgCyr             = orig_rhet_TgCyr          + (rollapply(orig_grid_output$rhet_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_rhet_lower_TgCyr       = orig_rhet_lower_TgCyr    + (rollapply(orig_grid_output$rhet_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_rhet_upper_TgCyr       = orig_rhet_upper_TgCyr    + (rollapply(orig_grid_output$rhet_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)   
             orig_nee_TgCyr              = orig_nee_TgCyr           + (rollapply(orig_grid_output$nee_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)   
             orig_nee_lower_TgCyr        = orig_nee_lower_TgCyr     + (rollapply(orig_grid_output$nee_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)   
             orig_nee_upper_TgCyr        = orig_nee_upper_TgCyr     + (rollapply(orig_grid_output$nee_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_nbe_TgCyr              = orig_nbe_TgCyr           + (rollapply(orig_grid_output$nbe_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_nbe_lower_TgCyr        = orig_nbe_lower_TgCyr     + (rollapply(orig_grid_output$nbe_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_nbe_upper_TgCyr        = orig_nbe_upper_TgCyr     + (rollapply(orig_grid_output$nbe_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_fire_TgCyr             = orig_fire_TgCyr          + (rollapply(orig_grid_output$fire_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_fire_lower_TgCyr       = orig_fire_lower_TgCyr    + (rollapply(orig_grid_output$fire_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_fire_upper_TgCyr       = orig_fire_upper_TgCyr    + (rollapply(orig_grid_output$fire_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_harvest_TgCyr          = orig_harvest_TgCyr       + (rollapply(orig_grid_output$harvest_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_harvest_lower_TgCyr    = orig_harvest_lower_TgCyr + (rollapply(orig_grid_output$harvest_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             orig_harvest_upper_TgCyr    = orig_harvest_upper_TgCyr + (rollapply(orig_grid_output$harvest_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             # Alternate
             alt_lai_grid[i_loc,j_loc,] = rollapply(alt_grid_output$lai_m2m2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
             alt_lai_m2m2               = alt_lai_m2m2            + alt_lai_grid[i_loc,j_loc,]
             alt_lai_lower_m2m2         = alt_lai_lower_m2m2      + rollapply(alt_grid_output$lai_m2m2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
             alt_lai_upper_m2m2         = alt_lai_upper_m2m2      + rollapply(alt_grid_output$lai_m2m2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
             alt_wood_TgC               = alt_wood_TgC            + rollapply(alt_grid_output$wood_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             alt_wood_lower_TgC         = alt_wood_lower_TgC      + rollapply(alt_grid_output$wood_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             alt_wood_upper_TgC         = alt_wood_upper_TgC      + rollapply(alt_grid_output$wood_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             alt_lit_TgC                = alt_lit_TgC             + rollapply(alt_grid_output$lit_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             alt_lit_lower_TgC          = alt_lit_lower_TgC       + rollapply(alt_grid_output$lit_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             alt_lit_upper_TgC          = alt_lit_upper_TgC       + rollapply(alt_grid_output$lit_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             if (length(which(names(alt_grid_output) == "litwood_gCm2")) > 0) {
                 alt_litwood_TgC        = alt_litwood_TgC         + rollapply(alt_grid_output$litwood_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
                 alt_litwood_lower_TgC  = alt_litwood_lower_TgC   + rollapply(alt_grid_output$litwood_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
                 alt_litwood_upper_TgC  = alt_litwood_upper_TgC   + rollapply(alt_grid_output$litwood_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             }
             alt_soil_TgC               = alt_soil_TgC            + rollapply(alt_grid_output$som_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             alt_soil_lower_TgC         = alt_soil_lower_TgC      + rollapply(alt_grid_output$som_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             alt_soil_upper_TgC         = alt_soil_upper_TgC      + rollapply(alt_grid_output$som_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             alt_gpp_TgCyr              = alt_gpp_TgCyr           + (rollapply(alt_grid_output$gpp_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_gpp_lower_TgCyr        = alt_gpp_lower_TgCyr     + (rollapply(alt_grid_output$gpp_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_gpp_upper_TgCyr        = alt_gpp_upper_TgCyr     + (rollapply(alt_grid_output$gpp_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_rauto_TgCyr            = alt_rauto_TgCyr         + (rollapply(alt_grid_output$rauto_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_rauto_lower_TgCyr      = alt_rauto_lower_TgCyr   + (rollapply(alt_grid_output$rauto_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_rauto_upper_TgCyr      = alt_rauto_upper_TgCyr   + (rollapply(alt_grid_output$rauto_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)         
             alt_rhet_TgCyr             = alt_rhet_TgCyr          + (rollapply(alt_grid_output$rhet_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_rhet_lower_TgCyr       = alt_rhet_lower_TgCyr    + (rollapply(alt_grid_output$rhet_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_rhet_upper_TgCyr       = alt_rhet_upper_TgCyr    + (rollapply(alt_grid_output$rhet_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)   
             alt_nee_TgCyr              = alt_nee_TgCyr           + (rollapply(alt_grid_output$nee_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)   
             alt_nee_lower_TgCyr        = alt_nee_lower_TgCyr     + (rollapply(alt_grid_output$nee_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)   
             alt_nee_upper_TgCyr        = alt_nee_upper_TgCyr     + (rollapply(alt_grid_output$nee_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_nbe_TgCyr              = alt_nbe_TgCyr           + (rollapply(alt_grid_output$nbe_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_nbe_lower_TgCyr        = alt_nbe_lower_TgCyr     + (rollapply(alt_grid_output$nbe_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_nbe_upper_TgCyr        = alt_nbe_upper_TgCyr     + (rollapply(alt_grid_output$nbe_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_fire_TgCyr             = alt_fire_TgCyr          + (rollapply(alt_grid_output$fire_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_fire_lower_TgCyr       = alt_fire_lower_TgCyr    + (rollapply(alt_grid_output$fire_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_fire_upper_TgCyr       = alt_fire_upper_TgCyr    + (rollapply(alt_grid_output$fire_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_harvest_TgCyr          = alt_harvest_TgCyr       + (rollapply(alt_grid_output$harvest_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_harvest_lower_TgCyr    = alt_harvest_lower_TgCyr + (rollapply(alt_grid_output$harvest_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
             alt_harvest_upper_TgCyr    = alt_harvest_upper_TgCyr + (rollapply(alt_grid_output$harvest_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
          
             ###
             ## Determining Drivers and trends   

             # Read in pixel driving data
             drivers = read_binary_file_format(paste(alt_PROJECT$datapath,alt_PROJECT$name,"_",alt_PROJECT$sites[n],".bin",sep=""))
             # Determine forest harvest intensity
             HarvestFraction[i_loc,j_loc] = sum(drivers$met[,8]) / nos_years
             # Determine mean annual fire intensity and frequency
             BurnedFraction[i_loc,j_loc] = sum(drivers$met[,9]) / nos_years
             FireFreq[i_loc,j_loc] = length(which(drivers$met[,9] > 0)) / nos_years
             # Load any priors
             SoilCPrior[i_loc,j_loc] = drivers$parpriors[23]
             # Clear missing data from and extract observed LAI
             drivers$obs[which(drivers$obs[,3] == -9999),3] = NA
             LAIobs[i_loc,j_loc,] = rollapply(drivers$obs[,3], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
             # If wood stock estimate available get that too
             tmp = which(drivers$obs[,13] > 0)
             if (length(tmp) > 0) {
                 for (t in seq(1, length(tmp))) {
                      # Observational constraint
                      WoodCobs[i_loc,j_loc,tmp[t]] = drivers$obs[tmp[t],13]
                      WoodCobs_CI[i_loc,j_loc,tmp[t]] = drivers$obs[tmp[t],14]
                      # Corresponding model output
                      orig_WoodC[i_loc,j_loc,tmp[t]] = orig_grid_output$wood_gCm2[n,mid_quant,tmp[t]]
                      alt_WoodC[i_loc,j_loc,tmp[t]] = alt_grid_output$wood_gCm2[n,mid_quant,tmp[t]]
                 } # loop time steps with obs
                 # Observed wood stock trends
                 obs_period_start = tmp[1] ; obs_period_end = tmp[length(tmp)] ; obs_period_years = length(c(obs_period_start:obs_period_end))      
                 WoodCobs_trend[n] = (coef(lm(WoodCobs[i_loc,j_loc,obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12) # *12 is month to yr adjustment
                 WoodCobs_trend_map[i_loc,j_loc] = WoodCobs_trend[n]
                 mean_obs_wood[n] = mean(WoodCobs[i_loc,j_loc,], na.rm=TRUE)
                 WoodCobs_mean_CI[n] = mean(drivers$obs[tmp,14], na.rm=TRUE)
                 # Original model wood stock trends
                 orig_wood_trend[n] = (coef(lm(orig_grid_output$wood_gCm2[n,mid_quant,obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12)
                 orig_wood_trend_map[i_loc,j_loc] = orig_wood_trend[n]
                 orig_mean_wood[n] = mean(orig_grid_output$mean_wood_gCm2[i_loc,j_loc,mid_quant])
                 # Alternate model wood stock trends
                 alt_wood_trend[n] = (coef(lm(alt_grid_output$wood_gCm2[n,mid_quant,obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12)
                 alt_wood_trend_map[i_loc,j_loc] = alt_wood_trend[n]
                 alt_mean_wood[n] = mean(alt_grid_output$mean_wood_gCm2[i_loc,j_loc,mid_quant])
             } # we have more than zero obs
         } # Is the location within the landfilter?
     } # Did this location run?
} # Site loop

## Original model output first
# LAI averaging
orig_lai_m2m2 = orig_lai_m2m2 / nos_sites_inc
orig_lai_lower_m2m2 = orig_lai_lower_m2m2 / nos_sites_inc
orig_lai_upper_m2m2 = orig_lai_upper_m2m2 / nos_sites_inc
# Now adjust units gC/yr -> TgC/yr
# Medians
orig_gpp_TgCyr     = orig_gpp_TgCyr * 1e-12
orig_rauto_TgCyr   = orig_rauto_TgCyr * 1e-12
orig_rhet_TgCyr    = orig_rhet_TgCyr * 1e-12
orig_nee_TgCyr     = orig_nee_TgCyr * 1e-12
orig_nbe_TgCyr     = orig_nbe_TgCyr * 1e-12
orig_fire_TgCyr    = orig_fire_TgCyr * 1e-12
orig_harvest_TgCyr = orig_harvest_TgCyr * 1e-12
orig_lit_TgC       = orig_lit_TgC * 1e-12
orig_litwood_TgC   = orig_litwood_TgC * 1e-12
orig_wood_TgC      = orig_wood_TgC * 1e-12
orig_soil_TgC      = orig_soil_TgC * 1e-12
# lower
orig_gpp_lower_TgCyr     = orig_gpp_lower_TgCyr * 1e-12
orig_rauto_lower_TgCyr   = orig_rauto_lower_TgCyr * 1e-12
orig_rhet_lower_TgCyr    = orig_rhet_lower_TgCyr * 1e-12
orig_nee_lower_TgCyr     = orig_nee_lower_TgCyr * 1e-12
orig_nbe_lower_TgCyr     = orig_nbe_lower_TgCyr * 1e-12
orig_fire_lower_TgCyr    = orig_fire_lower_TgCyr * 1e-12
orig_harvest_lower_TgCyr = orig_harvest_lower_TgCyr * 1e-12
orig_lit_lower_TgC       = orig_lit_lower_TgC * 1e-12
orig_litwood_lower_TgC   = orig_litwood_lower_TgC * 1e-12
orig_wood_lower_TgC      = orig_wood_lower_TgC * 1e-12
orig_soil_lower_TgC      = orig_soil_lower_TgC * 1e-12
# upper
orig_gpp_upper_TgCyr     = orig_gpp_upper_TgCyr * 1e-12
orig_rauto_upper_TgCyr   = orig_rauto_upper_TgCyr * 1e-12
orig_rhet_upper_TgCyr    = orig_rhet_upper_TgCyr * 1e-12
orig_nee_upper_TgCyr     = orig_nee_upper_TgCyr * 1e-12
orig_nbe_upper_TgCyr     = orig_nbe_upper_TgCyr * 1e-12
orig_fire_upper_TgCyr    = orig_fire_upper_TgCyr * 1e-12
orig_harvest_upper_TgCyr = orig_harvest_upper_TgCyr * 1e-12
orig_lit_upper_TgC       = orig_lit_upper_TgC * 1e-12
orig_litwood_upper_TgC   = orig_litwood_upper_TgC * 1e-12
orig_wood_upper_TgC      = orig_wood_upper_TgC * 1e-12
orig_soil_upper_TgC      = orig_soil_upper_TgC * 1e-12
## Alternate model output first
# LAI averaging
alt_lai_m2m2 = alt_lai_m2m2 / nos_sites_inc
alt_lai_lower_m2m2 = alt_lai_lower_m2m2 / nos_sites_inc
alt_lai_upper_m2m2 = alt_lai_upper_m2m2 / nos_sites_inc
# Now adjust units gC/yr -> TgC/yr
# Medians
alt_gpp_TgCyr     = alt_gpp_TgCyr * 1e-12
alt_rauto_TgCyr   = alt_rauto_TgCyr * 1e-12
alt_rhet_TgCyr    = alt_rhet_TgCyr * 1e-12
alt_nee_TgCyr     = alt_nee_TgCyr * 1e-12
alt_nbe_TgCyr     = alt_nbe_TgCyr * 1e-12
alt_fire_TgCyr    = alt_fire_TgCyr * 1e-12
alt_harvest_TgCyr = alt_harvest_TgCyr * 1e-12
alt_lit_TgC       = alt_lit_TgC * 1e-12
alt_litwood_TgC   = alt_litwood_TgC * 1e-12
alt_wood_TgC      = alt_wood_TgC * 1e-12
alt_soil_TgC      = alt_soil_TgC * 1e-12
# lower
alt_gpp_lower_TgCyr     = alt_gpp_lower_TgCyr * 1e-12
alt_rauto_lower_TgCyr   = alt_rauto_lower_TgCyr * 1e-12
alt_rhet_lower_TgCyr    = alt_rhet_lower_TgCyr * 1e-12
alt_nee_lower_TgCyr     = alt_nee_lower_TgCyr * 1e-12
alt_nbe_lower_TgCyr     = alt_nbe_lower_TgCyr * 1e-12
alt_fire_lower_TgCyr    = alt_fire_lower_TgCyr * 1e-12
alt_harvest_lower_TgCyr = alt_harvest_lower_TgCyr * 1e-12
alt_lit_lower_TgC       = alt_lit_lower_TgC * 1e-12
alt_litwood_lower_TgC   = alt_litwood_lower_TgC * 1e-12
alt_wood_lower_TgC      = alt_wood_lower_TgC * 1e-12
alt_soil_lower_TgC      = alt_soil_lower_TgC * 1e-12
# upper
alt_gpp_upper_TgCyr     = alt_gpp_upper_TgCyr * 1e-12
alt_rauto_upper_TgCyr   = alt_rauto_upper_TgCyr * 1e-12
alt_rhet_upper_TgCyr    = alt_rhet_upper_TgCyr * 1e-12
alt_nee_upper_TgCyr     = alt_nee_upper_TgCyr * 1e-12
alt_nbe_upper_TgCyr     = alt_nbe_upper_TgCyr * 1e-12
alt_fire_upper_TgCyr    = alt_fire_upper_TgCyr * 1e-12
alt_harvest_upper_TgCyr = alt_harvest_upper_TgCyr * 1e-12
alt_lit_upper_TgC       = alt_lit_upper_TgC * 1e-12
alt_litwood_upper_TgC   = alt_litwood_upper_TgC * 1e-12
alt_wood_upper_TgC      = alt_wood_upper_TgC * 1e-12
alt_soil_upper_TgC      = alt_soil_upper_TgC * 1e-12

###
## Plot overlayed PDFs for different data constraints and grid varying parameters and trait maps
###

### Maps of absolute harvest emission

# Then create new colours with high 'alpha', i.e. transparency
#c_colours = colorRampPalette(brewer.pal(8,"Accent"))
#c_colours = c_colours(8)
c_colours = c("blue","green")
nbins = 30 # desired number of catagories, you might not get this many

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_parameter_PDFs_by",outsuffix,".png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(6,6), mar = c(2,2,2,1))
# Loop parameters
for (p in seq(1, dim(orig_grid_parameters$parameters)[3]-1)) {
     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$parameters[,,p,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$parameters[,,p,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=paste("Parameter = ",p,sep=""), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with original
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) # Add next alternate
     if (p == 1) {legend("topleft",legend = c(orig_name, alt_name), col = c(c_colours[1],c_colours[2]), lty=1, lwd=1.6, bty = "n", cex=0.8)}
} # loop parameters
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_log_parameter_PDFs_by",outsuffix,".png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(6,6), mar = c(2,2,2,1))
# Loop parameters
for (p in seq(1, dim(orig_grid_parameters$parameters)[3]-1)) {
     # Set to local variables
     tmp1 = log(as.vector(orig_grid_parameters$parameters[,,p,mid_quant]))
     tmp2 = log(as.vector(alt_grid_parameters$parameters[,,p,mid_quant]))
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=paste("log(Parameter) = ",p,sep=""), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with original
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) # Add next alternate
     if (p == 1) {legend("topleft",legend = c(orig_name, alt_name), col = c(c_colours[1],c_colours[2]), lty=1, lwd=1.6, bty = "n", cex=0.8)}
} # loop parameters
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_traitv2_PDFs_by",outsuffix,".png",sep=""), width = 3000, height = 1800, res = 300)
par(mfrow=c(3,4), mar=c(2,2,2,1), omi=c(0.1,0.1,0.14,0.1))

     ## CUE

     # Set to local variables
     tmp1 = 1-as.vector(orig_grid_output$mean_rauto_gCm2day[,,mid_quant] / orig_grid_output$mean_gpp_gCm2day[,,mid_quant])
     tmp2 = 1-as.vector(alt_grid_output$mean_rauto_gCm2day[,,mid_quant] / alt_grid_output$mean_gpp_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("CUE (0-1)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) 
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     legend("topleft",legend = c(orig_name, alt_name), col = c(c_colours[1],c_colours[2]), lty=1, lwd=1.6, bty = "n", cex=0.8)

     ## NPP MgC/ha/yr to foliage

     # Set to local variables
     tmp1 = as.vector(orig_grid_output$mean_fnpp_gCm2day[,,mid_quant])*365.25*1e-2
     tmp2 = as.vector(alt_grid_output$mean_fnpp_gCm2day[,,mid_quant])*365.25*1e-2
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     xrange = range(x_axis) ; xrange[1] = min(0,min(xrange))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("NP",P[foliar]," (MgC h",a^-1,y^-1,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", xlim=xrange, ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 

     ## NPP MgC/ha/yr to fine roots

     # Set to local variables
     tmp1 = as.vector(orig_grid_output$mean_rnpp_gCm2day[,,mid_quant])*365.25*1e-2
     tmp2 = as.vector(alt_grid_output$mean_rnpp_gCm2day[,,mid_quant])*365.25*1e-2
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     xrange = range(x_axis) ; xrange[1] = min(0,min(xrange))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("NP",P[root]," (MgC h",a^-1,y^-1,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", xlim=xrange, ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## NPP MgC/ha/yr to wood

     # Set to local variables
     tmp1 = as.vector(orig_grid_output$mean_wnpp_gCm2day[,,mid_quant])*365.25*1e-2
     tmp2 = as.vector(alt_grid_output$mean_wnpp_gCm2day[,,mid_quant])*365.25*1e-2
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     xrange = range(x_axis) ; xrange[1] = min(0,min(xrange))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("NP",P[wood]," (MgC h",a^-1,y^-1,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", xlim=xrange, ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT foliage (years)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_foliar_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_foliar_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[foliar]," (y)",sep="")),
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT fine root (years)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_root_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_foliar_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[root]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
          
     ## MRT wood (years)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_wood_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_wood_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     # Now plot each of them
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[wood]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT DeadOrg (litter + CWD; if applicable)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_DeadOrg_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_DeadOrg_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[DeadOrg]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT soil

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_som_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_som_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[som]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## Leaf Carbon per unit leaf Area (gC/m2)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$parameters[,,17,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$parameters[,,17,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("LCA (gC",m^-2,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 

     ## Canopy Photosynthetic Efficiency (gC/m2/day)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$parameters[,,11,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$parameters[,,11,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("Ceff (gC",m^-2,d^-1,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 

dev.off()

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_trait_PDFs_by",outsuffix,".png",sep=""), width = 3000, height = 1800, res = 300)
par(mfrow=c(3,4), mar=c(2,2,2,1), omi=c(0.1,0.1,0.14,0.1))

     ## CUE

     # Set to local variables
     tmp1 = 1-as.vector(orig_grid_output$mean_rauto_gCm2day[,,mid_quant] / orig_grid_output$mean_gpp_gCm2day[,,mid_quant])
     tmp2 = 1-as.vector(alt_grid_output$mean_rauto_gCm2day[,,mid_quant] / alt_grid_output$mean_gpp_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("CUE (0-1)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) 
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     legend("topleft",legend = c(orig_name, alt_name), col = c(c_colours[1],c_colours[2]), lty=1, lwd=1.6, bty = "n", cex=0.8)

     ## NPP fraction to foliage

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$NPP_foliar_fraction[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$NPP_foliar_fraction[,,mid_quant])
     tmp1[which(tmp1 > 1)] = NA ; tmp2[which(tmp2 > 1)] = NA # prevents against precision error in codes not picking up on very small fl allocations but turn into large fractional ones
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     xrange = c(0,1)
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("NP",P[foliar]," (0-1)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", xlim=xrange, ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 

     ## NPP fraction to fine roots

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$NPP_root_fraction[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$NPP_root_fraction[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     xrange = c(0,1)
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("NP",P[root]," (0-1)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", xlim=xrange, ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## NPP fraction to wood

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$NPP_wood_fraction[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$NPP_wood_fraction[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     xrange = c(0,1)
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("NP",P[wood]," (0-1)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", xlim=xrange, ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT foliage (years)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_foliar_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_foliar_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[foliar]," (y)",sep="")),
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT fine root (years)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_root_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_foliar_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[root]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
          
     ## MRT wood (years)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_wood_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_wood_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     # Now plot each of them
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[wood]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT DeadOrg (litter + CWD; if applicable)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_DeadOrg_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_DeadOrg_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[DeadOrg]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## MRT soil

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$MTT_som_years[,,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$MTT_som_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("MR",T[som]," (y)",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 
     
     ## Leaf Carbon per unit leaf Area (gC/m2)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$parameters[,,17,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$parameters[,,17,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("LCA (gC",m^-2,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 

     ## Canopy Photosynthetic Efficiency (gC/m2/day)

     # Set to local variables
     tmp1 = as.vector(orig_grid_parameters$parameters[,,11,mid_quant])
     tmp2 = as.vector(alt_grid_parameters$parameters[,,11,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp1,tmp2), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp1,tmp2), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     tmp1 <- hist(tmp1, breaks = ax, plot = FALSE)
     tmp2 <- hist(tmp2, breaks = ax, plot = FALSE)
     x_axis = tmp1$mids
     tmp1 = tmp1$counts / sum(tmp1$counts) 
     tmp2 = tmp2$counts / sum(tmp2$counts)
     # Now plot them together
     ymax = max(c(tmp1,tmp2))
     plot(tmp1~x_axis, type="l", lwd=2, col = c_colours[1], main=expression(paste("Ceff (gC",m^-2,d^-1,")",sep="")), 
          xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     lines(tmp2~x_axis, col = c_colours[2], lwd=2) 

dev.off()

###
## Determine 1-posterior:prior ratio, i.e. how much have we learned?
###

# Extract parameter prior ranges from source code
prior_ranges = read_src_model_priors(orig_PROJECT)

# Create ratio array
orig_posterior_prior = array(NA, dim=c(dim(orig_grid_parameters$parameters)[1:2],length(prior_ranges$parmin)))
alt_posterior_prior = array(NA, dim=c(dim(orig_grid_parameters$parameters)[1:2],length(prior_ranges$parmin)))
for (n in seq(1, orig_PROJECT$nosites)) {

     # Check that location has run
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(alt_grid_output$i_location[n]) == FALSE & is.na(alt_grid_output$j_location[n]) == FALSE) {
         if (is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
             for (p in seq(1,length(prior_ranges$parmin))) {
                  # Original
                  tmp = orig_grid_parameters$parameters[orig_grid_output$i_location[n],orig_grid_output$j_location[n],p,high_quant] 
                  tmp = tmp - orig_grid_parameters$parameters[orig_grid_output$i_location[n],orig_grid_output$j_location[n],p,low_quant] 
                  orig_posterior_prior[orig_grid_output$i_location[n],orig_grid_output$j_location[n],p] = tmp / (prior_ranges$parmax[p]-prior_ranges$parmin[p])
                  # Alternate
                  tmp = alt_grid_parameters$parameters[orig_grid_output$i_location[n],orig_grid_output$j_location[n],p,high_quant] 
                  tmp = tmp - alt_grid_parameters$parameters[orig_grid_output$i_location[n],orig_grid_output$j_location[n],p,low_quant] 
                  alt_posterior_prior[orig_grid_output$i_location[n],orig_grid_output$j_location[n],p] = tmp / (prior_ranges$parmax[p]-prior_ranges$parmin[p])
             } # Loop parameters
         } # Inclued in land filter?
     } # Does site exist
} # Loop sites

# Generate some summary statistics
print("===Original all parameters===")
print(summary(apply(1-orig_posterior_prior, 3, mean, na.rm=TRUE)))
print("===Alternate all parameters===")
print(summary(apply(1-alt_posterior_prior, 3, mean, na.rm=TRUE)))
print("===Alternate-Original parameters===")
for (p in seq(1, dim(alt_posterior_prior)[3])) {
tmp1 = mean((1-alt_posterior_prior[,,p]), na.rm=TRUE) 
tmp2 = mean((1-orig_posterior_prior[,,p]), na.rm=TRUE)
     print(paste("p",p," mean = ",round(tmp1-tmp2,digits=3)," (",round(tmp1,digits=3),"-",round(tmp2,digits=3),")", sep=""))
}

# Generate some plots

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_posterior_prior_reductions",outsuffix,".png",sep=""), height = 700, width = 3000, res = 300)
par(mfrow=c(1,3), mar=c(0.01,1.5,0.2,7),omi=c(0.01,0.1,0.01,0.1))
var1 = apply(1-orig_posterior_prior,c(1,2),mean,na.rm=TRUE)
var1 = raster(vals = t(var1[,dim(var1)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = apply(1-alt_posterior_prior,c(1,2),mean,na.rm=TRUE)
var2 = raster(vals = t(var2[,dim(var2)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = var2-var1
plot(var1, zlim=c(0,1), col=colour_choices_gain, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main="")
mtext(orig_name, side = 3, cex = 1.2, padj = 1.3)
plot(landmask, add=TRUE)
mtext(expression('Mean posterior reduction (0-1)'), side = 2, cex = 0.9, padj = 0.0, adj = 0.5)
plot(var2, zlim=c(0,1), col=colour_choices_gain, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main="")
mtext(alt_name, side = 3, cex = 1.2, padj = 1.3)
plot(landmask, add=TRUE)
xrange = c(-1,1) * max(abs(range(values(var3), na.rm=TRUE)), na.rm=TRUE)
plot(var3, main="", zlim=xrange, col=colour_choices_sign, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side = 3, cex = 1.2, padj = 1.3)
dev.off()

###
## Loading and processing of independent observations
###

# The overall scheme here is to load multiple datasets for each observation type (NBE, GPP, Fire)
# combine them together to provide a mean estimate and an estimate of uncertainty based on the maximum and minimum values for each case

###
## Extract CarbonTracker Europe Inversion (NEE, NBE)

# Extract CarbonTrackerEurope (2000-2017)
CTE = nc_open("/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/CarbonTrackerEurope/flux1x1_all_years.nc")
cte_nee = ncvar_get(CTE,"bio_flux_opt")  # molC/m2/s (long,lat,ensemble,date) Fire not included
cte_fire = ncvar_get(CTE,"fire_flux_imp") # molC/m2/s 
cte_area = ncvar_get(CTE,"cell_area")    # m2
cte_days_since_2000 = ncvar_get(CTE,"date")
cte_lat = ncvar_get(CTE,"latitude")
cte_long = ncvar_get(CTE,"longitude")
nc_close(CTE)

# Convert cte_nee into nbe
cte_nbe = cte_nee + cte_fire

# Adjust units
cte_nee = cte_nee * 12 * 86400 * 365.25 # gC/m2/yr
cte_fire = cte_fire * 12 * 86400 * 365.25 # gC/m2/yr
cte_nbe = cte_nbe * 12 * 86400 * 365.25 # gC/m2/yr

# Search for africa locations and slot into africa only grid for matching
# Make into CARDAMOM paired masks.
# NOTE: filter 2000-2017 and 0.025,0.25,0.5,0.75,0.975 quantiles
cte_years = c(2000:2017) 
overlap_cte = intersect(cte_years,run_years)
overlap_start = which(cte_years == overlap_cte[1])
overlap_end = which(cte_years == overlap_cte[length(overlap_cte)])
# Create data arrays now
cte_nbe_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(overlap_cte)))
cte_nee_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(overlap_cte)))
cte_fire_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(overlap_cte)))
cte_m2 = array(NA, dim=dim(orig_grid_output$mean_lai_m2m2)[1:2])
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,cte_lat,cte_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         cte_nbe_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = cte_nbe[i1,j1,overlap_start:overlap_end]
         cte_nee_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = cte_nee[i1,j1,overlap_start:overlap_end]
         cte_fire_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = cte_fire[i1,j1,overlap_start:overlap_end]
         cte_m2[orig_grid_output$i_location[n],orig_grid_output$j_location[n]] = cte_area[i1,j1] 
     }
}

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,cte_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = cte_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - cte_years[length(cte_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(cte_nbe_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        cte_nbe_gCm2yr = abind(add_beginning,cte_nbe_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(cte_nbe_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        cte_nbe_gCm2yr = abind(cte_nbe_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Calculate domain wide mean net sink for overlap period
cte_domain_nbe_gCm2yr = mean(apply(cte_nbe_gCm2yr,c(1,2),mean, na.rm=TRUE), na.rm=TRUE) # gC/m2/yr
cte_domain_nbe_TgCyr =   sum(apply(cte_nbe_gCm2yr,c(1,2),mean, na.rm=TRUE) * cte_m2 * 1e-12, na.rm=TRUE) # TgC/yr

# Estimate the long term trend
func_lm <-function(var_in) { 
  tmp = length(which(is.na(var_in) == FALSE)) 
  if (tmp > 1) {
      return(coef(lm(var_in ~ c(1:length(var_in))))[2])
  } else {
      return(NA)
  } 
}
cte_domain_nbe_gCm2yr_trend = apply(cte_nbe_gCm2yr,c(1,2),func_lm) # gC/m2/yr

# Tidy up
rm(cte_nee,cte_fire,cte_nbe,cte_lat,cte_long)

###
## Extract CarbonTracker Ensembles
## Two ensembles estimates exist using either flask and oco2 observations

# Make a list of all available biosphere flux estimates
avail_files = list.files("/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/CarbonTrackerEurope/Gerbrand_ensemble/biofluxopt/", full.names=TRUE)

## First the flask based analysis (2009-2017)

# Select the analyses using flask data
flask_files = avail_files[grepl("flask", avail_files) == TRUE]
# Divide between the biosphere flux and fire flux
# \\ means don't consider . as a wildcard, $ means at the end of the string
flask_files_nee = flask_files[grepl("biofluxopt\\.nc$", flask_files) == TRUE]
flask_files_fire = flask_files[grepl("firefluximp\\.nc$", flask_files) == TRUE]
# Restrict the biosphere fluxes to those which we have the iposed fire emissions, 
# so that we can convert nee into nbe
check_fire_version = unlist(strsplit(x = flask_files_fire, split = "firefluximp.nc"))
tmp = 0
for (i in seq(1,length(check_fire_version))) {
     tmpp = which(grepl(check_fire_version[i],flask_files_nee))
     if (length(tmp) > 0) {
         tmp = append(tmp,tmpp)
     }
}
flask_files_nee = flask_files_nee[tmp]
rm(tmp,tmpp,check_fire_version)
if (length(flask_files_nee) != length(flask_files_fire)) {stop("Oh dear there is a problem with the number of fire and biosphere flux files for flask")}

# Read in the first NEE file to extract spatial and temporal information
flask = nc_open(flask_files_nee[1])
flask_lat = ncvar_get(flask, "latitude")
flask_long = ncvar_get(flask, "longitude")
flask_date = ncvar_get(flask, "date") 
# Read the first values (mol/m2/s)
flask_nee = ncvar_get(flask,"bio_flux_opt")
# Tidy
nc_close(flask)

# Read in the first Fire file to extract spatial information
flask = nc_open(flask_files_fire[1])
# Read the first values (mol/m2/s)
flask_fire = ncvar_get(flask,"fire_flux_imp")
# Tidy
nc_close(flask)

# Estimate step size
flask_step = abs(flask_date[1]-flask_date[2])
# Estimate the number of days in each year since 2000 (the reference point for)
create_years = c(0,2000:2020)
nos_days = 0 ; for (i in seq(2,length(create_years))) { nos_days = append(nos_days,nos_days_in_year(create_years[i])) }
# Convert all into decimal year
for (i in seq(1, length(flask_date))) {
     tmp = which(cumsum(nos_days) > flask_date[i])[1]
     flask_date[i] = create_years[tmp] + ((flask_date[i]-sum(nos_days[1:(tmp-1)]))/nos_days[tmp])
}
# Determine where we will clip the datasets in time
flask_years = floor(flask_date) ; flask_years_keep = 0
for (i in seq(1, length(unique(flask_years)))) {
     if (length(which(flask_years == unique(flask_years)[i])) > 48) {
         flask_years_keep = append(flask_years_keep, which(flask_years == unique(flask_years)[i]))
     } 
}
flask_years_keep = flask_years_keep[-1]
# Select time periods we want only
flask_years = flask_years[flask_years_keep]
overlap_flask = intersect(flask_years,run_years)
flask_nee = flask_nee[,,flask_years_keep]
flask_fire = flask_fire[,,flask_years_keep]

# Restructure to hold all ensemble members
flask_nee = array(flask_nee, dim=c(dim(flask_nee)[1:3],length(flask_files)))
flask_fire = array(flask_fire, dim=c(dim(flask_nee)[1:3],length(flask_files)))
# Loop through all files to get our full ensemble
for (i in seq(2, length(flask_files_nee))) {
     # Open the new file
     flask_bio = nc_open(flask_files_nee[i])
     flask_fir = nc_open(flask_files_fire[i])
     # Read NEE and fire (mol/m2/s)
     tmp_nee = ncvar_get(flask_bio, "bio_flux_opt")
     tmp_fire = ncvar_get(flask_fir, "fire_flux_imp")

     # Close file
     nc_close(flask_bio) ; nc_close(flask_fir)
     # Trim to desired time period
     flask_nee[,,,i] = tmp_nee[,,flask_years_keep]     
     flask_fire[,,,i] = tmp_fire[,,flask_years_keep]     
}

# Now apply units correction (mol/m2/s) to gC/m2/day
flask_nee = flask_nee * 12 * 86400
flask_fire = flask_fire * 12 * 86400
# Create NBE
flask_nbe = flask_nee + flask_fire

# Loop through each year to estimate the annual means
flask_nee_gCm2yr = array(NA, dim=c(dim(flask_nee)[1:2],length(unique(flask_years)),dim(flask_nee)[4]))
flask_nbe_gCm2yr = array(NA, dim=c(dim(flask_nee)[1:2],length(unique(flask_years)),dim(flask_nee)[4]))
flask_fire_gCm2yr = array(NA, dim=c(dim(flask_nee)[1:2],length(unique(flask_years)),dim(flask_nee)[4]))
for (i in seq(1, length(unique(flask_years)))) {
     tmp = which(flask_years == unique(flask_years)[i])
     # Average across each year and scale to annual total
     flask_nee_gCm2yr[,,i,] = apply(flask_nee[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     flask_nbe_gCm2yr[,,i,] = apply(flask_nbe[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     flask_fire_gCm2yr[,,i,] = apply(flask_fire[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
}
# Remove the existing output as not needed now
rm(flask_nee,flask_fire,flask_nbe)

# Update flask years to their annuals only
flask_years = unique(flask_years)

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
flask_cardamom_nee_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(flask_years),dim(flask_fire_gCm2yr)[4]))
flask_cardamom_nbe_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(flask_years),dim(flask_fire_gCm2yr)[4]))
flask_cardamom_fire_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(flask_years),dim(flask_fire_gCm2yr)[4]))
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,flask_lat,flask_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         flask_cardamom_nee_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],,] = flask_nee_gCm2yr[i1,j1,,]
         flask_cardamom_nbe_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],,] = flask_nbe_gCm2yr[i1,j1,,]
         flask_cardamom_fire_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],,] = flask_fire_gCm2yr[i1,j1,,]
     }
}

# Tidy up
rm(flask_nee_gCm2yr,flask_nbe_gCm2yr,flask_fire_gCm2yr)

## Second the OCO-2 based analysis (2015-2017)

# Select the analyses using oco2 data
oco2_files = avail_files[grepl("oco2", avail_files) == TRUE]
# Divide between the biosphere flux and fire flux
# \\ means don't consider . as a wildcard, $ means at the end of the string
oco2_files_nee = oco2_files[grepl("biofluxopt\\.nc$", oco2_files) == TRUE]
oco2_files_fire = oco2_files[grepl("firefluximp\\.nc$", oco2_files) == TRUE]
# Restrict the biosphere fluxes to those which we have the imposed fire emissions, 
# so that we can convert nee into nbe
check_fire_version = unlist(strsplit(x = oco2_files_fire, split = "firefluximp.nc"))
tmp = 0
for (i in seq(1,length(check_fire_version))) {
     tmpp = which(grepl(check_fire_version[i],oco2_files_nee))
     if (length(tmp) > 0) {
         tmp = append(tmp,tmpp)
     }
}
oco2_files_nee = oco2_files_nee[tmp]
rm(tmp,tmpp,check_fire_version)
if (length(oco2_files_nee) != length(oco2_files_fire)) { stop("Oh dear there is a problem with the number of fire and biosphere flux files for oco2") }

# Read in the first file to extract spatial and temporal information
oco2 = nc_open(oco2_files[1])
oco2_lat = ncvar_get(oco2, "latitude")
oco2_long = ncvar_get(oco2, "longitude")
oco2_date = ncvar_get(oco2, "date") 
# Read the first values (mol/m2/s)
oco2_nee = ncvar_get(oco2,"bio_flux_opt")
# Tidy
nc_close(oco2)

# Read in the first Fire file to extract spatial and temporal information
oco2 = nc_open(oco2_files_fire[1])
# Read the first values (mol/m2/s)
oco2_fire = ncvar_get(oco2,"fire_flux_imp")
# Tidy
nc_close(oco2)

# Estimate step size
oco2_step = abs(oco2_date[1]-oco2_date[2])
# Estimate the number of days in each year since 2000 (the reference point for )
create_years = c(0,2000:2020)
nos_days = 0 ; for (i in seq(2,length(create_years))) { nos_days = append(nos_days,nos_days_in_year(create_years[i]))}
# Convert all into decimal year
for (i in seq(1, length(oco2_date))) {
     tmp = which(cumsum(nos_days) > oco2_date[i])[1]
     oco2_date[i] = create_years[tmp] + ((oco2_date[i]-sum(nos_days[1:(tmp-1)]))/nos_days[tmp])
}
# Determine where we will clip the datasets in time
oco2_years = floor(oco2_date) ; oco2_years_keep = 0
for (i in seq(1, length(unique(oco2_years)))) {
     if (length(which(oco2_years == unique(oco2_years)[i])) > 48) {
         oco2_years_keep = append(oco2_years_keep, which(oco2_years == unique(oco2_years)[i]))
     } 
}
oco2_years_keep = oco2_years_keep[-1]
# Select time periods we want only
oco2_years = oco2_years[oco2_years_keep]
oco2_nee = oco2_nee[,,oco2_years_keep]
oco2_fire = oco2_fire[,,oco2_years_keep]
overlap_oco2 = intersect(oco2_years,run_years)

# Restructure to hold all ensemble members
oco2_nee = array(oco2_nee, dim=c(dim(oco2_nee)[1:3],length(oco2_files)))
oco2_fire = array(oco2_fire, dim=c(dim(oco2_nee)[1:3],length(oco2_files)))
# Loop through all files to get our full ensemble
for (i in seq(2, length(oco2_files_nee))) {
     # Open the new file
     oco2_bio = nc_open(oco2_files_nee[i])
     oco2_fir = nc_open(oco2_files_fire[i])
     # Read NEE (mol/m2/s)
     tmp_nee = ncvar_get(oco2_bio, "bio_flux_opt")
     tmp_fire = ncvar_get(oco2_fir, "fire_flux_imp")
     # Close file
     nc_close(oco2_bio) ; nc_close(oco2_fir)
     # Trim to desired time period
     oco2_nee[,,,i] = tmp_nee[,,oco2_years_keep]     
     oco2_fire[,,,i] = tmp_fire[,,oco2_years_keep]     
}

# Now apply units correction (mol/m2/s) to gC/m2/day
oco2_nee = oco2_nee * 12 * 86400
oco2_fire = oco2_fire * 12 * 86400
oco2_nbe = oco2_nee + oco2_fire

# Loop through each year to estimate the annual means
oco2_nee_gCm2yr = array(NA, dim=c(dim(oco2_nee)[1:2],length(unique(oco2_years)),dim(oco2_nee)[4]))
oco2_fire_gCm2yr = array(NA, dim=c(dim(oco2_nee)[1:2],length(unique(oco2_years)),dim(oco2_nee)[4]))
oco2_nbe_gCm2yr = array(NA, dim=c(dim(oco2_nee)[1:2],length(unique(oco2_years)),dim(oco2_nee)[4]))
for (i in seq(1, length(unique(oco2_years)))) {
     tmp = which(oco2_years == unique(oco2_years)[i])
     # Average across each year and scale to annual total
     oco2_nee_gCm2yr[,,i,] = apply(oco2_nee[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     oco2_fire_gCm2yr[,,i,] = apply(oco2_fire[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     oco2_nbe_gCm2yr[,,i,] = apply(oco2_nbe[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
}
# Remove the existing output as not needed now
rm(oco2_nee,oco2_fire,oco2_nbe)
# Now update oco2 years to their annuals only
oco2_years = unique(oco2_years)

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
oco2_cardamom_nee_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(oco2_years),dim(oco2_fire_gCm2yr)[4]))
oco2_cardamom_nbe_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(oco2_years),dim(oco2_fire_gCm2yr)[4]))
oco2_cardamom_fire_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(oco2_years),dim(oco2_fire_gCm2yr)[4]))
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,oco2_lat,oco2_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         oco2_cardamom_nee_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],,] = oco2_nee_gCm2yr[i1,j1,,]
         oco2_cardamom_nbe_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],,] = oco2_nbe_gCm2yr[i1,j1,,]
         oco2_cardamom_fire_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],,] = oco2_fire_gCm2yr[i1,j1,,]
     } # valid value exists
} # loop sites

## Combine the NBE estimates from our datasets

# How many unique years in total
obs_nbe_years = unique(c(flask_years,oco2_years))

# Define the combined timeseries datasets
obs_nbe_gCm2yr = array(NA, dim = c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(obs_nbe_years),sum(c(dim(oco2_cardamom_nbe_gCm2yr)[4],dim(flask_cardamom_nbe_gCm2yr)[4]))))
obs_nee_gCm2yr = array(NA, dim = c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(obs_nbe_years),sum(c(dim(oco2_cardamom_nbe_gCm2yr)[4],dim(flask_cardamom_nbe_gCm2yr)[4]))))
obs_fire_gCm2yr = array(NA, dim = c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(obs_nbe_years),sum(c(dim(oco2_cardamom_nbe_gCm2yr)[4],dim(flask_cardamom_nbe_gCm2yr)[4]))))
for (i in seq(1,length(obs_nbe_years))) {
     # determine whether the flask dataset has any values for this year
     tmp = which(flask_years == obs_nbe_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1 ; i_e = dim(flask_cardamom_nee_gCm2yr)[4]
         obs_nbe_gCm2yr[,,i,i_s:i_e] = flask_cardamom_nbe_gCm2yr[,,tmp,]
         obs_nee_gCm2yr[,,i,i_s:i_e] = flask_cardamom_nee_gCm2yr[,,tmp,]
         obs_fire_gCm2yr[,,i,i_s:i_e] = flask_cardamom_fire_gCm2yr[,,tmp,]
     }
     # determine whether the oco2 dataset has any values for this year
     tmp = which(oco2_years == obs_nbe_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+dim(flask_cardamom_nee_gCm2yr)[4] ; i_e = i_s - 1 + dim(oco2_cardamom_nee_gCm2yr)[4]
         obs_nbe_gCm2yr[,,i,i_s:i_e] = oco2_cardamom_nbe_gCm2yr[,,tmp,]
         obs_nee_gCm2yr[,,i,i_s:i_e] = oco2_cardamom_nee_gCm2yr[,,tmp,]
         obs_fire_gCm2yr[,,i,i_s:i_e] = oco2_cardamom_fire_gCm2yr[,,tmp,]
     }
} # loop years

# Extract the time step mean / min / max for each of these fluxes now
obs_nbe_mean_gCm2yr = apply(obs_nbe_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_nbe_min_gCm2yr = apply(obs_nbe_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_nbe_max_gCm2yr = apply(obs_nbe_gCm2yr,c(1,2,3),max, na.rm=TRUE)
obs_nee_mean_gCm2yr = apply(obs_nee_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_nee_min_gCm2yr = apply(obs_nee_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_nee_max_gCm2yr = apply(obs_nee_gCm2yr,c(1,2,3),max, na.rm=TRUE)
obs_fire_mean_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_fire_min_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_fire_max_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),max, na.rm=TRUE)
# Filter out the Inf values to NaN
obs_nbe_min_gCm2yr[is.infinite(obs_nbe_min_gCm2yr) == TRUE] = NA
obs_nbe_max_gCm2yr[is.infinite(obs_nbe_max_gCm2yr) == TRUE] = NA
obs_nee_min_gCm2yr[is.infinite(obs_nee_min_gCm2yr) == TRUE] = NA
obs_nee_max_gCm2yr[is.infinite(obs_nee_max_gCm2yr) == TRUE] = NA
obs_fire_min_gCm2yr[is.infinite(obs_fire_min_gCm2yr) == TRUE] = NA
obs_fire_max_gCm2yr[is.infinite(obs_fire_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,obs_nbe_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = obs_nbe_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - obs_nbe_years[length(obs_nbe_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_nbe_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_nbe_mean_gCm2yr  = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_nbe_min_gCm2yr   = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_nbe_max_gCm2yr   = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
        obs_nee_mean_gCm2yr  = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_nee_min_gCm2yr   = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_nee_max_gCm2yr   = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
        obs_fire_mean_gCm2yr = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_fire_min_gCm2yr  = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_fire_max_gCm2yr  = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_nbe_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_nbe_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_nbe_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_nbe_max_gCm2yr = abind(obs_nbe_max_gCm2yr,add_afterward, along=3)
        obs_nee_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_nee_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_nee_max_gCm2yr = abind(obs_nbe_max_gCm2yr, along=3)
        obs_fire_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_fire_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_fire_max_gCm2yr = abind(obs_nbe_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Generate aggregate values at the domain level
obs_nbe_mean_domain_TgCyr = apply(obs_nbe_mean_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_nbe_min_domain_TgCyr = apply(obs_nbe_min_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_nbe_max_domain_TgCyr = apply(obs_nbe_max_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

# tidy 
rm(oco2_lat,oco2_long,flask_years,oco2_years,obs_nbe_years,
   oco2_cardamom_nbe_gCm2yr,oco2_cardamom_nee_gCm2yr,oco2_cardamom_fire_gCm2yr,
   flask_cardamom_nbe_gCm2yr,flask_cardamom_nee_gCm2yr,flask_cardamom_fire_gCm2yr)

###
## Extract GPP estimates from Copernicus, FLUXCOM & MODIS

## Extract FLUXCOM GPP CRUJRA (2000-2017)

# Read first file to get additional information
fc_years = c(2000:2017)
fc_years = intersect(fc_years,run_years)

# Read in the first file to extract spatial and temporal information
FC = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/FLUXCOM/CRUJRA/GPP.RS_METEO.FP-ALL.MLM-ALL.METEO-CRUJRA_v1.720_360.monthly.",fc_years[1],".nc",sep=""))
fc_gpp = ncvar_get(FC,"GPP")  # gC/m2/d (long,lat,date) Fire not included
fc_lat = ncvar_get(FC,"lat")
fc_long = ncvar_get(FC,"lon")
nc_close(FC)

# Calculate mean annual from the monthly data
fc_gpp = apply(fc_gpp,c(1,2),mean)

# Update dimensions to accomodate all years about to be read in
fc_gpp = array(fc_gpp, dim=c(dim(fc_gpp)[1:2],length(fc_years)))

# Loop years now reading in each file in turn and adding to the output variable
for (t in seq(2,length(fc_years))) {
     # Read file
     FC = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/FLUXCOM/CRUJRA/GPP.RS_METEO.FP-ALL.MLM-ALL.METEO-CRUJRA_v1.720_360.monthly.",fc_years[t],".nc",sep=""))
     tmp = ncvar_get(FC,"GPP") 
     # Monthly to annual
     fc_gpp[,,t] = apply(tmp,c(1,2),mean)
     # tidy
     nc_close(FC) ; rm(tmp)
}

# Adjust units
fc_gpp = fc_gpp * 365.25 # gC/m2/yr

# Aggregate from 0.5 x 0.5 deg to 1x1
#fc_gpp_tmp = array(NA,dim=c(360,180,length(fc_years)))
#fc_gpp_mad_tmp = array(NA,dim=c(360,180,length(fc_years)))
## Adjust lat / long vectors
#fc_lat = rollapply(fc_lat,by = 2, width = 2, mean)
#fc_long = rollapply(fc_long,by = 2, width = 2, mean)
#x = 1 ; y = 1
#for (i in seq(1,dim(fc_gpp)[1],2)){
#    for (j in seq(1,dim(fc_gpp)[2],2)){
#         fc_gpp_tmp[x,y,] = apply(fc_gpp[i:(i+1),j:(j+1),],3,mean)
#         fc_gpp_mad_tmp[x,y,] = apply(fc_gpp_mad[i:(i+1),j:(j+1),],3,mean)
#         y = y + 1
#    }
#    x = x + 1 ; y = 1
#}
## Once aggregation done replace the original variables and tidy up
#fc_gpp = fc_gpp_tmp ; fc_gpp_mad = fc_gpp_mad_tmp 
#rm(fc_gpp_tmp,fc_gpp_mad_tmp)

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
fc_cardamom_gpp_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(fc_years)))
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,fc_lat,fc_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         fc_cardamom_gpp_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = fc_gpp[i1,j1,]
     } # valid value exists
} # loop sites

## Extract Copernicus GPP (2000-2019)

# "PointsOfChange"

# Read first file to get additional information
copernicus_years = c(2000:2019)
copernicus_years = intersect(copernicus_years,run_years)

# Make a list of all available files
copernicus_files = list.files("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/copernius/global_1deg/", full.names = TRUE)

# Loop years now reading in each file in turn and adding to the output variable
for (y in seq(1,length(copernicus_years))) {
     # Make a list of all files for this year
     infiles = copernicus_files[which(grepl(paste("c_gls_GDMP_",copernicus_years[y],sep=""),copernicus_files) == TRUE)]
     # Loop through all files in the year to store these too
     for (t in seq(1, length(infiles))) {
          # Open file
          copernicus = nc_open(infiles[t])
          # If this is the first year extract some spatial information
          if (y == 1 & t == 1) {
              # Read lat / long
              copernicus_lat = ncvar_get(copernicus,"lat")
              copernicus_long = ncvar_get(copernicus,"lon")         
              # Create output variable we will accumulate into 
              copernicus_gpp = array(NA, dim=c(length(copernicus_long),length(copernicus_lat),length(copernicus_years)))
          } # first year actions only
          # For first step of the year create the new years loading variable
          if (t == 1) { tmp = array(NA, dim=c(length(copernicus_long),length(copernicus_lat),length(infiles))) }
          # Read in variable
          tmp[,,t] = ncvar_get(copernicus, "GDMP")
          # Remove any missing data flags
          tmp[,,t][which(tmp[,,t] == -9999)] = NA
          # Tidy up
          nc_close(copernicus)
     } # loop time step in year
     # Monthly to annual
     copernicus_gpp[,,y] = apply(tmp,c(1,2),mean)
     # Tidy
     rm(tmp)
} # Loop copernicus years

# Adjust units
copernicus_gpp = copernicus_gpp * 365.25 #gC/m2/day -> gC/m2/yr

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
copernicus_cardamom_gpp_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(copernicus_years)))
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,copernicus_lat,copernicus_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         copernicus_cardamom_gpp_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = copernicus_gpp[i1,j1,]
     } # valid value exists
} # loop sites


## Extract FluxSat v2 GPP (2000-2019)

# "PointsOfChange"

# Read first file to get additional information
fluxsat_years = c(2000:2019)
fluxsat_years = intersect(fluxsat_years,run_years)

# Make a list of all available files
fluxsat_files = list.files("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/FluxSat/global_1deg_monthly/", full.names = TRUE)

# Loop years now reading in each file in turn and adding to the output variable
for (y in seq(1,length(fluxsat_years))) {
     # Make a list of all files for this year
     infiles = fluxsat_files[which(grepl(paste("GPP_FluxSat_daily_v2_",fluxsat_years[y],sep=""),fluxsat_files) == TRUE)]
     # Loop through all files in the year to store these too
     for (t in seq(1, length(infiles))) {
          # Open file
          copernicus = nc_open(infiles[t])
          # If this is the first year extract some spatial information
          if (y == 1 & t == 1) {
              # Read lat / long
              fluxsat_lat = ncvar_get(copernicus,"lat")
              fluxsat_long = ncvar_get(copernicus,"lon")         
              # Create output variable we will accumulate into 
              fluxsat_gpp = array(NA, dim=c(length(fluxsat_long),length(fluxsat_lat),length(fluxsat_years)))
          } # first year actions only
          # For first step of the year create the new years loading variable
          if (t == 1) { tmp = array(NA, dim=c(length(fluxsat_long),length(fluxsat_lat),length(infiles))) }
          # Read in variable
          tmp[,,t] = ncvar_get(copernicus, "GPP")
          # Remove any missing data flags
          tmp[,,t][which(tmp[,,t] == -9999)] = NA
          # Tidy up
          nc_close(copernicus)
     } # loop time step in year
     # Monthly to annual
     fluxsat_gpp[,,y] = apply(tmp,c(1,2),mean)
     # Tidy
     rm(tmp)
} # Loop copernicus years

# Adjust units
fluxsat_gpp = fluxsat_gpp * 365.25 #gC/m2/day -> gC/m2/yr

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
fluxsat_cardamom_gpp_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(fluxsat_years)))
fluxsat_cardamom_gpp_gCm2yr_trend = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2]))
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE &
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,fluxsat_lat,fluxsat_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         fluxsat_cardamom_gpp_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = fluxsat_gpp[i1,j1,]
         # Estimate the GPP trend at this time too
         if (length(which(is.na(fluxsat_gpp[i1,j1,]) == FALSE)) > 1) {
             fluxsat_cardamom_gpp_gCm2yr_trend[orig_grid_output$i_location[n],orig_grid_output$j_location[n]] = coef(lm(fluxsat_gpp[i1,j1,]~fluxsat_years))[2]
         } 
     } # valid value exists
} # loop sites

## Combine the GPP estimates from our datasets

# How many unique years in total
obs_gpp_years = unique(c(fc_years,copernicus_years,fluxsat_years))

# Define the combined timeseries datasets
nos_gpp_databases = 3
obs_gpp_gCm2yr = array(NA, dim = c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(obs_gpp_years),nos_gpp_databases))
for (i in seq(1,length(obs_gpp_years))) {
     # determine whether the fluxcom dataset has any values for this year
     tmp = which(fc_years == obs_gpp_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1 ; i_e = 1
         obs_gpp_gCm2yr[,,i,i_s:i_e] = fc_cardamom_gpp_gCm2yr[,,tmp]
     }
     # determine whether the copernicus dataset has any values for this year
     tmp = which(copernicus_years == obs_gpp_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+1 ; i_e = i_s - 1 + 1
         obs_gpp_gCm2yr[,,i,i_s:i_e] = copernicus_cardamom_gpp_gCm2yr[,,tmp]
     }
     # determine whether the fluxsatv2 dataset has any values for this year
     tmp = which(fluxsat_years == obs_gpp_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+1+1 ; i_e = i_s - 1 + 1
         obs_gpp_gCm2yr[,,i,i_s:i_e] = fluxsat_cardamom_gpp_gCm2yr[,,tmp]
     }
} # loop years

# Extract the time step mean / min / max for each of these fluxes now
obs_gpp_mean_gCm2yr = apply(obs_gpp_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_gpp_min_gCm2yr = apply(obs_gpp_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_gpp_max_gCm2yr = apply(obs_gpp_gCm2yr,c(1,2,3),max, na.rm=TRUE)
# Filter out the Inf values to NaN
obs_gpp_min_gCm2yr[is.infinite(obs_gpp_min_gCm2yr) == TRUE] = NA
obs_gpp_max_gCm2yr[is.infinite(obs_gpp_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,obs_gpp_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = obs_gpp_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - obs_gpp_years[length(obs_gpp_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_gpp_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_gpp_mean_gCm2yr = abind(add_beginning,obs_gpp_mean_gCm2yr, along=3)
        obs_gpp_min_gCm2yr  = abind(add_beginning,obs_gpp_min_gCm2yr, along=3)
        obs_gpp_max_gCm2yr  = abind(add_beginning,obs_gpp_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_gpp_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_gpp_mean_gCm2yr = abind(obs_gpp_mean_gCm2yr,add_afterward, along=3)
        obs_gpp_min_gCm2yr  = abind(obs_gpp_min_gCm2yr,add_afterward, along=3)
        obs_gpp_max_gCm2yr  = abind(obs_gpp_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

apply(copernicus_cardamom_gpp_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(copernicus_cardamom_gpp_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

# Generate aggregate values at the domain level
obs_gpp_mean_domain_TgCyr = apply(obs_gpp_mean_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_gpp_min_domain_TgCyr = apply(obs_gpp_min_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_gpp_max_domain_TgCyr = apply(obs_gpp_max_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

###
## Independent fire emissions estimate
## Two estimates available, gfed and gfas
###

## Extract GFEDv4.1s (2001-2016)

# Read first file to get additional information
gfed_years = c(2001:2016)
gfed_years = intersect(gfed_years,run_years)

gfed = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GFED4/1.0deg/monthly_emissions/GFED4.1s_",gfed_years[1],".nc",sep=""))
gfed_fire = ncvar_get(gfed,"emissions")  # gC/m2/month (long,lat,date) 
gfed_lat = ncvar_get(gfed,"lat")
gfed_long = ncvar_get(gfed,"lon")
nc_close(gfed)

# Aggregate to annual flux
gfed_fire = apply(gfed_fire,c(1,2),sum,na.rm=TRUE)
# correct dimension
gfed_fire = array(gfed_fire, dim=c(dim(gfed_fire)[1:2],length(gfed_years)))
for (t in seq(2,length(gfed_years))) {
     gfed = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GFED4/1.0deg/monthly_emissions/GFED4.1s_",gfed_years[t],".nc",sep=""))
     tmp = ncvar_get(gfed,"emissions")  # gC/m2/month (long,lat,date) 
     tmp = apply(tmp,c(1,2),sum,na.rm=TRUE)
     gfed_fire[,,t] = tmp
     nc_close(gfed)
}
# Search for africa locations and slot into africa only grid for matching
# Make into CARDAMOM paired masks.
gfed_cardamom_fire_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(gfed_years)))
gfed_m2 = array(NA, dim=dim(orig_grid_output$mean_lai_m2m2)[1:2])
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,gfed_lat,gfed_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         gfed_cardamom_fire_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = gfed_fire[i1,j1,]
     }
}

## Annual Emissions: (units: Teragrams (Tg) carbon)
# 2004-2017. GFAS
# Based on MODIS radiative energy emissions and converted to biomass based on PFT conversion factors

# Time period
gfas_years = c(2004:2017)
gfas_years = intersect(gfas_years,run_years)
gfas = nc_open(paste("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/cssp_brazil/GFAS_annual/cams_gfas_co2fire_",gfas_years[1],".nc",sep=""))
gfas_fire = ncvar_get(gfas, "AnnualFire")
gfas_lat = ncvar_get(gfas,"latitude")
gfas_long = ncvar_get(gfas,"longitude")

# Make space for all years
gfas_fire = array(gfas_fire, dim=c(dim(gfas_fire)[1:2],length(gfas_years)))
# Read in all the files and collect into single array
for (t in seq(2, length(gfas_years))) {
     gfas = nc_open(paste("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/cssp_brazil/GFAS_annual/cams_gfas_co2fire_",gfas_years[t],".nc",sep=""))
     tmp = ncvar_get(gfas, "AnnualFire")
     gfas_fire[,,t] = tmp
}
# Tidy
nc_close(gfas)

# Search for africa locations and slot into africa only grid for matching onto CARDAMOM
gfas_cardamom_fire_gCm2yr = array(NA, dim=c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(gfas_years)))
for (n in seq(1,orig_PROJECT$nosites)) {
     if (is.na(orig_grid_output$i_location[n]) == FALSE & is.na(orig_grid_output$j_location[n]) == FALSE & 
         is.na(landfilter[orig_grid_output$i_location[n],orig_grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,gfas_lat,gfas_long,grid_lat[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],grid_long[orig_grid_output$i_location[n],orig_grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         gfas_cardamom_fire_gCm2yr[orig_grid_output$i_location[n],orig_grid_output$j_location[n],] = gfas_fire[i1,j1,]
     }
}

## Combine the Fire estimates from our datasets

# How many unique years in total
obs_fire_years = unique(c(gfed_years,gfas_years))

# Define the combined timeseries datasets
nos_fire_databases = 2
obs_fire_gCm2yr = array(NA, dim = c(dim(orig_grid_output$mean_lai_m2m2)[1:2],length(obs_fire_years),nos_fire_databases))
for (i in seq(1,length(obs_fire_years))) {
     # determine whether the flask dataset has any values for this year
     tmp = which(gfed_years == obs_fire_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1 ; i_e = 1
         obs_fire_gCm2yr[,,i,i_s:i_e] = gfed_cardamom_fire_gCm2yr[,,tmp]
     }
     # determine whether the oco2 dataset has any values for this year
     tmp = which(gfas_years == obs_fire_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+1 ; i_e = i_s - 1 + 1
         obs_fire_gCm2yr[,,i,i_s:i_e] = gfas_cardamom_fire_gCm2yr[,,tmp]
     }
} # loop years

# Extract the time step mean / min / max for each of these fluxes now
obs_fire_mean_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_fire_min_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_fire_max_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),max, na.rm=TRUE)
# Filter out the Inf values to NaN
obs_fire_min_gCm2yr[is.infinite(obs_fire_min_gCm2yr) == TRUE] = NA
obs_fire_max_gCm2yr[is.infinite(obs_fire_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,obs_fire_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = obs_fire_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - obs_fire_years[length(obs_fire_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_fire_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_fire_mean_gCm2yr = abind(add_beginning,obs_fire_mean_gCm2yr, along=3)
        obs_fire_min_gCm2yr = abind(add_beginning,obs_fire_min_gCm2yr, along=3)
        obs_fire_max_gCm2yr = abind(add_beginning,obs_fire_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_fire_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_fire_mean_gCm2yr = abind(obs_fire_mean_gCm2yr,add_afterward, along=3)
        obs_fire_min_gCm2yr = abind(obs_fire_min_gCm2yr,add_afterward, along=3)
        obs_fire_max_gCm2yr = abind(obs_fire_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Generate aggregate values at the domain level
obs_fire_mean_domain_TgCyr = apply(obs_fire_mean_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_fire_min_domain_TgCyr = apply(obs_fire_min_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_fire_max_domain_TgCyr = apply(obs_fire_max_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

# Tidy
rm(gfas_cardamom_fire_gCm2yr,gfas_long,gfas_lat,gfas_years,gfas_fire,
   gfed_cardamom_fire_gCm2yr,gfed_long,gfed_lat,gfed_years,gfed_fire)

###
## Plot Observations

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_observed_modelled_wood_timeseries",outsuffix,".png",sep=""), height = 2000, width = 3500, res = 300)
par(mfrow=c(1,1), mar=c(4.2,4.7,2.8,2),omi=c(0.01,0.01,0.01,0.01))
## Plot modelled\ median and CI timeseries with corresponding observation and uncertainty, if available
var1 = NA ; var2 = NA ; var3 = NA ; var4 = NA ; var5 = NA
# Modelled wood
var1 = orig_wood_TgC ; var2 = orig_wood_lower_TgC ; var3 = orig_wood_upper_TgC
var6 = alt_wood_TgC ; var7 = alt_wood_lower_TgC ; var8 = alt_wood_upper_TgC
# Observed wood
var4 = apply(WoodCobs*array(landfilter*area,dim=dim(WoodCobs)),3,sum,na.rm=TRUE) ; var4[which(var4 == 0)] = NA
var4 = rollapply(var4, FUN = mean, by = 12, width = 12, na.rm=TRUE)*1e-12
var5 = apply(WoodCobs_CI**2*array(landfilter*area,dim=dim(WoodCobs)),3,sum,na.rm=TRUE)
var5 = sqrt(rollapply(var5, FUN = mean, by = 12, width = 12, na.rm=TRUE)*1e-12)
# Begin plotting
zrange = range(c(var1,var2,var3,var4,var5,var6,var7,var8), na.rm=TRUE)*c(0.8,1.3)
plot(var1~run_years, ylim=zrange, ylab="", xlab="", pch=16, cex.lab=2.3, cex.axis = 2.2, cex=2, col="white")
mtext(expression(paste('Year',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Wood stocks (TgC)',sep="")), side = 2, cex = 2.4, padj = -1.3)
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) #; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var3~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[2], lwd=3, lty = 1) #; points(var6~run_years, col=model_colours[2], pch=16)
lines(var7~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var7~run_years, col=model_colours[2], pch=16)
lines(var8~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var8~run_years, col=model_colours[2], pch=16)
if (length(which(is.na(var4) == FALSE)) > 0) {
    plotCI(x = run_years, y = var4, uiw = var5, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
           col="black", lwd=4, ylab="", xlab="", add=TRUE)
}
legend("topleft", legend = c("EO estimate",orig_name,alt_name), col = c("black",model_colours[1],model_colours[2]), lty = c(1,1,1), pch=c(NA,NA,NA), 
        horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_wood_litter_soil_timeseries_comparison_plusCI",outsuffix,".png",sep=""), height = 3800, width = 2500, res = 300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot wood stocks
var1 = orig_wood_TgC ; var2 = orig_wood_lower_TgC ; var3 = orig_wood_upper_TgC
var4 = alt_wood_TgC ; var5 = alt_wood_lower_TgC ; var6 = alt_wood_upper_TgC
# Observed wood
var7 = apply(WoodCobs*array(landfilter*area,dim=dim(WoodCobs)),3,sum,na.rm=TRUE) ; var7[which(var7 == 0)] = NA
var7 = rollapply(var7, FUN = mean, by = 12, width = 12, na.rm=TRUE)*1e-12
var8 = apply(WoodCobs_CI**2*array(landfilter*area,dim=dim(WoodCobs)),3,sum,na.rm=TRUE)
var8 = sqrt(rollapply(var8, FUN = mean, by = 12, width = 12, na.rm=TRUE)*1e-12)
zrange = range(c(var1,var2,var3), na.rm=TRUE)*c(0.8,1.3)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) #; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[2], lwd=3, lty = 1) #; points(var4~run_years, col=model_colours[2], pch=16)
lines(var5~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var5~run_years, col=model_colours[2], pch=16)
lines(var6~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var6~run_years, col=model_colours[2], pch=16)
if (length(which(is.na(var7) == FALSE)) > 0) {
    plotCI(x = run_years, y = var7, uiw = var8, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
           col="black", lwd=4, ylab="", xlab="", add=TRUE)
}
mtext(expression(paste("Wood (TgC)",sep="")), side=2, padj=-2.05,cex=1.5)
legend("topleft", legend = c(orig_name,alt_name), col = c(model_colours[1:2]), 
       lty = rep(1,2), pch=rep(NA,2), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
# Now plot foliar and fine root litter stocks
var1 = orig_lit_TgC ; var2 = orig_lit_lower_TgC ; var3 = orig_lit_upper_TgC
var4 = alt_lit_TgC ; var5 = alt_lit_lower_TgC ; var6 = alt_lit_upper_TgC
zrange = range(c(var1,var2,var3), na.rm=TRUE)*c(0.8,1.2)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) #; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[2], lwd=3, lty = 1) #; points(var4~run_years, col=model_colours[2], pch=16)
lines(var5~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var5~run_years, col=model_colours[2], pch=16)
lines(var6~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var6~run_years, col=model_colours[2], pch=16)
mtext(expression(paste("Litter (TgC)",sep="")), side=2, padj=-2.05, cex=1.5)
# Now plot soil stocks
var1 = orig_soil_TgC ; var2 = orig_soil_lower_TgC ; var3 = orig_soil_upper_TgC
var4 = alt_soil_TgC ; var5 = alt_soil_lower_TgC ; var6 = alt_soil_upper_TgC
zrange = range(c(var1,var2,var3), na.rm=TRUE)*c(0.8,1.2)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) #; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[2], lwd=3, lty = 1) #; points(var4~run_years, col=model_colours[2], pch=16)
lines(var5~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var5~run_years, col=model_colours[2], pch=16)
lines(var6~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var6~run_years, col=model_colours[2], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Soil (TgC)",sep="")), side=2, padj=-2.05,cex=1.5)
dev.off()
           
# Compare analyses against the observational constraints (LAI, Soil C prior, Cwood stock, potAGB)
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_compare_observation",outsuffix,".png",sep=""), height = 4000, width = 4500, res = 300)
par(mfrow=c(2,2), mar=c(3,4.2,3,2), omi = c(0.35,0.4,0.1,0.1))
# Plot LAI mean annual
var1 = as.vector(LAIobs) ; var1 = var1[which(is.na(var1) == FALSE)]
var2 = as.vector(orig_lai_grid) ; var2 = var2[which(is.na(var2) == FALSE)] 
plot(var2 , var1, col=model_colours[1],
     pch=1, cex = 1.6, cex.lab=2.4, cex.axis = 2.4, cex.main=2.0, ylab="", xlab="", main="")
var2 = as.vector(alt_lai_grid) ; var2 = var2[which(is.na(var2) == FALSE)] 
points(var2 , var1, pch=1, col=model_colours[2], cex = 1.6)
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Annual LAI (',m^2,'/',m^2,')',sep="")), side = 2, cex = 2.4, padj = -1.05)
legend("topleft", legend = c(orig_name,alt_name), col = model_colours[1:2], lty = 2, pch=rep(16,2), horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
abline(0,1, col="grey", lwd=3)
# Plot wood
plot(as.vector(1e-2*WoodCobs) ~ as.vector(1e-2*orig_WoodC), pch=1, cex = 1.6, cex.lab=2.0, cex.axis = 2.3, cex.main=2.0, ylab="", xlab="", main="", col=model_colours[1])
points(as.vector(1e-2*WoodCobs) ~ as.vector(1e-2*alt_WoodC), pch=1, cex = 1.6, col=model_colours[2])
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Wood stocks (MgC h',a^-1,')',sep="")), side = 2, cex = 2.4, padj = -1.00)
abline(0,1, col="grey", lwd=3)
# Now plot LAI time series
var3  = apply(LAIobs,3,mean,na.rm=TRUE)
var4  = orig_lai_m2m2   
var5  = alt_lai_m2m2 
zrange = range(c(var3,var4,var5), na.rm=TRUE) * c(0.8,1.2)
plot(var3~run_years, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
      col="black", type="l", lwd=4, ylab="", xlab="")
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[2], lwd=3, lty = 2) ; points(var5~run_years, col=model_colours[2], pch=16)
legend("topleft", legend = "Copernicus", col = "black", lty = 1, pch=NA, horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
mtext(expression(paste('Year',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Analysis-wide LAI (',m^2,'/',m^2,')',sep="")), side = 2, cex = 2.4, padj = -1.05)
abline(0,1, col="grey", lwd=3)
# Now plot initial soil
plot(as.vector(1e-2*orig_grid_parameters$parameters[,,23,mid_quant]) ~ as.vector(1e-2*SoilCPrior), pch=1, cex = 1.6, cex.lab=2.4, cex.axis = 2.4, cex.main=2.0, ylab="", xlab="", main="", col=model_colours[1])
points(as.vector(1e-2*alt_grid_parameters$parameters[,,23,mid_quant]) ~ as.vector(1e-2*SoilCPrior), pch=1, cex = 1.6, col=model_colours[2])
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Initial soil C (MgC h',a^-1,')',sep="")), side = 2, cex = 2.4, padj = -1.00)
abline(0,1, col="grey", lwd=3)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
model_flags=c(orig_name,alt_name)
obs_flags=c("CTE","FC/Copernicus/FluxSatv2","GFEDv4.1s / GFAS")
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NBE_GPP_Fire_timeseries_comparison_plusCI",outsuffix,".png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot NBE, annual time series TgC/yr
dims = dim(cte_nbe_gCm2yr)
var1  = c(apply(cte_nbe_gCm2yr * array(cte_m2, dim=dims),c(3),sum, na.rm=TRUE) * 1e-12)
var2  = cbind(cbind(c(obs_nbe_mean_domain_TgCyr),c(obs_nbe_min_domain_TgCyr)),c(obs_nbe_max_domain_TgCyr))
var3  = orig_nbe_TgCyr ; var4  = orig_nbe_lower_TgCyr ; var5  = orig_nbe_upper_TgCyr
var6  = alt_nbe_TgCyr ; var7  = alt_nbe_lower_TgCyr ; var8  = alt_nbe_upper_TgCyr
zrange = range(c(var1,var2,var3,var4,var5,var6,var7,var8), na.rm=TRUE)
zrange[2] = zrange[2] + 500
plot(var3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
plotconfidence(var2,run_years,2,obs_colours[1])
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 1) #; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=3, lty = 2) #; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[2], lwd=3, lty = 1) #; points(var6~run_years, col=model_colours[2], pch=16)
lines(var7~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var7~run_years, col=model_colours[2], pch=16)
lines(var8~run_years, col=model_colours[2], lwd=3, lty = 2) #; points(var8~run_years, col=model_colours[2], pch=16)
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(obs_flags,model_flags), col = c(obs_colours[1:3],model_colours), 
       lty = c(rep(1,length(obs_flags)),rep(1,length(model_flags))), pch=rep(NA,length(c(obs_flags,model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Net Biome Exchange (TgC y",r^-1,")",sep="")), side=2, padj=-2.65,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot GPP
var3  = cbind(cbind(c(obs_gpp_mean_domain_TgCyr),c(obs_gpp_min_domain_TgCyr)),c(obs_gpp_max_domain_TgCyr))
var4  = orig_gpp_TgCyr ; var5  = orig_gpp_lower_TgCyr ; var6  = orig_gpp_upper_TgCyr   
var7  = alt_gpp_TgCyr ; var8  = alt_gpp_lower_TgCyr ; var9  = alt_gpp_upper_TgCyr   
zrange = range(c(var3,var4,var5,var6,var7,var8,var9), na.rm=TRUE)*c(0.9,1.0)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd = 4, ylab="", xlab="", lty = 2)
plotconfidence(var3,run_years,2,obs_colours[2])
lines(var4~run_years, col=model_colours[1], lwd = 4, lty = 1) #; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd = 4, lty = 2) #; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[1], lwd = 4, lty = 2) #; points(var6~run_years, col=model_colours[1], pch=16)
lines(var7~run_years, col=model_colours[2], lwd = 4, lty = 1) #; points(var7~run_years, col=model_colours[2], pch=16)
lines(var8~run_years, col=model_colours[2], lwd = 4, lty = 2) #; points(var8~run_years, col=model_colours[2], pch=16)
lines(var9~run_years, col=model_colours[2], lwd = 4, lty = 2) #; points(var9~run_years, col=model_colours[2], pch=16)

#legend("bottomright", legend = c(obs_flags[-5],model_flags), col = c(obs_colours[1:4],model_colours), 
#       lty = c(rep(1,length(obs_flags[-5])),rep(2,length(model_flags))), pch=rep(NA,length(c(obs_flags[-5],model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
#mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Gross Primary Productivity (TgC y",r^-1,")",sep="")), side=2, padj=-2.65, cex=1.5)

# Now plot fire
var3  = cbind(cbind(c(obs_fire_mean_domain_TgCyr),c(obs_fire_min_domain_TgCyr)),c(obs_fire_max_domain_TgCyr))
var4  = orig_fire_TgCyr  ; var5  = orig_fire_lower_TgCyr ; var6  = orig_fire_upper_TgCyr
var7  = alt_fire_TgCyr   ; var8  = alt_fire_lower_TgCyr  ; var9  = alt_fire_upper_TgCyr
zrange = range(c(var3,var4,var5,var6,var7,var8,var9), na.rm=TRUE)*c(0.9,1.1)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, lty=2, ylab="", xlab="")
plotconfidence(var3,run_years,2,obs_colours[3])
lines(var4~run_years, col=model_colours[1], lwd=4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var6~run_years, col=model_colours[1], pch=16)
lines(var7~run_years, col=model_colours[2], lwd=4, lty = 1) ; points(var7~run_years, col=model_colours[2], pch=16)
lines(var8~run_years, col=model_colours[2], lwd=4, lty = 2) ; points(var8~run_years, col=model_colours[2], pch=16)
lines(var9~run_years, col=model_colours[2], lwd=4, lty = 2) ; points(var9~run_years, col=model_colours[2], pch=16)

mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Fire Emissions (TgC y",r^-1,")",sep="")), side=2, padj=-2.65,cex=1.5)
dev.off()

# Assimilated observations overlap
# Assign variables
var1  = orig_grid_output$lai_assim_data_overlap_fraction
var2  = orig_grid_output$gpp_assim_data_overlap_fraction
var3  = orig_grid_output$nbe_assim_data_overlap_fraction
var4  = orig_grid_output$wood_assim_data_overlap_fraction
var5  = orig_grid_output$soil_assim_data_overlap_fraction
var6  = alt_grid_output$lai_assim_data_overlap_fraction
var7  = alt_grid_output$gpp_assim_data_overlap_fraction
var8  = alt_grid_output$nbe_assim_data_overlap_fraction
var9  = alt_grid_output$wood_assim_data_overlap_fraction
var10 = alt_grid_output$soil_assim_data_overlap_fraction
var11 = var6-var1
var12 = var7-var2
var13 = var8-var3
var14 = var9-var4
var15 = var10-var5
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
var10[which(is.na(landfilter))] = NA
var11[which(is.na(landfilter))] = NA
var12[which(is.na(landfilter))] = NA
var13[which(is.na(landfilter))] = NA
var14[which(is.na(landfilter))] = NA
var15[which(is.na(landfilter))] = NA
# Convert to raster
var1  = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2  = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3  = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4  = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5  = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6  = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7  = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8  = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9  = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var10 = raster(vals = t((var10)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var11 = raster(vals = t((var11)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var12 = raster(vals = t((var12)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var13 = raster(vals = t((var13)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var14 = raster(vals = t((var14)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var15 = raster(vals = t((var15)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1  = c(0,1)
zrange2  = c(0,1)
zrange3  = c(0,1)
zrange4  = c(0,1)
zrange5  = c(0,1)
zrange6  = c(0,1)
zrange7  = c(0,1)
zrange8  = c(0,1)
zrange9  = c(0,1)
zrange10 = c(0,1)
zrange11 = c(-1,1)
zrange12 = c(-1,1)
zrange13 = c(-1,1)
zrange14 = c(-1,1)
zrange15 = c(-1,1)
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_LAI_GPP_NBE_wood_soil_assimilated_observations_fraction_overlap_change",outsuffix,".png",sep=""), height = 3800, width = 5000, res = 300)
par(mfrow=c(3,5), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("LAI (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex = 2.0, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("GPP (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NBE (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Wood (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var5, zlim=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Soil (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var6, zlim=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2, cex = 2.0, padj = -0.5)
plot(var7, zlim=zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var9, zlim=zrange9, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var10, zlim=zrange10, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var11, zlim = zrange11, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Difference (-1-1)",sep=""), side=2, cex = 2.0, padj = -0.5)
plot(var12, zlim = zrange12, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var13, zlim = zrange13, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var14, zlim = zrange14, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var15, zlim = zrange15, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) relative difference in assimilated LAI overlap (",alt_name,"-",orig_name,")  = ",round(mean(as.vector(var11),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in assimilated GPP overlap (",alt_name,"-",orig_name,")  = ",round(mean(as.vector(var12),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in assimilated NBE overlap (",alt_name,"-",orig_name,")  = ",round(mean(as.vector(var13),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in assimilated wood overlap (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var14),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in assimilated soil overlap (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var15),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Are CARDAMOM models consistent with their assimilated observations
# Assimilated observations overlap
# Assign variables
var1 = orig_grid_output$lai_assim_data_overlap_fraction
var2 = orig_grid_output$wood_assim_data_overlap_fraction
var3 = orig_grid_output$soil_assim_data_overlap_fraction
var4 = alt_grid_output$lai_assim_data_overlap_fraction
var5 = alt_grid_output$wood_assim_data_overlap_fraction
var6 = alt_grid_output$soil_assim_data_overlap_fraction
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)
zrange2 = c(0,1)
zrange3 = c(0,1)
zrange4 = c(-1,1)
zrange5 = c(-1,1)
zrange6 = c(-1,1)
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_LAI_wood_soil_assimilated_observations_fraction_overlap_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("LAI (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex = 2.0, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Wood (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Soil (0-1)", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2, cex = 2.0, padj = -0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Difference (-1-1)",sep=""), side=2, cex = 2.0, padj = -0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) relative difference in assimilated LAI overlap (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in assimilated wood overlap (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in assimilated soil overlap (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()


###
## Plot carbon fluxes

# C fluxes
# Assign variables
var1 = orig_grid_output$mean_nbe_gCm2day[,,mid_quant]*1e-2*365.25
var2 = orig_grid_output$mean_gpp_gCm2day[,,mid_quant]*1e-2*365.25
var3 = orig_grid_output$mean_reco_gCm2day[,,mid_quant]*1e-2*365.25 
var4 = orig_grid_output$mean_fire_gCm2day[,,mid_quant]*1e-2*365.25 
var5 = alt_grid_output$mean_nbe_gCm2day[,,mid_quant]*1e-2*365.25
var6 = alt_grid_output$mean_gpp_gCm2day[,,mid_quant]*1e-2*365.25
var7 = alt_grid_output$mean_reco_gCm2day[,,mid_quant]*1e-2*365.25 
var8 = alt_grid_output$mean_fire_gCm2day[,,mid_quant]*1e-2*365.25 
var9 = (alt_grid_output$mean_nbe_gCm2day[,,mid_quant]-orig_grid_output$mean_nbe_gCm2day[,,mid_quant])*1e-2*365.25
var10 = (alt_grid_output$mean_gpp_gCm2day[,,mid_quant]-orig_grid_output$mean_gpp_gCm2day[,,mid_quant])*1e-2*365.25
var11 = (alt_grid_output$mean_reco_gCm2day[,,mid_quant]-orig_grid_output$mean_reco_gCm2day[,,mid_quant])*1e-2*365.25
var12 = (alt_grid_output$mean_fire_gCm2day[,,mid_quant]-orig_grid_output$mean_fire_gCm2day[,,mid_quant])*1e-2*365.25
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
var10[which(is.na(landfilter))] = NA
var11[which(is.na(landfilter))] = NA
var12[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var10 = raster(vals = t((var10)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var11 = raster(vals = t((var11)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var12 = raster(vals = t((var12)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(-1,1)*max(abs(range(c(values(var1),values(var5)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var6)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var7)),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(c(values(var4),values(var8)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var10)),na.rm=TRUE)))
zrange7 = c(-1,1)*max(abs(range(c(values(var11)),na.rm=TRUE)))
zrange8 = c(-1,1)*max(abs(range(c(values(var12)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_C_fluxes",outsuffix,".png",sep=""), height = 3000, width = 4900, res = 300)
par(mfrow=c(3,4), mar=c(0.5,0.5,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NBE (MgC h",a^-1,"y",r^-1,")", sep="")), col=rev(colour_choices_default))
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=1.8, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("GPP (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Reco (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_default))
plot(landmask, add=TRUE)
mtext(alt_name, side=2, cex=1.8, padj = -0.5)
plot(var6, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var7, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Difference
plot(var9, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=1.8, padj = -0.5)
plot(var10, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var11, zlim = zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var12, zlim = zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
dev.off()

# Must be paired with the above figure to get the right variables
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_C_fluxes_rel_change",outsuffix,".png",sep=""), height = 3000, width = 4900, res = 300)
par(mfrow=c(3,4), mar=c(0.5,0.5,2.8,7),omi=c(0.13,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NBE (MgC h",a^-1,"y",r^-1,")", sep="")), col=rev(colour_choices_default))
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=1.8, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("GPP (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Reco (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_default))
plot(landmask, add=TRUE)
mtext(alt_name, side=2, cex=1.8, padj = -0.5)
plot(var6, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var7, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Difference
var9 = var9 / abs(var1) ; var10 = var10 / abs(var2) ; var11 = var11 / abs(var3) ; var12 = var12 / abs(var4)
var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 ; var10[var10 > 1] = 1 ; var10[var10 < -1] = -1
var11[var11 > 1] = 1 ; var11[var11 < -1] = -1 ; var12[var12 > 1] = 1 ; var12[var12 < -1] = -1
zrange5 = c(-1,1) * max(abs(range(c(values(var9),values(var10),values(var11),values(var12)), na.rm=TRUE)))
zrange6 = zrange5 ; zrange7 = zrange5 ; zrange8 = zrange5
plot(var9, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=1.8, padj = -0.5)
plot(var10, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var11, zlim = zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var12, zlim = zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
dev.off()

# C fluxes
# Assign variables
var1 = (orig_grid_output$mean_nbe_gCm2day[,,high_quant]-orig_grid_output$mean_nbe_gCm2day[,,low_quant])*1e-2*365.25
var2 = (orig_grid_output$mean_gpp_gCm2day[,,high_quant]-orig_grid_output$mean_gpp_gCm2day[,,low_quant])*1e-2*365.25
var3 = (orig_grid_output$mean_reco_gCm2day[,,high_quant]-orig_grid_output$mean_reco_gCm2day[,,low_quant])*1e-2*365.25
var4 = (orig_grid_output$mean_fire_gCm2day[,,high_quant]-orig_grid_output$mean_fire_gCm2day[,,low_quant])*1e-2*365.25
var5 = (alt_grid_output$mean_nbe_gCm2day[,,high_quant]-alt_grid_output$mean_nbe_gCm2day[,,low_quant])*1e-2*365.25
var6 = (alt_grid_output$mean_gpp_gCm2day[,,high_quant]-alt_grid_output$mean_gpp_gCm2day[,,low_quant])*1e-2*365.25
var7 = (alt_grid_output$mean_reco_gCm2day[,,high_quant]-alt_grid_output$mean_reco_gCm2day[,,low_quant])*1e-2*365.25
var8 = (alt_grid_output$mean_fire_gCm2day[,,high_quant]-alt_grid_output$mean_fire_gCm2day[,,low_quant])*1e-2*365.25
var9 = var5-var1
var10 = var6-var2
var11 = var7-var3
var12 = var8-var4
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
var10[which(is.na(landfilter))] = NA
var11[which(is.na(landfilter))] = NA
var12[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var10 = raster(vals = t((var10)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var11 = raster(vals = t((var11)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var12 = raster(vals = t((var12)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var2),values(var3),values(var4),values(var5),values(var6),values(var7),values(var8)),na.rm=TRUE)))
zrange2 = zrange1
zrange3 = zrange1
zrange4 = zrange1
zrange5 = c(-1,1)*max(abs(range(c(values(var9),values(var10),values(var11),values(var12)),na.rm=TRUE)))
zrange6 = zrange5
zrange7 = zrange5
zrange8 = zrange5
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_C_fluxes_CI",outsuffix,".png",sep=""), height = 3000, width = 4900, res = 300)
par(mfrow=c(3,4), mar=c(0.5,0.5,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NBE CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=1.8, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("GPP CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Reco CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=1.8, padj = -0.5)
plot(var6, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var7, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var9, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=1.8, padj = -0.5)
plot(var10, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var11, zlim = zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var12, zlim = zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

# Keep with the figure above to ensure that the correct variables are available and used
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_C_fluxes_CI_rel_change",outsuffix,".png",sep=""), height = 3000, width = 4900, res = 300)
par(mfrow=c(3,4), mar=c(0.5,0.5,2.8,7),omi=c(0.13,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NBE CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=1.8, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("GPP CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Reco CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=1.8, padj = -0.5)
plot(var6, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var7, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
var9 = var9 / abs(var1) ; var10 = var10 / abs(var2) ; var11 = var11 / abs(var3) ; var12 = var12 / abs(var4)
var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 ; var10[var10 > 1] = 1 ; var10[var10 < -1] = -1
var11[var11 > 1] = 1 ; var11[var11 < -1] = -1 ; var12[var12 > 1] = 1 ; var12[var12 < -1] = -1
zrange5 = c(-1,1)*max(abs(range(c(values(var9),values(var10),values(var11),values(var12)),na.rm=TRUE)))
zrange6 = zrange5 ; zrange7 = zrange5 ; zrange8 = zrange5 ; 
# Difference
plot(var9, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=1.8, padj = -0.5)
plot(var10, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var11, zlim = zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var12, zlim = zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

###
## Plot the final C stocks, change and uncertainty

# Final stocks
# Assign variables
var1 = orig_grid_output$final_totalC_gCm2[,,mid_quant]*1e-2
var2 = orig_grid_output$final_biomass_gCm2[,,mid_quant]*1e-2 
var3 = orig_grid_output$final_dom_gCm2[,,mid_quant]*1e-2 
var4 = alt_grid_output$final_totalC_gCm2[,,mid_quant]*1e-2
var5 = alt_grid_output$final_biomass_gCm2[,,mid_quant]*1e-2
var6 = alt_grid_output$final_dom_gCm2[,,mid_quant]*1e-2
var7 = (alt_grid_output$final_totalC_gCm2[,,mid_quant]-orig_grid_output$final_totalC_gCm2[,,mid_quant])*1e-2
var8 = (alt_grid_output$final_biomass_gCm2[,,mid_quant]-orig_grid_output$final_biomass_gCm2[,,mid_quant])*1e-2
var9 = (alt_grid_output$final_dom_gCm2[,,mid_quant]-orig_grid_output$final_dom_gCm2[,,mid_quant])*1e-2
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var8[var8 > 1] = 1 ; var9[var9 > 1] = 1
var7[var7 < -1] = -1 ; var8[var8 < -1] = -1 ; var9[var9 < -1] = -1
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total (MgC h",a^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex = 2.0, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass (MgC h",a^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM (MgC h",a^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2, cex = 2.0, padj = -0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex = 2.0, padj = -0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) relative difference in Total (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in Biomass (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in DOM (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Change stocks
# Assign variables
var1 = orig_grid_output$final_dCtotalC_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var2 = orig_grid_output$final_dCbio_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var3 = orig_grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var4 = alt_grid_output$final_dCtotalC_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var5 = alt_grid_output$final_dCbio_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var6 = alt_grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var7 = (alt_grid_output$final_dCtotalC_gCm2[,,mid_quant]-orig_grid_output$final_dCtotalC_gCm2[,,mid_quant])*1e-2*(1/nos_years)
var8 = (alt_grid_output$final_dCbio_gCm2[,,mid_quant]-orig_grid_output$final_dCbio_gCm2[,,mid_quant])*1e-2*(1/nos_years)
var9 = (alt_grid_output$final_dCdom_gCm2[,,mid_quant]-orig_grid_output$final_dCdom_gCm2[,,mid_quant])*1e-2*(1/nos_years)
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var8[var8 > 1] = 1 ; var9[var9 > 1] = 1
var7[var7 < -1] = -1 ; var8[var8 < -1] = -1 ; var9[var9 < -1] = -1
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(-1,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(-1,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(-1,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_change_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_default)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_default)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_default)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_default)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_default)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_default)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) relative difference in Total change (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in Biomass change (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in DOM change (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Final stocks CI
# Assign variables
var1 = (orig_grid_output$final_totalC_gCm2[,,high_quant]-orig_grid_output$final_totalC_gCm2[,,low_quant])*1e-2 
var2 = (orig_grid_output$final_biomass_gCm2[,,high_quant]-orig_grid_output$final_biomass_gCm2[,,low_quant])*1e-2  
var3 = (orig_grid_output$final_dom_gCm2[,,high_quant]-orig_grid_output$final_dom_gCm2[,,low_quant])*1e-2  
var4 = (alt_grid_output$final_totalC_gCm2[,,high_quant]-alt_grid_output$final_totalC_gCm2[,,low_quant])*1e-2 
var5 = (alt_grid_output$final_biomass_gCm2[,,high_quant]-alt_grid_output$final_biomass_gCm2[,,low_quant])*1e-2  
var6 = (alt_grid_output$final_dom_gCm2[,,high_quant]-alt_grid_output$final_dom_gCm2[,,low_quant])*1e-2  
var7 = (alt_grid_output$final_totalC_gCm2[,,high_quant] - alt_grid_output$final_totalC_gCm2[,,low_quant])
var7 = (var7 - (orig_grid_output$final_totalC_gCm2[,,high_quant] - orig_grid_output$final_totalC_gCm2[,,low_quant])) * 1e-2
var8 = (alt_grid_output$final_biomass_gCm2[,,high_quant] - alt_grid_output$final_biomass_gCm2[,,low_quant])
var8 = (var8 - (orig_grid_output$final_biomass_gCm2[,,high_quant] - orig_grid_output$final_biomass_gCm2[,,low_quant])) * 1e-2
var9 = (alt_grid_output$final_dom_gCm2[,,high_quant] - alt_grid_output$final_dom_gCm2[,,low_quant])
var9 = (var9 - (orig_grid_output$final_dom_gCm2[,,high_quant] - orig_grid_output$final_dom_gCm2[,,low_quant])) * 1e-2
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var8[var8 > 1] = 1 ; var9[var9 > 1] = 1
var7[var7 < -1] = -1 ; var8[var8 < -1] = -1 ; var9[var9 < -1] = -1
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var2),values(var3),values(var4),values(var5),values(var6)),na.rm=TRUE)))
zrange2 = zrange1
zrange3 = zrange1
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.8,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) relative difference in CI Total (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in CI Biomass (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in CI DOM (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Change stocks CI
# Assign variables
var1 = (orig_grid_output$final_dCtotalC_gCm2[,,high_quant]-orig_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var2 = (orig_grid_output$final_dCbio_gCm2[,,high_quant]-orig_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var3 = (orig_grid_output$final_dCdom_gCm2[,,high_quant]-orig_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
var4 = (alt_grid_output$final_dCtotalC_gCm2[,,high_quant]-alt_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var5 = (alt_grid_output$final_dCbio_gCm2[,,high_quant]-alt_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var6 = (alt_grid_output$final_dCdom_gCm2[,,high_quant]-alt_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
var7 = (alt_grid_output$final_dCtotalC_gCm2[,,high_quant]-alt_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var7 = var7 - ((orig_grid_output$final_dCtotalC_gCm2[,,high_quant]-orig_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years))
var8 = (alt_grid_output$final_dCbio_gCm2[,,high_quant]-alt_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var8 = var8 - ((orig_grid_output$final_dCbio_gCm2[,,high_quant]-orig_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years))
var9 = (alt_grid_output$final_dCdom_gCm2[,,high_quant]-alt_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
var9 = var9 - ((orig_grid_output$final_dCdom_gCm2[,,high_quant]-orig_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years))
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var8[var8 > 1] = 1 ; var9[var9 > 1] = 1
var7[var7 < -1] = -1 ; var8[var8 < -1] = -1 ; var9[var9 < -1] = -1
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_change_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.9,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) relative difference in CI Total change (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in CI Biomass change (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) relative difference in CI DOM change (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Final stocks
# Assign variables
var1 = orig_grid_output$final_totalC_gCm2[,,mid_quant]*1e-2
var2 = orig_grid_output$final_biomass_gCm2[,,mid_quant]*1e-2 
var3 = orig_grid_output$final_dom_gCm2[,,mid_quant]*1e-2 
var4 = alt_grid_output$final_totalC_gCm2[,,mid_quant]*1e-2
var5 = alt_grid_output$final_biomass_gCm2[,,mid_quant]*1e-2
var6 = alt_grid_output$final_dom_gCm2[,,mid_quant]*1e-2
var7 = (alt_grid_output$final_totalC_gCm2[,,mid_quant]-orig_grid_output$final_totalC_gCm2[,,mid_quant])*1e-2
var8 = (alt_grid_output$final_biomass_gCm2[,,mid_quant]-orig_grid_output$final_biomass_gCm2[,,mid_quant])*1e-2
var9 = (alt_grid_output$final_dom_gCm2[,,mid_quant]-orig_grid_output$final_dom_gCm2[,,mid_quant])*1e-2
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total (MgC h",a^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex = 2.0, padj = -0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass (MgC h",a^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM (MgC h",a^-1,")", sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2, cex = 2.0, padj = -0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex = 2.0, padj = -0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha) difference in Total (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in Biomass (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in DOM (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Change stocks
# Assign variables
var1 = orig_grid_output$final_dCtotalC_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var2 = orig_grid_output$final_dCbio_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var3 = orig_grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var4 = alt_grid_output$final_dCtotalC_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var5 = alt_grid_output$final_dCbio_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var6 = alt_grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var7 = (alt_grid_output$final_dCtotalC_gCm2[,,mid_quant]-orig_grid_output$final_dCtotalC_gCm2[,,mid_quant])*1e-2*(1/nos_years)
var8 = (alt_grid_output$final_dCbio_gCm2[,,mid_quant]-orig_grid_output$final_dCbio_gCm2[,,mid_quant])*1e-2*(1/nos_years)
var9 = (alt_grid_output$final_dCdom_gCm2[,,mid_quant]-orig_grid_output$final_dCdom_gCm2[,,mid_quant])*1e-2*(1/nos_years)
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(-1,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(-1,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(-1,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_default)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_default)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_default)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_default)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_default)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_default)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_sign)
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in Total change (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in Biomass change (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in DOM change (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Final stocks CI
# Assign variables
var1 = (orig_grid_output$final_totalC_gCm2[,,high_quant]-orig_grid_output$final_totalC_gCm2[,,low_quant])*1e-2 
var2 = (orig_grid_output$final_biomass_gCm2[,,high_quant]-orig_grid_output$final_biomass_gCm2[,,low_quant])*1e-2  
var3 = (orig_grid_output$final_dom_gCm2[,,high_quant]-orig_grid_output$final_dom_gCm2[,,low_quant])*1e-2  
var4 = (alt_grid_output$final_totalC_gCm2[,,high_quant]-alt_grid_output$final_totalC_gCm2[,,low_quant])*1e-2 
var5 = (alt_grid_output$final_biomass_gCm2[,,high_quant]-alt_grid_output$final_biomass_gCm2[,,low_quant])*1e-2  
var6 = (alt_grid_output$final_dom_gCm2[,,high_quant]-alt_grid_output$final_dom_gCm2[,,low_quant])*1e-2  
var7 = (alt_grid_output$final_totalC_gCm2[,,high_quant] - alt_grid_output$final_totalC_gCm2[,,low_quant])
var7 = (var7 - (orig_grid_output$final_totalC_gCm2[,,high_quant] - orig_grid_output$final_totalC_gCm2[,,low_quant])) * 1e-2
var8 = (alt_grid_output$final_biomass_gCm2[,,high_quant] - alt_grid_output$final_biomass_gCm2[,,low_quant])
var8 = (var8 - (orig_grid_output$final_biomass_gCm2[,,high_quant] - orig_grid_output$final_biomass_gCm2[,,low_quant])) * 1e-2
var9 = (alt_grid_output$final_dom_gCm2[,,high_quant] - alt_grid_output$final_dom_gCm2[,,low_quant])
var9 = (var9 - (orig_grid_output$final_dom_gCm2[,,high_quant] - orig_grid_output$final_dom_gCm2[,,low_quant])) * 1e-2
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var2),values(var3),values(var4),values(var5),values(var6)),na.rm=TRUE)))
zrange2 = zrange1
zrange3 = zrange1
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.8,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha) difference in CI Total (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in CI Biomass (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in CI DOM (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Change stocks CI
# Assign variables
var1 = (orig_grid_output$final_dCtotalC_gCm2[,,high_quant]-orig_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var2 = (orig_grid_output$final_dCbio_gCm2[,,high_quant]-orig_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var3 = (orig_grid_output$final_dCdom_gCm2[,,high_quant]-orig_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
var4 = (alt_grid_output$final_dCtotalC_gCm2[,,high_quant]-alt_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var5 = (alt_grid_output$final_dCbio_gCm2[,,high_quant]-alt_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var6 = (alt_grid_output$final_dCdom_gCm2[,,high_quant]-alt_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
var7 = (alt_grid_output$final_dCtotalC_gCm2[,,high_quant]-alt_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var7 = var7 - ((orig_grid_output$final_dCtotalC_gCm2[,,high_quant]-orig_grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years))
var8 = (alt_grid_output$final_dCbio_gCm2[,,high_quant]-alt_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var8 = var8 - ((orig_grid_output$final_dCbio_gCm2[,,high_quant]-orig_grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years))
var9 = (alt_grid_output$final_dCdom_gCm2[,,high_quant]-alt_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
var9 = var9 - ((orig_grid_output$final_dCdom_gCm2[,,high_quant]-orig_grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years))
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_final_stocks_change_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.9,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM CI (MgC h",a^-1,"y",r^-1,")", sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in CI Total change (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in CI Biomass change (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in CI DOM change (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()


###
## Plot the MRTwood and NPPwood

# Traits
# Assign variables
var1 = orig_grid_parameters$MTT_wood_years[,,mid_quant]
var2 = orig_grid_parameters$NPP_wood_fraction[,,mid_quant]
var3 = orig_grid_parameters$SS_wood_gCm2[,,mid_quant]*1e-2
var4 = alt_grid_parameters$MTT_wood_years[,,mid_quant]
var5 = alt_grid_parameters$NPP_wood_fraction[,,mid_quant]
var6 = alt_grid_parameters$SS_wood_gCm2[,,mid_quant]*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# print summary information to user
print(paste("Mean (years) difference in wMTT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in wNPP (",alt_name,"-",orig_name,")  = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in wSS (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPP_MRT_SS",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPP_MRT_SS_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
# print summary information to user
print(paste("Mean relative (-1-1) difference in wMTT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wSS (",alt_name,"-",orig_name,")  = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

# Traits
# Assign variables
var1 = orig_grid_parameters$MTT_wood_years[,,mid_quant]
var2 = orig_grid_output$mean_wnpp_gCm2day[,,mid_quant]*365.25*1e-2
var3 = orig_grid_parameters$SS_wood_gCm2[,,mid_quant]*1e-2
var4 = alt_grid_parameters$MTT_wood_years[,,mid_quant]
var5 = alt_grid_output$mean_wnpp_gCm2day[,,mid_quant]*365.25*1e-2
var6 = alt_grid_parameters$SS_wood_gCm2[,,mid_quant]*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# print summary information to user
print(paste("Mean (years) difference in median wMTT (",alt_name,"-",orig_name,")        = ",round(mean(as.vector((var7)), na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in median wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector((var8)),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in median wSS (",alt_name,"-",orig_name,")        = ",round(mean(as.vector((var9)),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPPflx_MRT_SS",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPPflx_MRT_SS_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in wMTT (",alt_name,"-",orig_name,")    = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wSS (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Traits CI
# Assign variables
var1 = orig_grid_parameters$MTT_wood_years[,,high_quant] - orig_grid_parameters$MTT_wood_years[,,low_quant]
var2 = orig_grid_parameters$NPP_wood_fraction[,,high_quant] - orig_grid_parameters$NPP_wood_fraction[,,low_quant]
var3 = (orig_grid_parameters$SS_wood_gCm2[,,high_quant]*1e-2) - (orig_grid_parameters$SS_wood_gCm2[,,low_quant]*1e-2)
var4 = alt_grid_parameters$MTT_wood_years[,,high_quant] - alt_grid_parameters$MTT_wood_years[,,low_quant]
var5 = alt_grid_parameters$NPP_wood_fraction[,,high_quant] - alt_grid_parameters$NPP_wood_fraction[,,low_quant]
var6 = (alt_grid_parameters$SS_wood_gCm2[,,high_quant]*1e-2) - (alt_grid_parameters$SS_wood_gCm2[,,low_quant]*1e-2)
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPP_MRT_SS_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (years) difference in CI wMTT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (0-1) difference in CI wNPP (",alt_name,"-",orig_name,")   = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in CI wSS (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPP_MRT_SS_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in CI wMTT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wSS (",alt_name,"-",orig_name,")  = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Traits CI
# Assign variables
var1 = orig_grid_parameters$MTT_wood_years[,,high_quant] - orig_grid_parameters$MTT_wood_years[,,low_quant]
var2 = (orig_grid_output$mean_wnpp_gCm2day[,,high_quant] - orig_grid_output$mean_wnpp_gCm2day[,,low_quant])*365.25*1e-2
var3 = (orig_grid_parameters$SS_wood_gCm2[,,high_quant]*1e-2) - (orig_grid_parameters$SS_wood_gCm2[,,low_quant]*1e-2)
var4 = alt_grid_parameters$MTT_wood_years[,,high_quant] - alt_grid_parameters$MTT_wood_years[,,low_quant]
var5 = (alt_grid_output$mean_wnpp_gCm2day[,,high_quant] - alt_grid_output$mean_wnpp_gCm2day[,,low_quant])*365.25*1e-2
var6 = (alt_grid_parameters$SS_wood_gCm2[,,high_quant]*1e-2) - (alt_grid_parameters$SS_wood_gCm2[,,low_quant]*1e-2)
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPPflx_MRT_SS_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (years) difference in CI wMTT (",alt_name,"-",orig_name,")        = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in CI wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha) difference in CI wSS (",alt_name,"-",orig_name,")        = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_NPPflx_MRT_SS_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("SS wood CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in CI wMTT (",alt_name,"-",orig_name,")    = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wSS (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

###
## Plot the Foliage, fine root and wood NPP, NPPflx and NPPwood

# Traits
# Assign variables
var1 = orig_grid_parameters$NPP_foliar_fraction[,,mid_quant]
var2 = orig_grid_parameters$NPP_root_fraction[,,mid_quant]
var3 = orig_grid_parameters$NPP_wood_fraction[,,mid_quant]
var4 = alt_grid_parameters$NPP_foliar_fraction[,,mid_quant]
var5 = alt_grid_parameters$NPP_root_fraction[,,mid_quant]
var6 = alt_grid_parameters$NPP_wood_fraction[,,mid_quant]
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# print summary information to user
print(paste("Mean (-1-1) difference in fNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in rNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in wNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPP",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPP_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (0-1)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
# print summary information to user
print(paste("Mean relative (-1-1) difference in fNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

# Traits
# Assign variables
var1 = orig_grid_output$mean_fnpp_gCm2day[,,mid_quant]*365.25*1e-2
var2 = orig_grid_output$mean_rnpp_gCm2day[,,mid_quant]*365.25*1e-2
var3 = orig_grid_output$mean_wnpp_gCm2day[,,mid_quant]*365.25*1e-2
var4 = alt_grid_output$mean_fnpp_gCm2day[,,mid_quant]*365.25*1e-2
var5 = alt_grid_output$mean_rnpp_gCm2day[,,mid_quant]*365.25*1e-2
var6 = alt_grid_output$mean_wnpp_gCm2day[,,mid_quant]*365.25*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in median fNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector((var7)), na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in median rNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector((var8)),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in median wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector((var9)),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPPflx",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPPflx_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in fNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Traits CI
# Assign variables
var1 = orig_grid_parameters$NPP_foliar_fraction[,,high_quant] - orig_grid_parameters$NPP_foliar_fraction[,,low_quant]
var2 = orig_grid_parameters$NPP_root_fraction[,,high_quant] - orig_grid_parameters$NPP_root_fraction[,,low_quant]
var3 = orig_grid_parameters$NPP_wood_fraction[,,high_quant] - orig_grid_parameters$NPP_wood_fraction[,,low_quant]
var4 = alt_grid_parameters$NPP_foliar_fraction[,,high_quant] - alt_grid_parameters$NPP_foliar_fraction[,,low_quant]
var5 = alt_grid_parameters$NPP_root_fraction[,,high_quant] - alt_grid_parameters$NPP_root_fraction[,,low_quant]
var6 = alt_grid_parameters$NPP_wood_fraction[,,high_quant] - alt_grid_parameters$NPP_wood_fraction[,,low_quant]
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPP_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) difference in CI fNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in CI rNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in CI wNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPP_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (0-1)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in CI fNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI rNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wNPP (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Traits CI
# Assign variables
var1 = (orig_grid_output$mean_fnpp_gCm2day[,,high_quant] - orig_grid_output$mean_fnpp_gCm2day[,,low_quant])*365.25*1e-2
var2 = (orig_grid_output$mean_rnpp_gCm2day[,,high_quant] - orig_grid_output$mean_rnpp_gCm2day[,,low_quant])*365.25*1e-2
var3 = (orig_grid_output$mean_wnpp_gCm2day[,,high_quant] - orig_grid_output$mean_wnpp_gCm2day[,,low_quant])*365.25*1e-2
var4 = (alt_grid_output$mean_fnpp_gCm2day[,,high_quant] - alt_grid_output$mean_fnpp_gCm2day[,,low_quant])*365.25*1e-2
var5 = (alt_grid_output$mean_rnpp_gCm2day[,,high_quant] - alt_grid_output$mean_rnpp_gCm2day[,,low_quant])*365.25*1e-2
var6 = (alt_grid_output$mean_wnpp_gCm2day[,,high_quant] - alt_grid_output$mean_wnpp_gCm2day[,,low_quant])*365.25*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPPflx_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in CI fNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in CI rNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in CI wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_NPPflx_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP foliar CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP root CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("NPP wood CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in CI fNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI rNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wNPPflx (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

###
## Plot the Foliage, fine root and wood wMRT 

# Traits
# Assign variables
var1 = orig_grid_parameters$MTT_foliar_years[,,mid_quant]
var2 = orig_grid_parameters$MTT_root_years[,,mid_quant]
var3 = orig_grid_parameters$MTT_wood_years[,,mid_quant]
var4 = alt_grid_parameters$MTT_foliar_years[,,mid_quant]
var5 = alt_grid_parameters$MTT_root_years[,,mid_quant]
var6 = alt_grid_parameters$MTT_wood_years[,,mid_quant]
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# print summary information to user
print(paste("Mean (-1-1) difference in fMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in rMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in wMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_MRT",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT foliar (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT root (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_MRT_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT foliar (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT root (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood (years)",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
# print summary information to user
print(paste("Mean relative (-1-1) difference in fMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
dev.off()

# Traits CI
# Assign variables
var1 = orig_grid_parameters$MTT_foliar_years[,,high_quant] - orig_grid_parameters$MTT_foliar_years[,,low_quant]
var2 = orig_grid_parameters$MTT_root_years[,,high_quant] - orig_grid_parameters$MTT_root_years[,,low_quant]
var3 = orig_grid_parameters$MTT_wood_years[,,high_quant] - orig_grid_parameters$MTT_wood_years[,,low_quant]
var4 = alt_grid_parameters$MTT_foliar_years[,,high_quant] - alt_grid_parameters$MTT_foliar_years[,,low_quant]
var5 = alt_grid_parameters$MTT_root_years[,,high_quant] - alt_grid_parameters$MTT_root_years[,,low_quant]
var6 = alt_grid_parameters$MTT_wood_years[,,high_quant] - alt_grid_parameters$MTT_wood_years[,,low_quant]
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_MRT_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT foliar CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT root CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (years) difference in CI fMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (years) difference in CI rMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (years) difference in CI wMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_fol_root_wood_MRT_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT foliar CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT root CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("MRT wood CI (years)",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in CI fMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI rMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

###
## Partition the importance of disturbance on residence times

# Estimate the proportion of turnover determined by natural
orig_turn = (orig_grid_parameters$MTT_wood_years[,,mid_quant]**-1)
alt_turn = (alt_grid_parameters$MTT_wood_years[,,mid_quant]**-1)
var1 = orig_grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1 
var2 = orig_grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var3 = orig_grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
var4 = alt_grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1 
var5 = alt_grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var6 = alt_grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
# Now make proportional
var1 = var1 / orig_turn ; var2 = var2 / orig_turn ; var3 = var3 / orig_turn
var4 = var4 / alt_turn ; var5 = var5 / alt_turn ; var6 = var6 / alt_turn
# Filter for the AGB map locations
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
# Further filtering
var1[which(landfilter == 1 & is.na(var1) == TRUE)] = 0
var2[which(landfilter == 1 & is.na(var2) == TRUE)] = 0
var3[which(landfilter == 1 & is.na(var3) == TRUE)] = 0
var4[which(landfilter == 1 & is.na(var4) == TRUE)] = 0
var5[which(landfilter == 1 & is.na(var5) == TRUE)] = 0
var6[which(landfilter == 1 & is.na(var6) == TRUE)] = 0
# print summary information to user
print(paste("Mean natural ",orig_name," wMTT (years)     = ",round(mean(as.vector(var1),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean fire ",orig_name," wMTT (years)        = ",round(mean(as.vector(var2),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean forest loss ",orig_name," wMTT (years) = ",round(mean(as.vector(var3),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean natural ",alt_name," wMTT (years)      = ",round(mean(as.vector(var4),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean fire ",alt_name," wMTT (years)         = ",round(mean(as.vector(var5),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean forest loss ",alt_name," wMTT (years)  = ",round(mean(as.vector(var6),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# specify ranges
zrange1 = c(0,1)
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_wood_turnover_contribution",outsuffix,".png",sep=""), height = 2500, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Natural MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj = -0.5)
plot(var2, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass removal MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj = -0.5)
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
dev.off()

# Fire driven tissue specific emissions of C
# Assign variables
var1 = orig_grid_output$mean_FIREemiss_foliar_gCm2yr[,,mid_quant]*1e-2
var2 = orig_grid_output$mean_FIREemiss_root_gCm2yr[,,mid_quant]*1e-2
var3 = orig_grid_output$mean_FIREemiss_wood_gCm2yr[,,mid_quant]*1e-2
var4 = alt_grid_output$mean_FIREemiss_foliar_gCm2yr[,,mid_quant]*1e-2
var5 = alt_grid_output$mean_FIREemiss_root_gCm2yr[,,mid_quant]*1e-2
var6 = alt_grid_output$mean_FIREemiss_wood_gCm2yr[,,mid_quant]*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireEmiss",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted foliage (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted root (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted wood (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in fFireEmiss (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in rFireEmiss (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in wFireEmiss (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireEmiss_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted foliage (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted root (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted wood (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in fFireEmiss (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rFireEmiss (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wFireEmiss (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Fire driven tissue specific litter of C
# Assign variables
var1 = orig_grid_output$mean_FIRElitter_foliar_gCm2yr[,,mid_quant]*1e-2
var2 = orig_grid_output$mean_FIRElitter_root_gCm2yr[,,mid_quant]*1e-2
var3 = orig_grid_output$mean_FIRElitter_wood_gCm2yr[,,mid_quant]*1e-2
var4 = alt_grid_output$mean_FIRElitter_foliar_gCm2yr[,,mid_quant]*1e-2
var5 = alt_grid_output$mean_FIRElitter_root_gCm2yr[,,mid_quant]*1e-2
var6 = alt_grid_output$mean_FIRElitter_wood_gCm2yr[,,mid_quant]*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireLitter",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Foliage fire mortality (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Root fire mortality (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Wood fire mortality (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in fFireLitter (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in rFireLitter (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in wFireLitter (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireLitter_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Foliage fire mortality (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Root fire mortality (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Wood fire mortality (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in fFireLitter (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rFireLitter (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wFireLitter (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Fire driven tissue specific emissions of C
# Assign variables
var1 = (orig_grid_output$mean_FIREemiss_foliar_gCm2yr[,,high_quant] - orig_grid_output$mean_FIREemiss_foliar_gCm2yr[,,low_quant])*1e-2
var2 = (orig_grid_output$mean_FIREemiss_root_gCm2yr[,,high_quant] - orig_grid_output$mean_FIREemiss_root_gCm2yr[,,low_quant])*1e-2
var3 = (orig_grid_output$mean_FIREemiss_wood_gCm2yr[,,high_quant] - orig_grid_output$mean_FIREemiss_wood_gCm2yr[,,low_quant])*1e-2
var4 = (alt_grid_output$mean_FIREemiss_foliar_gCm2yr[,,high_quant] - alt_grid_output$mean_FIREemiss_foliar_gCm2yr[,,low_quant])*1e-2
var5 = (alt_grid_output$mean_FIREemiss_root_gCm2yr[,,high_quant] - alt_grid_output$mean_FIREemiss_root_gCm2yr[,,low_quant])*1e-2
var6 = (alt_grid_output$mean_FIREemiss_wood_gCm2yr[,,high_quant] - alt_grid_output$mean_FIREemiss_wood_gCm2yr[,,low_quant])*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireEmiss_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted foliage CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted root CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted wood CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in fFireEmiss CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in rFireEmiss CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in wFireEmiss CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireEmiss_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted foliage CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted root CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Combusted wood CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in fFireEmiss CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rFireEmiss CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wFireEmiss CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Fire driven tissue specific litter of C
# Assign variables
var1 = (orig_grid_output$mean_FIRElitter_foliar_gCm2yr[,,high_quant] - orig_grid_output$mean_FIRElitter_foliar_gCm2yr[,,low_quant])*1e-2
var2 = (orig_grid_output$mean_FIRElitter_root_gCm2yr[,,high_quant] - orig_grid_output$mean_FIRElitter_root_gCm2yr[,,low_quant])*1e-2
var3 = (orig_grid_output$mean_FIRElitter_wood_gCm2yr[,,high_quant] - orig_grid_output$mean_FIRElitter_wood_gCm2yr[,,low_quant])*1e-2
var4 = (alt_grid_output$mean_FIRElitter_foliar_gCm2yr[,,high_quant] - alt_grid_output$mean_FIRElitter_foliar_gCm2yr[,,low_quant])*1e-2
var5 = (alt_grid_output$mean_FIRElitter_root_gCm2yr[,,high_quant] - alt_grid_output$mean_FIRElitter_root_gCm2yr[,,low_quant])*1e-2
var6 = (alt_grid_output$mean_FIRElitter_wood_gCm2yr[,,high_quant] - alt_grid_output$mean_FIRElitter_wood_gCm2yr[,,low_quant])*1e-2
var7 = var4-var1
var8 = var5-var2
var9 = var6-var3
# Apply filter
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var4)),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var5)),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(c(values(var3),values(var6)),na.rm=TRUE)))
zrange4 = c(-1,1)*max(abs(range(c(values(var7)),na.rm=TRUE)))
zrange5 = c(-1,1)*max(abs(range(c(values(var8)),na.rm=TRUE)))
zrange6 = c(-1,1)*max(abs(range(c(values(var9)),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireLitter_CI",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Foliage fire mortality CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Root fire mortality CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Wood fire mortality CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Difference
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=rev(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (MgC/ha/yr) difference in fFireLitter CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in rFireLitter CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (MgC/ha/yr) difference in wFireLitter CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

zrange4 = c(-1,1) ; zrange5 = zrange4 ; zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_FireLitter_CI_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.11,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Foliage fire mortality CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(orig_name, side=2, cex=2.0, padj=-0.5)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Root fire mortality CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Wood fire mortality CI (MgC h",a^-1,y^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_CI)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in fFireLitter CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in rFireLitter CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in wFireLitter CI (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

###
## Partition the importance of disturbance on residence times

# Estimate the proportion of turnover determined by natural
orig_turn = (orig_grid_parameters$MTT_wood_years[,,mid_quant]**-1)
alt_turn = (alt_grid_parameters$MTT_wood_years[,,mid_quant]**-1)
var1 = orig_grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1 
var2 = orig_grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var3 = orig_grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
var4 = alt_grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1 
var5 = alt_grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var6 = alt_grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
# Now make proportional
var1 = var1 / orig_turn ; var2 = var2 / orig_turn ; var3 = var3 / orig_turn
var4 = var4 / alt_turn ; var5 = var5 / alt_turn ; var6 = var6 / alt_turn
# Create difference
var7 = var4 - var1 ; var8 = var5 - var2 ; var9 = var6 - var3
# Filter for the AGB map locations
var1[which(is.na(landfilter))] = NA
var2[which(is.na(landfilter))] = NA
var3[which(is.na(landfilter))] = NA
var4[which(is.na(landfilter))] = NA
var5[which(is.na(landfilter))] = NA
var6[which(is.na(landfilter))] = NA
var7[which(is.na(landfilter))] = NA
var8[which(is.na(landfilter))] = NA
var9[which(is.na(landfilter))] = NA
# Further filtering
var1[which(landfilter == 1 & is.na(var1) == TRUE)] = 0
var2[which(landfilter == 1 & is.na(var2) == TRUE)] = 0
var3[which(landfilter == 1 & is.na(var3) == TRUE)] = 0
var4[which(landfilter == 1 & is.na(var4) == TRUE)] = 0
var5[which(landfilter == 1 & is.na(var5) == TRUE)] = 0
var6[which(landfilter == 1 & is.na(var6) == TRUE)] = 0
var7[which(landfilter == 1 & is.na(var7) == TRUE)] = 0
var8[which(landfilter == 1 & is.na(var8) == TRUE)] = 0
var9[which(landfilter == 1 & is.na(var9) == TRUE)] = 0
# print summary information to user
print(paste("Mean natural ",orig_name," wMTT (years)     = ",round(mean(as.vector(var1),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean fire ",orig_name," wMTT (years)        = ",round(mean(as.vector(var2),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean forest loss ",orig_name," wMTT (years) = ",round(mean(as.vector(var3),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean natural ",alt_name," wMTT (years)      = ",round(mean(as.vector(var4),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean fire ",alt_name," wMTT (years)         = ",round(mean(as.vector(var5),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean forest loss ",alt_name," wMTT (years)  = ",round(mean(as.vector(var6),na.rm=TRUE),digits=3),sep=""))
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var9 = raster(vals = t((var9)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# specify ranges
zrange1 = c(0,1)
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_wood_turnover_contribution_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Natural MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2,cex=2.0, padj=-0.5)
plot(var2, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass removal MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Difference
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Difference",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean (-1-1) difference in natMRT (",alt_name,"-",orig_name,")     = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in fireMRT (",alt_name,"-",orig_name,")    = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean (-1-1) difference in harvestMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

# Relative change version
png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_wood_turnover_contribution_rel_change",outsuffix,".png",sep=""), height = 4000, width = 4900, res = 300)
par(mfrow=c(3,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Original
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Natural MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(orig_name, side=2,cex=2.0, padj=-0.5)
plot(var2, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass removal MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Alternate
plot(var4, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
mtext(alt_name, side=2,cex=2.0, padj=-0.5)
plot(var5, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=colour_choices_loss)
plot(landmask, add=TRUE)
# Relative Difference
var7 = var7 / abs(var1) ; var8 = var8 / abs(var2) ; var9 = var9 / abs(var3)
var7[var7 > 1] = 1 ; var7[var7 < -1] = -1 ; var8[var8 > 1] = 1 ; var8[var8 < -1] = -1 ; var9[var9 > 1] = 1 ; var9[var9 < -1] = -1 
zrange4 = c(-1,1) * max(abs(range(c(values(var7),values(var8),values(var9)), na.rm=TRUE)))
zrange5 = zrange4 ; zrange6 = zrange5
plot(var7, zlim = zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
mtext(paste("Relative difference (-1-1)",sep=""), side=2, cex=2.0, padj=-0.5)
plot(var8, zlim = zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
plot(var9, zlim = zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "", col=(colour_choices_sign))
plot(landmask, add=TRUE)
# print summary information to user
print(paste("Mean relative (-1-1) difference in CI fMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var7),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI rMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var8),na.rm=TRUE),digits=3),sep=""))
print(paste("Mean relative (-1-1) difference in CI wMRT (",alt_name,"-",orig_name,") = ",round(mean(as.vector(var9),na.rm=TRUE),digits=3),sep=""))
dev.off()

###
## Extract pixels constaining the 200 sites

skip = TRUE
if (skip == FALSE) {
tmp = read.csv("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/SECO/miombo_site_information/extract_points2.csv", header=TRUE)
sites_cardamom=paste(names(tmp)[1],tmp$ID,sep="")
# lat/long of sites, if type = "grid"then these these are bottom left and top right corners
sites_cardamom_lat=tmp$lat
sites_cardamom_long=tmp$long

# Loop all 200 sites and find their locations
WoodCobs_trend = rep(NA, dim(tmp)[1])
orig_wood_trend = rep(NA, dim(tmp)[1])
alt_wood_trend = rep(NA, dim(tmp)[1])
orig_obs_wood_stock = rep(NA, dim(tmp)[1])
alt_obs_wood_stock = rep(NA, dim(tmp)[1])
orig_wood_natural_mrt = rep(NA, dim(tmp)[1])
alt_wood_natural_mrt = rep(NA, dim(tmp)[1])
orig_wood_mrt = rep(NA, dim(tmp)[1]) # years
alt_wood_mrt = rep(NA, dim(tmp)[1]) # years
orig_wood_npp = rep(NA, dim(tmp)[1]) # fraction
alt_wood_npp = rep(NA, dim(tmp)[1]) # fraction
obs_period_end = 5*12 ; obs_period_years = 5
for (n in seq(1, dim(tmp)[1])) {
     ij = unlist(closest2d(n,grid_lat,grid_long,tmp$lat,tmp$long,2))
     nn = which(as.numeric(alt_PROJECT$sites) == ((ij[2]-1) * alt_PROJECT$long_dim) + ij[1])
     ## Trends
     if (length(which(is.na(WoodCobs[ij[1],ij[2],]) == FALSE)) > 0) {
         WoodCobs_trend[n] = (coef(lm(WoodCobs[ij[1],ij[2],] ~ c(1:dim(WoodCobs)[3])))[2] * 12) # *12 is month to yr adjustment
         orig_wood_trend[n] = (coef(lm(orig_grid_output$wood_gCm2[nn,mid_quant,1:obs_period_end] ~ c(1:obs_period_end)))[2] * 12)
         alt_wood_trend[n] = (coef(lm(alt_grid_output$wood_gCm2[nn,mid_quant,1:obs_period_end] ~ c(1:obs_period_end)))[2] * 12)
         orig_obs_wood_stock[n] = (WoodCobs[ij[1],ij[2],which(WoodCobs[ij[1],ij[2],] > 0)[1]])
         alt_obs_wood_stock[n] = mean(WoodCobs[ij[1],ij[2],], na.rm=TRUE)
     }
     ## Static
     orig_wood_natural_mrt[n] = as.vector(orig_grid_parameters$parameters[ij[1],ij[2],6,mid_quant]*365.25)**-1
     alt_wood_natural_mrt[n] = as.vector(alt_grid_parameters$parameters[ij[1],ij[2],6,mid_quant]*365.25)**-1
     orig_wood_mrt[n] = as.vector(orig_grid_parameters$MTT_wood_years[ij[1],ij[2],mid_quant])
     alt_wood_mrt[n] = as.vector(alt_grid_parameters$MTT_wood_years[ij[1],ij[2],mid_quant])
     orig_wood_npp[n] = as.vector(orig_grid_parameters$NPP_wood_fraction[ij[1],ij[2],mid_quant])
     alt_wood_npp[n] = as.vector(alt_grid_parameters$NPP_wood_fraction[ij[1],ij[2],mid_quant])
}

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"wood_trend_199pixels",outsuffix,".png",sep=""), 
    height = 3500, width = 4000, res = 300)
xrange = range(c(orig_wood_trend,alt_wood_trend), na.rm=TRUE)
par(mfrow=c(2,2), mar=c(4.2,5.2,3.0,2),omi=c(0.01,0.01,0.01,0.01))
plot(WoodCobs_trend ~ orig_wood_trend, ylab = expression(paste("Obs wood trend (gC ",m^-2,"y",r^-1,")",sep="")), 
     main = paste(orig_name,sep=""), xlab = expression(paste("Model wood trend (gC ",m^-2,"y",r^-1,")",sep="")), 
     pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2, xlim=xrange)
abline(0,1,col="red", lwd=3) ; abline(0,0,col="grey", lwd=2) ; abline(v = 0,col="grey", lwd=2)
plot(WoodCobs_trend ~ alt_wood_trend, main=paste(alt_name,sep=""), ylab = "", 
     xlab = expression(paste("Model wood trend (gC ",m^-2,"y",r^-1,")",sep="")), 
     pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2, xlim=xrange)
abline(0,1,col="red", lwd=3) ; abline(0,0,col="grey", lwd=2) ; abline(v = 0,col="grey", lwd=2)
yrange = range(c((orig_wood_trend - WoodCobs_trend),(alt_wood_trend - WoodCobs_trend)), na.rm=TRUE)
plot((orig_wood_trend - WoodCobs_trend)  ~ orig_obs_wood_stock, 
     ylab = expression(paste("Model - Obs trend bias",sep="")), 
     xlab = expression(paste("Mean assimilated wood stock (gC ",m^-2,")",sep="")), 
     pch = 16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2, ylim=yrange)
abline(0,0,col="grey", lwd=2)
plot((alt_wood_trend - WoodCobs_trend)  ~ alt_obs_wood_stock, ylab=expression(paste("Model - Obs trend bias",sep="")), 
     xlab = expression(paste("Mean assimilated wood stock (gC ",m^-2,")",sep="")), pch = 16, cex = 2, 
     cex.main=2, cex.axis = 2.2, cex.lab=2.2, ylim=yrange)
abline(0,0,col="grey", lwd=2)
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_wood_MRT_199pixels",outsuffix,".png",sep=""), height = 4000, width = 2500, res = 300)
par(mfrow=c(4,2), mar=c(5.0,5.0,5.0,2),omi=c(0.1,0.1,0.1,0.1))
plot(orig_wood_natural_mrt ~ orig_wood_trend, main=orig_name, ylab="Wood natural MRT (years)", xlab="", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(alt_wood_natural_mrt ~ alt_wood_trend, main=alt_name, ylab="", xlab="", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(orig_wood_mrt ~ orig_wood_trend, main="", ylab="Wood MRT (years)", xlab="Model wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(alt_wood_mrt ~ alt_wood_trend, main="", ylab="", xlab="Model wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)

plot(orig_wood_natural_mrt ~ WoodCobs_trend, main="", ylab="Wood natural MRT (years)", xlab="", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(alt_wood_natural_mrt ~ WoodCobs_trend, main="", ylab="", xlab="", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(orig_wood_mrt ~ WoodCobs_trend, main="", ylab="Wood MRT (years)", xlab="Obs wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(alt_wood_mrt ~ WoodCobs_trend, main="", ylab="", xlab="Obs wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_wood_NPP_199pixels",outsuffix,".png",sep=""), height = 2500, width = 2500, res = 300)
par(mfrow=c(2,2), mar=c(5.0,5.0,5.0,2),omi=c(0.1,0.1,0.1,0.1))
plot(orig_wood_npp ~ orig_wood_trend, main=orig_name, ylab="Wood NPP (0-1)", xlab="Model wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(alt_wood_npp ~ alt_wood_trend, main=alt_name, ylab="", xlab="Model wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(orig_wood_npp ~ WoodCobs_trend, main="", ylab="Wood NPP (0-1)", xlab="Obs wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
plot(alt_wood_npp ~ WoodCobs_trend, main="", ylab="", xlab="Obs wood trend", pch=16, cex = 2, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
abline(v = 0, col="red", lwd=3)
dev.off()

###
## Compare against site level observations

###
## Read available site observations

# ID2
bob = read.csv("/home/lsmallma/gcel/miombo/CARDAMOM/FLX_ZA-Kru_FLUXNET2015_SUBSET_MM_2000-2010_1-3.csv", header=TRUE)
# Define obs
kru_gpp = rep(NA,12)
kru_reco = rep(NA,12)
kru_nee = rep(NA,12)
kru_nee_unc = rep(NA,12)

# Average to monthly means as the time periods do not overlap
for (i in seq(1,12)) {
     # Obs
     get = seq(i,dim(bob)[1],12)
     kru_gpp[i] = mean(bob$GPP_NT_VUT_REF[get])
     kru_reco[i] = mean(bob$RECO_NT_VUT_REF[get])
     kru_nee[i] = mean(bob$NEE_VUT_REF[get])
     kru_nee_unc[i] = mean(bob$NEE_VUT_REF_RANDUNC[get])

}
# Define function to estimate the average month
avg_month <-function(var) {
   out_var = rep(NA, 12) 
    for (n in seq(1,12)) {
         out_var[n] = mean(var[seq(n,length(var),12)], na.rm=TRUE)
    } 
    return(out_var) 
}

png(file=paste(out_dir,"/",gsub("%","_",orig_PROJECT$name),"_Comparison_of_gpp_nee_lai_ID2",outsuffix,".png",sep=""), width = 2100, height = 2000, res = 300)
par(mfrow=c(2,2))
# find location in grid for this site
n = 2
ij = unlist(closest2d(n,grid_lat,grid_long,tmp$lat,tmp$long,2))
nn = which(as.numeric(alt_PROJECT$sites) == ((ij[2]-1) * alt_PROJECT$long_dim) + ij[1])
# read assimilated observations
drivers = read_binary_file_format(paste(alt_PROJECT$datapath,alt_PROJECT$name,"_",alt_PROJECT$sites[nn],".bin",sep=""))
# GPP
yrange = range(c(orig_grid_output$gpp_gCm2day[nn,mid_quant,],alt_grid_output$gpp_gCm2day[nn,mid_quant,],kru_gpp))
plot(avg_month(orig_grid_output$gpp_gCm2day[nn,mid_quant,]), main="ID2", 
     ylab=expression(paste("Mean GPP (gC ",m^-2,d^-1,")",sep="")), xlab="Mean Month", type="l", col="red", lwd=3, ylim=yrange)
lines(avg_month(alt_grid_output$gpp_gCm2day[nn,mid_quant,]), lwd=3, col="blue")
points(kru_gpp, pch=16)
# NEE
yrange = range(c(orig_grid_output$nee_gCm2day[nn,mid_quant,],alt_grid_output$nee_gCm2day[nn,mid_quant,],kru_nee))
plot(avg_month(orig_grid_output$nee_gCm2day[nn,mid_quant,]), main="ID2", 
     ylab=expression(paste("Mean NEE (gC ",m^-2,d^-1,")",sep="")), 
     xlab="Mean Month", type="l", col="red", lwd=3, ylim=yrange)
lines(avg_month(alt_grid_output$nee_gCm2day[nn,mid_quant,]), lwd=3, col="blue")
points(kru_nee, pch=16)
legend("topleft", col=c("red","blue"), legend = c("One","All"), pch = c(NA,NA), lty = c(1,1))
# LAI
yrange = range(c(orig_grid_output$lai_m2m2[nn,mid_quant,],alt_grid_output$lai_m2m2[nn,mid_quant,],drivers$obs[which(drivers$obs[,3] > 0),3]))
plot(orig_grid_output$lai_m2m2[nn,mid_quant,], main="ID2", ylab=expression(paste("LAI (",m^2,m^-2,")",sep="")), 
     xlab="Simulation Month", type="l", col="red", lwd=3, ylim=yrange)
lines(alt_grid_output$lai_m2m2[nn,mid_quant,], lwd=3, col="blue")
points(drivers$obs[,3], pch=16)
# Wood
yrange = range(c(orig_grid_output$wood_gCm2[nn,mid_quant,],alt_grid_output$wood_gCm2[nn,mid_quant,],drivers$obs[which(drivers$obs[,13] > 0),13]))
plot(orig_grid_output$wood_gCm2[nn,mid_quant,], main="ID2", ylab=expression(paste("Wood (gC ",m^-2,")",sep="")), 
     xlab="Simulation Month", type="l", col="red", lwd=3, ylim=yrange)
lines(alt_grid_output$wood_gCm2[nn,mid_quant,], lwd=3, col="blue")
points(drivers$obs[,13], pch=16)
dev.off()

}
