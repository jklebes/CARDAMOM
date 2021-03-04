###
## Script to control the creation and operation of a CARDAMOM analysis
###

local({
## Prepare
## Compute
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")
## Print result
})

###
## Options

## Load libraries and internal functions used to manage CARDAMOM
source("./R_functions/load_all_cardamom_functions.r")

## use parallel functions?
use_parallel = FALSE
numWorkers = 6 # number of cores to assign to parallel job
## about you (for when using remote computer, e.g. Eddie)
username="lsmallma" # university uun
home_computer="ssh.geos.ed.ac.uk"
## projname
# Give a runid
projname="CARDAMOM_CAMP"

## Language
# i.e. "Fortran", "C"
language="Fortran" # Only Fortran code maintained

## Compiler options (Fortan only)
compiler="ifort" # "ifort" or "gfortan"
timing=FALSE     # 
debug=FALSE

## Model
# Currently operational "DALEC_CDEA_LU_FIRES", "DALEC_GSI_DFOL_CWD_FR", "DALEC_GSI_BUCKET"
model="DALEC_GSI_DFOL_CWD_FR"
pft_specific_parameters=FALSE # set TRUE if PFT map will be used to control some model options (e.g. use crop model)

## MDF method
# Currently operational "MHMCMC"
method="MHMCMC"

## Land cover map
# which land cover map to use, informs on land/sea mask
use_lcm="ECMWF" # choices are "CORINE2006", "LCM2007", "CORINE2006_1km", "ECMWF"
pft_wanted=FALSE # set TRUE if specific PFT values will be stored, must be true if pft_specific_parameters == TRUE

## Met paths must match with option selected from ## Data streams below unless "site_specific" has been set
path_to_met_source="/disk/scratch/local.2/crujra/"
path_to_lai="/disk/scratch/local.2/copernicus/LAI_1km_linked/"
path_to_hwsd_Csom="/exports/csce/datastore/geos/groups/gcel/HWSD/processed_file/global_0.25deg/"
path_to_forestry=" "
path_to_gleam="/exports/csce/datastore/geos/groups/gcel/GLEAM/v3.1a_soilmoisture_prior_TLS/"
path_to_site_obs="/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/example_files/"
met_interp=TRUE

## Data streams
# Drivers
met_source="site_specific" # "CRUJRA" or "ERA"
sand_clay_source="HWSD"    # "HWSD" or "site_specific"
burnt_area_source=" "      # " ", "site_specific", "GFED"
deforestation_source=" "   # " ", "site_specific" or "GFW"
# Obserations timeseries
lai_source="site_specific"     # "MODIS" or "COPERNICUS" or "site_specific"
soilwater_initial_source = " " # " ", "site_specific", "GLEAM"
Evap_source=" "                # " " or "site_specific"
woodinc_source=" "             # " " or "site_specific"
GPP_source=" "                 # " " or "site_specific"
Reco_source=" "                # " " or "site_specific"
NEE_source="site_specific"     # " " or "site_specific"
Cfol_stock_source=" "          # " " or "site_specific"
Cfolmax_stock_source=" "       # " " or "site_specific"
Cwood_stock_source=" "         # " " or "site_specific"
Cstem_stock_source=" "         # " " or "site_specific"
Cbranch_stock_source=" "       # " " or "site_specific"
Cagb_stock_source=" "          # " " or "site_specific"
Ccoarseroot_stock_source=" "   # " " or "site_specific"
Croots_stock_source=" "        # " " or "site_specific"
Clit_stock_source=" "          # " " or "site_specific"
Csom_stock_source=" "          # " " or "site_specific"
snow_source=" "
# Observations initial conditions
Cfol_initial_source=" "        # " " or "site_specific"
Cwood_initial_source=" "       # " ", "avitabile", "site_specific"
Croots_initial_source=" "      # " " or "site_specific"
Clit_initial_source=" "        # " " or "site_specific"
Csom_source="HWSD"             # " ", "HWSD", "site_specific"
crop_management_source=" "     # "_" or "site_specific" or "sacks_crop_calendar"

## sites for analysis
# start year
years_to_do=as.character(c(2001:2006)) # c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009")
# is this run "site" level or over a "grid"?
cardamom_type="site"
cardamom_grid_type=" "
# if type = "grid" then what resolution in m?
cardamom_resolution=1e5

# site names if specific locations e.g. "UKGri"
sites_cardamom="USHa1" #forest_Duke_site[6] # c("Cherokee") #tsb_site # run_sites
# lat/long of sites, if type = "grid"then these these are bottom left and top right corners
sites_cardamom_lat=c(42.54)#forest_Duke_lat[6] # 36.425914 # tsb_lat # run_sites_lat # c(52,58)
sites_cardamom_long=c(-72.17)#forest_Duke_long[6] # -81.953958 # tsb_long # run_sites_long # c(-3,0.5)
timestep_type="monthly"

## Stage 5: Driver modifications
# Define a proportional change to weather drivers.
airt_factor  = 1  # NOTE: this is also impact air temperature used in GSI model
swrad_factor = 1
co2_factor   = 1
rainfall_factor = 1
wind_spd_factor = 1
vpd_factor      = 1
# Define a proportional change on disturbance drivers.
# NOTE: at this point the only the intensity is impacted rather than frequency of events
deforestation_factor = 1
burnt_area_factor    = 1

## Stage
# stage -1 : Fix or create project first time (load source to eddie)
# stage  0 : Re-compile without re-creating project (you will be asked the stage 1 questions but these will not be saved)
# stage  1 : Create met / obs containing files for each site
# stage  2 : Submit the project to eddie / local machine
# stage  3 : (Copy back results and) process parameters into stocks and fluxes
# stage  4 : Create some standard plots
# stage  5 : Do single factor driver modification
stage=-1
repair=0 # to force (=1) re-run processed results or driver files if they already exist
grid_override=FALSE # force site specific files to be saved and figures to be generated when in "grid" operation

##
# Call CARDAMOM with specific stages
cardamom(projname,model,method,stage)
