
###
## Example script to control CARDAMOM for a single site run
###

# Set working directory in which the CARDAMOM code base can be found
setwd("<Path to where CARDAMOM has been placed>")

###
## Options

## Load needed libraries and internal functions
source("./R_functions/load_all_cardamom_functions.r")

## projname
# Give a runid
projname="FI-Hyy_example"

## Language
# i.e. "Fortran", "C"
language="Fortran"
## Compiler options (Fortan only)
compiler="ifort" #"ifort", "gfortran"
timing=FALSE
debug=FALSE

## about you (only valid if working on UoE remote server)
username=" " # put your Edinburgh uun here
home_computer=" "

## use parallel functions?
use_parallel = FALSE
numWorkers = 6 # number of cores to assign to parallel job

## Model - which DALEC 
# i.e. "DALEC_CDEA_LU_FIRES","DALEC_1005","DALEC_1005a","DALEC_CDEA_ACM2_BUCKET","DALEC_CDEA_ACM2_BUCKET_RmRg", "DALEC_CDEA_ACM2_BUCKET_CWD", "DALEC_GSI_DFOL_CWD_FR", "DALEC_GSI_BUCKET"
model="DALEC_CDEA_LU_FIRES"
pft_specific_parameters=FALSE # impacts crop model only

## MDF method
# i.e. MHMCMC or other. Note that while MHMCMC remains the only available choice this has been updated to an AP-MCMC
method="MHMCMC"

## Land cover map
# which land cover map to use. Used to determine the land / sea mask
use_lcm="ECMWF" # coded choices exist for other maps however only "ECMWF" map is provided with source code
pft_wanted=FALSE # Impacts crop model only

## Data paths - assumes that R code has been updated to deal with various gridded datasets
path_to_site_obs="./example_files/inputs/" # If not using gridded data then only the site obs path needs setting
path_to_met_source=" "
path_to_lai=" "
path_to_crop_management=" " # crop model only
path_to_sand_clay=" " # sand / clay fractions for ACM2 hydraulics model
path_to_Csom=" " # Currently just applies to soil C on initial conditions parameter
path_to_forestry=" " # i.e. forest harvest / management
path_to_Cwood=" " # Time specific wood stock information (above + below)
path_to_Cwood_initial=" " # Constraint on initial wood stock parameter
path_to_Cwood_potential=" " #  Constraint on the steady state attractor for wood stocks
path_to_gleam=" " # Evapotranspiration constraint
path_to_nbe = " " # Net biome exchange
path_to_burnt_area=" " # imposed as a fraction
path_to_lca = " " # leaf carbon per unit area (gC/m2)
met_interp=TRUE # apply linear interpolation to met drivers if the provided datasets are not correct

## Data streams - The currently coded data streams which can be used to drive or constrain the models
met_source="site_specific" # "trendy_v9" or other options coded into the R wrapper
lai_source="site_specific" # "COPERNICUS" or "site_specific"
Csom_source="site_specific" # "SoilGrids" or "site_specific"
soilwater_initial_source = " " # initial soil water fraction (m3/m3)
sand_clay_source="site_specific" # "SoilGrids" or "site_specific
Evap_source=" "        # " " or "site_specific"
woodinc_source=" " 	# " " or "site_specific"
GPP_source=" " 	# " " or "site_specific"
Reco_source=" " 	# " " or "site_specific"
NEE_source="site_specific" # "site_specific" 	# " " or "site_specific"
nbe_source = " "
# i.e. single value valid for beginning of simulation
Cfol_initial_source=" " #"site_specific" 	# " " or "site_specific"
Cwood_initial_source=" " #"site_specific" 	# " " or "site_specific"
Croots_initial_source=" " #"site_specific" 	# " " or "site_specific"
Clit_initial_source=" " #"site_specific"  	# " " or "site_specific"
# i.e. time series of stock estimates
Cfol_stock_source=" " 	# " " or "site_specific"
Cfolmax_stock_source=" " 	# " " or "site_specific"
Cwood_stock_source="site_specific" 	# " " or "site_specific"
Cstem_stock_source=" "      # " " or "site_specific"
Cagb_stock_source=" " 	# " " or "site_specific"
Ccoarseroot_stock_source=" " 	# " " or "site_specific"
Croots_stock_source=" " 	# " " or "site_specific"
Clit_stock_source=" "  	# " " or "site_specific"
Csom_stock_source=" "  	# " " or "site_specific"
lca_source = " " # "Butler" or " " or "site_specific"
# Steady state attractor
Cwood_potential_source = " " # "site_specific" or ""
# Management drivers
burnt_area_source=" "
deforestation_source=" " # " ", "site_specific" or "maryland_forestry_commission"
crop_management_source=" " # "_" or "site_specific" or "sacks_crop_calendar"
snow_source=" "

## sites for analysis
# start year
years_to_do=as.character(c(1999:2014)) 
# is this run "site" level or over a "grid"?
cardamom_type="site"
cardamom_grid_type=" " # "UK" or "wgs84", no value needed if site run
# if type = "grid" then at what spatial resolution (UK = m, wgs84 = degree)?
cardamom_resolution=1e5

# site names if specific locations e.g. "UKGri"
sites_cardamom=c("FI-Hyy")
# lat/long of sites, if type = "grid" then these these are bottom left and top right corners
sites_cardamom_lat=61.84741
sites_cardamom_long=24.29477
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

## Define the project setup
# NOTE: if these are not set CARDAMOM will ask you for them
request_nos_chains = 3         # Number of chains CARDAMOM should run for each location
request_nos_samples = 100e6    # Total number of parameter samples / iterations to be explored
request_nos_subsamples = 1e3   # Number of parameter sets to be sub-sampled from the chain
request_use_server = FALSE     # Use remote server? Currently coded for UoE Eddie.
request_runtime = 48           # How many hours of compute to request per job. Only applied for running on remote server
request_compile_server = FALSE # Copy and compile current source code on remote server
request_use_EDCs = TRUE        # Use EDCs

## Stage
# stage -1 : Create project first time (load source to eddie)
# stage  0 : Re-compile source code, if needed but without re-creating the PROJECT related infofile.RData
# stage  1 : Create met / obs containing files for the specifc project
# stage  2 : Submit the project to eddie
# stage  3 : Copy back results and process vectors
# stage  4 ; Do some standard result checking
# stage  5 ; Currently out of use
stage=-1
repair=1 # to force (=1) re-run processed results or driver files if they already exist
grid_override=FALSE # force site specific files to be saved and figures to be generated when in "grid" operation

##
# Call CARDAMOM with specific stages
cardamom(projname,model,method,stage)


