
###
## Example script to control CARDAMOM for a single site run
###

# Set working directory in which the CARDAMOM code base can be found
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

###
## Options

## Load needed libraries and internal functions
source("./R_functions/load_all_cardamom_functions.r")

## projname
# Give a runid
#projname="FI-Hyy_example_noflx"
projname="FI-Hyy_example"
#projname="FI-Hyy_example_O3"
#projname="FI-Hyy_example_Ofast"

## Language
# i.e. "Fortran", "C"
language="Fortran"
## Compiler options (Fortan only)
compiler="ifx" #"ifort", "gfortran"
compiler_optimisation = "-O2"
timing=FALSE
debug=FALSE

## about you
username="lsmallma"
home_computer="xrdp.geos.ed.ac.uk"
sshpass_key_home = "~/.ssh/id_rsa_geos.pub" # location of passkey for home server on remote server
sshpass_key_server = "~/.ssh/id_rsa_eddie.pub" # location of passkey for remote server on home server

## Model - which DALEC 
# see "MODEL_DESCRIPTIONS.md" for available models
model="DALEC.35."
pft_specific_parameters=FALSE # impacts crop model only

## MDF method
# i.e. MHMCMC or other. Note that while MHMCMC remains the only available choice this has been updated to an AP-MCMC
method="MHMCMC"

## Land cover map
# which land cover map to use. Used to determine the land / sea mask
use_lcm="ECMWF" # coded choices exist for other maps however only "ECMWF" map is provided with source code
pft_wanted=FALSE # Impacts crop model only
path_to_landsea="default" # If gridded analysis, any raster layer with >0 values will be taken as the mask area. To ignore = "default"

## Met paths
#path_to_met_source="/exports/geos.ed.ac.uk/gcel/spatial_datasets/meteorology/ERA5/global_0.1deg/"
path_to_met_source="/exports/geos.ed.ac.uk/gcel/spatial_datasets/meteorology/trendy/version_13/global_0.5deg/"
path_to_lai = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/leaf_area_index/MCD15A2H/collection_6.1/global_0.0625deg/"
path_to_fapar = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/fraction_absorbed_par/MCD15A2H/collection_6.1/global_0.0625deg/"
path_to_crop_management=" "
path_to_sand_clay="/exports/geos.ed.ac.uk/gcel/spatial_datasets/soil_texture/SoilGrids_v2/global_0.125deg/"
path_to_Csom = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/soil_carbon/SoilGrids_v2/global_0.0625deg/"
path_to_Cwood_inc = " "
path_to_Cwood_mortality = " "
#path_to_Cwood = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/total_woody_biomass/Xu2021/all/"
path_to_Cwood = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/total_woody_biomass/esa_cci/v6/global_0.5deg/"
path_to_Cwood_initial=" "
path_to_Cwood_potential=" "
path_to_soilwater=" "
path_to_nbe = "/home/lsmallma/gcel/AtmosphericInversions/OCO2v10_MIP/LNLG/"
path_to_gpp = "/home/lsmallma/gcel_ceph/spatial_datasets/gross_primary_production/MOD17A3HGF/collection_6.1/global_0.5deg/"
path_to_fire = "/exports/csce/datastore/geos/groups/gcel/FIRE_ESTIMATES/combined_fire/global_1deg_monthly/"
path_to_forestry="/exports/geos.ed.ac.uk/gcel/spatial_datasets/deforestation/GlobalForestWatch_FireRemoved/global_0.5deg/"
#path_to_burnt_area="/exports/geos.ed.ac.uk/gcel/spatial_datasets/burned_area/MCD64A1/collection_6.1/global_0.5deg/"
path_to_burnt_area="/exports/geos.ed.ac.uk/gcel/spatial_datasets/burned_area/GFED/v5.1/global_0.5deg/"
path_to_lca = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/LCA/Butler/global_0.5deg/"
path_to_landsea = "default"
path_to_co2 = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/meteorology/trendy/version_13/"
path_to_site_obs="./example_files/inputs/"
met_interp=FALSE

## Data streams
met_source="trendy" # "ERA" or "isimip3a" or "trendy"
lai_source="Gridded_nc" # "Gridded_nc" or "Gridded_tif" or "site_specific"
fapar_source="Gridded_nc" # "Gridded_nc" or "Gridded_tif" or "site_specific"
Csom_source="Gridded_tif" # "Gridded_nc" or "Gridded_tif" or "site_specific"
soilwater_source = " " # initial soil water fraction (m3/m3)
sand_clay_source="Gridded_tif" # "Gridded_nc" or "Gridded_tif" or "site_specific
et_source=" "
Cwood_inc_source = " " # " " or "site_specific" or "Rainfor"
Cwood_mortality_source = " " # "site_specific" or " " or "Rainfor"
gpp_source="Gridded_tif" 	# " " or "site_specific"
fire_source=" " # " " or "site_specific" or "Global_Combined"
Reco_source=" " 	# " " or "site_specific"
nee_source=" " # "Gridded_nc" or "Gridded_tif" or "site_specific"
nbe_source = " " # "Global_Combined" or "GEOSCHEM"
harvest_source = ""
foliage_to_litter_source = " " # " " or "site_specific"
# i.e. single value valid for beginning of simulation
Cfol_initial_source=" " #"site_specific" 	# " " or "site_specific"
Cwood_initial_source=" " #"site_specific" 	# " " or "site_specific"
Croots_initial_source=" " #"site_specific" 	# " " or "site_specific"
Clit_initial_source=" " #"site_specific"  	# " " or "site_specific"
# i.e. time series of stock estimates
Cfol_stock_source=" " 	# " " or "site_specific"
Cfolmax_stock_source=" " 	# " " or "site_specific"
Cwood_stock_source="Gridded_tif" 	# " " or "site_specific" or "Saatchi_2021" or "ESA_CCI_Biomass"
Cstem_stock_source=" "      # " " or "site_specific"
Cbranch_stock_source=" "      # " " or "site_specific"
Cagb_stock_source=" " 	# " " or "site_specific"
Ccoarseroot_stock_source=" " 	# " " or "site_specific"
Croots_stock_source=" " 	# " " or "site_specific"
Clit_stock_source=" "  	# " " or "site_specific"
Csom_stock_source=" "  	# " " or "site_specific"
# Parameter priors
lca_source = "Gridded_tif" # "Gridded_nc" or "Gridded_tif" or "site_specific"
frac_Cwood_coarse_root_source = "" # " " or "site_specific"
minLWP_source = "" # " " or "site_specific"
# Steady state attractor
Cwood_potential_source = " " # "site_specific" or ""
# Management drivers
burnt_area_source = "Gridded_nc" # "Gridded_nc" or "Gridded_tif" or "site_specific"
deforestation_source = "Gridded_tif" # "Gridded_nc" or "Gridded_tif" or "site_specific"
crop_management_source = " " # "_" or "site_specific" or "sacks_crop_calendar"
snow_source=" "

## sites for analysis
# start year
#years_to_do=as.character(c(1999:2014)) 
years_to_do=as.character(c(2003:2023)) 
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
select_country = FALSE # If gridded run and path_to_landsea = "default", 
                       # select country based on site_cardamom. Use function
                       # available_countries() for compatible country names.

## Define the project setup
# NOTE: if these are not set CARDAMOM will ask you for them
# Some interactive node options
use_parallel=FALSE             # use parallel functions or not
numWorkers=4                   # number of parallel tasks when using a interactive node
# Some Slurm server option
slurm_account = "geos_research"# Slurm research account, if using the slurm cluster
slurm_concurrent_cpus = 60     # maximum number of concurrent cpus for slurm, impacts stage 3
slurm_max_run_time = 12        # Number of hours per task to be requested in stage 3, if using slurm
# Control where to run
request_use_server = FALSE     # Use remote server? Currently coded for UoE Eddie.
request_use_local_slurm = FALSE # Only applies if request_use_server == FALSE
request_compile_server = FALSE # Copy and compile current source code on remote server
request_compile_local = TRUE   # Compile local executable even if not running on local
# Remote server options
request_runtime = 48           # How many hours of compute to request per job for stage 2. For Slurm and remote server. 
# MCMC specific options
request_nos_chains = 3         # Number of chains CARDAMOM should run for each location
request_nos_samples = 100e6    # Total number of parameter samples / iterations to be explored
request_nos_subsamples = 1e3   # Number of parameter sets to be sub-sampled from the chain
request_use_EDCs = TRUE        # Use EDCs
request_extended_mcmc = FALSE  # Extend the current MCMC by adding a further request_nos_extended_samples + request_nos_samples
request_nos_extended_samples = 40e6 # If request_extened_mcmc == TRUE then this is the number of additional proposals to be made
request_cost_function_scaling = 0 # 0 = Default, no normaliation of the likelihood score
                                  # 1 = Normaliation of the likelihood score by sample size
                                  # 2 = Normaliation of the likelihood score by sqrt(sample size)
                                  # 3 = Normaliation of the likelihood score by log(sample size)  

## Stage
# stage -1 : Create project first time (load source to eddie)
# stage  0 : Re-compile source code, if needed but without re-creating the PROJECT related infofile.RData
# stage  1 : Create met / obs containing files for the specifc project
# stage  2 : Submit the project to eddie
# stage  3 : Copy back results and process vectors
# stage  4 : Do some standard figure creation (and further processing for gridded analysis)
# stage  5 : Generic dump of RESULTS_PROCESSED files to netcdf
stage=3
repair=1 # to force (=1) re-run processed results or driver files if they already exist
grid_override=FALSE # force site specific files to be saved and figures to be generated when in "grid" operation

##
# Call CARDAMOM with specific stages
cardamom(projname,model,method,stage)


