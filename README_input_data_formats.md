## Readme: Input data formats

CARbon DAta MOdel fraMework (CARDAMOM) and the DALEC suite of models require large amounts of information,
across a large range of datasets. It is impractical to continuously develop the code base to deal with 
an ever growing number of potential data formats. Instead three file format situations have been decided on.

### Single site comma seperated (csv) files 

A combination of three comma separated (csv) files defining the time series of forcings (meteorology, disturbance),
time varying assimilatable observations (e.g. leaf area index), and time invariant information (e.g. initial conditions, soil texture).
This file format is expected to be used for individual site analyses, heavily drawing on in-situ observations rather than
from spatial datasets such as those from meterological reanalyses or satellite based Earth Observation.

The file containing forcings has the following naming convention: <site_name>_timeseries_met.csv
Where the <site_name> is to match the site name given in the control_*.r file used to control the CARDAMOM analysis.

Assumptions:
1) The number of time steps in the file matches that requested in the control file.
2) That all column should be provided, best practice, but those not needed for your 
   specific model version can be set to -9999 for each time step.
3) Missing values are not allowed in the meteorological variable.
4) Missing values in the disturbance related forcings uses -9999 flag.

Variable naming convention:
doy = Decimal Day Of Year (doy) of the analysis time step. If time step is greater than 1 day then the value is assumed to be
      the day of year at the end of the time step. e.g. the first value for a monthly time step analysis is 
      January = 31, February = 59 (i.e. 31+28).
mint_C = Minimum daily temperature (Celsius). If time step is greater than 1 day, then this is assumed to be to the average 
         of the minimum temperature of each day within the step.
maxt_C = Maximum daily temperature (Celcius). If time step is greater than 1 day, then this is assumed to be to the average 
         of the maximum temperature of each day within the step.
airt_C = Average daily temperature (Celsius). If time step is greater than 1 day, then this is assumed to be to the average 
         of the average temperature of each day within the step.
swrad_Wm2 = Daily average downwelling short wave radiation (W/m2), if time step is greater than 1 day, then this is assumed 
            to be the average of each day within the time step. Note only one of swrad_Wm2 or swrad_MJm2day is required.
swrad_MJm2day = Daily average downwelling short wave radiation (MJ/m2/day), if time step is greater than 1 day, then this is assumed 
                to be the average of each day within the time step. Note only one of swrad_Wm2 or swrad_MJm2day is required.            
co2_ppm = Atmospheric CO2 concentration (ppm), averaged over the time step.
precip_kgm2s = Liquid + ice precipitation (kgH2O.m-2.s-1), averaged over the time step.
vpd_kPa = Vapour pressure deficit (kPa), averaged over the time step. 
rh_fraction = Relative Humidity (0-1). This can be an alternative to the vpd_kPa but not recommended. RH is converted in the code
              based on the daily information provided. However, this will introduce an error compared to estimating from sub-daily.
wind_spd_ms = Wind speed (m/s) averaged over the time step.

burnt_area_fraction = Fraction of the surface burned by fire in the time step. Note that while a value of 1 
                      indicates the whole area was burned. It does not imply a severity, combustion completeness or mortality.
deforestation_fraction = The fraction of the surface cleared. Note that this does not mean that all biomass is removed, just that 
                        'management_type' has been carried out across the area.
management_type = Catagorical definition of what type of management activity was imposed over the deforestation_fraction area.
                  The default value is = 2. But for specific definitions see your specific model code.
lai_loss = Used for models .M2. only, directly forces a leaf area reduction (m2/m2). 

The file containing time varying observations to assimilate has the following naming convention: <site_name>_timeseries_obs.csv
Where the <site_name> is to match the site name given in the control_*.r file used to control the CARDAMOM analysis.

Assumptions:
1) The number of time steps in the file matches that requested in the control file.
2) That all column should be provided, best practice, but those not needed for your 
   specific model version can be set to -9999 for each time step.
3) Missing values uses -9999 flag.
4) The uncertainty value given is intended to be the best estimate of the observation uncertainty. 
   A model structural uncertainty is added in the extract_obs.r function which will be used in the actual analysis.
5) That the user has checked their specific model code (MODEL_LIKELIHOOD.f90) that their variable of interest is 
   currently coded into the log-likelihood cost function.

Variable naming convention:
NBE_gCm2day = Net Biome Exchange of CO2 (NBE, gC/m2/day), averaged across the model time step. Negative is a net sink of CO2.
              This flux is the net of GPP, Reco, and fire, but excludes lateral C transfer from mechanical disturbance (e.g. deforestation).  
NBE_unc_gCm2day = Uncertainty for Net Biome Exchange of CO2 (gC/m2/day).  
LAI_m2m2 = Leaf area index (LAI, m2/m2), averaged across the model time step. LAI is the 1-sided estimate of leaf area per unit ground area.  
LAI_unc_m2m2 = Uncertainty for the leaf area index (m2/m2).  
fAPAR_fraction = Fraction of Absorbed Photosynthetically Active Radiation (fAPAR, 0-1), averaged across the model time step.  
fAPAR_unc_fraction = Uncertainty for the fAPAR.  
Cfol_stock_gCm2 = Carbon stock in the foliage (gC/m2), averaged across the model time step.  
Cfol_stock_unc_gCm2 = Uncertainty for th Cfol_stock_gCm2.  
harvest_gCm2day = Estimate of carbon extracted as harvested biomass (gC/m2/day), averaged over the harvest_lag_step period. 
                  This could, subject to code check, be applied to deforestation, arable harvest and grassland management.  
harvest_uncertainty_gCm2day = Uncertainty for the harvest_gCm2day (gC/m2/day).  
harvest_lag_step = Number of model time steps over which to lag (average) the calculation of harvest_gCm2day.  
Cwood_increment_gCm2day = Estimate of carbon allocated to the wood pool (gC/m2/day), averaged of the Cwood_increment_lag_step period.  
Cwood_increment_uncertainty_gCm2day = Uncertainty for the Cwood_increment_gCm2day (gC/m2/day).  
Cwood_increment_lag_step = Number of model time steps over which to lag (average) the calculation of Cwood_mortality_gCm2day.  
Cwood_mortality_gCm2day = Estimate of carbon lost from the wood pool (gC/m2/day) from all sources, averaged of the Cwood_mortality_lag_step period.  
Cwood_mortality_uncertainty_gCm2day = Uncertainty for the Cwood_mortality_gCm2day (gC/m2/day).  
Cwood_mortality_lag_step = Number of model time steps over which to lag (average) the calculation of Cwood_mortality_gCm2day.  
foliage_to_litter_gCm2day = Estimate of carbon lost from the foliage pool to litter pool (gC/m2/day) from all sources, 
                            averaged of the foliage_to_litter_lag_step period.  
foliage_to_litter_unc_gCm2day = Uncertainty for the foliage_to_litter_gCm2day (gC/m2/day).  
foliage_to_litter_lag_step = Number of model time steps over which to lag (average) the calculation of foliage_to_litter_gCm2day.  
GPP_gCm2day = Gross Primary Production (GPP, gC/m2/day), i.e. canopy scale photosynthesis, averaged across the model time step. 
              Positive is CO2 uptake by the plant.  
GPP_unc_gCm2day = Uncertainty for gross primary production (gC/m2/day).  
fire_gCm2day = Carbon emissions to the atmosphere from fire (GPP, gC/m2/day). Note this variable does not discriminate between different carbon species and 
               does not include fire driven mortality. Emissions < 0.01 gC/m2/day are automatically removed to prevent biasing the log-likelihood calculation.
               Finally, this fix is dependent on there being a consistent burnt_area_fraction variable provided as there are currently no DALEC models with a fire generation sub-model.  
fire_unc_gCm2day = Uncertainty for carbon emissions to the atmosphere from fire (gC/m2/day).     
Evap_kgH2Om2day = Evapotranspiration (ET, kgH2O/m2/day), i.e. the combined total of transpiration,
                  soil evaporation, evaporation of canopy intercepted rainfall and where appropriate snow sublimation.
                  Averaged across the model time step. Positive is evaporation from the surface to the atmosphere.  
Evap_unc_kgH2Om2day = Uncertainty for Evapotranspiration (kgH2O/m2/day).  
Reco_gCm2day = Ecosystem respiration (Reco, gC/m2/day), i.e. the combination of plant and heterotrophic respiration, 
               averaged across the model time step. Check your specific model for details as to which terms are being 
               explicitly simulated. Positive is CO2 released from the ecosystem to the atmosphere.   
Reco_unc_gCm2day = Uncertainty for ecosystem respiration (gC/m2/day).  
NEE_gCm2day = Net Ecosystem Exchange of CO2 (NEE, gC/m2/day), averaged across the model time step. Negative is a net sink of CO2.
              This flux is the net of GPP and Reco, but excludes fire and lateral C transfer from mechanical disturbance (e.g. deforestation).  
NEE_unc_gCm2day = Uncertainty for Net Ecosystem Exchange of CO2 (gC/m2/day).  
Cwood_stock_gCm2 = Carbon stock in the wood (gC/m2), averaged across the model time step. Wood here is the combination of both above and below ground.  
Cwood_stock_unc_gCm2 = Uncertainty for Cwood_stock_gCm2.  
Cagb_stock_gCm2 = Carbon stock in the above ground wood biomass (gC/m2), averaged across the model time step. Note this is only available for some model versions
                  as it relies on having a parameter to partition the below from above ground wood which is only weakly constrained.  
Cagb_stock_unc_gCm2 = Uncertainty for the Cagb_stock_gCm2.  
Croots_stock_gCm2 = Carbon stock in the fine roots (gC/m2), averaged across the model time step. Note that coarse roots are assumed to be part of the wood stock.  
Croots_stock_unc_gCm2 = Uncertainty for Croots_stock_gCm2.  
Clit_stock_gCm2 = Carbon stock in the fine litter (foliage + fine root, gC/m2), averaged across the model time step. 
                  Note this is for foliage and fine root litter only. Other arrangements are made for wood.  
Clit_stock_unc_gCm2 = Uncertainty for Clit_stock_gCm2.  
Csom_stock_gCm2 = Carbon stock in the soil organic matter (gC/m2), averaged across the model time step.   
Csom_stock_unc_gCm2 = Uncertainty for Csom_stock_gCm2.  
Ccoarseroot_stock_gCm2 = Carbon stock in the coarse root (gC/m2), averaged across the model time step. Note this is only available for some model versions
                  as it relies on having a parameter to partition the below from above ground wood which is only weakly constrained.
                  Below ground wood is assumed to be coarse root.  
Ccoarseroot_stock_unc_gCm2 = Uncertainty for Ccoarseroot_stock_gCm2.  
Cfolmax_stock_gCm2 = Maximum carbon stock in the foliage pool each year (gC/m2). Matches the maximum value of foliage carbon in each year an observations is given.
                     This variable supports the use of seasonal maximum foliage information reported in literature but often without an explicit date. It also allows 
                     greater flexibility when mixing satellite based Earth Observation of leaf area with in-situ information which may have different temporal dynamics.   
Cfolmax_stock_unc_gCm2 = Uncertainty for Cfolmax_stock_gCm2.  
snow_water_kgH2Om2 = Water equivalent storage in snow (kgH2O/m2 or mm), averaged across the model time step. Available only in some models. 
                     Please check your specific code.  
snow_water_unc_kgH2Om2 = Uncertainty for snow_water_kgH2Om2.  

The file containing three types of information. First is value which map directly onto a model's parameters (e.g. LCA), the second is time-invarient observations 
of emergent ecosystem properties (WUE), third are time invarient forcings (e.g. soil texture). The file has the following naming convention: <site_name>_initial_obs.csv.
Where the <site_name> is to match the site name given in the control_*.r file used to control the CARDAMOM analysis.

Assumptions:
1) That all column should be provided, best practice, but those not needed for your 
   specific model version can be set to -9999 for each time step.
2) Missing values uses -9999 flag.
3) The uncertainty value given is intended to be the best estimate of the observation uncertainty. 
   A model structural uncertainty is added in the extract_obs.r function which will be used in the actual analysis.
4) That the user has checked their specific model code (MODEL_LIKELIHOOD.f90) that their variable of interest is 
   currently coded into the log-likelihood cost function.

Variable naming convention starting with assimilated variables, followed by forcings:
Csom_initial_gCm2 = Prior estimate on the initial soil carbon stock (gC/m2). Nominal assumption is this is carbon to 1 m soil depth. 
                    But is actually determined by the representative depth of the observations.  
Csom_initial_unc_gCm2 = Uncertainty for Csom_initial_gCm2 (gC/m2).  
Cfol_initial_gCm2 = Prior estimate on the initial foliage carbon stock (gC/m2).   
Cfol_initial_unc_gCm2 = Uncertainty for Cfol_initial_gCm2 (gC/m2).  
Cwood_initial_gCm2 = Prior estimate on the initial above and below wood carbon stock (gC/m2).   
Cwood_initial_unc_gCm2 = Uncertainty for Cwood_initial_gCm2 (gC/m2).  
Croots_initial_gCm2 = Prior estimate on the initial fine roots carbon stock (gC/m2).   
Croots_initial_unc_gCm2 = Uncertainty for Croots_initial_gCm2 (gC/m2).  
Clit_initial_gCm2 = Prior estimate on the initial foliage and fine root litter carbon stock (gC/m2).   
Clit_initial_unc_gCm2 = Uncertainty for Clit_initial_gCm2 (gC/m2).  
soil_water_fraction = Prior estimate on the initial soil water fraction (0-1) of the top soil layer (0-30 cm).
                      This applies only to some model, please check you appropriate source code before use.  
soil_water_unc_fraction = Uncertainty for soil_water_fraction (0-1).  
Cwood_potential_gCm2 = Prior estimate on the steady state estimate of above and below groun wood carbon stock (gC/m2). 
                       This applies only to some model, please check you appropriate source code before use.  
Cwood_potential_unc_gCm2 = Uncertainty for Cwood_potential_gCm2 (gC/m2).  
LCA_gCm2 = Prior estimate on the leaf carbon per unit area (gC/m2), the carbon equivalent of the leaf mass area.
           This applies only to some model, please check you appropriate source code before use.  
LCA_unc_gCm2 = Uncertainty for LCA_gCm2 (gC/m2).  
frac_Cwood_coarse_root_prior = Prior estimate on the fraction of wood stock which is below ground, i.e. assumed to be coarse root.
                               This applies only to some model, please check you appropriate source code before use.
                               Note that were a model uses this parameter and a value is not given here, a default value will be 
                               calculated following the allometric equation from
                               Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899  
frac_Cwood_coarse_root_prior_unc = Uncertainty for frac_Cwood_coarse_root_prior (0-1).  
minLWP_MPa = Prior estimate on the minimum tolerated leaf water potential (MPa).
             This applies only to some model, please check you appropriate source code before use.  
minLWP_unc_MPa = Uncertainty for minLWP_MPa (MPa).  
planting_doy_initial = Prior estimate on the day of planting (day of year). Intended for use with agriculture model version DALEC.14 and DALEC.15.  
planting_doy_unc_initial = Uncertainty for planting_doy_initial (gC/m2).  
growing_season_doy_initial = Prior estimate on the duration of the growing season after planting (days).
                             Intended for use with agriculture model version DALEC.14 and DALEC.15.  
growing_season_doy_unc_initial = Uncertainty for growing_season_doy_initial (gC/m2).  

top_sand_initial_percent = Percentage content of sand in the 0-30cm soil layer. Note not used in all models.  
bot_sand_initial_percent = Percentage content of sand in the 31-100cm soil layer. Note not used in all models.  
top_clay_initial_percent = Percentage content of clay in the 0-30cm soil layer. Note not used in all models.  
bot_clay_initial_percent = Percentage content of clay in the 31-100cm soil layer. Note not used in all models.  

### Meteorological datasets

Meteorological information is either available from in-situ information (thus provided in csv format, described above) or
drawn from large scale reanalyses (or a combination of spatial datasets). Meteorological information is mandatory to run 
CARDAMOM with information for every time step and location across at least 6 different variable. As there are a large number
of variables a specific section is given here describing the file format requirements.

In contrast to the assimilated observation, described in the 'Spatial datasets' section and subsections, the only permitted
file format is netCDF with strict file name and in file variable conditions. While there is some commonality with requirements 
of other files these are given here for clarity.

The file format...

THIS NEEDS TO BE DECIDED ON - THE CURRENTLY THERE ARE TWO CHOICES WHICH ARE IMPLICIT ON THEIR ASSUMPTIONS OR EITHER MONTHLY OR DAILY 
TIMESTEPS. NEITHER IS A GOOD FINAL CHOICE AND ANY CHOICE WILL RESULT IN REWRITING SIGNIFICANT NUMBER OF FILES FOR THIS TO WORK...

### Spatial datasets - general information

A large amount of information used in CARDAMOM is drawn from existing datasets potentially covering the whole globe for many years. 
Commonly used file formats are netCDF and geotifs, these formats are optimised for storing spatial data efficiently with trade-offs between
these two format choices. Therefore, we support both netCDF and geotif formats so long as the files are formatted to specific criterion. 
This will involve processing your source datasets into the formats outlind below. CARDAMOM is designed to read in data which is either 
temporally static or time varying, these two cases have slighty different formatting requirements. 

Below there will be three sections, the first outlines the format requirements for netcdf files, geotifs are described in the second section. 
The third section provides a list of available data types which can be read along with their corresponding specific file and variable naming differences.

Assumptions:
1) Missing data should be specified as NA or -9999.
2) Formatting must be carried out case perfect or CARDAMOM will not be able to read the files correctly.

#### Spatial datasets - using netcdf

NetCDF files have become a default format option for storing large complicated multi-dimensional datasets. 
NetCDF files can be used in CARDAMOM to provide both time varying and static information. To simplfy the CARDAMOM 
code we have strict formating requirements.

NetCDF files for time varying information have the following file name convention: <variable_name>_<unit>_YYYY.nc
NetCDF files for static information have the following file name convention: <variable_name>_<unit>.nc
In each case <variable_name> and <unit> are replaced by a combination currently coded for, see the below section which
provides a list of availble data types.

The files must each store information describing the latitude and longitude, and for time varying data a time variable.

Mandatory variables are:
doy = the day of year on which the data are valid (1-366). Using this information the code will assign each observation 
      to the correct model time step. Note not used for static information
lat = y-dimension information. The default expectation is that this will be the latitude (-90/90) using EPSG:4326.
lon = x-dimension information. The default expectation is that this will be the longitude (-180/180) using EPSG:4326.

Optional variables are (default assumptions will be made in their absence):
<variable>_lag_days = the number of days prior to the given doy that the observation is averaged over. 
                      The code will convert the number of days into the model time steps requested.
Description = Use the 'global attribute' option to add a description of the processing of the file and who is responsible. 
              The description should include the relevent EPSG projection code so that CARDAMOM knows for certain whether 
              it needs to reproject. Recommended that the code is given e.g 'epsg: 4326'.

Recommended variables are (useful for you):
year = in file variable giving the year of the observation.
cover_fraction = the fraction of the pixel the variable is representative of. This variable allows for the variable 
                 to be provided at per m2 but scalable to the pixel scale based on the pixel area * cover_fraction.
                 This is useful for land seas mask, or data covering specific land covers. Note that this value is assumed
                 to be temporally invarient, i.e. its dimension for any given year should vary in space only.

See also the below section for the list of required names for the in file variable names.

Assumptions:
1) That each netCDF file contains information for a single year only.
2) Unless explicitly given, that all observations represent their given time step, i.e. that their lag period = 1.
3) Each assimilated observation has a corresponding uncertainty variable.
4) Units cannot be changed from those specified below, these are expected by CARDAMOM.

#### Spatial datasets - using geotifs

Geotif files have become a default option for storing large spatial datasets. They have the advantage that GIS users will 
also have experience using them, which is less likely than netCDF. There are however trade-offs. It is less straightforward
to use geotifs to store time varying information and additional meta-data.  
To simplfy the CARDAMOM code we have strict formating requirements.

Geotif files for time varying information have the following file name convention: <variable_name>_<unit>_YYYY-DOY.tif
Geotif files for static information have the following file name convention: <variable_name>_<unit>.tif
In each case <variable_name> and <unit> are replaced by a combination currently coded for, see the below section which
provides a list of availble data types. YYYY is the calendar year (e.g. 2023) and DOY is the day of year (001-366).

For each of the data files there must be a corresponding uncertainty file with the following name conventions.
Geotif files for time varying information have the following file name convention: <variable_name>_<unit>_uncertainty_YYYY-DOY.tif
Geotif files for static information have the following file name convention: <variable_name>_<unit>_uncertainty.tif

See also the below section for the list of required names for the in file variable names.

Assumptions:
1) That each file contains information for a single time step only.
2) No lag period can be provided, so assumption is that the lag is = 1.
3) Each assimilated observation has a corresponding uncertainty in a uncertainty file.

#### Available data types and variable names

Some description for what NA means for forcing only information.

NEED TO ADD LAG VARIABLE IN HERE SSON AS THE CODE IS UPDATED TO REFLECT THIS...
COULD ALSO INCLUDE A THIRD ASSUMED GEOTIF FILE _LAGDAYS_?

|        <variable_name>        |   <unit>  | infile variable name | infile uncertainty name | Forcing / Assimilated |
|-------------------------------|-----------|----------------------|-------------------------|-----------------------|
| "forest_loss"                 | NA        | "forest_loss"        | NA                      | Forcing               |
| "BurnedFraction"              | NA        | "BurnedFraction"     | NA                      | Forcing               |
| "leaf_area_index"             | "m2m2"    | "LAI"                | "LAI_SD"                | Assimilated           |
| "fraction_absorbed_par"       | NA        | "fAPAR"              | "fAPAR_SD"              | Assimilated           |
| "net_biome_exchange"          | "gCm2day" | "NBE"                | "NBE_SD"                | Assimilated           |
| "gross_primary_production"    | "gCm2day" | "GPP"                | "GPP_SD"                | Assimilated           |
| "fire_carbon_emissions"       | "gCm2day" | "Fire"               | "Fire_SD"               | Assimilated           |
| "wood_stock"                  | "gCm2"    | "wood_stock"         | "wood_stock_SD"         | Assimilated           |
| "wood_stock_production"       | "gCm2day" | "wood_production"    | "wood_production_SD"    | Assimilated           |
| "wood_stock_mortality"        | "gCm2day" | "wood_mortality"     | "wood_mortality_SD"     | Assimilated           |
| "sand_percent_mean_0to30cm"   | NA        | "sand_content"       | "sand_content_unc"      | Forcing               |
| "sand_percent_mean_30to100cm" | NA        | "sand_content"       | "sand_content_unc"      | Forcing               |
| "clay_percent_mean_0to30cm"   | NA        | "clay_content"       | "clay_content_unc"      | Forcing               |
| "clay_percent_mean_30to100cm" | NA        | "clay_content"       | "clay_content_unc"      | Forcing               |
| "soil_stock"                  | "gCm2"    | "soil_stock"         | "soil_stock_SD"         | Assimilated           |
| "leaf_carbon_area"            | "gCm2"    | "leaf_carbon_area"   | "leaf_carbon_area_SD"   | Assimilated           |

  
