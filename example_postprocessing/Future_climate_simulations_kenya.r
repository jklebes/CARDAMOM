
###
## Load CARDAMOM analysis and run for future climate scenarios
###

# Until when do we want to loop our meteorology to?
new_end_year = 2100
# State the quantiles wanted in the output (don't be greedy)
quantiles_wanted = c(0.025,0.05,0.25,0.50,0.75,0.95,0.975)
mid_quant = 0.50 ; low_quant = 0.05 ; high_quant = 0.95
# State the climate scenarios wanted
climate_scenarios = c("ssp126","ssp434","ssp245","ssp370","ssp585")

# State the climate models wanted
ESM = c("MOHC")
# Path to CMIP6 files
cmip6dir = "/home/lsmallma/gcel/cmip6/TLS/ScenarioMIP/"
# Path to land use scenarios (LUH2)
luh2dir = "/home/lsmallma/gcel/LUH2/"
# Load libraries and CARDAMOM R code
setwd("~/WORK/GREENHOUSE/models/CARDAMOM/")
source("./R_functions/load_all_cardamom_functions.r")
library(RcppRoll)

##
# Load project information

output_prefix = "Kenya_" ; output_suffix = "climate_change_"

model_choice = "C1" ; potAGB = FALSE
if (model_choice == "C1") {

    if (potAGB) {

    } else { 
        # Load nopotAGB experiment
        load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/ODA_Kenya/infofile.RData")
        load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/ODA_Kenya/RESULTS_PROCESSED/ODA_Kenya_stock_flux.RData")
        load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/ODA_Kenya/RESULTS_PROCESSED/ODA_Kenya_parameter_maps.RData")
    }

} else if (model_choice == "C6") {


} else if (model_choice == "C7") {


} else if (model_choice == "C10") {


} else if (model_choice == "C11") {


} else if (model_choice == "G1") {


} else if (model_choice == "G2") {


} else if (model_choice == "G3") {


} else if (model_choice == "G4") {


} # model choice

# Store for later
original_end = as.numeric(PROJECT$end_year)

##
# Load or create information need for climate change meteorological drivers

if (PROJECT$model$timestep == "monthly") {
    print("...model will use calender monthly timestep ")
    all_years = as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)
    nos_days = nos_days_in_year(all_years[1])
    if (nos_days == 366) {timestep_days=c(31,29,31,30,31,30,31,31,30,31,30,31)} else {timestep_days=c(31,28,31,30,31,30,31,31,30,31,30,31)}
    for (y in seq(2, length(all_years))) {
         # calculate increment
         nos_days = nos_days_in_year(all_years[y])
         if (nos_days == 366) {
             timestep_days = append(timestep_days,c(31,29,31,30,31,30,31,31,30,31,30,31))
         } else {
             timestep_days = append(timestep_days,c(31,28,31,30,31,30,31,31,30,31,30,31))
         }
    } # loop through days
    # clean up
    rm(all_years)
} else if (PROJECT$model$timestep == "weekly") {
    print("...model will use weekly timestep ")
    all_years = as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)
    nos_days = nos_days_in_year(all_years[1])
    if (nos_days == 366) {timestep_days=c(rep(7,times=51),9)} else {timestep_days=c(rep(7,times=51),8)}
    for (y in seq(2, length(all_years))) {
         # calculate increment
         nos_days=nos_days_in_year(all_years[y])
         if (nos_days == 366) {
             timestep_days = append(timestep_days,c(rep(7,times=51),9))
         } else {
             timestep_days = append(timestep_days,c(rep(7,times=51),8))
         }
    } # loop through days
    # clean up
    rm(all_years)
} else if (PROJECT$model$timestep == "daily") {
    print("...model will use daily timestep")
    timestep_days = 1
} else {
    timestep_days = 1
    print("WARNING: user has not specified a time step option ('daily' or 'monthly'), default assumption is daily")
} # time step diagnosis loop

# generate the lat / long grid again
output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
latitude = output$lat ; longitude = output$long
# remove the values we don't want
if (length(remove) > 0) {latitude = latitude[-PROJECT$waterpixels] ; longitude = longitude[-PROJECT$waterpixels]}
# create lat / long axes, assumes WGS-84 grid
latitude_nc = seq(PROJECT$latitude[1]+(PROJECT$resolution*0.5),PROJECT$latitude[2]-(PROJECT$resolution*0.5), length.out = PROJECT$lat_dim)
longitude_nc = seq(PROJECT$longitude[1]+(PROJECT$resolution*0.5),PROJECT$longitude[2]-(PROJECT$resolution*0.5), length.out = PROJECT$long_dim)

# How many times do we need to loop our meteorology to reach the target year?
# -1 adjusts for the current meteorological dataset
nos_loops = ceiling(((new_end_year-as.numeric(PROJECT$start_year)) / length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))) - 1)

##
# Which quantiles do we want to extract information for...
nos_quantiles = length(quantiles_wanted)

# Loop through climate scenarios
for (s in seq(1,length(climate_scenarios))) {
     # Loop through Earth System Models providing the scenarios
     for (m in seq(1, length(ESM))) {
          
          # Create the static output variables for geotiff
          # ESM C stock comparisons
          ESM_bio_gCm2_initial = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))          
          ESM_bio_gCm2_2020 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2030 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2040 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2050 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2060 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2070 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2080 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2090 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_bio_gCm2_2100 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_initial = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_2020 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))             
          ESM_dom_gCm2_2030 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))          
          ESM_dom_gCm2_2040 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_2050 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_2060 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_2070 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))          
          ESM_dom_gCm2_2080 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_2090 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ESM_dom_gCm2_2100 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # Ecosystem property correlations
          # Biomass 2100
          bio_lca_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_cue_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_gpp_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_fMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_rMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_wMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_lMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_sMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_fNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_rNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_wNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_bioSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          bio_domSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Dead Organic Matter 2100
          dom_lca_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_cue_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_gpp_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_fMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_rMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_wMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_lMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))          
          dom_sMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_fNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_rNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_wNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_bioSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dom_domSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Biomass change (2015-2100)
          dbio_lca_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_cue_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_gpp_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_fMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_rMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_wMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_lMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_sMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_fNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_rNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_wNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_bioSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          dbio_domSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Dead Organic Matter change (2015-2100)
          ddom_lca_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_cue_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_gpp_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_fMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_rMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_wMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_lMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_sMTT_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_fNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_rNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_wNPP_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_bioSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          ddom_domSSprox_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # Find all the anomalies files for the current scenario / ESM climate anomaly information
          list_of_files = list.files(paste(cmip6dir,ESM[m],"/",climate_scenarios[s],"/month",sep=""), full.names=TRUE, recursive = TRUE)
          if (new_end_year > 2049) { 
              list_of_files = list_of_files[grepl("201501",list_of_files) == TRUE | grepl("205001",list_of_files) == TRUE]
              expected_nos_files = 2
          } else {
              list_of_files = list_of_files[grepl("201501",list_of_files)==TRUE]
              expected_nos_files = 1
          } 
          # Filter to consider only files which have anomalies in the title
          # Check that files exist for each of the variables we want
          if (length(list_of_files[grepl("/pr/",list_of_files)]) == expected_nos_files & 
              length(list_of_files[grepl("/rsds/",list_of_files)]) == expected_nos_files & 
              length(list_of_files[grepl("/sfcWind/",list_of_files)]) == expected_nos_files & 
              length(list_of_files[grepl("/tasmax/",list_of_files)]) == expected_nos_files &
              length(list_of_files[grepl("/tasmin/",list_of_files)]) == expected_nos_files & 
              length(list_of_files[grepl("/tas/",list_of_files)]) == expected_nos_files & 
              length(list_of_files[grepl("/vpd/",list_of_files)]) == expected_nos_files &
              length(list_of_files[grepl("/co2mass/",list_of_files)]) == expected_nos_files &
              length(list_of_files[grepl("/cLeaf/",list_of_files)]) == expected_nos_files &
              length(list_of_files[grepl("/cRoot/",list_of_files)]) == expected_nos_files &
              length(list_of_files[grepl("/cStem/",list_of_files)]) == expected_nos_files &
              length(list_of_files[grepl("/cSoil/",list_of_files)]) == expected_nos_files ) {

              # Open the climate variables we want
              sub_list_of_files = list_of_files[grepl("201501",list_of_files) == TRUE]
              rainfall_anomaly_file = nc_open(sub_list_of_files[grepl("/pr/",sub_list_of_files)])
              swrad_anomaly_file    = nc_open(sub_list_of_files[grepl("/rsds/",sub_list_of_files)])
              wind_spd_anomaly_file = nc_open(sub_list_of_files[grepl("/sfcWind/",sub_list_of_files)])
              maxt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmax/",sub_list_of_files)]) 
              mint_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmin/",sub_list_of_files)]) 
              airt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tas/",sub_list_of_files)])    
              vpd_anomaly_file      = nc_open(sub_list_of_files[grepl("/vpd/",sub_list_of_files)])
              co2_anomaly_file      = nc_open(sub_list_of_files[grepl("/co2mass/",sub_list_of_files)])
              # Open the C stock comparisons we want 
              Cfol_anomaly_file      = nc_open(sub_list_of_files[grepl("/cLeaf/",sub_list_of_files)])
              Croot_anomaly_file     = nc_open(sub_list_of_files[grepl("/cRoot/",sub_list_of_files)])
              Cstem_anomaly_file     = nc_open(sub_list_of_files[grepl("/cStem/",sub_list_of_files)])
              Csoil_anomaly_file     = nc_open(sub_list_of_files[grepl("/cSoil/",sub_list_of_files)])
              # Now read in each variable, applying relevant unit conversions
              # e.g. rainfall_anomaly_file$var$pr$units
              rainfall_kgm2s = ncvar_get(rainfall_anomaly_file,"pr")
              swrad_MJm2day  = ncvar_get(swrad_anomaly_file,"rsds") * 86400 * 1e-6 # W/m2 -> MJ/m2/day
              wind_spd_ms    = ncvar_get(wind_spd_anomaly_file,"sfcWind") 
              maxt_C         = ncvar_get(maxt_anomaly_file,"tasmax") - 273.15 # K->C
              mint_C         = ncvar_get(mint_anomaly_file,"tasmin") - 273.15 # K->C
              airt_C         = ncvar_get(airt_anomaly_file,"tas") - 273.15
              vpd_Pa         = ncvar_get(vpd_anomaly_file,"vpd")
              co2_ppm        = ncvar_get(co2_anomaly_file,"co2mass")/1e+12/7.804816 # kg -> ppm 
              # Where:
              # 1 ppmv of CO2= 2.13 Gt of C
              # 44.01 CO2 mass
              # 1 ppm CO2 = 2.13/12.0107*44.01 Gt = 7.804816
              # convert to Gt
    
              # Read in the ESM stock files (kgC/m2) convertion to gC/m2 later
              ESM_bio_gCm2 = ncvar_get(Cfol_anomaly_file,"cLeaf") 
              tmp = ncvar_get(Croot_anomaly_file,"cRoot") ; ESM_bio_gCm2 = ESM_bio_gCm2 + tmp
              tmp = ncvar_get(Cstem_anomaly_file,"cStem") ; ESM_bio_gCm2 = ESM_bio_gCm2 + tmp ; rm(tmp)
              ESM_dom_gCm2 = ncvar_get(Csoil_anomaly_file,"cSoil") 
              # Convert to gC/m2
              ESM_bio_gCm2 = ESM_bio_gCm2 * 1e3
              ESM_dom_gCm2 = ESM_dom_gCm2 * 1e3

              # Read in lat / long and time information from one of the files - it should be the same in all cases
              days_since_1850 = ncvar_get(vpd_anomaly_file, "time")
              lat_esm = ncvar_get(vpd_anomaly_file, "lat") # degrees north (-89.5 / 89.5)
              long_esm = ncvar_get(vpd_anomaly_file, "lon") # degrees east (-179.5 / 179.5)
              # Come out 0-360, convert to correct
              long_esm[which(long_esm > 180)] = long_esm[which(long_esm > 180)]-360 

              # Tidy up currently open files
              nc_close(rainfall_anomaly_file) ; nc_close(swrad_anomaly_file) ; nc_close(wind_spd_anomaly_file)
              nc_close(maxt_anomaly_file) ; nc_close(mint_anomaly_file) ; nc_close(airt_anomaly_file)
              nc_close(vpd_anomaly_file) ; nc_close(co2_anomaly_file) ; nc_close(Cfol_anomaly_file)
              nc_close(Croot_anomaly_file) ; nc_close(Cstem_anomaly_file) ; nc_close(Csoil_anomaly_file)

              # If we need to second half of the time series then read them in here
              if (expected_nos_files == 2) {
                  # Open the climate variables we want
                  sub_list_of_files = list_of_files[grepl("20500",list_of_files)==TRUE]
                  rainfall_anomaly_file = nc_open(sub_list_of_files[grepl("/pr/",sub_list_of_files)])
                  swrad_anomaly_file    = nc_open(sub_list_of_files[grepl("/rsds/",sub_list_of_files)])
                  wind_spd_anomaly_file = nc_open(sub_list_of_files[grepl("/sfcWind/",sub_list_of_files)])
                  maxt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmax/",sub_list_of_files)]) 
                  mint_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmin/",sub_list_of_files)]) 
                  airt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tas/",sub_list_of_files)])    
                  vpd_anomaly_file      = nc_open(sub_list_of_files[grepl("/vpd/",sub_list_of_files)])
                  co2_anomaly_file      = nc_open(sub_list_of_files[grepl("/co2mass/",sub_list_of_files)])
                  # Open the C stock comparisons we want 
                  Cfol_anomaly_file     = nc_open(sub_list_of_files[grepl("/cLeaf/",sub_list_of_files)])
                  Croot_anomaly_file    = nc_open(sub_list_of_files[grepl("/cRoot/",sub_list_of_files)])
                  Cstem_anomaly_file    = nc_open(sub_list_of_files[grepl("/cStem/",sub_list_of_files)])
                  Csoil_anomaly_file    = nc_open(sub_list_of_files[grepl("/cSoil/",sub_list_of_files)])
                  # Now read in each variable, applying relevant unit conversions
                  # e.g. rainfall_anomaly_file$var$pr$units
                  tmp_rainfall_kgm2s = ncvar_get(rainfall_anomaly_file,"pr")
                  tmp_swrad_MJm2day  = ncvar_get(swrad_anomaly_file,"rsds") * 86400 * 1e-6 # W/m2 -> MJ/m2/day
                  tmp_wind_spd_ms    = ncvar_get(wind_spd_anomaly_file,"sfcWind")
                  tmp_maxt_C         = ncvar_get(maxt_anomaly_file,"tasmax") - 273.15 # K->C
                  tmp_mint_C         = ncvar_get(mint_anomaly_file,"tasmin") - 273.15 # K->C
                  tmp_airt_C         = ncvar_get(airt_anomaly_file,"tas") - 273.15
                  tmp_vpd_Pa         = ncvar_get(vpd_anomaly_file,"vpd")
                  co2_ppm            = append(co2_ppm,ncvar_get(co2_anomaly_file,"co2mass")/1e+12/7.804816) # kg -> ppm 
                  # Where:
                  # 1 ppmv of CO2= 2.13 Gt of C
                  # 44.01 CO2 mass
                  # 1 ppm CO2 = 2.13/12.0107*44.01 Gt = 7.804816
                  # convert to Gt

                  # Read in the ESM stock files (kgC/m2) convertion to gC/m2 later
                  tmp_ESM_bio_gCm2 = ncvar_get(Cfol_anomaly_file,"cLeaf") 
                  tmp = ncvar_get(Croot_anomaly_file,"cRoot") ; tmp_ESM_bio_gCm2 = tmp_ESM_bio_gCm2 + tmp
                  tmp = ncvar_get(Cstem_anomaly_file,"cStem") ; tmp_ESM_bio_gCm2 = tmp_ESM_bio_gCm2 + tmp ; rm(tmp)
                  tmp_ESM_dom_gCm2 = ncvar_get(Csoil_anomaly_file,"cSoil") 
                  # Convert to gC/m2
                  tmp_ESM_bio_gCm2 = tmp_ESM_bio_gCm2 * 1e3
                  tmp_ESM_dom_gCm2 = tmp_ESM_dom_gCm2 * 1e3

                  # Read in lat / long and time information from one of the files - it should be the same in all cases
                  days_since_1850 = append(days_since_1850,ncvar_get(vpd_anomaly_file, "time"))

                  # Combine both datasets together
                  dims = dim(rainfall_kgm2s) ; dims[3] = length(days_since_1850)
                  rainfall_kgm2s = array(c(rainfall_kgm2s,tmp_rainfall_kgm2s), dim=dims)
                  swrad_MJm2day  = array(c(swrad_MJm2day,tmp_swrad_MJm2day), dim=dims)
                  wind_spd_ms    = array(c(wind_spd_ms,tmp_wind_spd_ms), dim=dims)
                  maxt_C         = array(c(maxt_C,tmp_maxt_C), dim=dims)
                  mint_C         = array(c(mint_C,tmp_mint_C), dim=dims)
                  airt_C         = array(c(airt_C,tmp_airt_C), dim=dims)
                  vpd_Pa         = array(c(vpd_Pa,tmp_vpd_Pa), dim=dims)
                  # Now C stocks combined
                  ESM_bio_gCm2 = array(c(ESM_bio_gCm2,tmp_ESM_bio_gCm2), dim=dims)
                  ESM_dom_gCm2 = array(c(ESM_dom_gCm2,tmp_ESM_dom_gCm2), dim=dims)

                  # Tidy up the tmp_* variables
                  rm(tmp_rainfall_kgm2s,tmp_swrad_MJm2day,tmp_wind_spd_ms,tmp_maxt_C,tmp_mint_C,tmp_airt_C,tmp_vpd_Pa,
                     tmp_ESM_bio_gCm2,tmp_ESM_dom_gCm2)

                  # Tidy up currently open files
                  nc_close(rainfall_anomaly_file) ; nc_close(swrad_anomaly_file) ; nc_close(wind_spd_anomaly_file)
                  nc_close(maxt_anomaly_file) ; nc_close(mint_anomaly_file) ; nc_close(airt_anomaly_file)
                  nc_close(vpd_anomaly_file) ; nc_close(co2_anomaly_file)  ; nc_close(Cfol_anomaly_file)
                  nc_close(Croot_anomaly_file) ; nc_close(Cstem_anomaly_file) ; nc_close(Csoil_anomaly_file)
 
              } # open second file if needed

              # Open harvest scenario files
              harvest_file = nc_open(paste(luh2dir,"/",climate_scenarios[s],"/LUH2_biomass_removal_scenario.nc",sep=""))
              harvest = ncvar_get(harvest_file, "harvest_fraction")
              lat_luh2 = ncvar_get(harvest_file, "lat") ; long_luh2 = ncvar_get(harvest_file, "lon")
              years_since_2015 = ncvar_get(harvest_file, "time")

              # Tidy up files
              nc_close(harvest_file)

              # Estimate division needed to the annual fraction of biomass harvest to that needed to be applied for each timestep
              harvest = harvest / (365.25 / mean(timestep_days))
              # Determine the start point in the management scenario, i.e. the year after the observed period
              begin_anomaly = which(years_since_2015 == (max(0,(original_end-2015))+1))
              harvest = harvest[,,begin_anomaly:length(years_since_2015)]
              rm(begin_anomaly)

              # Input is monthly but we want to know this information for estimating the anomaly points
              months_per_year = 12 ; days_per_month = 30 ; days_per_year = 360 
      
              # What year of the climate scenario is the anomaly from?
              begin_anomaly = max(0,(original_end - 2015))
              begin_anomaly = (begin_anomaly * months_per_year) + 1
              # Restrict to the desired time period
              end_anomaly = dim(rainfall_kgm2s)[3]
              rainfall_kgm2s = rainfall_kgm2s[,,begin_anomaly:end_anomaly]
              swrad_MJm2day = swrad_MJm2day[,,begin_anomaly:end_anomaly]
              wind_spd_ms = wind_spd_ms[,,begin_anomaly:end_anomaly]
              maxt_C = maxt_C[,,begin_anomaly:end_anomaly]
              mint_C = mint_C[,,begin_anomaly:end_anomaly]
              airt_C = airt_C[,,begin_anomaly:end_anomaly]
              vpd_Pa = vpd_Pa[,,begin_anomaly:end_anomaly]
              co2_ppm = co2_ppm[begin_anomaly:end_anomaly]
              ESM_bio_gCm2 = ESM_bio_gCm2[,,begin_anomaly:end_anomaly]
              ESM_dom_gCm2 = ESM_dom_gCm2[,,begin_anomaly:end_anomaly]

              # Calculate the anomaly for the whole grid
              for (x in seq(1,dim(rainfall_kgm2s)[1])) {
                   for (y in seq(1, dim(rainfall_kgm2s)[2])) {
                        rainfall_kgm2s[x,y,] = rainfall_kgm2s[x,y,] - rainfall_kgm2s[x,y,1:months_per_year]
                        swrad_MJm2day[x,y,] = swrad_MJm2day[x,y,] - swrad_MJm2day[x,y,1:months_per_year]
                        wind_spd_ms[x,y,] = wind_spd_ms[x,y,] - wind_spd_ms[x,y,1:months_per_year]
                        maxt_C[x,y,] = maxt_C[x,y,] - maxt_C[x,y,1:months_per_year]
                        mint_C[x,y,] = mint_C[x,y,] - mint_C[x,y,1:months_per_year]
                        airt_C[x,y,] = airt_C[x,y,] - airt_C[x,y,1:months_per_year]
                        vpd_Pa[x,y,] = vpd_Pa[x,y,] - vpd_Pa[x,y,1:months_per_year]
                   } # loop y axis
              } # loop x axis
              # CO2 is a global average value only - i.e. no grid
              co2_ppm = co2_ppm - co2_ppm[1:months_per_year]
              # Now remove the zeros from the beginning
              rainfall_kgm2s = rainfall_kgm2s[,,(months_per_year+1):dim(rainfall_kgm2s)[3]]
              swrad_MJm2day = swrad_MJm2day[,,(months_per_year+1):dim(swrad_MJm2day)[3]]
              wind_spd_ms = wind_spd_ms[,,(months_per_year+1):dim(wind_spd_ms)[3]]
              maxt_C = maxt_C[,,(months_per_year+1):dim(maxt_C)[3]]
              mint_C = mint_C[,,(months_per_year+1):dim(mint_C)[3]]
              airt_C = airt_C[,,(months_per_year+1):dim(airt_C)[3]]
              vpd_Pa = vpd_Pa[,,(months_per_year+1):dim(vpd_Pa)[3]]
              co2_ppm = co2_ppm[(months_per_year+1):length(co2_ppm)]

              ## 
              # Begin looping through each site 
              create_nc_vars = TRUE
              for (n in seq(1, PROJECT$nosites)) {

                   # Report to the user
                   if (n%%100 == 0) {print(paste("Site = ",n," of ",PROJECT$nosites,sep=""))}
  
                   # Construct the input file name...
                   infile = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
                   # ...and if the file exists begin running it
                   if (file.exists(infile)) {
    
                       # Load local CARDAMOM output file
                       load(infile)
  
                       ##
                       # Create site specific meteorological drivers

                       if (exists("steps_per_year") == FALSE) {
                           steps_per_year = dim(drivers$met)[1] / length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
                       }
    
                       # row bind the entire met array to expand by adding the times series on each time
                       new_drivers = drivers
                       for (r in seq(1,nos_loops)) {
                            new_drivers$met = rbind(new_drivers$met,drivers$met)
                            new_drivers$obs = rbind(new_drivers$obs,drivers$obs)
                       }
              
                       # Just looping the existing drivers doesn't work. 
                       # Clear obs from the post climate change period
                       new_drivers$obs[(dim(drivers$met)[1]+1):dim(new_drivers$met)[1],] = -9999
                       # The simulation day need to be updated so continuously increase.
                       # We can achieve this by cumulative sum of the day of year variable
                       new_drivers$met[,1] = timestep_days ; new_drivers$met[,1] = cumsum(new_drivers$met[,1])
                       # Remove biomass removals (i.e. management) - but leave fire
                       new_drivers$met[(dim(drivers$met)[1]+1):dim(new_drivers$met)[1],8] = 0

                       # Find site within the harvest scenario analysis grid
                       output = closest2d(1,lat_luh2,long_luh2,latitude[n],longitude[n],3)
                       i1 = unlist(output)[1] ; j1 = unlist(output)[2]
                       # Now impose LUH2 biomass removal scenarios
                       start = (dim(drivers$met)[1]+1) ; finish = dim(new_drivers$met)[1]
                       new_drivers$met[start:finish,8] = rep(harvest[i1,j1,], each = steps_per_year)[1:length(c(start:finish))]
                       new_drivers$met[is.na(new_drivers$met[,8]) == TRUE,8] = 0

                       # Adjust to make input array match the maximum length of the anomaly
                       new_drivers$met = new_drivers$met[1:(dim(drivers$met)[1]+length(co2_ppm)),]
                       new_drivers$obs = new_drivers$obs[1:(dim(drivers$met)[1]+length(co2_ppm)),]
                       # Determine the number of years in the simulation in total
                       nos_years = dim(new_drivers$met)[1] / steps_per_year
                       # Determine the number of years in the anomaly
                       nos_years_anomaly = nos_years - length(c(as.numeric(PROJECT$start_year):as.numeric(original_end)))
                       # Number years needed for output from dalec.so, make sure that the original PROJECT time period is updated
                       analysis_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
                       if (analysis_years != nos_years) {
                           PROJECT$end_year = as.numeric(PROJECT$end_year) + (nos_years - analysis_years)
                       }

                       # Now we have updated the PROJECT end year we will create, if we have not already done so, 
                       # the output variables needed for the NetCDF file. We assume fresh variables for each new grid.
                       if (create_nc_vars) {
                           # Update nc variable creation flag
                           create_nc_vars = FALSE
                           # DRIVERS
                           AIRT_MIN = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           AIRT_MAX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           SWRAD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           CO2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           DOY = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           PRECIP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           FLOSS_FRAC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1])) 
                           BURNT_FRAC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           WINDSPD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           VPD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           # OBSERVATIONS
                           LAI_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           LAI_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           WOOD_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           WOOD_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           SOIL_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           SOIL_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(new_drivers$met)[1]))
                           # STATES (gC/m2)
                           LAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           TOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           BIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           DOM = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           LAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           FOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           ROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           WOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           LIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           SOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           #WLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           # FLUXES (gC/m2/day)
                           GPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           RAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           RHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           NPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           FIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           HARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           RECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           NEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           NBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           NBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           # BIOPHYSICAL
                           CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           # NPP allocation (labile, foliar, fine root, wood, gC/m2/day)
                           fNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           wNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           rNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           # NPP allocation (labile, foliar, fine root, wood, fraction)
                           fNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           wNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           rNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           # MRT (foliar, wood, fine root, litter, soil; years)
                           fMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           wMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           rMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           lMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                           sMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,dim(new_drivers$met)[1]))
                       } # create_nc_vars

                       # Find site within the climate scenario analysis grid
                       output = closest2d(1,lat_esm,long_esm,latitude[n],longitude[n],3)
                       i1 = unlist(output)[1] ; j1 = unlist(output)[2]

                       # Estimate the start point for applying the anomaly
                       start = dim(drivers$met)[1]+1 ; finish = dim(new_drivers$met)[1]
                       # Strictly speaking the anomaly is applied to the final year of the original time series...
                       s_final_yr = dim(drivers$met)[1]-(steps_per_year-1) ; f_final_yr = dim(drivers$met)[1]
                       new_drivers$met[start:finish,2]  = new_drivers$met[s_final_yr:f_final_yr,2]  + mint_C[i1,j1,]
                       new_drivers$met[start:finish,3]  = new_drivers$met[s_final_yr:f_final_yr,3]  + maxt_C[i1,j1,]
                       new_drivers$met[start:finish,4]  = new_drivers$met[s_final_yr:f_final_yr,4]  + swrad_MJm2day[i1,j1,]
                       new_drivers$met[start:finish,5]  = new_drivers$met[s_final_yr:f_final_yr,5]  + co2_ppm
                       new_drivers$met[start:finish,7]  = new_drivers$met[s_final_yr:f_final_yr,7]  + rainfall_kgm2s[i1,j1,]
                       new_drivers$met[start:finish,10] = new_drivers$met[s_final_yr:f_final_yr,10] + mint_C[i1,j1,]
                       new_drivers$met[start:finish,12] = new_drivers$met[s_final_yr:f_final_yr,12] + vpd_Pa[i1,j1,]
                       new_drivers$met[start:finish,14] = new_drivers$met[s_final_yr:f_final_yr,14] + airt_C[i1,j1,]
                       new_drivers$met[start:finish,15] = new_drivers$met[s_final_yr:f_final_yr,15] + wind_spd_ms[i1,j1,]
                       new_drivers$met[start:finish,16] = new_drivers$met[s_final_yr:f_final_yr,16] + vpd_Pa[i1,j1,]

                       # Ensure realism is obeyed
                       new_drivers$met[new_drivers$met[,4] < 0,4] = 0   # Can't have less than zero radiation
                       new_drivers$met[new_drivers$met[,7] < 0,7] = 0   # Can't have less than zero rainfall
                       new_drivers$met[new_drivers$met[,12] < 0,12] = 0 # Can't have less than zero VPD
                       new_drivers$met[new_drivers$met[,15] < 0,15] = 0 # Can't have less than zero wind speed
                       new_drivers$met[new_drivers$met[,16] < 0,16] = 0 # Can't have less than zero VPD

                       ##
                       # Now run the actual model

                       # run subsample of parameters for full results / propogation
                       soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
                       states_all = simulate_all(n,PROJECT,PROJECT$model$name,new_drivers$met,parameters[1:PROJECT$model$nopars[n],,],
                                                 drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                                                 PROJECT$exepath,soil_info)

                       # Very near zero values can be returned as NaN, we want to set those as zero here
                       states_all$gpp_gCm2day[is.na(states_all$gpp_gCm2day) == TRUE] = 0
                       states_all$rauto_gCm2day[is.na(states_all$rauto_gCm2day) == TRUE] = 0
                       states_all$rhet_gCm2day[is.na(states_all$rhet_gCm2day) == TRUE] = 0
                       states_all$reco_gCm2day[is.na(states_all$reco_gCm2day) == TRUE] = 0
                       states_all$nee_gCm2day[is.na(states_all$nee_gCm2day) == TRUE] = 0
                       states_all$fire_gCm2day[is.na(states_all$fire_gCm2day) == TRUE] = 0
                       states_all$harvest_C_gCm2day[is.na(states_all$harvest_C_gCm2day) == TRUE] = 0
                       states_all$CiCa[is.na(states_all$CiCa) == TRUE] = 0
                       states_all$lai_m2m2[is.na(states_all$lai_m2m2) == TRUE] = 0
                       states_all$lab_gCm2[is.na(states_all$lab_gCm2) == TRUE] = 0
                       states_all$fol_gCm2[is.na(states_all$fol_gCm2) == TRUE] = 0
                       states_all$root_gCm2[is.na(states_all$root_gCm2) == TRUE] = 0
                       states_all$wood_gCm2[is.na(states_all$wood_gCm2) == TRUE] = 0
                       states_all$lit_gCm2[is.na(states_all$lit_gCm2) == TRUE] = 0
                       if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
                           states_all$litwood_gCm2[is.na(states_all$litwood_gCm2) == TRUE] = 0
                       }
                       states_all$som_gCm2[is.na(states_all$som_gCm2) == TRUE] = 0
                       states_all$bio_gCm2[is.na(states_all$bio_gCm2) == TRUE] = 0
                       states_all$MTT[is.na(states_all$MTT) == TRUE] = 0
                       states_all$aMTT[is.na(states_all$aMTT) == TRUE] = 0
                       states_all$aNPP[is.na(states_all$aNPP) == TRUE] = 0
                       states_all$SS[is.na(states_all$SS) == TRUE] = 0
                       # Create some derived values
                       states_all$npp_gCm2day = states_all$gpp_gCm2day - states_all$rauto_gCm2day
                       states_all$fnpp_gCm2day = states_all$npp_gCm2day * states_all$aNPP[,1]
                       states_all$rnpp_gCm2day = states_all$npp_gCm2day * states_all$aNPP[,2]
                       states_all$wnpp_gCm2day = states_all$npp_gCm2day * states_all$aNPP[,3]
                       states_all$nbe_gCm2day = states_all$nee_gCm2day + states_all$fire_gCm2day
                       states_all$nbp_gCm2day = -states_all$nbe_gCm2day - states_all$harvest_C_gCm2day
                       states_all$totalC_gCm2 = states_all$bio_gCm2 + states_all$som_gCm2 + states_all$lit_gCm2
                       states_all$dom_gCm2 = states_all$som_gCm2 + states_all$lit_gCm2
                       if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
                           states_all$totalC_gCm2 = states_all$totalC_gCm2 + states_all$litwood_gCm2
                           states_all$dom_gCm2 = states_all$dom_gCm2 + states_all$litwood_gCm2
                       }

                       ###
                       # Extract information for NetCDF file

                       # DRIVERS
                       AIRT_MIN[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,2] # mint C
                       AIRT_MAX[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,3] # maxt C
                       SWRAD[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,4] # SWRAD MJ/m2/day
                       CO2[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,5] # CO2 ppm
                       DOY[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,6] # Julian day of year
                       PRECIP[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,7] # precipitation kgH2O/m2/s
                       FLOSS_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,8] # Forest loss fraction
                       BURNT_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,9] # Burned fraction
                       WINDSPD[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,15] # Wind speed m/s
                       VPD[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$met[,16] # Vapour pressure deficit Pa
                       # OBSERVATIONS
                       LAI_OBS[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$obs[,3] # LAI m2/m2
                       LAI_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = new_drivers$obs[,4] # LAI UNC m2/m2
                       # ...first assign priors to the 1st time step to simplify our storage
                       WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = new_drivers$parpriors[21]
                       WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = new_drivers$parpriorunc[21]
                       SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = new_drivers$parpriors[23]
                       SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = new_drivers$parpriorunc[23]
                       # ...second assign time series information if any exists
                       WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],2:dim(new_drivers$met)[1]] = new_drivers$obs[2:dim(new_drivers$met)[1],13]
                       WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:dim(new_drivers$met)[1]] = new_drivers$obs[2:dim(new_drivers$met)[1],14]
                       SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],2:dim(new_drivers$met)[1]] = new_drivers$obs[2:dim(new_drivers$met)[1],19]
                       SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:dim(new_drivers$met)[1]] = new_drivers$obs[2:dim(new_drivers$met)[1],20]

                       # STATES
                       LAI[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$lai_m2m2,2,quantile, prob=quantiles_wanted)
                       TOT[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$totalC_gCm2,2,quantile, prob=quantiles_wanted)
                       BIO[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$bio_gCm2,2,quantile, prob=quantiles_wanted)
                       DOM[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$dom_gCm2,2,quantile, prob=quantiles_wanted)
                       LAB[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$lab_gCm2,2,quantile, prob=quantiles_wanted)
                       FOL[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$fol_gCm2,2,quantile, prob=quantiles_wanted)
                       ROOT[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$root_gCm2,2,quantile, prob=quantiles_wanted)
                       WOOD[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$wood_gCm2,2,quantile, prob=quantiles_wanted)
                       LIT[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$lit_gCm2,2,quantile, prob=quantiles_wanted)
                       SOIL[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$som_gCm2,2,quantile, prob=quantiles_wanted)
#                       WLIT[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$litwood_gCm2,2,quantile, prob=quantiles_wanted)
                       # FLUXES
                       GPP[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$gpp_gCm2day,2,quantile, prob=quantiles_wanted)
                       RAU[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$rauto_gCm2day,2,quantile, prob=quantiles_wanted)
                       RHE[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$rhet_gCm2day,2,quantile, prob=quantiles_wanted)
                       NPP[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$npp_gCm2day,2,quantile, prob=quantiles_wanted)
                       FIR[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$fire_gCm2day,2,quantile, prob=quantiles_wanted)
                       HARV[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$harvest_C_gCm2day,2,quantile, prob=quantiles_wanted)
                       RECO[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$reco_gCm2day,2,quantile, prob=quantiles_wanted)
                       NEE[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$nee_gCm2day,2,quantile, prob=quantiles_wanted)
                       NBE[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$nbe_gCm2day,2,quantile, prob=quantiles_wanted)
                       NBP[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$nbp_gCm2day,2,quantile, prob=quantiles_wanted)
                       # BIOPHYSICAL
                       CiCa[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$CiCa,2,quantile, prob=quantiles_wanted)
                       # NPP (foliar, root, wood; gC/m2/day)
                       fNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$fnpp_gCm2day,2,quantile, prob=quantiles_wanted)
                       rNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$rnpp_gCm2day,2,quantile, prob=quantiles_wanted)
                       wNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = apply(states_all$wnpp_gCm2day,2,quantile, prob=quantiles_wanted)
        
                       # NPP (fraction) and MRT years are requested to have same number of time steps as stocks and fluxes
                       # This is awkward as no easy way to repeat specific elements without loop for variables which have no meaningful value at sub-annual timescales
                       # (and are therefore calculated as annuals)
                       for (q in seq(1, nos_quantiles)) {
                            # MRT (years)
                            fMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(apply(states_all$aMTT[,1,],2,quantile,prob=quantiles_wanted[q]), each = steps_per_year)
                            rMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(apply(states_all$aMTT[,2,],2,quantile,prob=quantiles_wanted[q]), each = steps_per_year)
                            wMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(apply(states_all$aMTT[,3,],2,quantile,prob=quantiles_wanted[q]), each = steps_per_year)
                            lMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(apply(states_all$aMTT[,4,],2,quantile,prob=quantiles_wanted[q]), each = steps_per_year)
                            sMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(apply(states_all$aMTT[,5,],2,quantile,prob=quantiles_wanted[q]), each = steps_per_year)
                            # NPP (fraction)
                            fNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(quantile(states_all$aNPP[,1], prob=quantiles_wanted[q]), each = dim(new_drivers$met)[1])
                            rNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(quantile(states_all$aNPP[,2], prob=quantiles_wanted[q]), each = dim(new_drivers$met)[1])
                            wNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(quantile(states_all$aNPP[,3], prob=quantiles_wanted[q]), each = dim(new_drivers$met)[1])
                       } # loop quantiles

                       ### 
                       # Extract static information for use in geotiffs

                       # Post-hoc calculation of parameter correlations with key C-cycle variables              
                       end_time = dim(states_all$nee_gCm2day)[2]
                       dom = states_all$bio_gCm2[,end_time]+states_all$som_gCm2[,end_time]+states_all$lit_gCm2[,end_time]  
                       if (length(which(names(states_all) == "litwood_gCm2")) > 0) {dom = dom + states_all$litwood_gCm2[,end_time]}
                       ddom = dom-(states_all$bio_gCm2[,1]+states_all$som_gCm2[,1]+states_all$lit_gCm2[,1])
                       if (length(which(names(states_all) == "litwood_gCm2")) > 0) {ddom = ddom - states_all$litwood_gCm2[,1]}
                       cue = pmax(0,pmin(1,apply(1-(states_all$rauto_gCm2day / states_all$gpp_gCm2day),1,mean, na.rm=TRUE)))
                       npp = apply(states_all$gpp_gCm2day - states_all$rauto_gCm2day,1,mean, na.rm=TRUE)
                       bioSSprox = apply(states_all$SS[,1:3],1,sum, na.rm=TRUE) - (states_all$bio_gCm2[,end_time] - states_all$bio_gCm2[,1])
                       domSSprox = apply(states_all$SS[,4:5],1,sum, na.rm=TRUE) - ddom
                       # Biomass 2100
                       bio_lca_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(as.vector(parameters[17,,]),states_all$bio_gCm2[,end_time])
                       bio_cue_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(cue,states_all$bio_gCm2[,end_time])
                       bio_gpp_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(apply(states_all$gpp_gCm2,1,mean),states_all$bio_gCm2[,end_time])
                       bio_fMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,1],states_all$bio_gCm2[,end_time])
                       bio_rMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,2],states_all$bio_gCm2[,end_time])
                       bio_wMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,3],states_all$bio_gCm2[,end_time])
                       bio_lMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,4],states_all$bio_gCm2[,end_time])
                       bio_sMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,5],states_all$bio_gCm2[,end_time])
#                       bio_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,1],states_all$bio_gCm2[,end_time])
#                       bio_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,2],states_all$bio_gCm2[,end_time])
#                       bio_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,3],states_all$bio_gCm2[,end_time])
                       bio_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,1]*npp),states_all$bio_gCm2[,end_time])
                       bio_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,2]*npp),states_all$bio_gCm2[,end_time])
                       bio_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,3]*npp),states_all$bio_gCm2[,end_time])
                       bio_bioSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(bioSSprox,states_all$bio_gCm2[,end_time])
                       bio_domSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(domSSprox,states_all$bio_gCm2[,end_time])
                       # Dead Organic Matter 2100
                       dom_lca_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(as.vector(parameters[17,,]),dom)
                       dom_cue_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(cue,dom)
                       dom_gpp_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(apply(states_all$gpp_gCm2,1,mean),dom)
                       dom_fMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,1],dom)
                       dom_rMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,2],dom)
                       dom_wMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,3],dom)
                       dom_lMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,4],dom)
                       dom_sMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,5],dom)
#                       dom_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,1],dom)
#                       dom_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,2],dom)
#                       dom_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,3],dom)
                       dom_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,1]*npp),dom)
                       dom_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,2]*npp),dom)
                       dom_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,3]*npp),dom)
                       dom_bioSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(bioSSprox,dom)
                       dom_domSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(domSSprox,dom)
                       # Biomass change (2015-2100)
                       dbio_lca_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(as.vector(parameters[17,,]),(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_cue_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(cue,(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_gpp_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(apply(states_all$gpp_gCm2,1,mean),(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_fMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,1],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_rMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,2],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_wMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,3],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_lMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,4],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_sMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,5],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
#                       dbio_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,1],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
#                       dbio_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,2],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
#                       dbio_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,3],(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,1]*npp),(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,2]*npp),(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,3]*npp),(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_bioSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(bioSSprox,(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       dbio_domSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(domSSprox,(states_all$bio_gCm2[,end_time]-states_all$bio_gCm2[,1]))
                       # Dead Organic Matter change (2015-2100)
                       ddom_lca_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(as.vector(parameters[17,,]),ddom)
                       ddom_cue_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(cue,ddom)
                       ddom_gpp_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(apply(states_all$gpp_gCm2,1,mean),ddom)
                       ddom_fMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,1],ddom)
                       ddom_rMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,2],ddom)
                       ddom_wMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,3],ddom)
                       ddom_lMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,4],ddom)
                       ddom_sMTT_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$MTT[,5],ddom)
#                       ddom_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,1],ddom)
#                       ddom_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,2],ddom)
#                       ddom_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(states_all$aNPP[,3],ddom)
                       ddom_fNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,1]*npp),ddom)
                       ddom_rNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,2]*npp),ddom)
                       ddom_wNPP_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor((states_all$aNPP[,3]*npp),ddom)
                       ddom_bioSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(bioSSprox,ddom)
                       ddom_domSSprox_cor[grid_output$i_location[n],grid_output$j_location[n]] = cor(domSSprox,ddom)

                       ###
                       # Extract static information for ESM Climate & C stock for geotiffs

                       # Initial
                       ESM_bio_gCm2_initial[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,1]
                       ESM_dom_gCm2_initial[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,1]

                       # 2020
                       year_of_run = length(as.numeric(PROJECT$start_year):2020)
                       timestep_wanted = (steps_per_year * year_of_run) # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2020[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2020[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2030
                       year_of_run = length(as.numeric(PROJECT$start_year):2030)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2030[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2030[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2040
                       year_of_run = length(as.numeric(PROJECT$start_year):2040)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2040[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2040[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2050
                       year_of_run = length(as.numeric(PROJECT$start_year):2049)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2050[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2050[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2060
                       year_of_run = length(as.numeric(PROJECT$start_year):2059)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2060[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2060[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2070
                       year_of_run = length(as.numeric(PROJECT$start_year):2069)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2070[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2070[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2080
                       year_of_run = length(as.numeric(PROJECT$start_year):2079)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2080[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2080[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2090
                       year_of_run = length(as.numeric(PROJECT$start_year):2089)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2090[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2090[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
                       # 2100
                       year_of_run = length(as.numeric(PROJECT$start_year):2100)
                       timestep_wanted = (steps_per_year * year_of_run)-f_final_yr # offset for difference in time series of DALEC run vs anomaly period
                       conv_to_rate = year_of_run-(f_final_yr/steps_per_year)
                       ESM_bio_gCm2_2100[grid_output$i_location[n],grid_output$j_location[n]] = ESM_bio_gCm2[i1,j1,timestep_wanted]
                       ESM_dom_gCm2_2100[grid_output$i_location[n],grid_output$j_location[n]] = ESM_dom_gCm2[i1,j1,timestep_wanted]
    
                  } # does the current site file exist?
    
              } # loop through sites
   
              ###
              # Write to netcdf file
      
              ## define dimension
              lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude_nc )
              long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude_nc )
              time_dimen <- ncdim_def( "time", units="", 1:dim(new_drivers$met)[1])
              quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
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
              var11 = ncvar_def("cTotal",            unit="gC.m-2", longname = "Carbon live and dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

              # create the empty file
              output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,PROJECT$start_year,"_",PROJECT$end_year,".nc",sep="")
              new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                                 var1,var2,var3,var4,var5,var6,
                                                                #var7,
                                                                 var8,var9,var10,var11),
                                                                 force_v4 = TRUE)

              ## Load data into output variable

              ## TIMING
              ncvar_put(new_file, var0, new_drivers$met[,1])
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
              var11 = ncvar_def("CiCa_Ensemble",  unit="1", longname = "Internal:Ambiant CO2 ratio - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

              # create the empty file
              output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,PROJECT$start_year,"_",PROJECT$end_year,".nc",sep="")
              new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                                 var1,var2,var3,var4,var5,var6,var7,        
                                                                 var8,var9,var10,var11),
                                                                 force_v4 = TRUE)

              ## Load data into output variable

              ## TIMING
              ncvar_put(new_file, var0, new_drivers$met[,1])

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
              output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_MRT_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,PROJECT$start_year,"_",PROJECT$end_year,".nc",sep="")
              new_file=nc_create(filename=output_name, vars=list(var0,                                     
                                                                 var1,var2,var3,var4,var5,var6,var7,       
                                                                 var8,var9,var10,var11),
                                                                 force_v4 = TRUE)

              ## Load data into output variable
              
              ## TIMING
              ncvar_put(new_file, var0, new_drivers$met[,1])

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
              var1 = ncvar_def("tas_min", unit="C", longname = "Mean daily minimum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
              # Maximum air temperature
              var2 = ncvar_def("tas_max", unit="C", longname = "Mean daily maximum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
              # Short Wave radiation
              var3 = ncvar_def("rsds", unit="MJ.m-2.d-1", longname = "Mean downwelling short wave radiation", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
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
   
              # create the empty file
              output_name = paste(PROJECT$results_processedpath,output_prefix,"DRIVERS_OBS_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,PROJECT$start_year,"_",PROJECT$end_year,".nc",sep="")
              new_file=nc_create(filename=output_name, vars=list(var0,                                     
                                                                 var1,var2,var3,var4,var5,var6,var7,       
                                                                 var8,var9,var10,var11,var12,var13,var14,var15),
                                                                 force_v4 = TRUE)
              ## Load data into output variable

              ## TIMING
              ncvar_put(new_file, var0, new_drivers$met[,1])

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

              ###
              # Output to geoTIFFs
  
              # Create raster template with the correct spatial information
              # This function makes use of the WGS-84 coordinate system
              pt = raster(resolution = PROJECT$resolution, xmn = PROJECT$longitude[1], xmx = PROJECT$longitude[2], 
                                                           ymn = PROJECT$latitude[1],  ymx = PROJECT$latitude[2], crs = "+init=epsg:4326")

              # Convert arrays into raster format
              lat_dim = dim(ESM_bio_gCm2_initial)[2]
              # ESM C stocks
              ESM_bio_gCm2_initial = raster(t(ESM_bio_gCm2_initial[,lat_dim:1]))
              ESM_dom_gCm2_initial = raster(t(ESM_dom_gCm2_initial[,lat_dim:1]))
              ESM_bio_gCm2_2020 = raster(t(ESM_bio_gCm2_2020[,lat_dim:1]))
              ESM_dom_gCm2_2020 = raster(t(ESM_dom_gCm2_2020[,lat_dim:1]))
              ESM_bio_gCm2_2030 = raster(t(ESM_bio_gCm2_2030[,lat_dim:1]))
              ESM_dom_gCm2_2030 = raster(t(ESM_dom_gCm2_2030[,lat_dim:1]))
              ESM_bio_gCm2_2040 = raster(t(ESM_bio_gCm2_2040[,lat_dim:1]))
              ESM_dom_gCm2_2040 = raster(t(ESM_dom_gCm2_2040[,lat_dim:1]))
              ESM_bio_gCm2_2050 = raster(t(ESM_bio_gCm2_2050[,lat_dim:1]))
              ESM_dom_gCm2_2050 = raster(t(ESM_dom_gCm2_2050[,lat_dim:1]))
              ESM_bio_gCm2_2060 = raster(t(ESM_bio_gCm2_2060[,lat_dim:1]))
              ESM_dom_gCm2_2060 = raster(t(ESM_dom_gCm2_2060[,lat_dim:1]))
              ESM_bio_gCm2_2070 = raster(t(ESM_bio_gCm2_2070[,lat_dim:1]))
              ESM_dom_gCm2_2070 = raster(t(ESM_dom_gCm2_2070[,lat_dim:1]))
              ESM_bio_gCm2_2080 = raster(t(ESM_bio_gCm2_2080[,lat_dim:1]))
              ESM_dom_gCm2_2080 = raster(t(ESM_dom_gCm2_2080[,lat_dim:1]))
              ESM_bio_gCm2_2090 = raster(t(ESM_bio_gCm2_2090[,lat_dim:1]))
              ESM_dom_gCm2_2090 = raster(t(ESM_dom_gCm2_2090[,lat_dim:1]))
              ESM_bio_gCm2_2100 = raster(t(ESM_bio_gCm2_2100[,lat_dim:1]))
              ESM_dom_gCm2_2100 = raster(t(ESM_dom_gCm2_2100[,lat_dim:1]))

              # Ecosystem property correlations
              # Biomass 2100
              bio_lca_cor = raster(t(bio_lca_cor[,lat_dim:1]))
              bio_cue_cor = raster(t(bio_cue_cor[,lat_dim:1]))
              bio_gpp_cor = raster(t(bio_gpp_cor[,lat_dim:1]))
              bio_fMTT_cor = raster(t(bio_fMTT_cor[,lat_dim:1]))
              bio_rMTT_cor = raster(t(bio_rMTT_cor[,lat_dim:1]))
              bio_wMTT_cor = raster(t(bio_wMTT_cor[,lat_dim:1]))
              bio_lMTT_cor = raster(t(bio_lMTT_cor[,lat_dim:1]))
              bio_sMTT_cor = raster(t(bio_sMTT_cor[,lat_dim:1]))
              bio_fNPP_cor = raster(t(bio_fNPP_cor[,lat_dim:1]))
              bio_rNPP_cor = raster(t(bio_rNPP_cor[,lat_dim:1]))
              bio_wNPP_cor = raster(t(bio_wNPP_cor[,lat_dim:1]))
              bio_bioSSprox_cor = raster(t(bio_bioSSprox_cor[,lat_dim:1]))
              bio_domSSprox_cor = raster(t(bio_domSSprox_cor[,lat_dim:1]))
              # Dead Organic Matter 2100
              dom_lca_cor = raster(t(dom_lca_cor[,lat_dim:1]))
              dom_cue_cor = raster(t(dom_cue_cor[,lat_dim:1]))
              dom_gpp_cor = raster(t(dom_gpp_cor[,lat_dim:1]))
              dom_fMTT_cor = raster(t(dom_fMTT_cor[,lat_dim:1]))
              dom_rMTT_cor = raster(t(dom_rMTT_cor[,lat_dim:1]))
              dom_wMTT_cor = raster(t(dom_wMTT_cor[,lat_dim:1]))
              dom_lMTT_cor = raster(t(dom_lMTT_cor[,lat_dim:1]))
              dom_sMTT_cor = raster(t(dom_sMTT_cor[,lat_dim:1]))
              dom_fNPP_cor = raster(t(dom_fNPP_cor[,lat_dim:1]))
              dom_rNPP_cor = raster(t(dom_rNPP_cor[,lat_dim:1]))
              dom_wNPP_cor = raster(t(dom_wNPP_cor[,lat_dim:1]))
              dom_bioSSprox_cor = raster(t(dom_bioSSprox_cor[,lat_dim:1]))
              dom_domSSprox_cor = raster(t(dom_domSSprox_cor[,lat_dim:1]))
              # Biomass change (2015-2100)
              dbio_lca_cor = raster(t(dbio_lca_cor[,lat_dim:1]))
              dbio_cue_cor = raster(t(dbio_cue_cor[,lat_dim:1]))
              dbio_gpp_cor = raster(t(dbio_gpp_cor[,lat_dim:1]))
              dbio_fMTT_cor = raster(t(dbio_fMTT_cor[,lat_dim:1]))
              dbio_rMTT_cor = raster(t(dbio_rMTT_cor[,lat_dim:1]))
              dbio_wMTT_cor = raster(t(dbio_wMTT_cor[,lat_dim:1]))
              dbio_lMTT_cor = raster(t(dbio_lMTT_cor[,lat_dim:1]))              
              dbio_sMTT_cor = raster(t(dbio_sMTT_cor[,lat_dim:1]))
              dbio_fNPP_cor = raster(t(dbio_fNPP_cor[,lat_dim:1]))
              dbio_rNPP_cor = raster(t(dbio_rNPP_cor[,lat_dim:1]))
              dbio_wNPP_cor = raster(t(dbio_wNPP_cor[,lat_dim:1]))
              dbio_bioSSprox_cor = raster(t(dbio_bioSSprox_cor[,lat_dim:1]))
              dbio_domSSprox_cor = raster(t(dbio_domSSprox_cor[,lat_dim:1]))
              # Dead Organic Matter change (2015-2100)
              ddom_lca_cor = raster(t(ddom_lca_cor[,lat_dim:1]))
              ddom_cue_cor = raster(t(ddom_cue_cor[,lat_dim:1]))
              ddom_gpp_cor = raster(t(ddom_gpp_cor[,lat_dim:1]))
              ddom_fMTT_cor = raster(t(ddom_fMTT_cor[,lat_dim:1]))
              ddom_rMTT_cor = raster(t(ddom_rMTT_cor[,lat_dim:1]))
              ddom_wMTT_cor = raster(t(ddom_wMTT_cor[,lat_dim:1]))
              ddom_lMTT_cor = raster(t(ddom_lMTT_cor[,lat_dim:1]))              
              ddom_sMTT_cor = raster(t(ddom_sMTT_cor[,lat_dim:1]))
              ddom_fNPP_cor = raster(t(ddom_fNPP_cor[,lat_dim:1]))
              ddom_rNPP_cor = raster(t(ddom_rNPP_cor[,lat_dim:1]))
              ddom_wNPP_cor = raster(t(ddom_wNPP_cor[,lat_dim:1]))
              ddom_bioSSprox_cor = raster(t(ddom_bioSSprox_cor[,lat_dim:1]))
              ddom_domSSprox_cor = raster(t(ddom_domSSprox_cor[,lat_dim:1]))

              # Impose correct lat / long
              # ESM C stocks
              extent(ESM_bio_gCm2_initial) <- extent(pt)
              extent(ESM_dom_gCm2_initial) <- extent(pt)
              extent(ESM_bio_gCm2_2020) <- extent(pt)
              extent(ESM_dom_gCm2_2020) <- extent(pt)
              extent(ESM_bio_gCm2_2030) <- extent(pt)
              extent(ESM_dom_gCm2_2030) <- extent(pt)
              extent(ESM_bio_gCm2_2040) <- extent(pt)
              extent(ESM_dom_gCm2_2040) <- extent(pt)
              extent(ESM_bio_gCm2_2050) <- extent(pt)
              extent(ESM_dom_gCm2_2050) <- extent(pt)
              extent(ESM_bio_gCm2_2060) <- extent(pt)
              extent(ESM_dom_gCm2_2060) <- extent(pt)
              extent(ESM_bio_gCm2_2070) <- extent(pt)
              extent(ESM_dom_gCm2_2070) <- extent(pt)
              extent(ESM_bio_gCm2_2080) <- extent(pt)
              extent(ESM_dom_gCm2_2080) <- extent(pt)
              extent(ESM_bio_gCm2_2090) <- extent(pt)
              extent(ESM_dom_gCm2_2090) <- extent(pt)
              extent(ESM_bio_gCm2_2100) <- extent(pt)
              extent(ESM_dom_gCm2_2100) <- extent(pt)

              # Ecosystem property correlations
              # Biomass 2100
              extent(bio_lca_cor) = extent(pt)
              extent(bio_cue_cor) = extent(pt)
              extent(bio_gpp_cor) = extent(pt)
              extent(bio_fMTT_cor) = extent(pt)
              extent(bio_rMTT_cor) = extent(pt)
              extent(bio_wMTT_cor) = extent(pt)
              extent(bio_lMTT_cor) = extent(pt)
              extent(bio_sMTT_cor) = extent(pt)
              extent(bio_fNPP_cor) = extent(pt)
              extent(bio_rNPP_cor) = extent(pt)
              extent(bio_wNPP_cor) = extent(pt)
              extent(bio_bioSSprox_cor) = extent(pt)
              extent(bio_domSSprox_cor) = extent(pt)
              # Dead Organic Matter 2100
              extent(dom_lca_cor) = extent(pt)
              extent(dom_cue_cor) = extent(pt)
              extent(dom_gpp_cor) = extent(pt)
              extent(dom_fMTT_cor) = extent(pt)
              extent(dom_rMTT_cor) = extent(pt)
              extent(dom_wMTT_cor) = extent(pt)
              extent(dom_lMTT_cor) = extent(pt)
              extent(dom_sMTT_cor) = extent(pt)
              extent(dom_fNPP_cor) = extent(pt)
              extent(dom_rNPP_cor) = extent(pt)
              extent(dom_wNPP_cor) = extent(pt)
              extent(dom_bioSSprox_cor) = extent(pt)
              extent(dom_domSSprox_cor) = extent(pt)
              # Biomass change (2015-2100)
              extent(dbio_lca_cor) = extent(pt)
              extent(dbio_cue_cor) = extent(pt)
              extent(dbio_gpp_cor) = extent(pt)
              extent(dbio_fMTT_cor) = extent(pt)
              extent(dbio_rMTT_cor) = extent(pt)
              extent(dbio_wMTT_cor) = extent(pt)
              extent(dbio_lMTT_cor) = extent(pt)
              extent(dbio_sMTT_cor) = extent(pt)
              extent(dbio_fNPP_cor) = extent(pt)
              extent(dbio_rNPP_cor) = extent(pt)
              extent(dbio_wNPP_cor) = extent(pt)
              extent(dbio_bioSSprox_cor) = extent(pt)
              extent(dbio_domSSprox_cor) = extent(pt)
              # Dead Organic Matter change (2015-2100)
              extent(ddom_lca_cor) = extent(pt)
              extent(ddom_cue_cor) = extent(pt)
              extent(ddom_gpp_cor) = extent(pt)
              extent(ddom_fMTT_cor) = extent(pt)
              extent(ddom_rMTT_cor) = extent(pt)
              extent(ddom_wMTT_cor) = extent(pt)
              extent(ddom_lMTT_cor) = extent(pt)
              extent(ddom_sMTT_cor) = extent(pt)
              extent(ddom_fNPP_cor) = extent(pt)
              extent(ddom_rNPP_cor) = extent(pt)
              extent(ddom_wNPP_cor) = extent(pt)
              extent(ddom_bioSSprox_cor) = extent(pt)
              extent(ddom_domSSprox_cor) = extent(pt)

              # ESM Biomass Stock
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_initial_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_initial,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2020_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2020,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2030_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2030,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2040_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2040,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2050_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2050,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2060_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2060,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2070_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2070,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2080_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2080,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2090_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2090,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_bio_gCm2_2100_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_bio_gCm2_2100,file = outfile, format = "GTiff",overwrite=TRUE)
              # ESM Dead Organic Matter
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_initial_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_initial,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2020_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2020,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2030_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2030,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2040_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2040,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2050_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2050,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2060_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2060,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2070_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2070,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2080_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2080,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2090_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2090,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"ESM_dom_gCm2_2100_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ESM_dom_gCm2_2100,file = outfile, format = "GTiff",overwrite=TRUE)

              # Ecosystem property correlations
              # Biomass 2100
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_lca_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_lca_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_cue_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_cue_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_gpp_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_gpp_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_fMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_fMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_rMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_rMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_wMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_wMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_lMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_lMTT_cor,file = outfile, format = "GTiff",ovewrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_sMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_sMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_fNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_fNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_rNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_rNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_wNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_wNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_bioSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_bioSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_domSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(bio_domSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              # Dead Organic Matter 2100
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_lca_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_lca_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_cue_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_cue_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_gpp_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_gpp_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_fMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_fMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_rMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_rMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_wMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_wMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_lMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_lMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_sMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_sMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_fNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_fNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_rNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_rNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_wNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_wNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_bioSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_bioSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_domSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dom_domSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              # Biomass change (2015-2100)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_lca_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_lca_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_cue_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_cue_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_gpp_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_gpp_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_fMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_fMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_rMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_rMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_wMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_wMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_lMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_lMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_sMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_sMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_fNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_fNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_rNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_rNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_wNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_wNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_bioSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_bioSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"bio_change_domSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(dbio_domSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)

              # Dead Organic Matter change (2015-2100)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_lca_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_lca_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_cue_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_cue_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_gpp_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_gpp_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_fMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_fMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_rMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_rMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_wMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_wMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_lMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_lMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_sMTT_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_sMTT_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_fNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_fNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_rNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_rNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_wNPP_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_wNPP_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_bioSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_bioSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)
              outfile = paste(PROJECT$results_processedpath,output_prefix,"dom_change_domSSprox_cor_",ESM[m],"_",climate_scenarios[s],"_",output_suffix,".tif",sep="")
              writeRaster(ddom_domSSprox_cor,file = outfile, format = "GTiff",overwrite=TRUE)

          } # files exist for this scenario
     } # ESM
} # climate scenarios

