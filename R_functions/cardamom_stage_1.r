#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Function to carry out stage 1 processes, i.e. creation of CARDAMOM input files
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

# Create function needed to process the site specific creation
write_bin_files<-function(n,PROJECT,latlon,timestep_days,met_all
                         ,lai_all,Csom_all,forest_all
                         ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                         ,sand_clay_all,crop_man_all,burnt_all,soilwater_all
                         ,nbe_all, lca_all,gpp_all,Cwood_inc_all,Cwood_mortality_all,fire_all
                         ,fapar_all) {

   # create the file name for the met/obs binary
   filename = paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")

   # Determine local latitude value, ensure it is in wgs-84 -90/90 regardless of grid projection
   if (PROJECT$grid_type != "epsg:4326") {
       # The required grid for calculations in this function does not match, 
       # do the required conversions
       lat_degrees = vect(cbind(latlon[n,2], latlon[n,1]), crs=PROJECT$grid_type) 
       lat_degrees = project(lat_degrees, "epsg:4326")
       lat_degrees = crds(lat_degrees,df=TRUE) # extract latitude, i.e. y-dimension only
       lat_degrees = as.vector(lat_degrees$y)
   } else {

       # The required grid is a match for that provided here, assign to local variable and move on
       lat_degrees = latlon[n,1]
   } # lat in degrees or not?

   # All CARDAMOM read gridded datasets now map onto the same projection, extent and resolution.
   # This means that we can extract the location of the current site within any grid just the
   # once and pass around the solution to all the extraction functions.
   # NOTE: that met_all is an exception to this as the for the analysis to work we must always have
   # meteorology driving the model, so we assume a nearest neighbour approach of valid
   # locations, rather than accepting data gaps as done in observations / disturbance drivers.
   output = closest2d_tif(cardamom_ext,latlon[n,1],latlon[n,2]) 
   grid_long_loc = output$i_loc ; grid_lat_loc = dim(cardamom_ext) [1] - output$j_loc + 1
   grid_n = output$n_loc
   rm(output)

   # For selecting the right meteorological information
   if (PROJECT$spatial_type == "grid") {
       # Determine whether we have a valid meteorology variable, 
       # and the correct wheat_from_chaff number for the location.
       wheat_n = which(met_all$wheat == grid_n)
   } else {     
       # Assume site analysis, therefore wheat_n == n
       wheat_n = n
   } # grid or site?

   if (length(wheat_n) == 1) {

       # Assuming we have not already created the file or we wish to force recreation
       if (file.exists(filename) == FALSE | repair == 1){
            # Extract meteorology
            met = extract_met_drivers(wheat_n,timestep_days,PROJECT$start_year,PROJECT$end_year,
                                      lat_degrees,met_all,met_source,PROJECT$sites[n],PROJECT$grid_type)
#            # Load met drivers for ACM or other models
#            if (PROJECT$model$name != "ACM") {
#                met = extract_met_drivers(n,timestep_days,PROJECT$start_year,PROJECT$end_year,latlon[n,],met_all,met_source,PROJECT$sites[n])
#            } else { # if (PROJECT$model$name != "ACM")
#                # assume ACM special case
#                met = extract_acm_met_drivers(PROJECT,latlon[n,],PROJECT$sites[n])
#            }
            # Load observations
            obs = extract_obs(grid_long_loc,grid_lat_loc,latlon[n,],lai_all,Csom_all,forest_all
                             ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                             ,sand_clay_all,crop_man_all,burnt_all,soilwater_all
                             ,nbe_all, lca_all, gpp_all,Cwood_inc_all,Cwood_mortality_all, fire_all
                             ,fapar_all
                             ,PROJECT$ctessel_pft[n],PROJECT$sites[n],PROJECT$start_year,PROJECT$end_year
                             ,timestep_days,PROJECT$spatial_type,PROJECT$resolution,PROJECT$grid_type,PROJECT$model$name)

            # update ctessel pft in the project and potentially the model information
            PROJECT$ctessel_pft[n] = obs$ctessel_pft
            # Load additional model information
            PROJECT$model = cardamom_model_details(PROJECT$model$name,pft_specific_parameters,PROJECT$ctessel_pft)
            # write out the relevant binary files
            binary_data(met,obs,filename,PROJECT$edc,lat_degrees,PROJECT$ctessel_pft[n],
                        PROJECT$model$name,PROJECT$parameter_type,PROJECT$model$nopars[n],noyears)
       } # if (file.exists(filename) == FALSE | repair == 1)
    
   } # do we have a valid meteorological location

} # end function

cardamom_stage_1<-function(PROJECT) {

   # Check for existing binary input files
   print("Check for existance of binary input files")

   # Check whether there are any files which still need creating
   check = FALSE
   if (repair == 0) {
       n = 0 ; check = TRUE ; missing = 0
       while (n < PROJECT$nosites) {
            # Increment
            n = n + 1
            # Check whether the current file exists
            if (file.exists(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")) == FALSE) {
                # Escape loop now we know some sites need making
                check = FALSE ; missing = append(missing,n)
            }
       }
       # Say how many are missing
       print(paste("There are ",length(missing)-1," missing files of ",PROJECT$nosites," sites",sep=""))
   } # repair == 0

   if (check == FALSE) {
       # Update the user
       print("Beginning creation of binary input files")

       # flag for met drivers load
       met_all = 0 ; lai_all = 0 ; Csom_all = 0 ; forest_all = 0 ; Cwood_all = 0
       # load from PROJECT time step information
       timestep_days = PROJECT$model$timestep_days
       noyears = length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
       # Determine location information
       if (cardamom_type == "grid") {
           print("Determining number / locations of grid points for this run ...")
           output = determine_lat_long_needed(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution,PROJECT$grid_type,PROJECT$waterpixels)
           print("Have now determined grid point locations")
           # Bind together the latitude / longitudes for the grid, extract grid information needed to aid further processing (cardamom_ext)
           latlon = cbind(output$lat,output$long) ; cardamom_ext = output$cardamom_ext
           # Extract the lat / long grids designed for determining the extraction locations for the gridded observation datasets
           obs_long_grid = output$obs_long_grid ; obs_lat_grid = output$obs_lat_grid
           # Tidy up
           rm(output) ; gc(reset=TRUE,verbose=FALSE)
       } else if (cardamom_type != "grid") {
           print("Determining number / locations of grid points for this run ...")
           # Combine the latitude / longitude from the site list
           latlon = cbind(PROJECT$latitude,PROJECT$longitude)
           # However we still need a reduced area domain for extracting the site level analyses. +c(-0.5,0.5) allows buffer
           # Assume that all single sites are based on WGS-84 grid, EPSG:4326
           output = determine_lat_long_needed(lat = range(PROJECT$latitude)+c(-0.5,0.5), long = range(PROJECT$longitude)+c(-0.5,0.5)
                                             ,resolution = 0.125*0.5, grid_type = "epsg:4326", remove = NULL)
           # Extract grid information needed to aid further processing (cardamom_ext)
           cardamom_ext = output$cardamom_ext
           # Extract the lat / long grids designed for determining the extraction locations for the gridded observation datasets
           obs_long_grid = output$obs_long_grid ; obs_lat_grid = output$obs_lat_grid
           # Tidy up
           rm(output) ; gc(reset=TRUE,verbose=FALSE)
           print("Have now determined grid point locations")
       }

       # Load available data from gridded datasets?
       if (PROJECT$model$name != "ACM") {
           ## Load all time varying spatial forcings
           # Meteorological forcings
           met_all = load_met_fields_for_extraction(latlon,met_source,PROJECT$model$name,PROJECT$start_year,PROJECT$end_year,PROJECT$spatial_type,cardamom_ext)
           # Mechanical disturbance (e.g., deforestation, 0-1)
           forestry_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             deforestation_source,path_to_forestry,prefix = "forest_loss_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "forest_loss",
                                                             unc_var_name_in = "",
                                                             lag_var_name_in = "", 
                                                             est_var_name_out = "loss_fraction",
                                                             unc_var_name_out = "",
                                                             lag_var_name_out = "")
           # Burned area (0-1)
           burnt_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             burnt_area_source,path_to_burnt_area,prefix = "BurnedFraction_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "BurnedFraction",
                                                             unc_var_name_in = "",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "burnt_area",
                                                             unc_var_name_out = "",
                                                             lag_var_name_out = "")

           ## Load all time varying spatial observations
           # Leaf area index (m2/m2)
           lai_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             lai_source,path_to_lai,prefix = "leaf_area_index_m2m2_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "LAI",
                                                             unc_var_name_in = "LAI_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "lai_m2m2",
                                                             unc_var_name_out = "lai_unc_m2m2",
                                                             lag_var_name_out = "")      
           # fraction of Absorbed Photosynthetically Active Radation (0-1)
           fapar_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             fapar_source,path_to_fapar,prefix = "fraction_absorbed_par_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "fAPAR",
                                                             unc_var_name_in = "fAPAR_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "fapar",
                                                             unc_var_name_out = "fapar_unc",
                                                             lag_var_name_out = "")                   
           # Net Biome Exchange (gC/m2/day)
           nbe_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             nbe_source,path_to_nbe,prefix = "net_biome_exchange_gCm2day_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "NBE",
                                                             unc_var_name_in = "NBE_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "nbe_gCm2day",
                                                             unc_var_name_out = "nbe_unc_gCm2day",
                                                             lag_var_name_out = "")
           # Gross Primary Production (gC/m2/day)
           gpp_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             gpp_source,path_to_gpp,prefix = "gross_primary_production_gCm2day_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "GPP",
                                                             unc_var_name_in = "GPP_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "gpp_gCm2day",
                                                             unc_var_name_out = "gpp_unc_gCm2day",
                                                             lag_var_name_out = "")           
           # Fire carbon emissions (gC/m2/day)
           fire_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             fire_source,path_to_fire,prefix = "fire_carbon_emissions_gCm2day_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "Fire",
                                                             unc_var_name_in = "Fire_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "fire_gCm2day",
                                                             unc_var_name_out = "fire_unc_gCm2day",
                                                             lag_var_name_out = "")           
           # Wood stock (gC/m2)
           Cwood_stock_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             Cwood_stock_source,path_to_Cwood,prefix = "wood_stock_gCm2_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "wood_stock",
                                                             unc_var_name_in = "wood_stock_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "biomass_gCm2",
                                                             unc_var_name_out = "biomass_uncertainty_gCm2",
                                                             lag_var_name_out = "")           
           # Wood stock production (gC/m2/day)
           Cwood_inc_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             Cwood_inc_source,path_to_Cwood_inc,prefix = "wood_stock_production_gCm2day_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "wood_production",
                                                             unc_var_name_in = "wood_production_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "Cwood_increment_gCm2day",
                                                             unc_var_name_out = "Cwood_increment_uncertainty_gCm2day",
                                                             lag_var_name_out = "Cwood_increment_lag")                      
           # Wood stock mortality (gC/m2/day)
           Cwood_mortality_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             Cwood_mortality_source,path_to_Cwood_mortality,prefix = "wood_stock_mortality_gCm2day_",
                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
                                                             est_var_name_in = "wood_mortality",
                                                             unc_var_name_in = "wood_mortality_SD",
                                                             lag_var_name_in = "",
                                                             est_var_name_out = "Cwood_mortality_gCm2day",
                                                             unc_var_name_out = "Cwood_mortality_uncertainty_gCm2day",
                                                             lag_var_name_out = "Cwood_mortality_lag")                      
           # Surface soil water content (m3/m3)
#           soilwater_all = load_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
#                                                             soilwater_initial_source,path_to_soil_water,prefix = "soil_water_m3m3_",
#                                                             as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),
#                                                             est_var_name_in = "soil_moisture",
#                                                             unc_var_name_in = "soil_moisture_SD",
#                                                             lag_var_name_in = "",
#                                                             est_var_name_out = "soil_moisture_m3m3",
#                                                             unc_var_name_out = "soil_moisture_unc_m3m3",
#                                                             lag_var_name_out = "")       

           ## Load all static spatial forcings
           # Sand / Clay (%)
           sand_clay_all = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             sand_clay_source,path_to_sand_clay,prefix = "sand_percent_mean_0to30cm",
                                                             est_var_name_in = "sand_content",
                                                             unc_var_name_in = "sand_content_unc",
                                                             est_var_name_out = "top_sand",
                                                             unc_var_name_out = "")
           tmp = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             sand_clay_source,path_to_sand_clay,prefix = "sand_percent_mean_30to100cm",
                                                             est_var_name_in = "sand_content",
                                                             unc_var_name_in = "sand_content_unc",
                                                             est_var_name_out = "bot_sand",
                                                             unc_var_name_out = "")
           sand_clay_all$bot_sand = tmp$bot_sand ; rm(tmp) # update the list object with the next variable
           tmp = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             sand_clay_source,path_to_sand_clay,prefix = "clay_percent_mean_0to30cm",
                                                             est_var_name_in = "clay_content",
                                                             unc_var_name_in = "clay_content_unc",
                                                             est_var_name_out = "top_clay",
                                                             unc_var_name_out = "")
           sand_clay_all$top_clay = tmp$top_clay ; rm(tmp) # update the list object with the next variable
           tmp = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             sand_clay_source,path_to_sand_clay,prefix = "clay_percent_mean_30to100cm",
                                                             est_var_name_in = "clay_content",
                                                             unc_var_name_in = "clay_content_unc",
                                                             est_var_name_out = "bot_clay",
                                                             unc_var_name_out = "")
           sand_clay_all$bot_clay = tmp$bot_clay ; rm(tmp) # update the list object with the next variable

           ## Load all static spatial observations

           # Soil C stocks (gC/m2)
           Csom_all = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             Csom_source,path_to_Csom,prefix = "soil_stock_gCm2_",
                                                             est_var_name_in = "soil_stock",
                                                             unc_var_name_in = "soil_stock_SD",
                                                             est_var_name_out = "Csom",
                                                             unc_var_name_out = "Csom_unc")
           # Initial wood stocks (gC/m2)
           Cwood_initial_all = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             Cwood_initial_source,path_to_Cwood_initial,prefix = "wood_stock_gCm2_",
                                                             est_var_name_in = "wood_stock",
                                                             unc_var_name_in = "wood_stock_SD",
                                                             est_var_name_out = "biomass_gCm2",
                                                             unc_var_name_out = "biomass_uncertainty_gCm2")           
           # Potential wood stocks (gC/m2)
           Cwood_potential_all = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             Cwood_potential_source,path_to_Cwood_potential,prefix = "wood_stock_potential_gCm2_",
                                                             est_var_name_in = "wood_stock_potential",
                                                             unc_var_name_in = "wood_stock_potential_SD",
                                                             est_var_name_out = "biomass_gCm2",
                                                             unc_var_name_out = "biomass_uncertainty_gCm2")                
           # Leaf Carbon per unit leaf Area (LCA, gC/m2)
           lca_all = load_static_observation_dataset_for_extraction(latlon,cardamom_ext,PROJECT$spatial_type,
                                                             lca_source,path_to_lca,prefix = "leaf_carbon_area_gCm2_",
                                                             est_var_name_in = "leaf_carbon_area",
                                                             unc_var_name_in = "leaf_carbon_area_SD",
                                                             est_var_name_out = "lca_gCm2",
                                                             unc_var_name_out = "lca_uncertainty_gCm2")               

           #lai_all = load_lai_fields_for_extraction(latlon,lai_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           #fapar_all = load_fapar_fields_for_extraction(latlon,fapar_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           #nbe_all = load_nbe_fields_for_extraction(latlon,nbe_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           #gpp_all = load_gpp_fields_for_extraction(latlon,GPP_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
           #fire_all = load_fire_emission_fields_for_extraction(latlon,fire_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
           #Csom_all = load_Csom_fields_for_extraction(latlon,Csom_source,cardamom_ext,PROJECT$spatial_type)
           #sand_clay_all = load_sand_clay_fields_for_extraction(latlon,sand_clay_source,cardamom_ext,PROJECT$spatial_type)
           #forest_all = load_forestry_fields_for_extraction(latlon,deforestation_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           #Cwood_initial_all = load_initial_biomass_maps_for_extraction(latlon,Cwood_initial_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
           #Cwood_stock_all = load_biomass_stocks_maps_for_extraction(latlon,Cwood_stock_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
           #Cwood_potential_all = load_potential_biomass_maps_for_extraction(latlon,Cwood_potential_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
           #burnt_all = load_burnt_area_fields_for_extraction(latlon,burnt_area_source,path_to_burnt_area,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
           #soilwater_all = load_soilwater_fields_for_extraction(latlon,soilwater_initial_source)
           #lca_all = load_lca_maps_for_extraction(latlon,lca_source,cardamom_ext,PROJECT$spatial_type)
           #Cwood_inc_all = load_wood_productivity_maps_for_extraction(Cwood_inc_source,cardamom_ext,PROJECT$spatial_type,latlon,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days)
           #Cwood_mortality_all = load_wood_mortality_maps_for_extraction(Cwood_mortality_source,cardamom_ext,PROJECT$spatial_type,latlon,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days)

       } # if (PROJECT$model$name != "ACM")

       # Update user
       print("Loading completed, beginning writing out file write out")

       if (use_parallel) {

          # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
          # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only
          cl <- min(PROJECT$nosites,numWorkers)
          dummy = mclapply(c(1:PROJECT$nosites), FUN = write_bin_files, mc.cores = cl,
                           PROJECT = PROJECT,latlon = latlon,timestep_days = timestep_days,
                           met_all = met_all, lai_all = lai_all, Csom_all = Csom_all,
                           forest_all = forest_all, Cwood_initial_all = Cwood_initial_all,
                           Cwood_stock_all = Cwood_stock_all, Cwood_potential_all = Cwood_potential_all,
                           sand_clay_all = sand_clay_all, crop_man_all = crop_man_all,
                           burnt_all = burnt_all, soilwater_all = soilwater_all, nbe_all = nbe_all, 
                           lca_all = lca_all, gpp_all = gpp_all, Cwood_inc_all = Cwood_inc_all,
                           Cwood_mortality_all = Cwood_mortality_all, fire_all = fire_all, 
                           fapar_all = fapar_all)

      } else { # use parallel

          # start looping through sites to create site specific files of obs and met
          for (n in seq(1, PROJECT$nosites)) {

               # Inform user
               print(paste("Site ",n," of ",PROJECT$nosites," ",Sys.time(),sep=""))
               write_bin_files(n,PROJECT,latlon,timestep_days,met_all
                              ,lai_all,Csom_all,forest_all
                              ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                              ,sand_clay_all,crop_man_all,burnt_all,soilwater_all
                              ,nbe_all, lca_all,gpp_all,Cwood_inc_all,Cwood_mortality_all,fire_all
                              ,fapar_all)    

          } # site loop

      } # use_parallel

      # clean up to remove large drains on computer memeory
      rm(met_all,lai_all) ; gc(reset=TRUE,verbose=FALSE)

   } # Files were needed to be created

   # copy files to eddie?
   if (PROJECT$ecdf) {
       failed = TRUE ; copy_to = "y" ; failed = FALSE
       while(failed) {
          copy_to = readline("Copy binary driver files to cluster? (y/n)")
           if (copy_to != "y" & copy_to != "n") {failed = TRUE} else {failed = FALSE}
       }
       if (copy_to == "y") {
          #home_computer=Sys.info()["nodename"]
           # Check whether an existing cardamom_inputs.zip exists
           if (file.exists(paste(PROJECT$datapath,"cardamom_inputs.zip",sep=""))) {
               # Delete this file before creating a new one of the latest input files
               system(paste("rm ",PROJECT$datapath,"cardamom_inputs.zip", sep=""))
           }
           # Update the user
           print("Begin compressing *.bin files into zip directory for transfer to remote server...")
           # Compress all input files into zip directory
           system(paste("zip -j -r -q ",PROJECT$datapath,"cardamom_inputs.zip ",PROJECT$datapath," -i '*.bin'",sep=""))
           # Update the user
           print("...zip directory creation completed")
           # Copy the zip directory to the remote server
           command = paste("scp -r -q ",username,"@",home_computer,":",PROJECT$datapath,"cardamom_inputs.zip ",PROJECT$edatapath,sep="")
           # Unzip on remote server
           command = c(command,paste("unzip -o -qq ",PROJECT$edatapath,"cardamom_inputs.zip -d ",PROJECT$edatapath, sep=""))
           # Remove the zip directory on remote server
           command = c(command,paste("rm ",PROJECT$edatapath,"cardamom_inputs.zip" ,sep=""))
           #command = paste("scp -r ",username,"@",home_computer,":",PROJECT$datapath,"* ",PROJECT$edatapath,sep="")
           #print(command)
           # Update the user
           print("About to issue commands to remote server to transfer files")
           # Execute command on remote server
           ecdf_execute(command,PROJECT$paths$cardamom_cluster)
           # Delete local copy of the zip directory
           system(paste("rm ",PROJECT$datapath,"cardamom_inputs.zip", sep=""))
       }
   } # copy to Eddie

   return(paste("CARDAMOM Report: 1 completed", sep=""))

} # end function cardamom_stage_1

## Use byte compile
cardamom_stage_1<-cmpfun(cardamom_stage_1)
