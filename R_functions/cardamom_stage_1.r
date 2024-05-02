###
## Function to carry out stage 1 processes,
## i.e. creation of CARDAMOM input files
###

# Author: T. Luke Smallman (02/05/2024)

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
            if (exists(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")) == FALSE) {
                # Escape loop now we know some sites need making
                check = FALSE ; missing = append(missing,n)
            }
       }
       # If all files are present we can stop the current proess
       if (check) {return(paste("CARDAMOM Report: ",stage," completed", sep=""))}
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
           output = determine_lat_long_needed(lat = range(PROJECT$latitude)+c(-0.5,0.5), long = range(PROJECT$longitude)+c(-0.5,0.5)
                                             ,resolution = 0.125*0.5, grid_type = "wgs84", remove = NULL)
           cardamom_ext = output$cardamom_ext
           # Extract the lat / long grids designed for determining the extraction locations for the gridded observation datasets
           obs_long_grid = output$obs_long_grid ; obs_lat_grid = output$obs_lat_grid
           # Tidy up
           rm(output) ; gc(reset=TRUE,verbose=FALSE)
           print("Have now determined grid point locations")
       }

       # Load available data from gridded datasets?
       if (PROJECT$model$name != "ACM") {
           # if this is the first time of all creating new met files this time round load the whole dataset for rapid access
           met_all = load_met_fields_for_extraction(latlon,met_source,PROJECT$model$name,PROJECT$start_year,PROJECT$end_year,PROJECT$spatial_type,cardamom_ext)
           lai_all = load_lai_fields_for_extraction(latlon,lai_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           fapar_all = load_fapar_fields_for_extraction(latlon,fapar_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           nbe_all = load_nbe_fields_for_extraction(latlon,nbe_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           gpp_all = load_gpp_fields_for_extraction(latlon,GPP_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
           fire_all = load_fire_emission_fields_for_extraction(latlon,fire_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
           Csom_all = load_Csom_fields_for_extraction(latlon,Csom_source,cardamom_ext,PROJECT$spatial_type)
           crop_man_all = load_sacks_calendar_fields_for_extraction(latlon,crop_management_source)
           sand_clay_all = load_sand_clay_fields_for_extraction(latlon,sand_clay_source,cardamom_ext,PROJECT$spatial_type)
           forest_all = load_forestry_fields_for_extraction(latlon,deforestation_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
           Cwood_initial_all = load_initial_biomass_maps_for_extraction(latlon,Cwood_initial_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
           Cwood_stock_all = load_biomass_stocks_maps_for_extraction(latlon,Cwood_stock_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
           Cwood_potential_all = load_potential_biomass_maps_for_extraction(latlon,Cwood_potential_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
           burnt_all = load_burnt_area_fields_for_extraction(latlon,burnt_area_source,path_to_burnt_area,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
           soilwater_all = load_soilwater_fields_for_extraction(latlon,soilwater_initial_source)
           lca_all = load_lca_maps_for_extraction(latlon,lca_source,cardamom_ext,PROJECT$spatial_type)
           Cwood_inc_all = load_wood_productivity_maps_for_extraction(Cwood_inc_source,cardamom_ext,PROJECT$spatial_type,latlon,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days)
           Cwood_mortality_all = load_wood_mortality_maps_for_extraction(Cwood_mortality_source,cardamom_ext,PROJECT$spatial_type,latlon,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days)
       } # if (PROJECT$model$name != "ACM")

       # Update user
       print("Loading completed, beginning writing out file write out")

       if (use_parallel) {
           # Create function needed to process the site specific creation
           write_bin_files<-function(n) {

               # create the file name for the met/obs binary
               filename = paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")

               # All CARDAMOM read gridded datasets now map onto the same projection, extent and resolution.
               # This means that we can extract the location of the current site within any grid just the
               # once and pass around the solution to all the extraction functions.
               # NOTE: that met_all is an exception to this as the for the analysis to work we must always have
               # meteorology driving the model, so we assume a nearest neighbour approach of valid
               # locations, rather than accepting data gaps as done in observations / disturbance drivers.
               output = closest2d_2(1,obs_lat_grid,obs_long_grid,latlon[n,1],latlon[n,2])
               grid_long_loc = unlist(output, use.names=FALSE)[1] ; grid_lat_loc = unlist(output, use.names=FALSE)[2]
               rm(output)

               # load met drivers and obs; convert to daily if needed
               # these data are site specific so
               if (file.exists(filename) == FALSE | repair == 1){
                   # Load met drivers for ACM or other models
                   if (PROJECT$model$name != "ACM") {
                       met = extract_met_drivers(n,timestep_days,PROJECT$start_year,PROJECT$end_year,latlon[n,],met_all,met_source,PROJECT$sites[n])
                   } else { # if (PROJECT$model$name != "ACM")
                       # assume ACM special case
                       met = extract_acm_met_drivers(PROJECT,latlon[n,],PROJECT$sites[n])
                   }
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
                   binary_data(met,obs,filename,PROJECT$edc,latlon[n,],PROJECT$ctessel_pft[n],
                               PROJECT$model$name,PROJECT$parameter_type,PROJECT$model$nopars[n],noyears)
               } # if (file.exists(filename) == FALSE | repair == 1)
           } # end function

          # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
          # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only
          cl <- min(PROJECT$nosites,numWorkers)
          dummy = mclapply(c(1:PROJECT$nosites), FUN = write_bin_files, mc.cores = cl)
      } else { # use parallel
          # start looping through sites to create site specific files of obs and met
          for (n in seq(1, PROJECT$nosites)) {

               print(paste("Site ",n," of ",PROJECT$nosites," ",Sys.time(),sep=""))
               # create the file name for the met/obs binary
               filename=paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")

               # All CARDAMOM read gridded datasets now map onto the same projection, extent and resolution.
               # This means that we can extract the location of the current site within any grid just the
               # once and pass around the solution to all the extraction functions.
               # NOTE: that met_all is an exception to this as the for the analysis to work we must always have
               # meteorology driving the model, so we assume a nearest neighbour approach of valid
               # locations, rather than accepting data gaps as done in observations / disturbance drivers.
               output = closest2d_2(1,obs_lat_grid,obs_long_grid,latlon[n,1],latlon[n,2])
               grid_long_loc = unlist(output, use.names=FALSE)[1] ; grid_lat_loc = unlist(output, use.names=FALSE)[2]
               rm(output)

               # load met drivers and obs from CTESSEL and convert to daily if needed
               # these data are site specific so
               if (file.exists(filename) == FALSE | repair == 1){
                   # Load met drivers for ACM or other models
                   if (PROJECT$model$name != "ACM") {
                       met = extract_met_drivers(n,timestep_days,PROJECT$start_year,PROJECT$end_year,latlon[n,],met_all,met_source,PROJECT$sites[n])
                   } else { # if (PROJECT$model$name != "ACM")
                       # assume ACM special case
                       met = extract_acm_met_drivers(PROJECT,latlon[n,],PROJECT$sites[n])
                   }
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
                   binary_data(met,obs,filename,PROJECT$edc,latlon[n,],PROJECT$ctessel_pft[n],
                               PROJECT$model$name,PROJECT$parameter_type,PROJECT$model$nopars[n],noyears)
               } # if (file.exists(filename) == FALSE | repair == 1)
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
           # Compress all input files into zip directory
           system(paste("zip -j -r -q ",PROJECT$datapath,"cardamom_inputs.zip ",PROJECT$datapath," -i '*.bin'",sep=""))
           # Copy the zip directory to the remote server
           command = paste("scp -r -q ",username,"@",home_computer,":",PROJECT$datapath,"cardamom_inputs.zip ",PROJECT$edatapath,sep="")
           # Unzip on remote server
           command = c(command,paste("unzip -o -qq ",PROJECT$edatapath,"cardamom_inputs.zip -d ",PROJECT$edatapath, sep=""))
           # Remove the zip directory on remote server
           command = c(command,paste("rm ",PROJECT$edatapath,"cardamom_inputs.zip" ,sep=""))
           #command = paste("scp -r ",username,"@",home_computer,":",PROJECT$datapath,"* ",PROJECT$edatapath,sep="")
           print(command)
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
