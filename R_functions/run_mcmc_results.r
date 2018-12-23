
###
## Function to run CARDAMOM parameters via the chosen model
###

run_each_site<-function(n,PROJECT,stage,repair,grid_override,stage5modifiers) {

	if (stage == 5) {
		stage5name = as.vector(unlist(stage5modifiers))
		stage5name = paste(stage5name[1],stage5name[2],stage5name[3],stage5name[4],stage5name[5],stage5name[6],stage5name[7],stage5name[8],sep="_")
		outfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_",stage5name,".RData",sep="")
		outfile1 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_",stage5name,"_parameters.RData",sep="")
	} else {
		outfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")
		outfile1 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
		outfile2 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")
	}

	if (file.exists(outfile1) == FALSE | repair == 1) {
		# load only the desired latter fraction of the parameter vectors
		# output is order dimensions(npar+1,iter,chain)
		parameters=read_parameter_chains(PROJECT,n,3)

		# ok so if we actually have some parameters we will
		if (parameters[1] != -9999) {
			# test for convergence and whether or not there is any single chain which can be removed in they do not converge
			if (dim(parameters)[3] > 2) {
				converged=have_chains_converged(parameters)
				# if log-likelihood has passed then we are not interested
				if (converged[length(converged)] == "FAIL") {
					notconv = TRUE ; i = 1 ; max_likelihood = rep(NA, length.out=dim(parameters)[3]) ; CI90 = rep(NA,length.out=c(2))
					while (notconv){
						max_likelihood[i] = max(parameters[dim(parameters)[1],,i])
						converged = have_chains_converged(parameters[,,-i]) ; i = i + 1
						# if removing one of the chains get convergence then great
						if (converged[length(converged)] == "PASS") {
							# but we need to check for the possibility that the chain we have removed is actually better than the others
							CI90[1] = quantile(parameters[dim(parameters)[1],,(i-1)], prob=c(0.10)) ; CI90[2] = quantile(parameters[dim(parameters)[1],,-(i-1)], prob=c(0.90))
							# if the rejected chain is significantly better (at 90 % CI) than the converged chains then we have a problem
							if (CI90[1] > CI90[2]) {
								# rejected chain (while others converge) is actually better and the others have gotten stuck in a local minima.
								# we will now assume that we use the single good chain instead...
								parameters = array(parameters[,,(i-1)],dim=c(dim(parameters)[1:2],2))
								notconv = FALSE ; i = (i-1) * -1
							} else {
								# if the non-converged chain is worse or just the same in likelihood terms as the others then we will ditch it
								notconv = FALSE ; i = i-1 # converged now?
							}
						}
						# if we have tried removing each chain and we still have not converged then give up anyway
						#if (i > dim(parameters)[3] & notconv) {notconv=FALSE ; i=-9999}
						# or actually just remove the lowest average likelihood chain
						if (i > dim(parameters)[3] & notconv) {notconv=FALSE ; i=which(max_likelihood == min(max_likelihood)) ; i=i[1]}
					} # while to removing chains
					# if we successfully found only chain to remove then remove it from the rest of the analysis.
					if (i > 0) {parameters = parameters[,,-i] ; print(paste("chain rejected = ",i,sep=""))}
					if (i < 0) {print(paste("chain ",i*-1," only has been accepted",sep=""))}
				} # if likelihood not converged
			} # if more than 2 chains

			# load the met data for each site
			drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
			# if Stage 5 sensitivity analysis is enacted - apply factor scaled adjustment to the drivers
			if (stage == 5) {

				# do weather drivers first
				drivers$met[,2]  = drivers$met[,2]  * stage5modifiers$airt_factor # min temperature
				drivers$met[,3]  = drivers$met[,3]  * stage5modifiers$airt_factor # max temperature
				drivers$met[,4]  = drivers$met[,4]  * stage5modifiers$swrad_factor
				drivers$met[,5]  = drivers$met[,5]  * stage5modifiers$co2_factor
				drivers$met[,7]  = drivers$met[,7]  * stage5modifiers$rainfall_factor
				drivers$met[,15] = drivers$met[,15] * stage5modifiers$wind_spd_factor
				drivers$met[,16] = drivers$met[,16] * stage5modifiers$vpd_factor
				# do disturbance drivers next
				disturbed_locations = which(drivers$met[,8] > 0)
				drivers$met[disturbed_locations,8] = drivers$met[disturbed_locations,8] * stage5modifiers$deforestation_factor
				disturbed_locations = which(drivers$met[,9] > 0)
				drivers$met[disturbed_locations,9] = drivers$met[disturbed_locations,9] * stage5modifiers$burnt_area_factor
				# do GSI related parameters next
				drivers$met[,10] = drivers$met[,10] * stage5modifiers$airt_factor # 21 day rolling average daily minimum temperature

			}
			# now run the parameters to generate state variables
			# only use 500 for full propogation of states and fluxes
			if (dim(parameters)[2] > 500) {
				print("Note only 500 parameter sets PER CHAIN will be used in full propogation of parameter / flux / state uncertainties")
				sub_parameter = parameters[,sample(dim(parameters)[2],500,replace=FALSE),]
			} else {
				sub_parameter = parameters
			}
			# run subsample of parameters for full results / propogation
			soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
			states_all = simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,sub_parameter[1:PROJECT$model$nopars[n],,],
																drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
																PROJECT$exepath,soil_info)
			# pass to local variable for saving
			site_ctessel_pft = PROJECT$ctessel_pft[n]
			aNPP = states_all$aNPP
			# store the results now in binary file
			save(parameters,drivers,sub_parameter,site_ctessel_pft,aNPP,file=outfile1)#, compress="gzip")
			# determine whether this is a gridded run (or one with the override in place)
			if (PROJECT$spatial_type == "site" | grid_override == TRUE) {
				# ...if this is a site run save the full ensemble and everything else...
				save(parameters,drivers,states_all,sub_parameter,site_ctessel_pft,aNPP,file=outfile)#, compress="gzip")
			} else {
				# ...otherwise this is a grid and we want straight forward reduced dataset of common stocks and fluxes
				num_quantiles = c(0.025,0.25,0.5,0.75,0.975)
				# Stocks first
				site_output = list(labile_gCm2 = apply(states_all$lab,2,quantile,prob=num_quantiles),
				foliage_gCm2 = apply(states_all$fol,2,quantile,prob=num_quantiles),
				roots_gCm2 = apply(states_all$root,2,quantile,prob=num_quantiles),
				wood_gCm2 = apply(states_all$wood,2,quantile,prob=num_quantiles),
				lit_gCm2 = apply(states_all$lit,2,quantile,prob=num_quantiles),
				som_gCm2 = apply(states_all$som,2,quantile,prob=num_quantiles),
				cwd_gCm2 = apply(states_all$litwood,2,quantile,prob=num_quantiles),
				# Fluxes second
				gpp_gCm2day = apply(states_all$gpp,2,quantile,prob=num_quantiles),
				rauto_gCm2day = apply(states_all$rauto,2,quantile,prob=num_quantiles),
				rhet_gCm2day = apply(states_all$rhet,2,quantile,prob=num_quantiles),
				harvest_gCm2day = apply(states_all$harvest_C,2,quantile,prob=num_quantiles),
				fire_gCm2day = apply(states_all$fire,2,quantile,prob=num_quantiles))
				# Finally water cycle specific if available
				if (length(which(names(states_all) == "evap")) > 0) {
						# currently water in the soil surface layer (0-10 cm)
						site_output$SurfWater_kgH2Om2 = apply(states_all$rootwater,2,quantile,prob=num_quantiles)
						# plant apparent soil water potential (MPa)
						site_output$wSWP_kgH2Om2 = apply(states_all$wSWP,2,quantile,prob=num_quantiles)
						# evapotranspiration (Etrans + Esoil + Ewetcanopy)
						site_output$evap_kgH2Om2day = apply(states_all$evap,2,quantile,prob=num_quantiles)
				}
				# save to pixel specific file for the moment... in "run_mcmc_results" these will be combined into a single grid
				save(site_output,file=outfile2)
			}
			dummy = 0

		} else {

			dummy = -1

		} # parameters[1] != -9999

	} else {

		dummy = -1
		print('Already extracted result vectors (set repair = 1 if re-run is needed)')

	} # *parameters.RData already exists

	return(dummy)

} # end of run_each_site
## Use byte compile
run_each_site<-cmpfun(run_each_site)

run_mcmc_results <- function (PROJECT,stage,repair,grid_override) {

	print('Welcome to RUN_MCMC_RESULTS!!')

	# bundle needed functions down the chain
	functions_list=c("read_parameter_chains","read_binary_file_format","simulate_all",
	"read_binary_response_surface","crop_development_parameters","have_chains_converged","psrf")
	# start marker
	stime = proc.time()["elapsed"]

	# how many plots in total do we have
	nos_plots = 1:PROJECT$nosites

	# now check which ones we need to calculate, but only if override not in play
	if (repair != 1) {
		print("...beginning filterings for sites we have already processed")
		keep_list = 0
		for (i in seq(1, length(nos_plots))) {
			outfile1 = paste(PROJECT$results_processedpath,PROJECT$sites[i],"_parameters.RData",sep="")
			if (file.exists(outfile1) == FALSE) {keep_list=append(keep_list,i)}
		}
		# filter out the sites we already have then
		keep_list = keep_list[-1] ; print(paste("......removing ",length(nos_plots)-length(keep_list)," sites out of ",length(nos_plots)," from the analysis",sep=""))
		nos_plots = nos_plots[keep_list]
	}

	# Combine driver modifiers into a list to be passed to the run_each_site()
	stage5modifiers=list(airt_factor = airt_factor, swrad_factor = swrad_factor, co2_factor = co2_factor,
		rainfall_factor = rainfall_factor, wind_spd_factor = wind_spd_factor, vpd_factor = vpd_factor,
		deforestation_factor = deforestation_factor, burnt_area_factor = burnt_area_factor)

		# now request the creation of the plots
		if (use_parallel & length(nos_plots) > 1) {
			print("...beginning parallel operations")
			cl <- makeCluster(min(length(nos_plots),numWorkers), type = "PSOCK")
			clusterExport(cl,functions_list)
			# load R libraries in cluster
			clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
			dummy=parLapply(cl,nos_plots,fun=run_each_site,PROJECT=PROJECT,stage=stage,repair=repair,grid_override=grid_override,stage5modifiers=stage5modifiers)
			stopCluster(cl)
		} else {
			print("...beginning serial operations")
			# or use serial
			dummy=lapply(nos_plots,FUN=run_each_site,PROJECT=PROJECT,stage=stage,repair=repair,grid_override=grid_override,stage5modifiers=stage5modifiers)
		} # parallel option

		# now if this is a gridded run we want to take out individual site specific summary files and combine them into a single file
		if (PROJECT$spatial_type == "grid" & grid_override == FALSE) {

				# update the user
				print("...combing pixel level stock / flux output into single array")

				# determine what the output file name is here, so that we can check if one already exists
				outfile_grid = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
print(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
				if (file.exists(outfile_grid) == FALSE) {
print("arrays")
						# make a list of all the files we will be reading in
						to_do = list.files(PROJECT$results_processedpath, full.names=TRUE)
						to_do = to_do[grepl("_stock_fluxes",to_do)]
						# read in the first file so that we can then set up all the grids for combined summary output
						load(to_do[1])

						#
						# Generate the summary information
						# Time invariant but contain uncertainty, shaped into the full spatial grid
						# i.e. including areas not part of the analysis but within the spatail domain
						#
print(PROJECT$long_dim) ; print(PROJECT$lat_dim) ; print(dim(site_output$labile_gCm2)[1])
						# Stocks first
						grid_output = list(mean_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1])))
						grid_output$mean_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_lit_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_cwd_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						# Fluxes second
						grid_output$mean_gpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_rauto_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_rhet_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_harvest_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						grid_output$mean_fire_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))

						# Finally water cycle specific if available
						if (length(which(names(site_output) == "evap")) > 0) {
								# currently water in the soil surface layer (0-10 cm)
								grid_output$mean_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
								# plant apparent soil water potential (MPa)
								grid_output$mean_wSWP_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
								# evapotranspiration (Etrans + Esoil + Ewetcanopy)
								grid_output$mean_evap_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
						}

						#
						# Generate the variables needed to contain just the analysis locations
						# Fully time series provided along with uncertainty information but for the analysis pixels only
						# Each pixel will have a corresponding i and j location which related to its lat/long within the overall domain grid
						#

						# Stocks first
						grid_output$labile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$foliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$roots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$wood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$lit_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$som_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$cwd_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						# Fluxes second
						grid_output$gpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$rauto_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$rhet_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$harvest_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						grid_output$fire_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))

						# Finally water cycle specific if available
						if (length(which(names(site_output) == "evap")) > 0) {
								# currently water in the soil surface layer (0-10 cm)
								grid_output$SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
								# plant apparent soil water potential (MPa)
								grid_output$wSWP_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
								# evapotranspiration (Etrans + Esoil + Ewetcanopy)
								grid_output$evap_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
						}

						#
						# Generate the spatial information needed to relate i,j within the grid to long / lat and the summary vs detailed output
						#

						# vector into which the i,j position information within the grid will be stored
						grid_output$i_location = rep(NA, length.out = PROJECT$nosites)
						grid_output$j_location = rep(NA, length.out = PROJECT$nosites)

						# determine the latitude / longitude grid (dimension long = i, lat = j)
						grid_output$lat = seq(PROJECT$latitude[1],PROJECT$latitude[2],length.out=PROJECT$lat_dim)
						grid_output$long = seq(PROJECT$longitude[1],PROJECT$longitude[2],length.out=PROJECT$long_dim)
						grid_output$lat = array(grid_output$lat, dim=c(length(grid_output$lat),length(grid_output$long)))
						grid_output$lat = t(grid_output$lat)
						grid_output$long = array(grid_output$long, dim=dim(grid_output$lat))

				} else {

						# if file already exists load it up
						load(outfile_grid)

				} # have we already an output file

				for (n in seq(1,PROJECT$nosites)) {
						# load the current file
						outfile2 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")
						if (file.exists(outfile2)) {
								# load the site
								load(outfile2)
								# determine the lat / long location within the grid
					      slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
					      slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
					      if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
								# save for later
								grid_output$i_location[n] = slot_i ; grid_output$j_location[n] = slot_j

								# now assign to correct location in array
								# Stocks first
								grid_output$labile_gCm2[n,,] = site_output$labile_gCm2
								grid_output$foliage_gCm2[n,,] = site_output$foliage_gCm2
								grid_output$roots_gCm2[n,,] = site_output$roots_gCm2
								grid_output$wood_gCm2[n,,] = site_output$wood_gCm2
								grid_output$lit_gCm2[n,,] = site_output$lit_gCm2
								grid_output$som_gCm2[n,,] = site_output$som_gCm2
								grid_output$cwd_gCm2[n,,] = site_output$cwd_gCm2
								# Fluxes second
								grid_output$gpp_gCm2day[n,,] = site_output$gpp_gCm2day
								grid_output$rauto_gCm2day[n,,] = site_output$rauto_gCm2day
								grid_output$rhet_gCm2day[n,,] = site_output$rhet_gCm2day
								grid_output$harvest_gCm2day[n,,] = site_output$harvest_gCm2day
								grid_output$fire_gCm2day[n,,] = site_output$fire_gCm2day
								# Finally water cycle specific if available
								if (length(which(names(site_output) == "evap")) > 0) {
										# currently water in the soil surface layer (0-10 cm)
										grid_output$SurfWater_kgH2Om2[n,,] = site_output$SurfWater_kgH2Om2
										# plant apparent soil water potential (MPa)
										grid_output$wSWP_kgH2Om2[n,,] = site_output$wSWP_kgH2Om2
										# evapotranspiration (Etrans + Esoil + Ewetcanopy)
										grid_output$evap_kgH2Om2day[n,,] = site_output$evap_kgH2Om2day
								}

								# now assign to correct location in array
								# Stocks first
								grid_output$mean_labile_gCm2[slot_i,slot_j,] = apply(site_output$labile_gCm2,1,mean)
								grid_output$mean_foliage_gCm2[slot_i,slot_j,] = apply(site_output$foliage_gCm2,1,mean)
								grid_output$mean_roots_gCm2[slot_i,slot_j,] = apply(site_output$roots_gCm2,1,mean)
								grid_output$mean_wood_gCm2[slot_i,slot_j,] = apply(site_output$wood_gCm2,1,mean)
								grid_output$mean_lit_gCm2[slot_i,slot_j,] = apply(site_output$lit_gCm2,1,mean)
								grid_output$mean_som_gCm2[slot_i,slot_j,] = apply(site_output$som_gCm2,1,mean)
								grid_output$mean_cwd_gCm2[slot_i,slot_j,] = apply(site_output$cwd_gCm2,1,mean)
								# Fluxes second
								grid_output$mean_gpp_gCm2day[slot_i,slot_j,] = apply(site_output$gpp_gCm2day,1,mean)
								grid_output$mean_rauto_gCm2day[slot_i,slot_j,] = apply(site_output$rauto_gCm2day,1,mean)
								grid_output$mean_rhet_gCm2day[slot_i,slot_j,] = apply(site_output$rhet_gCm2day,1,mean)
								grid_output$mean_harvest_gCm2day[slot_i,slot_j,] = apply(site_output$harvest_gCm2day,1,mean)
								grid_output$mean_fire_gCm2day[slot_i,slot_j,] = apply(site_output$fire_gCm2day,1,mean)
								# Finally water cycle specific if available
								if (length(which(names(site_output) == "evap")) > 0) {
										# currently water in the soil surface layer (0-10 cm)
										grid_output$mean_SurfWater_kgH2Om2[slot_i,slot_j,] = apply(site_output$SurfWater_kgH2Om2,1,mean)
										# plant apparent soil water potential (MPa)
										grid_output$mean_wSWP_kgH2Om2[slot_i,slot_j,] = apply(site_output$wSWP_kgH2Om2,1,mean)
										# evapotranspiration (Etrans + Esoil + Ewetcanopy)
										grid_output$mean_evap_kgH2Om2day[slot_i,slot_j,] = apply(site_output$evap_kgH2Om2day,1,mean)
								}

								# now tidy away the file
								file.remove(outfile2)

						} # if site has been done and can now be made part of the overall...

				} # looping through sites

				# now save the combined grid file
				save(grid_output, file=outfile_grid)

		} # gridded run?

		# tell me whats happening
		print(paste("...time to process ",round((proc.time()["elapsed"]-stime)/60,1)," minutes",sep=""))

	} # end function run_mcmc_results
	## Use byte compile
	run_mcmc_results<-cmpfun(run_mcmc_results)
