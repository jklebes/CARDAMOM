
ensemble_within_range<-function(target,proposal) {

   # Determine what proportion of a proposed PDF is within a target range
   # Returned value 0-1

   t_range = range(target, na.rm = TRUE)
   in_range = length(which(proposal >= t_range[1] & proposal <= t_range[2]))
   return(in_range / length(proposal))

} # ensemble_within_range

rollapply_mean_annual_max<-function(var_in, step) {

   # Determine the mean of a maximums extracted from a rolling descrete windows
   # function specific command turned into a function to allow
   # for use within apply()
   return(mean(rollapply(var_in, by = step, width = step, FUN=max), na.rm = TRUE))

} # rollapply_mean_annual_max

rollapply_mean_annual<-function(var_in, step) {

   # Determine the mean of a maximums extracted from a rolling descrete windows
   # function specific command turned into a function to allow
   # for use within apply()
   return(rollapply(var_in, by = step, width = step, FUN=mean))

} # rollapply_mean_annual

###
## Function to run CARDAMOM parameters via the chosen model
###

run_each_site<-function(n,PROJECT,stage,repair,grid_override) {

  # Define the output file names
  outfile_site         = paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")
  outfile_parameters   = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
  outfile_stock_fluxes = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")

  # Set dummy output variable, the value may be changes by the code below
  dummy = -1

  if (file.exists(outfile_parameters) == FALSE | repair == 1) {
      # load only the desired latter fraction of the parameter vectors
      # output is order dimensions(npar+1,iter,chain)
      parameters = read_parameter_chains(PROJECT,n)
      parameter_covariance = read_parameter_covariance(PROJECT,n)

      # determine whether we have any actual completed chains and whether they include EDC consistent value only
      error_check = FALSE
      if (parameters[1] == -9999) {
          error_check = TRUE
      } else {
          if (length(which(as.vector(is.na(parameters)))) > 0 ) {
              error_check = TRUE
              print("NA found in likelihood score")
           } else if (min(as.vector(parameters)) == -Inf) {
              error_check = TRUE
              print("Inf found in likelihood score")
          } # NaN / Inf check
      } # error check

      # ok so if we actually have some parameters we will
      if (error_check == FALSE) {

          # test for convergence and whether or not there is any single chain which can be removed in they do not converge
          notconv = TRUE ; converged = rep("TRUE", times = max(PROJECT$model$nopars))
          while (dim(parameters)[3] > 2 & notconv) {
              print("begin convergence checking")
              converged = have_chains_converged(parameters)
              # if log-likelihood has passed then we are not interested
              if (converged[length(converged)] == "FAIL") {
                  print("...not converged begin removing potential parameter vectors")
                  i = 1 ; max_likelihood = rep(NA, length.out=dim(parameters)[3]) ; CI90 = rep(NA,length.out=c(2))
                  while (notconv){
                     print(paste("......trying removal of parameter vector ",i,sep=""))
                     # Track the maximum likelihood across each chain.
                     max_likelihood[i] = max(parameters[dim(parameters)[1],,i])
                     converged = have_chains_converged(parameters[,,-i]) ; i = i + 1
                     # if removing one of the chains get convergence then great
                     if (converged[length(converged)] == "PASS") {
                         print(".........convergence found on chain removal")
                         # likelihoods converge now but we need to check for the possibility that the chain we have removed is actually better than the others
                         CI90[1] = quantile(parameters[dim(parameters)[1],,(i-1)], prob=c(0.10)) ; CI90[2] = quantile(parameters[dim(parameters)[1],,-(i-1)], prob=c(0.90))
                         # if the rejected chain is significantly better (at 90 % CI) than the converged chains then we have a problem
                         if (CI90[1] > CI90[2]) {
                             # rejected chain (while others converge) is actually better and the others have gotten stuck in a local minima.
                             # we will now assume that we use the single good chain instead...
                             parameters = array(parameters[,,(i-1)],dim=c(dim(parameters)[1:2],2))
                             notconv = FALSE ; i = (i-1) * -1
                             print(paste("............chain ",i*-1," only has been accepted",sep=""))
                         } else {
                             # if the non-converged chain is worse or just the same in likelihood terms as the others then we will ditch it
                             notconv = FALSE ; i = i-1 # converged now?
                             parameters = parameters[,,-i]
                             print(paste("............chain rejected = ",i,sep=""))
                         }
                     }
                     # If removing one chain does not lead to convergence then lowest average likelihood chain could be removed
                     # NOTE: this should be made optional, such that non-converging locations are excluded from the analysis instead
                     if (i > dim(parameters)[3] & notconv) {
                         # Which is lowest likelihood
                         i = which(max_likelihood == min(max_likelihood))
                         # Remove from array
                         parameters = parameters[,,-i]
                         # Update the maximum likelihood vector also
                         max_likelihood = max_likelihood[-i]
                         # Update the user
                         print(paste(".........single chain removal couldn't find convergence; removing lowest likelihood chain = ",i,sep=""))
                         # reset counter
                         i = 1
                         # If we have removed chains down to 2 (or kept just one) we need to abort
                         if (dim(parameters)[3] < 3 & notconv) {
                             notconv = FALSE ; print(".........have removed all low likelihood chains without successful convergence")
                         }
                     }
                  } # while to removing chains
              } else {
                  print("All chains converge")
                  # we have conveged
                  notconv = FALSE
              } # if likelihood not converged
          } # if more than 2 chains

          # load the met data for each site
          drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))

          # run parameters for full results / propogation
          soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
          print("running model ensemble")
          states_all = simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,parameters[1:PROJECT$model$nopars[n],,],
                                    drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                                    PROJECT$exepath,soil_info)

          # Avoid running with ACM basically where not all fluxes exist
          if (grepl("DALEC",PROJECT$model$name)) {

              ###
              # Derive stocks and fluxes used in the calculation of gridded aggregates
              # These are variables which for a site analysis would be easy to calcuate
              # from the ensembles but difficult if not determined here and now before aggregation
              ###

              # If a combined ecosystem heterotrophic respiration flux does not
              # exist we shall calculate it
              if (exists(x = "rhet_gCm2day", where = states_all) == FALSE) {
                  if (exists(x = "rhet_dom_gCm2day", where = states_all)) {
                      states_all$rhet_gCm2day = states_all$rhet_dom_gCm2day
                  } else {
                      # Calculate the combined ecosystem heterotrophic respiration.
                      # All models have a som pool, so start with that
                      states_all$rhet_gCm2day = states_all$rhet_som_gCm2day
                      # If the model has a litter pool (foliar + fine root) add this
                      if (exists(x = "rhet_litter_gCm2day", where = states_all)) {
                         states_all$rhet_gCm2day = states_all$rhet_gCm2day + states_all$rhet_litter_gCm2day
                      }
                      # If the model has a wood litter pool add this
                      if (exists(x = "rhet_woodlitter_gCm2day", where = states_all)) {
                         states_all$rhet_gCm2day = states_all$rhet_gCm2day + states_all$rhet_woodlitter_gCm2day
                      }
                  } # does rhet_dom_gCm2day exist?
              } # does rhet_gCm2day exist?

              # Combine autotrophic and heterotrophic respiration into ecosystem respiration
              states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
              # Calculate the net ecosystem exchange of CO2
              states_all$nee_gCm2day = states_all$reco_gCm2day - states_all$gpp_gCm2day
              # Calculate net primary productivity
              states_all$npp_gCm2day = states_all$gpp_gCm2day - states_all$rauto_gCm2day
              # Begin calculation of net biome exchange and net biome productivity
              states_all$nbe_gCm2day = states_all$nee_gCm2day  # negative = sink
              states_all$nbp_gCm2day = -states_all$nee_gCm2day # positive = sink
              # If fire exists then update the NBE and NBP accordingly
              if (exists(x = "fire_gCm2day", where = states_all)) {
                  states_all$nbe_gCm2day = states_all$nbe_gCm2day + states_all$fire_gCm2day
                  states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$fire_gCm2day
              }
              # If a harvest flux exists update the NBP. NOTE: that this harvest flux
              # specifically accouts for C removed, there may be mortality due to harvest
              # but remains in system as residues.
              if (exists(x = "harvest_gCm2day", where = states_all)) {
                  states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$harvest_gCm2day
              }

              ###
              ## Post-hoc calculation of parameter correlations with key C-cycle variables

              tmp = t(array(as.vector(parameters[1:PROJECT$model$nopars[n],,]),dim=c(PROJECT$model$nopars[n],prod(dim(parameters)[2:3]))))
              states_all$nee_par_cor = cor(tmp,apply(states_all$nee_gCm2day,1,mean))
              states_all$gpp_par_cor = cor(tmp,apply(states_all$gpp_gCm2day,1,mean))
              states_all$rauto_par_cor = cor(tmp,apply(states_all$rauto_gCm2day,1,mean))
              states_all$rhet_par_cor = cor(tmp,apply(states_all$rhet_gCm2day,1,mean))
              # Avoid error flag when no fire
              if (max(as.vector(states_all$fire_gCm2day), na.rm = TRUE) > 0) {
                  states_all$fire_par_cor = cor(tmp,apply(states_all$fire_gCm2day,1,mean))
              } else {
                  states_all$fire_par_cor = array(0, dim = c(PROJECT$model$nopars[n],1))
              }

              ###
              ## Comparison with assimilated observation - to what extent does the ensemble overlap?

              ## GPP (gC/m2/day)
              obs_id = 1 ; unc_id = obs_id+1
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$gpp_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$gpp_gCm2day[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$gpp_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## LAI (m2/m2)
              obs_id = 3 ; unc_id = obs_id+1
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$lai_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$lai_m2m2[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$lai_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## NEE (gC/m2/day)
              obs_id = 5 ; unc_id = obs_id+1
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$nee_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$nee_gCm2day[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$nee_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## Wood (gC/m2)
              obs_id = 13 ; unc_id = obs_id+1
              # If there is a prior assign it to the first timestep of the observation timeseries
              if (drivers$parpriors[21] > 0) { drivers$obs[1,obs_id] = drivers$parpriors[21] ; drivers$obs[1,unc_id] = drivers$parpriorunc[21] }
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$wood_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$wood_gCm2[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$wood_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## Soil (gC/m2)
              obs_id = 15 ; unc_id = obs_id+1
              # If there is a prior assign it to the first timestep of the observation timeseries
              if (drivers$parpriors[23] > 0) { drivers$obs[1,obs_id] = drivers$parpriors[23] ; drivers$obs[1,unc_id] = drivers$parpriorunc[23] }
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$soil_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$soil_gCm2[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$soil_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## ET (kgH2O/m2/day)
              obs_id = 31 ; unc_id = obs_id+1
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$et_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$ET_kgH2Om2day[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$et_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## NBE (gC/m2/day)
              obs_id = 35 ; unc_id = obs_id+1 ; states_all$nbe_gCm2day = states_all$nee_gCm2day + states_all$fire_gCm2day
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$nbe_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$nbe_gCm2day[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$nbe_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

              ## Fire (gC/m2/day)
              obs_id = 7 ; unc_id = obs_id+1
              if (length(which(drivers$obs[,obs_id] != -9999)) > 0 & length(which(drivers$obs[,unc_id] != -9999)) > 0) {
                  # Loop through time to assess model overlap with observations
                  nobs = 0 ; states_all$fire_assim_data_overlap_fraction = 0
                  for (t in seq(1, length(drivers$met[,1]))) {
                       if (drivers$obs[t,obs_id] != -9999) {
                           # Estimate the min / max values for the observations
                           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
                           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_min,obs_max), m = states_all$fire_gCm2day[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           states_all$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction + tmp2
                           nobs = nobs + 1
                       }  # != -9999
                   } # time loop
                   # Average the overlap
                   if (nobs > 0) {
                       states_all$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction / nobs
                   } else {
                       states_all$fire_assim_data_overlap_fraction = 0
                   }
              } # was the obs assimilated?

          } # DALEC model or not?

          # pass to local variable for saving
          site_ctessel_pft = PROJECT$ctessel_pft[n]
          NPP_fraction = list(NPP_foliage_fraction = states_all$NPP_foliage_fraction)
          if (exists(x = "NPP_roots_fraction", where = states_all)) {NPP_fraction$NPP_roots_fraction = states_all$NPP_roots_fraction}
          if (exists(x = "NPP_wood_fraction", where = states_all)) {NPP_fraction$NPP_wood_fraction = states_all$NPP_wood_fraction}
          MTT_years = list(MTT_foliage_years = states_all$MTT_foliage_years)
          if (exists(x = "MTT_labile_years", where = states_all)) {MTT_years$MTT_labile_years = states_all$MTT_labile_years}
          if (exists(x = "MTT_roots_years", where = states_all)) {MTT_years$MTT_roots_years = states_all$MTT_roots_years}
          if (exists(x = "MTT_wood_years", where = states_all)) {MTT_years$MTT_wood_years = states_all$MTT_wood_years}
          if (exists(x = "MTT_litter_years", where = states_all)) {MTT_years$MTT_litter_years = states_all$MTT_litter_years}
          if (exists(x = "MTT_woodlitter_years", where = states_all)) {MTT_years$MTT_woodlitter_years = states_all$MTT_woodlitter_years}
          if (exists(x = "MTT_som_years", where = states_all)) {MTT_years$MTT_som_years = states_all$MTT_som_years}
          SS_gCm2 = list(SS_foliage_gCm2 = states_all$SS_foliage_gCm2)
          if (exists(x = "SS_labile_gCm2", where = states_all)) {SS_gCm2$SS_labile_gCm2 = states_all$SS_labile_gCm2}
          if (exists(x = "SS_roots_gCm2", where = states_all)) {SS_gCm2$SS_roots_gCm2 = states_all$SS_roots_gCm2}
          if (exists(x = "SS_wood_gCm2", where = states_all)) {SS_gCm2$SS_wood_gCm2 = states_all$SS_wood_gCm2}
          if (exists(x = "SS_litter_gCm2", where = states_all)) {SS_gCm2$SS_litter_gCm2 = states_all$SS_litter_gCm2}
          if (exists(x = "SS_woodlitter_gCm2", where = states_all)) {SS_gCm2$SS_woodlitter_gCm2 = states_all$SS_woodlitter_gCm2}
          if (exists(x = "SS_som_gCm2", where = states_all)) {SS_gCm2$SS_som_gCm2 = states_all$SS_som_gCm2}

          # Sanity check
          if (length(which(is.na(as.vector(NPP_fraction))) == TRUE) > 0) {
              print(paste("NA value found in NPP for site ",PROJECT$site[n],sep="")) ; dummy = -1 ; return(dummy)
          }
          print("processing and storing ensemble output")
          # store the results now in binary file
          save(parameter_covariance,parameters,drivers,site_ctessel_pft,NPP_fraction,MTT_years,SS_gCm2,
               file=outfile_parameters, compress="gzip", compression_level = 9)
          # determine whether this is a gridded run (or one with the override in place)
          if (PROJECT$spatial_type == "site" | grid_override == TRUE) {
              # ...if this is a site run save the full ensemble and everything else...
              save(parameters,drivers,states_all,site_ctessel_pft,file=outfile_site, compress="gzip", compression_level = 9)
          } else {
              # ...otherwise this is a grid and we want straight forward reduced dataset of common stocks and fluxes
              num_quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975) ; num_quantiles_agg = seq(0.0,1, length = 100)
              na_flag = TRUE

              # Determine some useful information for the analysis below
              nos_years = (as.numeric(PROJECT$end_year) - as.numeric(PROJECT$start_year))+1
              steps_per_year = floor(dim(drivers$met)[1] / nos_years)

              # Declare the site level output list object
              site_output = list(num_quantiles = num_quantiles)

              ###
              # Derive stocks and fluxes used in the calculation of gridded aggregates
              # These are variables which for a site analysis would be easy to calcuate
              # from the ensembles but difficult if not determined here and now before aggregation
              ###

              # The total allocation of C to the foliage pool can, depending on model,
              # be the combined total of direct allocation and that via a labile pool.
              # For many comparison we will need their combined total.
              if (exists(x = "alloc_foliage_gCm2day", where = states_all) &
                  exists(x = "labile_to_foliage_gCm2day", where = states_all)) {
                  states_all$combined_alloc_foliage_gCm2day = states_all$alloc_foliage_gCm2day + states_all$labile_to_foliage_gCm2day
              } else {
                  if (exists(x = "labile_to_foliage_gCm2day", where = states_all)) {
                      states_all$combined_alloc_foliage_gCm2day = states_all$labile_to_foliage_gCm2day
                  } else if (exists(x = "alloc_foliage_gCm2day", where = states_all)) {
                      states_all$combined_alloc_foliage_gCm2day = states_all$alloc_foliage_gCm2day
                  } else {
                      stop("Error, CARDAMOM cannnot determine where C allocation foliage has come from.")
                  }
              }

              # Determine the total biomass within the system
              # All models have a foliage pool
              states_all$biomass_gCm2 = states_all$foliage_gCm2
              states_all$biomass_to_litter_gCm2day = states_all$foliage_to_litter_gCm2day
              # Now look for accumulating optional disturbance drivers
              if (exists(x = "FIREemiss_foliage_gCm2day", where = states_all)) {
                  states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_foliage_gCm2day
              }
              if (exists(x = "FIRElitter_foliage_gCm2day", where = states_all)) {
                  states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_foliage_gCm2day
              }
              if (exists(x = "HARVESTextracted_foliage_gCm2day", where = states_all)) {
                  states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_foliage_gCm2day
                  states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_foliage_gCm2day
              }
              if (exists(x = "labile_gCm2", where = states_all)) {
                  states_all$biomass_gCm2 = states_all$biomass_gCm2 + states_all$labile_gCm2
                  # Now look for accumulating optional disturbance drivers
                  if (exists(x = "FIREemiss_labile_gCm2day", where = states_all)) {
                      states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_biomass_gCm2day + states_all$FIREemiss_labile_gCm2day
                  }
                  if (exists(x = "FIRElitter_labile_gCm2day", where = states_all)) {
                      states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_biomass_gCm2day + states_all$FIRElitter_labile_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_labile_gCm2day", where = states_all)) {
                      states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTextracted_labile_gCm2day
                      states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_biomass_gCm2day + states_all$HARVESTlitter_labile_gCm2day
                  }
              }
              if (exists(x = "roots_gCm2", where = states_all)) {
                  states_all$biomass_gCm2 = states_all$biomass_gCm2 + states_all$roots_gCm2
                  states_all$biomass_to_litter_gCm2day = states_all$biomass_to_litter_gCm2day + states_all$roots_to_litter_gCm2day
                  # Now look for accumulating optional disturbance drivers
                  if (exists(x = "FIREemiss_roots_gCm2day", where = states_all)) {
                      states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_biomass_gCm2day + states_all$FIREemiss_roots_gCm2day
                  }
                  if (exists(x = "FIRElitter_roots_gCm2day", where = states_all)) {
                      states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_biomass_gCm2day + states_all$FIRElitter_roots_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_roots_gCm2day", where = states_all)) {
                      states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTextracted_roots_gCm2day
                      states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_biomass_gCm2day + states_all$HARVESTlitter_roots_gCm2day
                  }

              }
              if (exists(x = "wood_gCm2", where = states_all)) {
                  states_all$biomass_gCm2 = states_all$biomass_gCm2 + states_all$wood_gCm2
                  states_all$biomass_to_litter_gCm2day = states_all$biomass_to_litter_gCm2day + states_all$wood_to_litter_gCm2day
                  # Now look for accumulating optional disturbance drivers
                  if (exists(x = "FIREemiss_wood_gCm2day", where = states_all)) {
                      states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_biomass_gCm2day + states_all$FIREemiss_wood_gCm2day
                  }
                  if (exists(x = "FIRElitter_wood_gCm2day", where = states_all)) {
                      states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_biomass_gCm2day + states_all$FIRElitter_wood_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_wood_gCm2day", where = states_all)) {
                      states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTextracted_wood_gCm2day
                      states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_biomass_gCm2day + states_all$HARVESTlitter_wood_gCm2day
                  }
              } # biomass

              # Determine the total dead organic matter within the system
              # plus their associated loss terms.
              # All models have a som pool so start there
              if (exists(x = "dom_gCm2", where = states_all) == FALSE) {
                  states_all$dom_gCm2 = states_all$som_gCm2
                  states_all$rhet_dom_gCm2day = states_all$rhet_som_gCm2day
                  # Create corresponding variables for the disturbance drivers.
                  # NOTE: that we do not need to account for disturbance induced
                  # litter fluxes as these always remain within the DOM pool
                  if (exists(x = "FIREemiss_som_gCm2day", where = states_all)) {
                      states_all$FIREemiss_dom_gCm2day = states_all$FIREemiss_som_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_som_gCm2day", where = states_all)) {
                      states_all$HARVESTextracted_dom_gCm2day = states_all$HARVESTextracted_som_gCm2day
                  }
                  # Check for presence of litter or wood litter pools
                  if (exists(x = "litter_gCm2", where = states_all)) {
                      states_all$dom_gCm2 = states_all$dom_gCm2 + states_all$litter_gCm2
                      states_all$rhet_dom_gCm2day = states_all$rhet_dom_gCm2day + states_all$rhet_litter_gCm2day
                      if (exists(x = "FIREemiss_litter_gCm2day", where = states_all)) {
                          states_all$FIREemiss_dom_gCm2day = states_all$FIREemiss_dom_gCm2day + states_all$FIREemiss_litter_gCm2day
                      }
                      if (exists(x = "HARVESTextracted_litter_gCm2day", where = states_all)) {
                          states_all$HARVESTextracted_dom_gCm2day = states_all$HARVESTextracted_dom_gCm2day + states_all$HARVESTextracted_litter_gCm2day
                      }
                  } # litter
                  if (exists(x = "woodlitter_gCm2", where = states_all)) {
                      states_all$dom_gCm2 = states_all$dom_gCm2 + states_all$woodlitter_gCm2
                      states_all$rhet_dom_gCm2day = states_all$rhet_dom_gCm2day + states_all$rhet_woodlitter_gCm2day
                      if (exists(x = "FIREemiss_woodlitter_gCm2day", where = states_all)) {
                          states_all$FIREemiss_dom_gCm2day = states_all$FIREemiss_dom_gCm2day + states_all$FIREemiss_woodlitter_gCm2day
                      }
                      if (exists(x = "HARVESTextracted_woodlitter_gCm2day", where = states_all)) {
                          states_all$HARVESTextracted_dom_gCm2day = states_all$HARVESTextracted_dom_gCm2day + states_all$HARVESTextracted_woodlitter_gCm2day
                      }
                  } # wood litter
              } # DOM pool exists

              # Determine the total C within the system
              # Biomass and dead organic matter have already been determined,
              # so combine.
              states_all$Ctotal_gCm2 = states_all$biomass_gCm2 + states_all$dom_gCm2

              ###
              # Extract likelihood / parameter / driver information
              ###

              # loop through parameters + likelihood
              site_output$parameters = array(NA, dim=c(dim(parameters)[1],length(num_quantiles)))
              for (p in seq(1, dim(parameters)[1])) {
                   site_output$parameters[p,] = quantile(as.vector(parameters[p,,]), prob=num_quantiles)
              }
              # track which parameters have converged + likelihood
              site_output$parameters_converged = rep(0, dim(parameters)[1])
              site_output$parameters_converged[which(converged == "PASS")] = 1
              # Mean meteorological conditions
              site_output$mean_temperature_C = mean((drivers$met[,3]+drivers$met[,2])*0.5)
              site_output$mean_radiation_MJm2day = mean(drivers$met[,4])
              site_output$mean_vpd_Pa = mean(drivers$met[,16])
              site_output$mean_precipitation_kgH2Om2yr = mean(drivers$met[,7]*86400*365.25)
              # Assimilated LAI information
              if (max(drivers$obs[,3]) > 0) {
                  site_output$assimilated_lai_max_m2m2 = max(drivers$obs[,3])
                  site_output$assimilated_lai_mean_m2m2 = mean(drivers$obs[which(drivers$obs[,3] != -9999),3])
                  site_output$assimilated_lai_sd_m2m2 = sd(drivers$obs[which(drivers$obs[,3] != -9999),3])
                  site_output$assimilated_lai_unc_m2m2 = mean(drivers$obs[which(drivers$obs[,4] != -9999),4])
              } else {
                  site_output$assimilated_lai_max_m2m2 = NA
                  site_output$assimilated_lai_mean_m2m2 = NA
                  site_output$assimilated_lai_sd_m2m2 = NA
                  site_output$assimilated_lai_unc_m2m2 = NA
              }
              # Assimilated wood stock / prior information
              # Do we have one or both wood stock prior and time series inforamtion
              if (drivers$parpriors[21] > 0 |
                  length(which(drivers$obs[,13] > -9999)) > 0) {
                  # Initialise the variable with zero to allow averaging
                  site_output$assimilated_wood_mean_gCm2 = 0
                  site_output$assimilated_wood_mean_unc_gCm2 = 0
                  # Do we have time series information?
                  if (length(which(drivers$obs[,13] > -9999)) > 0) { # Does time series info exist
                      site_output$assimilated_wood_mean_gCm2 = append(site_output$assimilated_wood_mean_gCm2,
                                                                      drivers$obs[which(drivers$obs[,13] != -9999),13])
                      site_output$assimilated_wood_mean_unc_gCm2 = append(site_output$assimilated_wood_mean_unc_gCm2,
                                                                          drivers$obs[which(drivers$obs[,14] != -9999),14])
                  }
                  # Do we have a prior estimate on initial wood stocks?
                  if (drivers$parpriors[21] > 0) { # Does a prior exist
                      site_output$assimilated_wood_mean_gCm2 = append(site_output$assimilated_wood_mean_gCm2,
                                                                      drivers$parpriors[21])
                      site_output$assimilated_wood_mean_unc_gCm2 = append(site_output$assimilated_wood_mean_unc_gCm2,
                                                                          drivers$parpriorunc[21])
                  }
                  # Average available observation and uncertainty information
                  # remembering to remove the initial value
                  site_output$assimilated_wood_mean_gCm2 = mean(site_output$assimilated_wood_mean_gCm2[-1], na.rm = na_flag)
                  site_output$assimilated_wood_mean_unc_gCm2 = mean(site_output$assimilated_wood_mean_unc_gCm2[-1], na.rm = na_flag)
              } else {
                 # assign missing value flag for consistency
                 site_output$assimilated_wood_mean_gCm2 = NA
                 site_output$assimilated_wood_mean_unc_gCm2 = NA
              } # wood stock or prior information was assimilated
              # Assimilated som stock / prior information
              # Do we have one or both som stock prior and time series inforamtion
              if (drivers$parpriors[23] > 0 |
                  length(which(drivers$obs[,19] > -9999)) > 0) {
                  # Initialise the variable with zero to allow averaging
                  site_output$assimilated_som_mean_gCm2 = 0
                  site_output$assimilated_som_mean_unc_gCm2 = 0
                  # Do we have time series information?
                  if (length(which(drivers$obs[,19] > -9999)) > 0) { # Does time series info exist
                      site_output$assimilated_som_mean_gCm2 = append(site_output$assimilated_som_mean_gCm2,
                                                                     drivers$obs[which(drivers$obs[,19] != -9999),19])
                      site_output$assimilated_som_mean_unc_gCm2 = append(site_output$assimilated_som_mean_unc_gCm2,
                                                                         drivers$obs[which(drivers$obs[,20] != -9999),20])
                  }
                  # Do we have a prior estimate on initial wood stocks?
                  if (drivers$parpriors[23] > 0) { # Does a prior exist
                      site_output$assimilated_som_mean_gCm2 = append(site_output$assimilated_som_mean_gCm2,
                                                                     drivers$parpriors[23])
                      site_output$assimilated_som_mean_unc_gCm2 = append(site_output$assimilated_som_mean_unc_gCm2,
                                                                         drivers$parpriorunc[23])
                  }
                  # Average available observation and uncertainty information
                  # remembering to remove the initial value
                  site_output$assimilated_som_mean_gCm2 = mean(site_output$assimilated_som_mean_gCm2[-1], na.rm = na_flag)
                  site_output$assimilated_som_mean_unc_gCm2 = mean(site_output$assimilated_som_mean_unc_gCm2[-1], na.rm = na_flag)
              } else {
                 # assign missing value flag for consistency
                 site_output$assimilated_som_mean_gCm2 = NA
                 site_output$assimilated_som_mean_unc_gCm2 = NA
              } # som stock or prior information was assimilated

              ###
              # Begin aggregating the states_all information into
              # quantile based site_output
              ###

              # Net primary production allocation fractions
              if (exists(x = "NPP_foliage_fraction", where = states_all)) {site_output$NPP_foliage_fraction = quantile(states_all$NPP_foliage_fraction, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "NPP_roots_fraction", where = states_all)) {site_output$NPP_roots_fraction = quantile(states_all$NPP_roots_fraction, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "NPP_wood_fraction", where = states_all)) {site_output$NPP_wood_fraction = quantile(states_all$NPP_wood_fraction, prob=num_quantiles, na.rm = na_flag)}
              # Analysis mean transit (residence) times (years)
              if (exists(x = "MTT_labile_years", where = states_all)) {site_output$MTT_labile_years = quantile(states_all$MTT_labile_years, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "MTT_foliage_years", where = states_all)) {site_output$MTT_foliage_years = quantile(states_all$MTT_foliage_years, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "MTT_roots_years", where = states_all)) {site_output$MTT_roots_years = quantile(states_all$MTT_roots_years, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "MTT_wood_years", where = states_all)) {site_output$MTT_wood_years = quantile(states_all$MTT_wood_years, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "MTT_litter_years", where = states_all)) {site_output$MTT_litter_years = quantile(states_all$MTT_litter_years, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "MTT_woodlitter_years", where = states_all)) {site_output$MTT_woodlitter_years = quantile(states_all$MTT_woodlitter_years, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "MTT_som_years", where = states_all)) {site_output$MTT_som_years = quantile(states_all$MTT_som_years, prob=num_quantiles, na.rm = na_flag)}
              # Steady state C stock estimates (gC/m2)
              if (exists(x = "SS_labile_gCm2", where = states_all)) {site_output$SS_labile_gCm2 = quantile(states_all$SS_labile_gCm2, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "SS_foliage_gCm2", where = states_all)) {site_output$SS_foliage_gCm2 = quantile(states_all$SS_foliage_gCm2, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "SS_roots_gCm2", where = states_all)) {site_output$SS_roots_gCm2 = quantile(states_all$SS_roots_gCm2, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "SS_wood_gCm2", where = states_all)) {site_output$SS_wood_gCm2 = quantile(states_all$SS_wood_gCm2, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "SS_litter_gCm2", where = states_all)) {site_output$SS_litter_gCm2 = quantile(states_all$SS_litter_gCm2, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "SS_woodlitter_gCm2", where = states_all)) {site_output$SS_woodlitter_gCm2 = quantile(states_all$SS_woodlitter_gCm2, prob=num_quantiles, na.rm = na_flag)}
              if (exists(x = "SS_som_gCm2", where = states_all)) {site_output$SS_som_gCm2 = quantile(states_all$SS_som_gCm2, prob=num_quantiles, na.rm = na_flag)}

              # State variables - NOTE: extraction of pixel specific means done here to account for different ensemble trajectories, i.e. correlation in time.
              site_output$lai_m2m2         = apply(states_all$lai_m2m2,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_lai_m2m2    = quantile(apply(states_all$lai_m2m2,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$Ctotal_gCm2      = apply(states_all$Ctotal_gCm2,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_Ctotal_gCm2 = quantile(apply(states_all$Ctotal_gCm2,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              # Fluxes second
              site_output$gpp_gCm2day          = apply(states_all$gpp_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_gpp_gCm2day     = quantile(apply(states_all$gpp_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$rauto_gCm2day        = apply(states_all$rauto_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_rauto_gCm2day   = quantile(apply(states_all$rauto_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$rhet_gCm2day         = apply(states_all$rhet_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_rhet_gCm2day    = quantile(apply(states_all$rhet_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$nee_gCm2day          = apply(states_all$nee_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_nee_gCm2day     = quantile(apply(states_all$nee_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$npp_gCm2day          = apply(states_all$npp_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_npp_gCm2day     = quantile(apply(states_all$npp_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$nbe_gCm2day          = apply(states_all$nbe_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_nbe_gCm2day     = quantile(apply(states_all$nbe_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$nbp_gCm2day          = apply(states_all$nbp_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_nbp_gCm2day     = quantile(apply(states_all$nbp_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$reco_gCm2day         = apply(states_all$reco_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_reco_gCm2day    = quantile(apply(states_all$reco_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$harvest_gCm2day      = apply(states_all$harvest_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_harvest_gCm2day = quantile(apply(states_all$harvest_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
              site_output$fire_gCm2day         = apply(states_all$fire_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
              site_output$mean_fire_gCm2day    = quantile(apply(states_all$fire_gCm2day,1,mean,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)

              ###
              # Track net pool change over time
              ###

              # Start with the common variables
              dCbio = states_all$Ctotal_gCm2 - states_all$Ctotal_gCm2[,1] # difference in total C from initial
              site_output$dCtotal_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

              dCbio = states_all$lai_m2m2 - states_all$lai_m2m2[,1] #  difference in lai from initial
              site_output$dlai_m2m2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

              # Then do model specific pool combinations
              # Include the pool, its net change, the allocation of C from GPP / NPP direct
              # and its output flow(s).

              # Biomass related pool, change and output variables
              if (exists(x = "biomass_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$biomass_gCm2 = apply(states_all$biomass_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_biomass_gCm2 = quantile(apply(states_all$biomass_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$biomass_gCm2 - states_all$biomass_gCm2[,1] # difference in labile from initial
                  site_output$dCbiomass_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Track natural biomass losses
                  site_output$biomass_to_litter_gCm2day = apply(states_all$biomass_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_biomass_to_litter_gCm2day = quantile(apply(states_all$biomass_to_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Begin accumulating the total output
                  site_output$outflux_biomass_gCm2day = states_all$biomass_to_litter_gCm2day
                  site_output$NaturalFractionOfTurnover_biomass = states_all$biomass_to_litter_gCm2day
                  # Other natural flux pathways should really go here before disturbance related
                  if (exists(x = "FIREemiss_biomass_gCm2day", where = states_all)) {
                      site_output$FIREemiss_biomass_gCm2day = apply(states_all$FIREemiss_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_biomass_gCm2day = quantile(apply(states_all$FIREemiss_biomass_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_biomass_gCm2day = site_output$outflux_biomass_gCm2day + states_all$FIREemiss_biomass_gCm2day
                      site_output$FireFractionOfTurnover_biomass = states_all$FIREemiss_biomass_gCm2day
                  }
                  if (exists(x = "FIRElitter_biomass_gCm2day", where = states_all)) {
                      site_output$FIRElitter_biomass_gCm2day = apply(states_all$FIRElitter_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_biomass_gCm2day = quantile(apply(states_all$FIRElitter_biomass_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_biomass_gCm2day = site_output$outflux_biomass_gCm2day + states_all$FIRElitter_biomass_gCm2day
                      site_output$FireFractionOfTurnover_biomass = site_output$FireFractionOfTurnover_biomass + states_all$FIRElitter_biomass_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_biomass_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_biomass_gCm2day = apply(states_all$HARVESTextracted_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_biomass_gCm2day = quantile(apply(states_all$HARVESTextracted_biomass_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$HARVESTlitter_biomass_gCm2day = apply(states_all$HARVESTlitter_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTlitter_biomass_gCm2day = quantile(apply(states_all$HARVESTlitter_biomass_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_biomass_gCm2day = site_output$outflux_biomass_gCm2day + states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTlitter_biomass_gCm2day
                      site_output$HarvestFractionOfTurnover_biomass = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTlitter_biomass_gCm2day
                  }
                  # Estimate the ecosystem mean transit (residence) times as a function of natural and disturbance processes
                  site_output$MTT_annual_biomass_years = apply(site_output$outflux_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_biomass_years = apply(states_all$biomass_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_biomass_years * 365.25)
                  site_output$MTT_annual_biomass_years = t(site_output$MTT_annual_biomass_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_biomass_years = apply(site_output$MTT_annual_biomass_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_biomass = quantile(apply(site_output$NaturalFractionOfTurnover_biomass,1,mean, na.rm = na_flag) /
                                                                           apply(site_output$outflux_biomass_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_biomass = quantile(apply(site_output$FireFractionOfTurnover_biomass,1,mean, na.rm = na_flag) /
                                                                        apply(site_output$outflux_biomass_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_biomass = quantile(apply(site_output$HarvestFractionOfTurnover_biomass,1,mean, na.rm = na_flag) /
                                                                           apply(site_output$outflux_biomass_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from biomass
                  site_output$mean_outflux_biomass_gCm2day = quantile(apply(site_output$outflux_biomass_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_biomass_gCm2day = apply(site_output$outflux_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }
              # Labile related pool, change, input and output variables
              if (exists(x = "labile_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$labile_gCm2 = apply(states_all$labile_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_labile_gCm2 = quantile(apply(states_all$labile_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$labile_gCm2 - states_all$labile_gCm2[,1] # difference in labile from initial
                  site_output$dClabile_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Determine the allocation to labile - in all cases this must be a direct variable
                  site_output$alloc_labile_gCm2day = apply(states_all$alloc_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_alloc_labile_gCm2day = quantile(apply(states_all$alloc_labile_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Check for the possible loss pathways
                  if (exists(x = "labile_to_foliage_gCm2day", where = states_all)) {
                      site_output$labile_to_foliage_gCm2day = apply(states_all$labile_to_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_labile_to_foliage_gCm2day = quantile(apply(states_all$labile_to_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$NaturalFractionOfTurnover_labile = states_all$labile_to_foliage_gCm2day
                      # Begin accumulating the total output
                      site_output$outflux_labile_gCm2day = states_all$labile_to_foliage_gCm2day
                  }
                  # Other natural flux pathways should really go here before disturbance related
                  if (exists(x = "FIREemiss_labile_gCm2day", where = states_all)) {
                      site_output$FIREemiss_labile_gCm2day = apply(states_all$FIREemiss_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_labile_gCm2day = quantile(apply(states_all$FIREemiss_labile_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_labile = states_all$FIREemiss_labile_gCm2day
                      site_output$outflux_labile_gCm2day = site_output$outflux_labile_gCm2day + states_all$FIREemiss_labile_gCm2day
                  }
                  if (exists(x = "FIRElitter_labile_gCm2day", where = states_all)) {
                      site_output$FIRElitter_labile_gCm2day = apply(states_all$FIRElitter_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_labile_gCm2day = quantile(apply(states_all$FIRElitter_labile_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_labile_gCm2day = site_output$outflux_labile_gCm2day + states_all$FIRElitter_labile_gCm2day
                      site_output$FireFractionOfTurnover_labile = site_output$FireFractionOfTurnover_labile + states_all$FIRElitter_labile_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_labile_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_labile_gCm2day = apply(states_all$HARVESTextracted_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_labile_gCm2day = quantile(apply(states_all$HARVESTextracted_labile_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$HARVESTlitter_labile_gCm2day = apply(states_all$HARVESTlitter_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTlitter_labile_gCm2day = quantile(apply(states_all$HARVESTlitter_labile_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_labile_gCm2day = site_output$outflux_labile_gCm2day + states_all$HARVESTextracted_labile_gCm2day + states_all$HARVESTlitter_labile_gCm2day
                      site_output$HarvestFractionOfTurnover_labile = states_all$HARVESTextracted_labile_gCm2day + states_all$HARVESTlitter_labile_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_labile_years = apply(site_output$outflux_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_labile_years = apply(states_all$labile_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_labile_years * 365.25)
                  site_output$MTT_annual_labile_years = t(site_output$MTT_annual_labile_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_labile_years = apply(site_output$MTT_annual_labile_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_labile = quantile(apply(site_output$NaturalFractionOfTurnover_labile,1,mean, na.rm = na_flag) /
                                                                          apply(site_output$outflux_labile_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_labile = quantile(apply(site_output$FireFractionOfTurnover_labile,1,mean, na.rm = na_flag) /
                                                                       apply(site_output$outflux_labile_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_labile = quantile(apply(site_output$HarvestFractionOfTurnover_labile,1,mean, na.rm = na_flag) /
                                                                          apply(site_output$outflux_labile_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from labile
                  site_output$mean_outflux_labile_gCm2day = quantile(apply(site_output$outflux_labile_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_labile_gCm2day = apply(site_output$outflux_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Foliage related pool, change, input and output variables
              if (exists(x = "foliage_gCm2", where = states_all)) {
                  # A combined total of C to foliage must always exist
                  site_output$combined_alloc_foliage_gCm2day = apply(states_all$combined_alloc_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_combined_alloc_foliage_gCm2day = quantile(apply(states_all$foliage_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Assign pool to site_output
                  site_output$foliage_gCm2 = apply(states_all$foliage_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_foliage_gCm2 = quantile(apply(states_all$foliage_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$foliage_gCm2 - states_all$foliage_gCm2[,1] # difference in root from initial
                  site_output$dCfoliage_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Check for the possible pathways
                  if (exists(x = "alloc_foliage_gCm2day", where = states_all)) {
                      site_output$alloc_foliage_gCm2day = apply(states_all$alloc_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_alloc_foliage_gCm2day = quantile(apply(states_all$alloc_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  }
                  if (exists(x = "foliage_to_litter_gCm2day", where = states_all)) {
                      site_output$foliage_to_litter_gCm2day = apply(states_all$foliage_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_foliage_to_litter_gCm2day = quantile(apply(states_all$foliage_to_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$NaturalFractionOfTurnover_foliage = states_all$foliage_to_litter_gCm2day
                      # Begin accumulating the total output fluxes here
                      site_output$outflux_foliage_gCm2day = states_all$foliage_to_litter_gCm2day
                  }
                  if (exists(x = "FIREemiss_foliage_gCm2day", where = states_all)) {
                      site_output$FIREemiss_foliage_gCm2day = apply(states_all$FIREemiss_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_foliage_gCm2day = quantile(apply(states_all$FIREemiss_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_foliage = states_all$FIREemiss_foliage_gCm2day
                      site_output$outflux_foliage_gCm2day = site_output$outflux_foliage_gCm2day + states_all$FIREemiss_foliage_gCm2day
                  }
                  if (exists(x = "FIRElitter_foliage_gCm2day", where = states_all)) {
                      site_output$FIRElitter_foliage_gCm2day = apply(states_all$FIRElitter_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_foliage_gCm2day = quantile(apply(states_all$FIRElitter_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_foliage_gCm2day = site_output$outflux_foliage_gCm2day + states_all$FIRElitter_foliage_gCm2day
                      site_output$FireFractionOfTurnover_foliage = site_output$FireFractionOfTurnover_foliage + states_all$FIRElitter_foliage_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_foliage_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_foliage_gCm2day = apply(states_all$HARVESTextracted_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_foliage_gCm2day = quantile(apply(states_all$HARVESTextracted_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$HARVESTlitter_foliage_gCm2day = apply(states_all$HARVESTlitter_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTlitter_foliage_gCm2day = quantile(apply(states_all$HARVESTlitter_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_foliage_gCm2day = site_output$outflux_foliage_gCm2day + states_all$HARVESTextracted_foliage_gCm2day + states_all$HARVESTlitter_foliage_gCm2day
                      site_output$HarvestFractionOfTurnover_foliage = states_all$HARVESTextracted_foliage_gCm2day + states_all$HARVESTlitter_foliage_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_foliage_years = apply(site_output$outflux_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_foliage_years = apply(states_all$foliage_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_foliage_years * 365.25)
                  site_output$MTT_annual_foliage_years = t(site_output$MTT_annual_foliage_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_foliage_years = apply(site_output$MTT_annual_foliage_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_foliage = quantile(apply(site_output$NaturalFractionOfTurnover_foliage,1,mean, na.rm = na_flag) /
                                                                           apply(site_output$outflux_foliage_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_foliage = quantile(apply(site_output$FireFractionOfTurnover_foliage,1,mean, na.rm = na_flag) /
                                                                        apply(site_output$outflux_foliage_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_foliage = quantile(apply(site_output$HarvestFractionOfTurnover_foliage,1,mean, na.rm = na_flag) /
                                                                           apply(site_output$outflux_foliage_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from foliage
                  site_output$mean_outflux_foliage_gCm2day = quantile(apply(site_output$outflux_foliage_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_foliage_gCm2day = apply(site_output$outflux_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Fine roots related pool, change, input and output variables
              if (exists(x = "roots_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$roots_gCm2 = apply(states_all$roots_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_roots_gCm2 = quantile(apply(states_all$roots_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$roots_gCm2 - states_all$roots_gCm2[,1] # difference in root from initial
                  site_output$dCroots_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Calculate the mean annual maximums
                  dCbio = apply(states_all$roots_gCm2, 1, rollapply_mean_annual_max, step = steps_per_year)
                  site_output$annual_max_roots_gCm2 = quantile(dCbio, prob=num_quantiles, na.rm = na_flag)
                  # Is rooting depth calculate (m) by this model?
                  if (exists(x = "RootDepth_m", where = states_all)) {
                      site_output$RootDepth_m = apply(states_all$RootDepth_m,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_RootDepth_m = quantile(apply(states_all$RootDepth_m,1,mean, na.rm = na_flag), prob=num_quantiles)
                      dCbio = states_all$RootDepth_m - states_all$RootDepth_m[,1] # difference in root from initial
                      site_output$dRootDepth_m = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  }
                  # Check for the possible pathways
                  if (exists(x = "alloc_roots_gCm2day", where = states_all)) {
                      site_output$alloc_roots_gCm2day = apply(states_all$alloc_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_alloc_roots_gCm2day = quantile(apply(states_all$alloc_roots_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  }
                  if (exists(x = "roots_to_litter_gCm2day", where = states_all)) {
                      site_output$roots_to_litter_gCm2day = apply(states_all$roots_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_roots_to_litter_gCm2day = quantile(apply(states_all$roots_to_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$NaturalFractionOfTurnover_roots = states_all$roots_to_litter_gCm2day
                      # Begin accumulation of output fluxes
                      site_output$outflux_roots_gCm2day = states_all$roots_to_litter_gCm2day
                  }
                  # If this one exists then maybe some other fluxes do
                  if (exists(x = "FIREemiss_roots_gCm2day", where = states_all)) {
                      site_output$FIREemiss_roots_gCm2day = apply(states_all$FIREemiss_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_roots_gCm2day = quantile(apply(states_all$FIREemiss_roots_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_roots = states_all$FIREemiss_roots_gCm2day
                      site_output$outflux_roots_gCm2day = site_output$outflux_roots_gCm2day + states_all$FIREemiss_roots_gCm2day
                  }
                  if (exists(x = "FIRElitter_roots_gCm2day", where = states_all)) {
                      site_output$FIRElitter_roots_gCm2day = apply(states_all$FIRElitter_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_roots_gCm2day = quantile(apply(states_all$FIRElitter_roots_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_roots_gCm2day = site_output$outflux_roots_gCm2day + states_all$FIRElitter_roots_gCm2day
                      site_output$FireFractionOfTurnover_roots = site_output$FireFractionOfTurnover_roots + states_all$FIRElitter_roots_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_roots_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_roots_gCm2day = apply(states_all$HARVESTextracted_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_roots_gCm2day = quantile(apply(states_all$HARVESTextracted_roots_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$HARVESTlitter_roots_gCm2day = apply(states_all$HARVESTlitter_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTlitter_roots_gCm2day = quantile(apply(states_all$HARVESTlitter_roots_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_roots_gCm2day = site_output$outflux_roots_gCm2day + states_all$HARVESTextracted_roots_gCm2day + states_all$HARVESTlitter_roots_gCm2day
                      site_output$HarvestFractionOfTurnover_roots = states_all$HARVESTextracted_roots_gCm2day + states_all$HARVESTlitter_roots_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_roots_years = apply(site_output$outflux_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_roots_years = apply(states_all$roots_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_roots_years * 365.25)
                  site_output$MTT_annual_roots_years = t(site_output$MTT_annual_roots_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_roots_years = apply(site_output$MTT_annual_roots_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_roots = quantile(apply(site_output$NaturalFractionOfTurnover_roots,1,mean, na.rm = na_flag) /
                                                                         apply(site_output$outflux_roots_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_roots = quantile(apply(site_output$FireFractionOfTurnover_roots,1,mean, na.rm = na_flag) /
                                                                      apply(site_output$outflux_roots_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_roots = quantile(apply(site_output$HarvestFractionOfTurnover_roots,1,mean, na.rm = na_flag) /
                                                                         apply(site_output$outflux_roots_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from roots
                  site_output$mean_outflux_roots_gCm2day = quantile(apply(site_output$outflux_roots_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_roots_gCm2day = apply(site_output$outflux_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Wood related pool, change, input and output variables
              if (exists(x = "wood_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$wood_gCm2 = apply(states_all$wood_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_wood_gCm2 = quantile(apply(states_all$wood_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$wood_gCm2 - states_all$wood_gCm2[,1] # difference in wood from initial
                  site_output$dCwood_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Calculate the mean annual maximums
                  dCbio = apply(states_all$wood_gCm2, 1, rollapply_mean_annual_max, step = steps_per_year)
                  site_output$annual_max_wood_gCm2 = quantile(dCbio, prob=num_quantiles, na.rm = na_flag)
                  # Check for the possible pathways
                  if (exists(x = "alloc_wood_gCm2day", where = states_all)) {
                      site_output$alloc_wood_gCm2day = apply(states_all$alloc_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_alloc_wood_gCm2day = quantile(apply(states_all$alloc_wood_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  }
                  if (exists(x = "wood_to_litter_gCm2day", where = states_all)) {
                      site_output$wood_to_litter_gCm2day = apply(states_all$wood_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_wood_to_litter_gCm2day = quantile(apply(states_all$wood_to_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$NaturalFractionOfTurnover_wood = states_all$wood_to_litter_gCm2day
                      # Begin accumulating output fluxes
                      site_output$outflux_wood_gCm2day = states_all$wood_to_litter_gCm2day
                  }
                  # If this one exists then maybe some other fluxes do
                  if (exists(x = "FIREemiss_wood_gCm2day", where = states_all)) {
                      site_output$FIREemiss_wood_gCm2day = apply(states_all$FIREemiss_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_wood_gCm2day = quantile(apply(states_all$FIREemiss_wood_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_wood = states_all$FIREemiss_wood_gCm2day
                      site_output$outflux_wood_gCm2day = site_output$outflux_wood_gCm2day + states_all$FIREemiss_wood_gCm2day
                  }
                  if (exists(x = "FIRElitter_wood_gCm2day", where = states_all)) {
                      site_output$FIRElitter_wood_gCm2day = apply(states_all$FIRElitter_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_wood_gCm2day = quantile(apply(states_all$FIRElitter_wood_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_wood_gCm2day = site_output$outflux_wood_gCm2day + states_all$FIRElitter_wood_gCm2day
                      site_output$FireFractionOfTurnover_wood = site_output$FireFractionOfTurnover_wood + states_all$FIRElitter_wood_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_wood_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_wood_gCm2day = apply(states_all$HARVESTextracted_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_wood_gCm2day = quantile(apply(states_all$HARVESTextracted_wood_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$HARVESTlitter_wood_gCm2day = apply(states_all$HARVESTlitter_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTlitter_wood_gCm2day = quantile(apply(states_all$HARVESTlitter_wood_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_wood_gCm2day = site_output$outflux_wood_gCm2day + states_all$HARVESTextracted_wood_gCm2day + states_all$HARVESTlitter_wood_gCm2day
                      site_output$HarvestFractionOfTurnover_wood = states_all$HARVESTextracted_wood_gCm2day + states_all$HARVESTlitter_wood_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_wood_years = apply(site_output$outflux_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_wood_years = apply(states_all$wood_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_wood_years * 365.25)
                  site_output$MTT_annual_wood_years = t(site_output$MTT_annual_wood_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_wood_years = apply(site_output$MTT_annual_wood_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_wood = quantile(apply(site_output$NaturalFractionOfTurnover_wood,1,mean, na.rm = na_flag) /
                                                                        apply(site_output$outflux_wood_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_wood = quantile(apply(site_output$FireFractionOfTurnover_wood,1,mean, na.rm = na_flag) /
                                                                     apply(site_output$outflux_wood_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_wood = quantile(apply(site_output$HarvestFractionOfTurnover_wood,1,mean, na.rm = na_flag) /
                                                                        apply(site_output$outflux_wood_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from wood
                  site_output$mean_outflux_wood_gCm2day = quantile(apply(site_output$outflux_wood_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_wood_gCm2day = apply(site_output$outflux_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Fine litter pool, change and output variables
              if (exists(x = "litter_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$litter_gCm2 = apply(states_all$litter_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_litter_gCm2 = quantile(apply(states_all$litter_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$litter_gCm2 - states_all$litter_gCm2[,1] # difference in litter from initial
                  site_output$dClitter_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Acumulate total output
                  site_output$outflux_litter_gCm2day = states_all$rhet_litter_gCm2day + states_all$litter_to_som_gCm2day
                  site_output$NaturalFractionOfTurnover_litter = site_output$outflux_litter_gCm2day
                  # Heterotrophic respiration
                  site_output$rhet_litter_gCm2day = apply(states_all$rhet_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_rhet_litter_gCm2day = quantile(apply(states_all$rhet_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Decomposition
                  site_output$litter_to_som_gCm2day = apply(states_all$litter_to_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_litter_to_som_gCm2day = quantile(apply(states_all$litter_to_som_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # If this one exists then maybe some other fluxes do
                  if (exists(x = "FIREemiss_litter_gCm2day", where = states_all)) {
                      site_output$FIREemiss_litter_gCm2day = apply(states_all$FIREemiss_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_litter_gCm2day = quantile(apply(states_all$FIREemiss_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_litter = states_all$FIREemiss_litter_gCm2day
                      site_output$outflux_litter_gCm2day = site_output$outflux_litter_gCm2day + states_all$FIREemiss_litter_gCm2day
                  }
                  if (exists(x = "FIRElitter_litter_gCm2day", where = states_all)) {
                      site_output$FIRElitter_litter_gCm2day = apply(states_all$FIRElitter_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_litter_gCm2day = quantile(apply(states_all$FIRElitter_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_litter_gCm2day = site_output$outflux_litter_gCm2day + states_all$FIRElitter_litter_gCm2day
                      site_output$FireFractionOfTurnover_litter = site_output$FireFractionOfTurnover_litter + states_all$FIRElitter_litter_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_litter_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_litter_gCm2day = apply(states_all$HARVESTextracted_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_litter_gCm2day = quantile(apply(states_all$HARVESTextracted_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_litter_gCm2day = site_output$outflux_litter_gCm2day + states_all$HARVESTextracted_litter_gCm2day
                      site_output$HarvestFractionOfTurnover_litter = states_all$HARVESTextracted_litter_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_litter_years = apply(site_output$outflux_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_litter_years = apply(states_all$litter_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_litter_years * 365.25)
                  site_output$MTT_annual_litter_years = t(site_output$MTT_annual_litter_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_litter_years = apply(site_output$MTT_annual_litter_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_litter = quantile(apply(site_output$NaturalFractionOfTurnover_litter,1,mean, na.rm = na_flag) /
                                                                          apply(site_output$outflux_litter_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_litter = quantile(apply(site_output$FireFractionOfTurnover_litter,1,mean, na.rm = na_flag) /
                                                                       apply(site_output$outflux_litter_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_litter = quantile(apply(site_output$HarvestFractionOfTurnover_litter,1,mean, na.rm = na_flag) /
                                                                          apply(site_output$outflux_litter_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from litter
                  site_output$mean_outflux_litter_gCm2day = quantile(apply(site_output$outflux_litter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_litter_gCm2day = apply(site_output$outflux_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Wood litter pool, change and output variables
              if (exists(x = "woodlitter_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$woodlitter_gCm2 = apply(states_all$woodlitter_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_woodlitter_gCm2 = quantile(apply(states_all$woodlitter_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$woodlitter_gCm2 - states_all$woodlitter_gCm2[,1] # difference in wood litter from initial
                  site_output$dCwoodlitter_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Begin accumulating losses from wood litter
                  site_output$outflux_woodlitter_gCm2day = states_all$rhet_woodlitter_gCm2day + states_all$woodlitter_to_som_gCm2day
                  site_output$NaturalFractionOfTurnover_woodlitter = site_output$outflux_woodlitter_gCm2day
                  # Heterotrophic respiration
                  site_output$rhet_woodlitter_gCm2day = apply(states_all$rhet_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_rhet_woodlitter_gCm2day = quantile(apply(states_all$rhet_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Decomposition
                  site_output$woodlitter_to_som_gCm2day = apply(states_all$woodlitter_to_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_woodlitter_to_som_gCm2day = quantile(apply(states_all$woodlitter_to_som_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # If this one exists then maybe some other fluxes do
                  if (exists(x = "FIREemiss_woodlitter_gCm2day", where = states_all)) {
                      site_output$FIREemiss_woodlitter_gCm2day = apply(states_all$FIREemiss_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_woodlitter_gCm2day = quantile(apply(states_all$FIREemiss_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_woodlitter = states_all$FIREemiss_woodlitter_gCm2day
                      site_output$outflux_woodlitter_gCm2day = site_output$outflux_woodlitter_gCm2day + states_all$FIREemiss_woodlitter_gCm2day
                  }
                  if (exists(x = "FIRElitter_woodlitter_gCm2day", where = states_all)) {
                      site_output$FIRElitter_woodlitter_gCm2day = apply(states_all$FIRElitter_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIRElitter_woodlitter_gCm2day = quantile(apply(states_all$FIRElitter_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_woodlitter_gCm2day = site_output$outflux_woodlitter_gCm2day + states_all$FIRElitter_woodlitter_gCm2day
                      site_output$FireFractionOfTurnover_woodlitter = site_output$FireFractionOfTurnover_woodlitter + states_all$FIRElitter_woodlitter_gCm2day
                  }
                  if (exists(x = "HARVESTextracted_woodlitter_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_woodlitter_gCm2day = apply(states_all$HARVESTextracted_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_woodlitter_gCm2day = quantile(apply(states_all$HARVESTextracted_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_woodlitter_gCm2day = site_output$outflux_woodlitter_gCm2day + states_all$HARVESTextracted_woodlitter_gCm2day
                      site_output$HarvestFractionOfTurnover_woodlitter = states_all$HARVESTextracted_woodlitter_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_woodlitter_years = apply(site_output$outflux_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_woodlitter_years = apply(states_all$woodlitter_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_woodlitter_years * 365.25)
                  site_output$MTT_annual_woodlitter_years = t(site_output$MTT_annual_woodlitter_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_woodlitter_years = apply(site_output$MTT_annual_woodlitter_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_woodlitter = quantile(apply(site_output$NaturalFractionOfTurnover_woodlitter,1,mean, na.rm = na_flag) /
                                                                              apply(site_output$outflux_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_woodlitter = quantile(apply(site_output$FireFractionOfTurnover_woodlitter,1,mean, na.rm = na_flag) /
                                                                           apply(site_output$outflux_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_woodlitter = quantile(apply(site_output$HarvestFractionOfTurnover_woodlitter,1,mean, na.rm = na_flag) /
                                                                              apply(site_output$outflux_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from wood
                  site_output$mean_outflux_woodlitter_gCm2day = quantile(apply(site_output$outflux_woodlitter_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_woodlitter_gCm2day = apply(site_output$outflux_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Soil organic matter pool, change and output variables
              if (exists(x = "som_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$som_gCm2 = apply(states_all$som_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_som_gCm2 = quantile(apply(states_all$som_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$som_gCm2 - states_all$som_gCm2[,1] # difference in som from initial
                  site_output$dCsom_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Begin accumulating losses from som
                  site_output$outflux_som_gCm2day = states_all$rhet_som_gCm2day
                  site_output$NaturalFractionOfTurnover_som = site_output$outflux_som_gCm2day
                  # Heterotrphic respiration
                  site_output$rhet_som_gCm2day = apply(states_all$rhet_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_rhet_som_gCm2day = quantile(apply(states_all$rhet_som_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # If this one exists then maybe some other fluxes do
                  if (exists(x = "FIREemiss_som_gCm2day", where = states_all)) {
                      site_output$FIREemiss_som_gCm2day = apply(states_all$FIREemiss_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_som_gCm2day = quantile(apply(states_all$FIREemiss_som_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_som = states_all$FIREemiss_som_gCm2day
                      site_output$outflux_som_gCm2day = site_output$outflux_som_gCm2day + states_all$FIREemiss_som_gCm2day
                  }
                  # NOTE: FIRElitter does not exist as there is not litter which leaves the som pool
                  if (exists(x = "HARVESTextracted_som_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_som_gCm2day = apply(states_all$HARVESTextracted_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_som_gCm2day = quantile(apply(states_all$HARVESTextracted_som_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_som_gCm2day = site_output$outflux_som_gCm2day + states_all$HARVESTextracted_som_gCm2day
                      site_output$HarvestFractionOfTurnover_som = states_all$HARVESTextracted_som_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_som_years = apply(site_output$outflux_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_som_years = apply(states_all$som_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_som_years * 365.25)
                  site_output$MTT_annual_som_years = t(site_output$MTT_annual_som_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_som_years = apply(site_output$MTT_annual_som_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_som = quantile(apply(site_output$NaturalFractionOfTurnover_som,1,mean, na.rm = na_flag) /
                                                                       apply(site_output$outflux_som_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_som = quantile(apply(site_output$FireFractionOfTurnover_som,1,mean, na.rm = na_flag) /
                                                                    apply(site_output$outflux_som_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_som = quantile(apply(site_output$HarvestFractionOfTurnover_som,1,mean, na.rm = na_flag) /
                                                                       apply(site_output$outflux_som_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from som
                  site_output$mean_outflux_som_gCm2day = quantile(apply(site_output$outflux_som_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_som_gCm2day = apply(site_output$outflux_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              # Dead organic matter pool, change and output variables
              # NOTE: this is potentially the combination of som, litter and wood litter
              if (exists(x = "dom_gCm2", where = states_all)) {
                  # Assign pool to site_output
                  site_output$dom_gCm2 = apply(states_all$dom_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_dom_gCm2 = quantile(apply(states_all$dom_gCm2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Determine net pool change over time
                  dCbio = states_all$dom_gCm2 - states_all$dom_gCm2[,1] # difference in dom from initial
                  site_output$dCdom_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  # Accumulating output fluxes from dead organic matter
                  site_output$outflux_dom_gCm2day = states_all$rhet_dom_gCm2day
                  site_output$NaturalFractionOfTurnover_dom = site_output$outflux_dom_gCm2day
                  # Heterotrophic respiration
                  site_output$rhet_dom_gCm2day = apply(states_all$rhet_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_rhet_dom_gCm2day = quantile(apply(states_all$rhet_dom_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # If this one exists then maybe some other fluxes do
                  if (exists(x = "FIREemiss_dom_gCm2day", where = states_all)) {
                      site_output$FIREemiss_dom_gCm2day = apply(states_all$FIREemiss_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_FIREemiss_dom_gCm2day = quantile(apply(states_all$FIREemiss_dom_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$FireFractionOfTurnover_dom = states_all$FIREemiss_dom_gCm2day
                      site_output$outflux_dom_gCm2day = site_output$outflux_dom_gCm2day + states_all$FIREemiss_dom_gCm2day
                  }
                  # NOTE: FIRElitter does not exist as there is not litter which leaves the dom pool
                  if (exists(x = "HARVESTextracted_dom_gCm2day", where = states_all)) {
                      site_output$HARVESTextracted_dom_gCm2day = apply(states_all$HARVESTextracted_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
                      site_output$mean_HARVESTextracted_dom_gCm2day = quantile(apply(states_all$HARVESTextracted_dom_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                      site_output$outflux_dom_gCm2day = site_output$outflux_dom_gCm2day + states_all$HARVESTextracted_som_gCm2day
                      site_output$HarvestFractionOfTurnover_dom = states_all$HARVESTextracted_dom_gCm2day
                  }
                  # Use this information to determine the mean residence times as it evolves over time
                  site_output$MTT_annual_dom_years = apply(site_output$outflux_dom_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                  site_output$MTT_annual_dom_years = apply(states_all$dom_gCm2,1, rollapply_mean_annual, step = steps_per_year) / (site_output$MTT_annual_dom_years * 365.25)
                  site_output$MTT_annual_dom_years = t(site_output$MTT_annual_dom_years) # rollapply inverts the dimensions from that wanted
                  site_output$MTT_annual_dom_years = apply(site_output$MTT_annual_dom_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
                  # Convert to fractions
                  site_output$NaturalFractionOfTurnover_dom = quantile(apply(site_output$NaturalFractionOfTurnover_dom,1,mean, na.rm = na_flag) /
                                                                       apply(site_output$outflux_dom_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$FireFractionOfTurnover_dom = quantile(apply(site_output$FireFractionOfTurnover_dom,1,mean, na.rm = na_flag) /
                                                                    apply(site_output$outflux_dom_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  site_output$HarvestFractionOfTurnover_dom = quantile(apply(site_output$HarvestFractionOfTurnover_dom,1,mean, na.rm = na_flag) /
                                                                       apply(site_output$outflux_dom_gCm2day,1,mean, na.rm = na_flag), prob = num_quantiles, na.rm = na_flag)
                  # Mean outflux from dom
                  site_output$mean_outflux_dom_gCm2day = quantile(apply(site_output$outflux_dom_gCm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Now aggregate across quantiles
                  site_output$outflux_dom_gCm2day = apply(site_output$outflux_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
              }

              ###
              # Aggregate Water cycle related fluxes
              ###

              # Water cycle specific if available
              if (exists(x = "ET_kgH2Om2day", where = states_all)) {
                  ## current water in the soil surface layer (0-30 cm)
                  site_output$SurfWater_kgH2Om2 = apply(states_all$SurfWater_kgH2Om2,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  site_output$mean_SurfWater_kgH2Om2 = quantile(apply(states_all$SurfWater_kgH2Om2,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$SurfWater_kgH2Om2 - states_all$SurfWater_kgH2Om2[,1] # difference in surface water from initial
                  site_output$dSurfWater_kgH2Om2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

                  # plant apparent soil water potential (MPa)
                  site_output$wSWP_MPa = apply(states_all$wSWP_MPa,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  site_output$mean_wSWP_MPa = quantile(apply(states_all$wSWP_MPa,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$wSWP_MPa - states_all$wSWP_MPa[,1] # difference from initial
                  site_output$dwSWP_MPa = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

                  # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                  site_output$ET_kgH2Om2day = apply(states_all$ET_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  site_output$mean_ET_kgH2Om2day = quantile(apply(states_all$ET_kgH2Om2day,1,mean, na.rm = na_flag), prob=num_quantiles)
              }

              # Snow related
              if (exists(x = "snow_kgH2Om2", where = states_all)) {
                  ## Snow on soil surface
                  site_output$snow_kgH2Om2 = apply(states_all$snow_kgH2Om2,2,quantile,prob=num_quantiles,na.rm = na_flag)
                  site_output$mean_snow_kgH2Om2 = quantile(apply(states_all$snow_kgH2Om2,1,mean, na.rm = na_flag), prob=num_quantiles)
              }

              ###
              # Aggregate ACM diagnositic information
              ###

              # Canopy process information if available
              if (exists(x = "APAR_MJm2day", where = states_all)) {
                  # Extract the absorbed photosynthetically active radiation by the canopy
                  site_output$APAR_MJm2day = apply(states_all$APAR_MJm2day,2,quantile, prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_APAR_MJm2day = quantile(apply(states_all$APAR_MJm2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$APAR_MJm2day - states_all$APAR_MJm2day[,1] # difference in dom from initial
                  site_output$dAPAR_MJm2day = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
              }
              if (exists(x = "CiCa", where = states_all)) {
                  # Extract the internal vs ambient CO2 ratio
                  site_output$CiCa = apply(states_all$CiCa,2,quantile, prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_CiCa = quantile(apply(states_all$CiCa,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$CiCa - states_all$CiCa[,1] # difference in dom from initial
                  site_output$dCiCa = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
              }
              if (exists(x = "gs_demand_supply_ratio", where = states_all)) {
                  # Extract the ratio of stomatal conductance relative to its maximum value,
                  # this metric provides information on the demand vs supply constrains on stomatal conductance
                  site_output$gs_demand_supply_ratio = apply(states_all$gs_demand_supply_ratio,2,quantile, prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_gs_demand_supply_ratio = quantile(apply(states_all$gs_demand_supply_ratio,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$gs_demand_supply_ratio - states_all$gs_demand_supply_ratio[,1] # difference in dom from initial
                  site_output$dgs_demand_supply_ratio = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
              }
              if (exists(x = "gs_mmolH2Om2day", where = states_all)) {
                  # Extract the canopy stomatal conductance
                  site_output$gs_mmolH2Om2day = apply(states_all$gs_mmolH2Om2day,2,quantile, prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_gs_mmolH2Om2day = quantile(apply(states_all$gs_mmolH2Om2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$gs_mmolH2Om2day - states_all$gs_mmolH2Om2day[,1] # difference in dom from initial
                  site_output$dgs_mmolH2Om2day = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
              }
              if (exists(x = "gb_mmolH2Om2day", where = states_all)) {
                  # Extract the canopy boundary layer conductance
                  site_output$gb_mmolH2Om2day = apply(states_all$gb_mmolH2Om2day,2,quantile, prob=num_quantiles, na.rm = na_flag)
                  site_output$mean_gb_mmolH2Om2day = quantile(apply(states_all$gb_mmolH2Om2day,1,mean, na.rm = na_flag), prob=num_quantiles)
                  # Calculate change over time
                  dCbio = states_all$gb_mmolH2Om2day - states_all$gb_mmolH2Om2day[,1] # difference in dom from initial
                  site_output$dgb_mmolH2Om2day = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
              }

              ###
              # Aggregate model ensemble - observation uncertainty consistency
              # information
              ###

              # Any time series assimilated data overlaps?
              if (exists(x = "gpp_assim_data_overlap_fraction", where = states_all)) {
                  site_output$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction
              }
              if (exists(x = "lai_assim_data_overlap_fraction", where = states_all)) {
                  site_output$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction
              }
              if (exists(x = "nee_assim_data_overlap_fraction", where = states_all)) {
                  site_output$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction
              }
              if (exists(x = "wood_assim_data_overlap_fraction", where = states_all)) {
                  site_output$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction
              }
              if (exists(x = "soil_assim_data_overlap_fraction", where = states_all)) {
                  site_output$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction
              }
              if (exists(x = "et_assim_data_overlap_fraction", where = states_all)) {
                  site_output$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction
              }
              if (exists(x = "nbe_assim_data_overlap_fraction", where = states_all)) {
                  site_output$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction
              }
              if (exists(x = "fire_assim_data_overlap_fraction", where = states_all)) {
                  site_output$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction
              }

              # C-cycle flux correlation with parameters
              site_output$nee_par_cor = states_all$nee_par_cor
              site_output$gpp_par_cor = states_all$gpp_par_cor
              site_output$rauto_par_cor = states_all$rauto_par_cor
              site_output$rhet_par_cor = states_all$rhet_par_cor
              site_output$fire_par_cor = states_all$fire_par_cor

              # save to pixel specific file for the moment... in "run_mcmc_results" these will be combined into a single grid
              save(site_output,file=outfile_stock_fluxes, compress = "gzip", compression_level = 9)
          }

          dummy = 0

      } else { # error_check == FALSE

          dummy = -1

      } # error_check == FALSE

  } else { # *parameters.RData already exists

      dummy = -1
      print('Already extracted result vectors (set repair = 1 if re-run is needed)')

  } # *parameters.RData already exists

  # Return
  return(dummy)

} # end of run_each_site

## Use byte compile
run_each_site<-cmpfun(run_each_site)

###
## Function to control the site level running of the re-processing by DALEC
###

run_mcmc_results <- function (PROJECT,stage,repair,grid_override) {

  print('Welcome to RUN_MCMC_RESULTS!!')

  # bundle needed functions down the chain
  functions_list=c("read_parameter_chains","read_binary_file_format","simulate_all",
                   "read_binary_response_surface","crop_development_parameters",
                   "have_chains_converged","psrf","read_parameter_covariance",
                   "ensemble_within_range","rollapply_mean_annual_max",
                   "rollapply_mean_annual")
  # start marker
  stime = proc.time()["elapsed"]

  # how many plots in total do we have
  nos_plots = 1:PROJECT$nosites

  # now check which ones we need to calculate, but only if override not in play
  if (repair != 1) {
      print("...beginning filterings for sites we have already processed")
      keep_list = 0
      for (i in seq(1, length(nos_plots))) {
           outfile_parameters = paste(PROJECT$results_processedpath,PROJECT$sites[i],"_parameters.RData",sep="")
           if (file.exists(outfile_parameters) == FALSE) {keep_list=append(keep_list,i)}
      }
      # filter out the sites we already have then
      keep_list = keep_list[-1] ; print(paste("......removing ",length(nos_plots)-length(keep_list)," sites out of ",length(nos_plots)," from the analysis",sep=""))
      nos_plots = nos_plots[keep_list]
  }

  # now request the creation of the plots
  if (use_parallel & length(nos_plots) > 1) {
      print("...beginning parallel operations")
      cl <- makeCluster(min(length(nos_plots),numWorkers), type = "PSOCK")
      clusterExport(cl,functions_list)
      # load R libraries in cluster
      clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
      dummy = parLapply(cl,nos_plots,fun=run_each_site,PROJECT=PROJECT,stage=stage,
                        repair=repair,grid_override=grid_override)
#print(dummy)
      stopCluster(cl)
      print("...finished parallel operations")
  } else {
      print("...beginning serial operations")
      # or use serial
      dummy = lapply(nos_plots,FUN=run_each_site,PROJECT=PROJECT,stage=stage,
                   repair=repair,grid_override=grid_override)
      print("...finished serial operations")
  } # parallel option

  # now if this is a gridded run we want to take out individual site specific summary files and combine them into a single file
  if (PROJECT$spatial_type == "grid" & grid_override == FALSE) {

      # update the user
      print("...combining pixel level stock / flux output into single array")

      # determine what the output file name is here, so that we can check if one already exists
      outfile_grid = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
      # read in example input files for using in some calculation
      driver_file = paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[1],".bin",sep="")
      if (exists(driver_file)) {
          drivers = read_binary_file_format(driver_file)
      } else {
          read_driver = FALSE ; nn = 2
          while (read_driver == FALSE) {
             driver_file = paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[nn],".bin",sep="")
             if (file.exists(driver_file) == TRUE) { read_driver = TRUE } else { nn = nn + 1}
          }
          drivers = read_binary_file_format(driver_file)
      }
      # Determine some useful information for the analysis below
      nos_years = (as.numeric(PROJECT$end_year) - as.numeric(PROJECT$start_year))+1
      steps_per_year = floor(dim(drivers$met)[1] / nos_years)

      # Begin creation of all variables for the output gridded dataset
      # otherwise load the existing but incomplete version from file
      if (file.exists(outfile_grid) == FALSE | repair == 1) {

          # make a list of all the files we will be reading in
          to_do = list.files(PROJECT$results_processedpath, full.names=TRUE)
          to_do = to_do[grepl("_stock_fluxes",to_do)]
          # read in the first file so that we can then set up all the grids for combined summary output
          load(to_do[1])

          #
          # Generate the summary information
          # Time invariant but contain uncertainty, shaped into the full spatial grid
          # i.e. including areas not part of the analysis but within the spatial domain
          #

          # Define grid_output object and add some basic information
          grid_output = list(readme = c("For gridded variables (e.g. mean_*, or final_*) the default dimesions are long/lat/quantile",
                                        "For time varying variables (e.g. wood_gCm2) default dimensions are site/quantile/time"),
                             num_quantiles = site_output$num_quantiles,
                             long_dim = PROJECT$long_dim, lat_dim = PROJECT$lat_dim, time_dim = dim(site_output$labile_gCm2)[2],
                             start_year = PROJECT$start_year, end_year = PROJECT$end_year,
                             steps_per_year = steps_per_year, nos_years = nos_years)

          ###
          # Define likelihood / parameter / driver information
          ###

          # loop through parameters + likelihood
          grid_output$parameters = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)+1,dim(site_output$labile_gCm2)[1]))
          # track which parameters have converged + likelihood
          grid_output$parameters_converged = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)+1))
          # Mean meteorological conditions
          grid_output$mean_temperature_C = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$mean_radiation_MJm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$mean_vpd_Pa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$mean_precipitation_kgH2Om2yr = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Assimilated leaf area index information
          grid_output$assimilated_lai_max_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_lai_mean_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_lai_sd_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_lai_unc_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Assimilated wood stock / prior information
          grid_output$assimilated_wood_mean_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_wood_mean_unc_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Assimilated som stock / prior information
          grid_output$assimilated_som_mean_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_som_mean_unc_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # Net primary production allocation fractions
          if (exists(x = "NPP_foliage_fraction", where = site_output)) {grid_output$NPP_foliage_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "NPP_roots_fraction", where = site_output)) {grid_output$NPP_roots_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "NPP_wood_fraction", where = site_output)) {grid_output$NPP_wood_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          # Analysis mean transit (residence) times (years)
          if (exists(x = "MTT_labile_years", where = site_output)) {grid_output$MTT_labile_years =array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "MTT_foliage_years", where = site_output)) {grid_output$MTT_foliage_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "MTT_roots_years", where = site_output)) {grid_output$MTT_roots_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "MTT_wood_years", where = site_output)) {grid_output$MTT_wood_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "MTT_litter_years", where = site_output)) {grid_output$MTT_litter_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "MTT_woodlitter_years", where = site_output)) {grid_output$MTT_woodlitter_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "MTT_som_years", where = site_output)) {grid_output$MTT_som_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          # Steady state C stock estimates (gC/m2)
          if (exists(x = "SS_labile_gCm2", where = site_output)) {grid_output$SS_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "SS_foliage_gCm2", where = site_output)) {grid_output$SS_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "SS_roots_gCm2", where = site_output)) {grid_output$SS_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "SS_wood_gCm2", where = site_output)) {grid_output$SS_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "SS_litter_gCm2", where = site_output)) {grid_output$SS_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "SS_woodlitter_gCm2", where = site_output)) {grid_output$SS_woodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (exists(x = "SS_som_gCm2", where = site_output)) {grid_output$SS_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}

          # Total ecosystem C is always present
          grid_output$mean_Ctotal_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_Ctotal_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCtotal_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Ecosystem leaf area index is always present
          grid_output$mean_lai_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dlai_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Mean ecosystem bulk fluxes are always present
          grid_output$mean_nee_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_gpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rauto_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rhet_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_reco_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_npp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_harvest_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_fire_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_nbe_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_nbp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Time varying pixel based (i.e. not within the grid) values with quantile based uncertainty for...
          # States
          grid_output$Ctotal_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$lai_m2m2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          # Fluxes
          grid_output$nee_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$gpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rauto_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rhet_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$reco_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$npp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$harvest_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$fire_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$nbe_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$nbp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))

          # Based on the presence of each pool define the grids for the mean and final values.
          # Also, create the time varying but quantile based values and time

          # Gridded biomass information
          if (exists(x = "biomass_gCm2", where = site_output)) {
              # Grid mean / finals for globally available variables
              grid_output$mean_biomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_biomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCbiomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_biomass_to_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$biomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCbiomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_biomass_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (exists(x = "FIREemiss_biomass_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_biomass_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_biomass_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$HARVESTlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded labile information
          if (exists(x = "labile_gCm2", where = site_output)) {
              # Grid averages
              grid_output$mean_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dClabile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_labile_to_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_alloc_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$labile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dClabile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$labile_to_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$alloc_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_labile_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (exists(x = "FIREemiss_labile_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_labile_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_labile_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$HARVESTlitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded foliage information
          if (exists(x = "foliage_gCm2", where = site_output)) {
              # Grided mean / final
              grid_output$mean_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCfoliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$foliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCfoliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_foliage_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (exists(x = "alloc_foliage_gCm2day", where = site_output)) {
                  grid_output$alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIREemiss_foliage_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_foliage_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_foliage_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$HARVESTlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded roots information
          if (exists(x = "roots_gCm2", where = site_output)) {
              # Gridded mean / final
              grid_output$mean_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCroots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_alloc_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$annual_max_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$roots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCroots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$alloc_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_roots_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Is rooting depth calculate (m) by this model?
              if (exists(x = "RootDepth_m", where = site_output)) {
                  grid_output$RootDepth_m = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$dRootDepth_m = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_RootDepth_m = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$final_dRootDepth_m = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIREemiss_roots_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_roots_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_roots_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$HARVESTlitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded wood information
          if (exists(x = "wood_gCm2", where = site_output)) {
              # Gridded mean / final variables
              grid_output$mean_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_alloc_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$annual_max_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$wood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCwood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$alloc_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_wood_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (exists(x = "FIREemiss_wood_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_wood_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_wood_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$HARVESTlitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded litter information
          if (exists(x = "litter_gCm2", where = site_output)) {
              # Gridded mean / final
              grid_output$mean_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dClitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_litter_to_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$litter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dClitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$litter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_litter_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (exists(x = "FIREemiss_litter_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_litter_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_litter_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded wood litter information
          if (exists(x = "woodlitter_gCm2", where = site_output)) {
              # Gridded mean / final
              grid_output$mean_woodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_woodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCwoodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$woodlitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCwoodlitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_woodlitter_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (exists(x = "FIREemiss_woodlitter_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "FIRElitter_woodlitter_gCm2day", where = site_output)) {
                  grid_output$FIRElitter_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_woodlitter_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded som information
          if (exists(x = "som_gCm2", where = site_output)) {
              # Gridded mean / finals
              grid_output$mean_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCsom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$som_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCsom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_som_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying
              if (exists(x = "FIREemiss_som_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_som_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Gridded dom information
          if (exists(x = "dom_gCm2", where = site_output)) {
              # Gridded mean / finals
              grid_output$mean_dom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCdom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$dom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCdom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$MTT_annual_dom_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of tunover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional
              if (exists(x = "FIREemiss_dom_gCm2day", where = site_output)) {
                  grid_output$FIREemiss_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
              if (exists(x = "HARVESTextracted_dom_gCm2day", where = site_output)) {
                  grid_output$HARVESTextracted_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }

          # Water cycle specific variables
          if (exists(x = "ET_kgH2Om2day", where = site_output)) {
              # currently water in the soil surface layer (0-30 cm)
              grid_output$mean_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dSurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dSurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # plant apparent soil water potential (MPa)
              grid_output$mean_wSWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_wSWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$wSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dwSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # evapotranspiration (Etrans + Esoil + Ewetcanopy)
              grid_output$mean_ET_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$ET_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          # Snow specific
          if (exists(x = "snow_kgH2Om2", where = site_output)) {
              ## snow on soil surface
              grid_output$snow_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$mean_snow_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          }
          # Canopy process variables
          if (exists(x = "APAR_MJm2day", where = site_output)) {
              # Absorbed photosynthetically active radation
              grid_output$mean_APAR_MJm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$APAR_MJm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (exists(x = "CiCa", where = site_output)) {
              # Canopy Ci:Ca
              grid_output$mean_CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$CiCa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (exists(x = "gs_demand_supply_ratio", where = site_output)) {
              # Ratio of stomatal conductance relative to its maximum value,
              # this metric provides information on the demand vs supply constrains on stomatal conductance
              grid_output$mean_gs_demand_supply_ratio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$gs_demand_supply_ratio = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (exists(x = "gs_mmolH2Om2day", where = site_output)) {
              # Canopy stomatal conductance
              grid_output$mean_gs_mmolH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$gs_mmolH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (exists(x = "gb_mmolH2Om2day", where = site_output)) {
              # Canopy boundary layer conductance
              grid_output$mean_gb_mmolH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$gb_mmolH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }

          # Create overlap statistics variables - may not always get filled in the end
#          if (exists(x = "gpp_assim_data_overlap_fraction", where = site_output)) {
              grid_output$gpp_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "lai_assim_data_overlap_fraction", where = site_output)) {
              grid_output$lai_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "nee_assim_data_overlap_fraction", where = site_output)) {
              grid_output$nee_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "wood_assim_data_overlap_fraction", where = site_output)) {
              grid_output$wood_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "soil_assim_data_overlap_fraction", where = site_output)) {
              grid_output$soil_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "et_assim_data_overlap_fraction", where = site_output)) {
              grid_output$et_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "nbe_assim_data_overlap_fraction", where = site_output)) {
              grid_output$nbe_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }
#          if (exists(x = "fire_assim_data_overlap_fraction", where = site_output)) {
              grid_output$fire_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
#          }

          # Time and uncertainty invarient information,
          # this is the correlation between ensemble members for parameter and C-cycle flux variables
          grid_output$nee_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$gpp_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$rauto_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$rhet_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$fire_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))

          #
          # Generate the spatial information needed to relate i,j within the grid to long / lat and the summary vs detailed output
          #

          # vector into which the i,j position information within the grid will be stored
          grid_output$i_location = rep(NA, length.out = PROJECT$nosites)
          grid_output$j_location = rep(NA, length.out = PROJECT$nosites)

          # generate the lat / long grid again
          output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
          grid_output$lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # Determine grid area (m2)
          grid_output$area = calc_pixel_area(grid_output$lat,grid_output$long,PROJECT$resolution)
          # this output is in vector form and we need matching array shapes so...
          grid_output$area = array(grid_output$area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # calculate land mask
          grid_output$landmask=array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

      } else {

          # if file already exists load it up
          load(outfile_grid)

      } # have we already an output file

      # Loop through all sites now in turn and load their specific information
      # into the overall gridded variables
      for (n in seq(1,PROJECT$nosites)) {
           # load the current file
           outfile_stock_fluxes = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")
           if (file.exists(outfile_stock_fluxes)) {
               # load the site
               load(outfile_stock_fluxes)
               # determine the lat / long location within the grid
               slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
               slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
               if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
               # save for later
               grid_output$i_location[n] = slot_i ; grid_output$j_location[n] = slot_j

               # loop through parameters + likelihood
               grid_output$parameters[slot_i,slot_j,,] = site_output$parameters
               # track which parameters have converged + likelihood
               grid_output$parameters_converged[slot_i,slot_j,] = site_output$parameters_converged
               # Mean meteorological conditions
               grid_output$mean_temperature_C[slot_i,slot_j] = site_output$mean_temperature_C
               grid_output$mean_radiation_MJm2day[slot_i,slot_j] = site_output$mean_radiation_MJm2day
               grid_output$mean_vpd_Pa[slot_i,slot_j] = site_output$mean_vpd_Pa
               grid_output$mean_precipitation_kgH2Om2yr[slot_i,slot_j] = site_output$mean_precipitation_kgH2Om2yr
               # Assimilated leaf area index information
               grid_output$assimilated_lai_max_m2m2[slot_i,slot_j] = site_output$assimilated_lai_max_m2m2
               grid_output$assimilated_lai_mean_m2m2[slot_i,slot_j] = site_output$assimilated_lai_mean_m2m2
               grid_output$assimilated_lai_sd_m2m2[slot_i,slot_j] = site_output$assimilated_lai_sd_m2m2
               grid_output$assimilated_lai_unc_m2m2[slot_i,slot_j] = site_output$assimilated_lai_unc_m2m2
               # Assimilated wood stock / prior information
               grid_output$assimilated_wood_mean_gCm2[slot_i,slot_j] = site_output$assimilated_wood_mean_gCm2
               grid_output$assimilated_wood_mean_unc_gCm2[slot_i,slot_j] = site_output$assimilated_wood_mean_unc_gCm2
               # Assimilated som stock / prior information
               grid_output$assimilated_som_mean_gCm2[slot_i,slot_j] = site_output$assimilated_som_mean_gCm2
               grid_output$assimilated_som_mean_unc_gCm2[slot_i,slot_j] = site_output$assimilated_som_mean_unc_gCm2

               # Net primary production allocation fractions
               if (exists(x = "NPP_foliage_fraction", where = site_output)) {grid_output$NPP_foliage_fraction[slot_i,slot_j,] = site_output$NPP_foliage_fraction}
               if (exists(x = "NPP_roots_fraction", where = site_output)) {grid_output$NPP_roots_fraction[slot_i,slot_j,] = site_output$NPP_roots_fraction}
               if (exists(x = "NPP_wood_fraction", where = site_output)) {grid_output$NPP_wood_fraction[slot_i,slot_j,] = site_output$NPP_wood_fraction}
               # Analysis mean transit (residence) times (years)
               if (exists(x = "MTT_labile_years", where = site_output)) {grid_output$MTT_labile_years[slot_i,slot_j,] = site_output$MTT_labile_years}
               if (exists(x = "MTT_foliage_years", where = site_output)) {grid_output$MTT_foliage_years[slot_i,slot_j,] = site_output$MTT_foliage_years}
               if (exists(x = "MTT_roots_years", where = site_output)) {grid_output$MTT_roots_years[slot_i,slot_j,] = site_output$MTT_roots_years}
               if (exists(x = "MTT_wood_years", where = site_output)) {grid_output$MTT_wood_years[slot_i,slot_j,] = site_output$MTT_wood_years}
               if (exists(x = "MTT_litter_years", where = site_output)) {grid_output$MTT_litter_years[slot_i,slot_j,] = site_output$MTT_litter_years}
               if (exists(x = "MTT_woodlitter_years", where = site_output)) {grid_output$MTT_woodlitter_years[slot_i,slot_j,] = site_output$MTT_woodlitter_years}
               if (exists(x = "MTT_som_years", where = site_output)) {grid_output$MTT_som_years[slot_i,slot_j,] = site_output$MTT_som_years}
               # Steady state C stock estimates (gC/m2)
               if (exists(x = "SS_labile_gCm2", where = site_output)) {grid_output$SS_labile_gCm2[slot_i,slot_j,] = site_output$SS_labile_gCm2}
               if (exists(x = "SS_foliage_gCm2", where = site_output)) {grid_output$SS_foliage_gCm2[slot_i,slot_j,] = site_output$SS_foliage_gCm2}
               if (exists(x = "SS_roots_gCm2", where = site_output)) {grid_output$SS_roots_gCm2[slot_i,slot_j,] = site_output$SS_roots_gCm2}
               if (exists(x = "SS_wood_gCm2", where = site_output)) {grid_output$SS_wood_gCm2[slot_i,slot_j,] = site_output$SS_wood_gCm2}
               if (exists(x = "SS_litter_gCm2", where = site_output)) {grid_output$SS_litter_gCm2[slot_i,slot_j,] = site_output$SS_litter_gCm2}
               if (exists(x = "SS_woodlitter_gCm2", where = site_output)) {grid_output$SS_woodlitter_gCm2[slot_i,slot_j,] = site_output$SS_woodlitter_gCm2}
               if (exists(x = "SS_som_gCm2", where = site_output)) {grid_output$SS_som_gCm2[slot_i,slot_j,] = site_output$SS_som_gCm2}

               # Total ecosystem C is always present
               grid_output$mean_Ctotal_gCm2[slot_i,slot_j,] = site_output$mean_Ctotal_gCm2
               grid_output$final_Ctotal_gCm2[slot_i,slot_j,] = site_output$Ctotal_gCm2[,grid_output$time_dim]
               grid_output$final_dCtotal_gCm2[slot_i,slot_j,] = site_output$dCtotal_gCm2[,grid_output$time_dim]
               # Ecosystem leaf area index is always present
               grid_output$mean_lai_m2m2[slot_i,slot_j,] = site_output$mean_lai_m2m2
               grid_output$final_dlai_m2m2[slot_i,slot_j,] = site_output$dlai_m2m2[,grid_output$time_dim]
               # Mean ecosystem bulk fluxes are always present
               grid_output$mean_nee_gCm2day[slot_i,slot_j,] = site_output$mean_nee_gCm2day
               grid_output$mean_gpp_gCm2day[slot_i,slot_j,] = site_output$mean_gpp_gCm2day
               grid_output$mean_rauto_gCm2day[slot_i,slot_j,] = site_output$mean_rauto_gCm2day
               grid_output$mean_rhet_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_gCm2day
               grid_output$mean_reco_gCm2day[slot_i,slot_j,] = site_output$mean_reco_gCm2day
               grid_output$mean_npp_gCm2day[slot_i,slot_j,] = site_output$mean_npp_gCm2day
               grid_output$mean_harvest_gCm2day[slot_i,slot_j,] = site_output$mean_harvest_gCm2day
               grid_output$mean_fire_gCm2day[slot_i,slot_j,] = site_output$mean_fire_gCm2day
               grid_output$mean_nbe_gCm2day[slot_i,slot_j,] = site_output$mean_nbe_gCm2day
               grid_output$mean_nbp_gCm2day[slot_i,slot_j,] = site_output$mean_nbp_gCm2day
               # Time varying pixel based (i.e. not within the grid) values with quantile based uncertainty for...
               # States
               grid_output$Ctotal_gCm2[n,,] = site_output$Ctotal_gCm2
               grid_output$lai_m2m2[n,,]    = site_output$lai_m2m2
               # Fluxes
               grid_output$nee_gCm2day[n,,]     = site_output$nee_gCm2day
               grid_output$gpp_gCm2day[n,,]     = site_output$gpp_gCm2day
               grid_output$rauto_gCm2day[n,,]   = site_output$rauto_gCm2day
               grid_output$rhet_gCm2day[n,,]    = site_output$rhet_gCm2day
               grid_output$reco_gCm2day[n,,]    = site_output$reco_gCm2day
               grid_output$npp_gCm2day[n,,]     = site_output$npp_gCm2day
               grid_output$harvest_gCm2day[n,,] = site_output$harvest_gCm2day
               grid_output$fire_gCm2day[n,,]    = site_output$fire_gCm2day
               grid_output$nbe_gCm2day[n,,]     = site_output$nbe_gCm2day
               grid_output$nbp_gCm2day[n,,]     = site_output$nbp_gCm2day

               # Based on the presence of each pool define the grids for the mean and final values.
               # Also, create the time varying but quantile based values and time

               # Gridded biomass information
               if (exists(x = "biomass_gCm2", where = site_output)) {
                   # Grid mean / finals for globally available variables
                   grid_output$mean_biomass_gCm2[slot_i,slot_j,] = site_output$mean_biomass_gCm2
                   grid_output$final_biomass_gCm2[slot_i,slot_j,] = site_output$biomass_gCm2[,grid_output$time_dim]
                   grid_output$final_dCbiomass_gCm2[slot_i,slot_j,] = site_output$dCbiomass_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_biomass_gCm2day
                   grid_output$mean_biomass_to_litter_gCm2[slot_i,slot_j,] = site_output$mean_biomass_to_litter_gCm2
                   # Time varying pixel specific with quantiles
                   grid_output$biomass_gCm2[n,,] = site_output$biomass_gCm2
                   grid_output$dCbiomass_gCm2[n,,] = site_output$dCbiomass_gCm2
                   grid_output$outflux_biomass_gCm2day[n,,] = site_output$outflux_biomass_gCm2day
                   grid_output$biomass_to_litter_gCm2day[n,,] = site_output$biomass_to_litter_gCm2day
                   # Annual information
                   grid_output$MTT_annual_biomass_years[n,,] = site_output$MTT_annual_biomass_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_biomass
                   grid_output$FireFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$FireFractionOfTurnover_biomass
                   grid_output$HarvestFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_biomass
                   # Conditional variables
                   if (exists(x = "FIREemiss_biomass_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_biomass_gCm2day[n,,] = site_output$FIREemiss_biomass_gCm2day
                       grid_output$mean_FIREemiss_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_biomass_gCm2day
                   }
                   if (exists(x = "FIRElitter_biomass_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_biomass_gCm2day[n,,] = site_output$FIRElitter_biomass_gCm2day
                       grid_output$mean_FIRElitter_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_biomass_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_biomass_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_biomass_gCm2day[n,,] = site_output$HARVESTextracted_biomass_gCm2day
                       grid_output$mean_HARVESTextracted_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_biomass_gCm2day
                       grid_output$HARVESTlitter_biomass_gCm2day[n,,] = site_output$HARVESTlitter_biomass_gCm2day
                       grid_output$mean_HARVESTlitter_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_biomass_gCm2day
                   }
               }
               # Gridded labile information
               if (exists(x = "labile_gCm2", where = site_output)) {
                   # Grid averages
                   grid_output$mean_labile_gCm2[slot_i,slot_j,] = site_output$mean_labile_gCm2
                   grid_output$final_labile_gCm2[slot_i,slot_j,] = site_output$labile_gCm2[,grid_output$time_dim]
                   grid_output$final_dClabile_gCm2[slot_i,slot_j,] = site_output$dClabile_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_labile_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_labile_gCm2day
                   grid_output$mean_labile_to_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_labile_to_foliage_gCm2day
                   grid_output$mean_alloc_labile_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_labile_gCm2day
                   # Time varying
                   grid_output$labile_gCm2[n,,] = site_output$labile_gCm2
                   grid_output$dClabile_gCm2[n,,] = site_output$dClabile_gCm2
                   grid_output$outflux_labile_gCm2day[n,,] = site_output$outflux_labile_gCm2day
                   grid_output$labile_to_foliage_gCm2day[n,,] = site_output$labile_to_foliage_gCm2day
                   grid_output$alloc_labile_gCm2day[n,,] = site_output$alloc_labile_gCm2day
                   # Annual information
                   grid_output$MTT_annual_labile_years[n,,] = site_output$MTT_annual_labile_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_labile[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_labile
                   grid_output$FireFractionOfTurnover_labile[slot_i,slot_j,] = site_output$FireFractionOfTurnover_labile
                   grid_output$HarvestFractionOfTurnover_labile[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_labile
                   # Conditional variables
                   if (exists(x = "FIREemiss_labile_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_labile_gCm2day[n,,] = site_output$FIREemiss_labile_gCm2day
                       grid_output$mean_FIREemiss_labile_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_labile_gCm2day
                   }
                   if (exists(x = "FIRElitter_labile_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_labile_gCm2day[n,,] = site_output$FIRElitter_labile_gCm2day
                       grid_output$mean_FIRElitter_labile_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_labile_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_labile_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_labile_gCm2day[n,,] = site_output$HARVESTextracted_labile_gCm2day
                       grid_output$mean_HARVESTextracted_labile_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_labile_gCm2day
                       grid_output$HARVESTlitter_labile_gCm2day[n,,] = site_output$HARVESTlitter_labile_gCm2day
                       grid_output$mean_HARVESTlitter_labile_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_labile_gCm2day
                   }
               }
               # Gridded foliage information
               if (exists(x = "foliage_gCm2", where = site_output)) {
                   # Grided mean / final
                   grid_output$mean_foliage_gCm2[slot_i,slot_j,] = site_output$mean_foliage_gCm2
                   grid_output$final_foliage_gCm2[slot_i,slot_j,] = site_output$foliage_gCm2[,grid_output$time_dim]
                   grid_output$final_dCfoliage_gCm2[slot_i,slot_j,] = site_output$dCfoliage_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_foliage_gCm2day
                   grid_output$mean_foliage_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_foliage_to_litter_gCm2day
                   grid_output$mean_combined_alloc_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_combined_alloc_foliage_gCm2day
                   # Time varying pixel specific
                   grid_output$foliage_gCm2[n,,] = site_output$foliage_gCm2
                   grid_output$dCfoliage_gCm2[n,,] = site_output$dCfoliage_gCm2
                   grid_output$outflux_foliage_gCm2day[n,,] = site_output$outflux_foliage_gCm2day
                   grid_output$foliage_to_litter_gCm2day[n,,] = site_output$foliage_to_litter_gCm2day
                   grid_output$combined_alloc_foliage_gCm2day[n,,] = site_output$combined_alloc_foliage_gCm2day
                   # Annual information
                   grid_output$MTT_annual_foliage_years[n,,] = site_output$MTT_annual_foliage_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_foliage
                   grid_output$FireFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$FireFractionOfTurnover_foliage
                   grid_output$HarvestFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_foliage
                   # Conditional variables
                   if (exists(x = "alloc_foliage_gCm2day", where = site_output)) {
                       grid_output$alloc_foliage_gCm2day[n,,] = site_output$alloc_foliage_gCm2day
                       grid_output$mean_alloc_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_foliage_gCm2day
                   }
                   if (exists(x = "FIREemiss_foliage_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_foliage_gCm2day[n,,] = site_output$FIREemiss_foliage_gCm2day
                       grid_output$mean_FIREemiss_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_foliage_gCm2day
                   }
                   if (exists(x = "FIRElitter_foliage_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_foliage_gCm2day[n,,] = site_output$FIRElitter_foliage_gCm2day
                       grid_output$mean_FIRElitter_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_foliage_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_foliage_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_foliage_gCm2day[n,,] = site_output$HARVESTextracted_foliage_gCm2day
                       grid_output$mean_HARVESTextracted_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_foliage_gCm2day
                       grid_output$HARVESTlitter_foliage_gCm2day[n,,] = site_output$HARVESTlitter_foliage_gCm2day
                       grid_output$mean_HARVESTlitter_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_foliage_gCm2day
                   }
               }
               # Gridded roots information
               if (exists(x = "roots_gCm2", where = site_output)) {
                   # Gridded mean / final
                   grid_output$mean_roots_gCm2[slot_i,slot_j,] = site_output$mean_roots_gCm2
                   grid_output$final_roots_gCm2[slot_i,slot_j,] = site_output$roots_gCm2[,grid_output$time_dim]
                   grid_output$final_dCroots_gCm2[slot_i,slot_j,] = site_output$dCroots_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_roots_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_roots_gCm2day
                   grid_output$mean_roots_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_roots_to_litter_gCm2day
                   grid_output$mean_alloc_roots_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_roots_gCm2day
                   grid_output$annual_max_roots_gCm2[slot_i,slot_j,] = site_output$annual_max_roots_gCm2
                   # Pixel specific time varying
                   grid_output$roots_gCm2[n,,] = site_output$roots_gCm2
                   grid_output$dCroots_gCm2[n,,] = site_output$dCroots_gCm2
                   grid_output$outflux_roots_gCm2day[n,,] = site_output$outflux_roots_gCm2day
                   grid_output$roots_to_litter_gCm2day[n,,] = site_output$roots_to_litter_gCm2day
                   grid_output$alloc_roots_gCm2day[n,,] = site_output$alloc_roots_gCm2day
                   # Annual information
                   grid_output$MTT_annual_roots_years[n,,] = site_output$MTT_annual_roots_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_roots[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_roots
                   grid_output$FireFractionOfTurnover_roots[slot_i,slot_j,] = site_output$FireFractionOfTurnover_roots
                   grid_output$HarvestFractionOfTurnover_roots[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_roots
                   # Is rooting depth calculate (m) by this model?
                   if (exists(x = "RootDepth_m", where = site_output)) {
                       grid_output$RootDepth_m[n,,] = site_output$RootDepth_m
                       grid_output$dRootDepth_m[n,,] = site_output$dRootDepth_m
                       grid_output$mean_RootDepth_m[slot_i,slot_j,] = site_output$mean_RootDepth_m
                       grid_output$final_dRootDepth_m[slot_i,slot_j,] = site_output$dRootDepth_m[,grid_output$time_dim]
                   }
                   # Conditional variables
                   if (exists(x = "FIREemiss_roots_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_roots_gCm2day[n,,] = site_output$FIREemiss_roots_gCm2day
                       grid_output$mean_FIREemiss_roots_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_roots_gCm2day
                   }
                   if (exists(x = "FIRElitter_roots_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_roots_gCm2day[n,,] = site_output$FIRElitter_roots_gCm2day
                       grid_output$mean_FIRElitter_roots_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_roots_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_roots_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_roots_gCm2day[n,,] = site_output$HARVESTextracted_roots_gCm2day
                       grid_output$mean_HARVESTextracted_roots_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_roots_gCm2day
                       grid_output$HARVESTlitter_roots_gCm2day[n,,] = site_output$HARVESTlitter_roots_gCm2day
                       grid_output$mean_HARVESTlitter_roots_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_roots_gCm2day
                   }
               }
               # Gridded wood information
               if (exists(x = "wood_gCm2", where = site_output)) {
                   # Gridded mean / final variables
                   grid_output$mean_wood_gCm2[slot_i,slot_j,] = site_output$mean_wood_gCm2
                   grid_output$final_wood_gCm2[slot_i,slot_j,] = site_output$wood_gCm2[,grid_output$time_dim]
                   grid_output$final_dCwood_gCm2[slot_i,slot_j,] = site_output$dCwood_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_wood_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_wood_gCm2day
                   grid_output$mean_wood_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_wood_to_litter_gCm2day
                   grid_output$mean_alloc_wood_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_wood_gCm2day
                   grid_output$annual_max_wood_gCm2[slot_i,slot_j,] = site_output$annual_max_wood_gCm2
                   # Pixel specific time varying
                   grid_output$wood_gCm2[n,,] = site_output$wood_gCm2
                   grid_output$dCwood_gCm2[n,,] = site_output$dCwood_gCm2
                   grid_output$outflux_wood_gCm2day[n,,] = site_output$outflux_wood_gCm2day
                   grid_output$wood_to_litter_gCm2day[n,,] = site_output$wood_to_litter_gCm2day
                   grid_output$alloc_wood_gCm2day[n,,] = site_output$alloc_wood_gCm2day
                   # Annual information
                   grid_output$MTT_annual_wood_years[n,,] = site_output$MTT_annual_wood_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_wood[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_wood
                   grid_output$FireFractionOfTurnover_wood[slot_i,slot_j,] = site_output$FireFractionOfTurnover_wood
                   grid_output$HarvestFractionOfTurnover_wood[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_wood
                   # Conditional variables
                   if (exists(x = "FIREemiss_wood_gCm2day", where = site_output)) {
                        grid_output$FIREemiss_wood_gCm2day[n,,] = site_output$FIREemiss_wood_gCm2day
                        grid_output$mean_FIREemiss_wood_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_wood_gCm2day
                   }
                   if (exists(x = "FIRElitter_wood_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_wood_gCm2day[n,,] = site_output$FIRElitter_wood_gCm2day
                       grid_output$mean_FIRElitter_wood_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_wood_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_wood_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_wood_gCm2day[n,,] = site_output$HARVESTextracted_wood_gCm2day
                       grid_output$mean_HARVESTextracted_wood_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_wood_gCm2day
                       grid_output$HARVESTlitter_wood_gCm2day[n,,] = site_output$HARVESTlitter_wood_gCm2day
                       grid_output$mean_HARVESTlitter_wood_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_wood_gCm2day
                   }
               }
               # Gridded litter information
               if (exists(x = "litter_gCm2", where = site_output)) {
                   # Gridded mean / final
                   grid_output$mean_litter_gCm2[slot_i,slot_j,] = site_output$mean_litter_gCm2
                   grid_output$final_litter_gCm2[slot_i,slot_j,] = site_output$litter_gCm2[,grid_output$time_dim]
                   grid_output$final_dClitter_gCm2[slot_i,slot_j,] = site_output$dClitter_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_litter_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_litter_gCm2day
                   grid_output$mean_litter_to_som_gCm2day[slot_i,slot_j,] = site_output$mean_litter_to_som_gCm2day
                   grid_output$mean_rhet_litter_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_litter_gCm2day
                   # Pixel specific time varying
                   grid_output$litter_gCm2[n,,] = site_output$litter_gCm2
                   grid_output$dClitter_gCm2[n,,] = site_output$dClitter_gCm2
                   grid_output$outflux_litter_gCm2day[n,,] = site_output$outflux_litter_gCm2day
                   grid_output$litter_to_som_gCm2day[n,,] = site_output$litter_to_som_gCm2day
                   grid_output$rhet_litter_gCm2day[n,,] = site_output$rhet_litter_gCm2day
                   # Annual information
                   grid_output$MTT_annual_litter_years[n,,] = site_output$MTT_annual_litter_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_litter[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_litter
                   grid_output$FireFractionOfTurnover_litter[slot_i,slot_j,] = site_output$FireFractionOfTurnover_litter
                   grid_output$HarvestFractionOfTurnover_litter[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_litter
                   # Conditional variables
                   if (exists(x = "FIREemiss_litter_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_litter_gCm2day[n,,] = site_output$FIREemiss_litter_gCm2day
                       grid_output$mean_FIREemiss_litter_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_litter_gCm2day
                   }
                   if (exists(x = "FIRElitter_litter_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_litter_gCm2day[n,,] = site_output$FIRElitter_litter_gCm2day
                       grid_output$mean_FIRElitter_litter_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_litter_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_litter_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_litter_gCm2day[n,,] = site_output$HARVESTextracted_litter_gCm2day
                       grid_output$mean_HARVESTextracted_litter_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_litter_gCm2day
                   }
               }
               # Gridded wood litter information
               if (exists(x = "woodlitter_gCm2", where = site_output)) {
                   # Gridded mean / final
                   grid_output$mean_woodlitter_gCm2[slot_i,slot_j,] = site_output$mean_woodlitter_gCm2
                   grid_output$final_woodlitter_gCm2[slot_i,slot_j,] = site_output$woodlitter_gCm2[,grid_output$time_dim]
                   grid_output$final_dCwoodlitter_gCm2[slot_i,slot_j,] = site_output$dCwoodlitter_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_woodlitter_gCm2day
                   grid_output$mean_woodlitter_to_som_gCm2day[slot_i,slot_j,] = site_output$mean_woodlitter_to_som_gCm2day
                   grid_output$mean_rhet_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_woodlitter_gCm2day
                   # Pixel specific time varying
                   grid_output$woodlitter_gCm2[n,,] = site_output$woodlitter_gCm2
                   grid_output$dCwoodlitter_gCm2[n,,] = site_output$dCwoodlitter_gCm2
                   grid_output$outflux_woodlitter_gCm2day[n,,] = site_output$outflux_woodlitter_gCm2day
                   grid_output$woodlitter_to_som_gCm2day[n,,] = site_output$woodlitter_to_som_gCm2day
                   grid_output$rhet_woodlitter_gCm2day[n,,] = site_output$rhet_woodlitter_gCm2day
                   # Annual information
                   grid_output$MTT_annual_woodlitter_years[n,,] = site_output$MTT_annual_woodlitter_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_woodlitter
                   grid_output$FireFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$FireFractionOfTurnover_woodlitter
                   grid_output$HarvestFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_woodlitter
                   # Conditional variables
                   if (exists(x = "FIREemiss_woodlitter_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_woodlitter_gCm2day[n,,] = site_output$FIREemiss_woodlitter_gCm2day
                       grid_output$mean_FIREemiss_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_woodlitter_gCm2day
                   }
                   if (exists(x = "FIRElitter_woodlitter_gCm2day", where = site_output)) {
                       grid_output$FIRElitter_woodlitter_gCm2day[n,,] = site_output$FIRElitter_woodlitter_gCm2day
                       grid_output$mean_FIRElitter_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_woodlitter_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_woodlitter_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_woodlitter_gCm2day[n,,] = site_output$HARVESTextracted_woodlitter_gCm2day
                       grid_output$mean_HARVESTextracted_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_woodlitter_gCm2day
                   }
               }
               # Gridded som information
               if (exists(x = "som_gCm2", where = site_output)) {
                   # Gridded mean / finals
                   grid_output$mean_som_gCm2[slot_i,slot_j,] = site_output$mean_som_gCm2
                   grid_output$final_som_gCm2[slot_i,slot_j,] = site_output$som_gCm2[,grid_output$time_dim]
                   grid_output$final_dCsom_gCm2[slot_i,slot_j,] = site_output$dCsom_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_som_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_som_gCm2day
                   grid_output$mean_rhet_som_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_som_gCm2day
                   # Pixel specific time varying
                   grid_output$som_gCm2[n,,] = site_output$som_gCm2
                   grid_output$dCsom_gCm2[n,,] = site_output$dCsom_gCm2
                   grid_output$outflux_som_gCm2day[n,,] = site_output$outflux_som_gCm2day
                   grid_output$rhet_som_gCm2day[n,,] = site_output$rhet_som_gCm2day
                   # Annual information
                   grid_output$MTT_annual_som_years[n,,] = site_output$MTT_annual_som_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_som[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_som
                   grid_output$FireFractionOfTurnover_som[slot_i,slot_j,] = site_output$FireFractionOfTurnover_som
                   grid_output$HarvestFractionOfTurnover_som[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_som
                   # Time varying
                   if (exists(x = "FIREemiss_som_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_som_gCm2day[n,,] = site_output$FIREemiss_som_gCm2day
                       grid_output$mean_FIREemiss_som_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_som_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_som_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_som_gCm2day[n,,] = site_output$HARVESTextracted_som_gCm2day
                       grid_output$mean_HARVESTextracted_som_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_som_gCm2day
                   }
               }
               # Gridded dom information
               if (exists(x = "dom_gCm2", where = site_output)) {
                   # Gridded mean / finals
                   grid_output$mean_dom_gCm2[slot_i,slot_j,] = site_output$mean_dom_gCm2
                   grid_output$final_dom_gCm2[slot_i,slot_j,] = site_output$dom_gCm2[,grid_output$time_dim]
                   grid_output$final_dCdom_gCm2[slot_i,slot_j,] = site_output$dCdom_gCm2[,grid_output$time_dim]
                   grid_output$mean_outflux_dom_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_dom_gCm2day
                   grid_output$mean_rhet_dom_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_dom_gCm2day
                   # Pixel specific time varying
                   grid_output$dom_gCm2[n,,] = site_output$dom_gCm2
                   grid_output$dCdom_gCm2[n,,] = site_output$dCdom_gCm2
                   grid_output$outflux_dom_gCm2day[n,,] = site_output$outflux_dom_gCm2day
                   grid_output$rhet_dom_gCm2day[n,,] = site_output$rhet_dom_gCm2day
                   # Annual information
                   grid_output$MTT_annual_dom_years[n,,] = site_output$MTT_annual_dom_years
                   # Fractional partitioning of tunover to different drivers - should they exist
                   grid_output$NaturalFractionOfTurnover_dom[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_dom
                   grid_output$FireFractionOfTurnover_dom[slot_i,slot_j,] = site_output$FireFractionOfTurnover_dom
                   grid_output$HarvestFractionOfTurnover_dom[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_dom
                   # Conditional
                   if (exists(x = "FIREemiss_dom_gCm2day", where = site_output)) {
                       grid_output$FIREemiss_dom_gCm2day[n,,] = site_output$FIREemiss_dom_gCm2day
                       grid_output$mean_FIREemiss_dom_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_dom_gCm2day
                   }
                   if (exists(x = "HARVESTextracted_dom_gCm2day", where = site_output)) {
                       grid_output$HARVESTextracted_dom_gCm2day[n,,] = site_output$HARVESTextracted_dom_gCm2day
                       grid_output$mean_HARVESTextracted_dom_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_dom_gCm2day
                   }
               }

               # Water cycle specific variables
               if (exists(x = "ET_kgH2Om2day", where = site_output)) {
                   # currently water in the soil surface layer (0-30 cm)
                   grid_output$mean_SurfWater_kgH2Om2[slot_i,slot_j,] = site_output$mean_SurfWater_kgH2Om2
                   grid_output$final_SurfWater_kgH2Om2[slot_i,slot_j,] = site_output$SurfWater_kgH2Om2[,grid_output$time_dim]
                   grid_output$final_dSurfWater_kgH2Om2[slot_i,slot_j,] = site_output$dSurfWater_kgH2Om2[,grid_output$time_dim]
                   grid_output$SurfWater_kgH2Om2[n,,] = site_output$SurfWater_kgH2Om2
                   grid_output$dSurfWater_kgH2Om2[n,,] = site_output$dSurfWater_kgH2Om2
                   # plant apparent soil water potential (MPa)
                   grid_output$mean_wSWP_MPa[slot_i,slot_j,] = site_output$mean_wSWP_MPa
                   grid_output$final_wSWP_MPa[slot_i,slot_j,] = site_output$wSWP_MPa[,grid_output$time_dim]
                   grid_output$wSWP_MPa[n,,] = site_output$wSWP_MPa
                   grid_output$dwSWP_MPa[n,,] = site_output$dwSWP_MPa
                   # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                   grid_output$mean_ET_kgH2Om2day[slot_i,slot_j,] = site_output$mean_ET_kgH2Om2day
                   grid_output$ET_kgH2Om2day[n,,] = site_output$ET_kgH2Om2day
               }
               # Snow specific
               if (exists(x = "snow_kgH2Om2", where = site_output)) {
                   ## snow on soil surface
                   grid_output$mean_snow_kgH2Om2[slot_i,slot_j,] = site_output$mean_snow_kgH2Om2
                   grid_output$snow_kgH2Om2[n,,] = site_output$snow_kgH2Om2
               }
               # Canopy process variables
               if (exists(x = "APAR_MJm2day", where = site_output)) {
                   # Absorbed photosynthetically active radation
                   grid_output$mean_APAR_MJm2day[slot_i,slot_j,] = site_output$mean_APAR_MJm2day
                   grid_output$APAR_MJm2day[n,,] = site_output$APAR_MJm2day
               }
               if (exists(x = "CiCa", where = site_output)) {
                   # Canopy Ci:Ca
                   grid_output$mean_CiCa[slot_i,slot_j,] = site_output$mean_CiCa
                   grid_output$CiCa[n,,] = site_output$CiCa
               }
               if (exists(x = "gs_demand_supply_ratio", where = site_output)) {
                   # Ratio of stomatal conductance relative to its maximum value,
                   # this metric provides information on the demand vs supply constrains on stomatal conductance
                   grid_output$mean_gs_demand_supply_ratio[slot_i,slot_j,] = site_output$mean_gs_demand_supply_ratio
                   grid_output$gs_demand_supply_ratio[n,,] = site_output$gs_demand_supply_ratio
               }
               if (exists(x = "gs_mmolH2Om2day", where = site_output)) {
                   # Canopy stomatal conductance
                   grid_output$mean_gs_mmolH2Om2day[slot_i,slot_j,] = site_output$mean_gs_mmolH2Om2day
                   grid_output$gs_mmolH2Om2day[n,,] = site_output$gs_mmolH2Om2day
               }
               if (exists(x = "gb_mmolH2Om2day", where = site_output)) {
                   # Canopy boundary layer conductance
                   grid_output$mean_gb_mmolH2Om2day[slot_i,slot_j,] = site_output$mean_gb_mmolH2Om2day
                   grid_output$gb_mmolH2Om2day[n,,] = site_output$gb_mmolH2Om2day
               }

               # Any time series assimilated data overlaps?
               if (exists(x = "gpp_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$gpp_assim_data_overlap_fraction[slot_i,slot_j] = site_output$gpp_assim_data_overlap_fraction
               }
               if (exists(x = "lai_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$lai_assim_data_overlap_fraction[slot_i,slot_j] = site_output$lai_assim_data_overlap_fraction
               }
               if (exists(x = "nee_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$nee_assim_data_overlap_fraction[slot_i,slot_j] = statsite_outputes_all$nee_assim_data_overlap_fraction
               }
               if (exists(x = "wood_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$wood_assim_data_overlap_fraction[slot_i,slot_j] = site_output$wood_assim_data_overlap_fraction
               }
               if (exists(x = "soil_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$soil_assim_data_overlap_fraction[slot_i,slot_j] = site_output$soil_assim_data_overlap_fraction
               }
               if (exists(x = "et_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$et_assim_data_overlap_fraction[slot_i,slot_j] = site_output$et_assim_data_overlap_fraction
               }
               if (exists(x = "nbe_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$nbe_assim_data_overlap_fraction[slot_i,slot_j] = site_output$nbe_assim_data_overlap_fraction
               }
               if (exists(x = "fire_assim_data_overlap_fraction", where = site_output)) {
                   grid_output$fire_assim_data_overlap_fraction[slot_i,slot_j] = site_output$fire_assim_data_overlap_fraction
               }
               # Parameter vs C-cycle flux correlation across ensemble member
               grid_output$nee_par_cor[slot_i,slot_j,] = site_output$nee_par_cor
               grid_output$gpp_par_cor[slot_i,slot_j,] = site_output$gpp_par_cor
               grid_output$rauto_par_cor[slot_i,slot_j,] = site_output$rauto_par_cor
               grid_output$rhet_par_cor[slot_i,slot_j,] = site_output$rhet_par_cor
               grid_output$fire_par_cor[slot_i,slot_j,] = site_output$fire_par_cor

               # now tidy away the file
               file.remove(outfile_stock_fluxes)

          } # if site has been done and can now be made part of the overall...

      } # looping through sites

      # now save the combined grid file
      save(grid_output, file=outfile_grid, compress = "gzip", compression_level = 9)

  } # gridded run?

  # tell me whats happening
  print(paste("...time to process ",round((proc.time()["elapsed"]-stime)/60,1)," minutes",sep=""))

} # end function run_mcmc_results

## Use byte compile
run_mcmc_results<-cmpfun(run_mcmc_results)
