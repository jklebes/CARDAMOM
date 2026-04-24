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
# Functions to derive stocks and fluxes used in the calculation of gridded aggregates
# These are variables which for a site analysis would be easy to calculate
# from the ensembles but difficult if not determined here and now before aggregation
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

post_process_dalec<-function(states_all,parameters,drivers,PROJECT,n) {

  # Determine some useful information for the analysis below
  nos_years = PROJECT$nos_years
  steps_per_year = floor(dim(drivers$met)[1] / nos_years)

  # Check the list variables in states_all which we will be searching
  check_list = names(states_all)

  # If a combined ecosystem heterotrophic respiration flux does not
  # exist we shall calculate it
  if (any(check_list == "rhet_gCm2day") == FALSE) {
      if (any(check_list == "rhet_dom_gCm2day")) {
          states_all$rhet_gCm2day = states_all$rhet_dom_gCm2day
      } else {
          # Calculate the combined ecosystem heterotrophic respiration.
          # All models have a som pool, so start with that
          states_all$rhet_gCm2day = states_all$rhet_som_gCm2day
          states_all$mean_rhet_gCm2day = states_all$mean_rhet_som_gCm2day
          states_all$mean_annual_rhet_gCm2day = states_all$mean_annual_rhet_som_gCm2day
          # If the model has a litter pool (foliar + fine root) add this
          if (any(check_list == "rhet_litter_gCm2day")) {
              states_all$rhet_gCm2day = states_all$rhet_gCm2day + states_all$rhet_litter_gCm2day
              states_all$mean_rhet_gCm2day = states_all$mean_rhet_gCm2day + states_all$mean_rhet_litter_gCm2day
              states_all$mean_annual_rhet_gCm2day = states_all$mean_annual_rhet_gCm2day + states_all$mean_annual_rhet_litter_gCm2day
          }
          # If the model has a wood litter pool add this
          if (any(check_list == "rhet_woodlitter_gCm2day")) {
              states_all$rhet_gCm2day = states_all$rhet_gCm2day + states_all$rhet_woodlitter_gCm2day
              states_all$mean_rhet_gCm2day = states_all$mean_rhet_gCm2day + states_all$mean_rhet_woodlitter_gCm2day
              states_all$mean_annual_rhet_gCm2day = states_all$mean_annual_rhet_gCm2day + states_all$mean_annual_rhet_woodlitter_gCm2day
          }
      } # does rhet_dom_gCm2day exist?
  } # does rhet_gCm2day exist?
  # Combine autotrophic and heterotrophic respiration into ecosystem respiration
  states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  states_all$mean_reco_gCm2day = states_all$mean_rauto_gCm2day + states_all$mean_rhet_gCm2day
  states_all$mean_annual_reco_gCm2day = states_all$mean_annual_rauto_gCm2day + states_all$mean_annual_rhet_gCm2day
  # Calculate the net ecosystem exchange of CO2
  states_all$nee_gCm2day = states_all$reco_gCm2day - states_all$gpp_gCm2day
  states_all$mean_nee_gCm2day = states_all$mean_reco_gCm2day - states_all$mean_gpp_gCm2day
  states_all$mean_annual_nee_gCm2day = states_all$mean_annual_reco_gCm2day - states_all$mean_annual_gpp_gCm2day
  # Calculate net primary productivity
  states_all$npp_gCm2day = states_all$gpp_gCm2day - states_all$rauto_gCm2day
  states_all$mean_npp_gCm2day = states_all$mean_gpp_gCm2day - states_all$mean_rauto_gCm2day
  states_all$mean_annual_npp_gCm2day = states_all$mean_annual_gpp_gCm2day - states_all$mean_annual_rauto_gCm2day
  # Begin calculation of net biome exchange and net biome productivity
  states_all$nbe_gCm2day = states_all$nee_gCm2day  # negative = sink
  states_all$mean_nbe_gCm2day = states_all$mean_nee_gCm2day  # negative = sink
  states_all$mean_annual_nbe_gCm2day = states_all$mean_annual_nee_gCm2day  # negative = sink
  states_all$nbp_gCm2day = -states_all$nee_gCm2day # positive = sink
  states_all$mean_nbp_gCm2day = -states_all$mean_nee_gCm2day # positive = sink
  states_all$mean_annual_nbp_gCm2day = -states_all$mean_annual_nee_gCm2day # positive = sink
  # If fire exists then update the NBE and NBP accordingly
  if (any(check_list == "fire_gCm2day")) {
      states_all$nbe_gCm2day = states_all$nbe_gCm2day + states_all$fire_gCm2day
      states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$fire_gCm2day
      states_all$mean_nbe_gCm2day = states_all$mean_nbe_gCm2day + states_all$mean_fire_gCm2day
      states_all$mean_nbp_gCm2day = states_all$mean_nbp_gCm2day - states_all$mean_fire_gCm2day
      states_all$mean_annual_nbe_gCm2day = states_all$mean_annual_nbe_gCm2day + states_all$mean_annual_fire_gCm2day
      states_all$mean_annual_nbp_gCm2day = states_all$mean_annual_nbp_gCm2day - states_all$mean_annual_fire_gCm2day
  }
  # If a harvest flux exists update the NBP. NOTE: that this harvest flux
  # specifically accouts for C removed, there may be mortality due to harvest
  # but remains in system as residues.
  if (any(check_list == "harvest_gCm2day")) {
      states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$harvest_gCm2day
      states_all$mean_nbp_gCm2day = states_all$mean_nbp_gCm2day - states_all$mean_harvest_gCm2day
      states_all$mean_annual_nbp_gCm2day = states_all$mean_annual_nbp_gCm2day - states_all$mean_annual_harvest_gCm2day
  }
  # In managed grassland systems (and possible, in the future others) part of the NBP is extracted via animal grazing
  if (any(check_list == "grazing_gCm2day")) {
      states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$grazing_gCm2day
      states_all$mean_nbp_gCm2day = states_all$mean_nbp_gCm2day - states_all$mean_grazing_gCm2day
      states_all$mean_annual_nbp_gCm2day = states_all$mean_annual_nbp_gCm2day - states_all$mean_annual_grazing_gCm2day
  }
  # In some cases, e.g. crop models, the harvest variable includes just the harvested yield. But the NBP requires
  # tracking of extracted non-yield C
  if (any(check_list == "extracted_residue_gCm2day")) {
      states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$extracted_residue_gCm2day
      states_all$mean_nbp_gCm2day = states_all$mean_nbp_gCm2day - states_all$mean_extracted_residue_gCm2day
      states_all$mean_annual_nbp_gCm2day = states_all$mean_annual_nbp_gCm2day - states_all$mean_annual_extracted_residue_gCm2day
  }

  # Now calculate the mean annual carbon use efficiency (NPP:GPP) as some models do now have a parameter for this
  states_all$mean_annual_cue = states_all$mean_annual_npp_gCm2day / states_all$mean_annual_gpp_gCm2day
  
  ###
  ## Post-hoc calculation of parameter correlations with key C-cycle variables

  # Construct and rearrange array of parameter suitable for correlation determination
  tmp = t(array(as.vector(parameters[1:PROJECT$model$nopars[n],,]),dim=c(PROJECT$model$nopars[n],prod(dim(parameters)[2:3]))))
  # Determine the correlation matrix between all parameters
  states_all$absolute_mean_parameter_correlation = cor(tmp)
  # Determine the mean of the absolute correlations from the matrix
  states_all$absolute_mean_parameter_correlation = mean(abs(states_all$absolute_mean_parameter_correlation[lower.tri(states_all$absolute_mean_parameter_correlation,diag=FALSE)]))

  # Determine multi-use variables for the correlation analysis below
  ensLAI = rowMeans(states_all$lai_m2m2)
  ensNBP = rowMeans(states_all$nbp_gCm2day)
  ensNEE = rowMeans(states_all$nee_gCm2day)
  ensGPP = rowMeans(states_all$gpp_gCm2day)
  ensRauto = rowMeans(states_all$rauto_gCm2day)
  ensRhet = rowMeans(states_all$rhet_gCm2day)

  # Determine correlations between parameter values and various state variables
  states_all$lai_parameter_correlation = cor(tmp,ensLAI)
  states_all$nbp_parameter_correlation = cor(tmp,ensNBP)
  states_all$nee_parameter_correlation = cor(tmp,ensNEE)
  states_all$gpp_parameter_correlation = cor(tmp,ensGPP)
  states_all$rauto_parameter_correlation = cor(tmp,ensRauto)
  states_all$rhet_parameter_correlation = cor(tmp,ensRhet)
  states_all$CiCa_parameter_correlation = cor(tmp,rowMeans(states_all$CiCa))
  # Avoid error flag when no fire
  if (any(check_list == "fire_gCm2day")) {
      if (max(as.vector(states_all$fire_gCm2day)) > 0) {
          states_all$fire_parameter_correlation = cor(tmp,rowMeans(states_all$fire_gCm2day))
      } else {
          states_all$fire_parameter_correlation = array(0, dim = c(PROJECT$model$nopars[n],1))
      }
  }
  # Avoid error flag when no LWP
  if (any(check_list == "LWP_MPa")) {
      states_all$LWP_parameter_correlation = cor(tmp,rowMeans(states_all$LWP_MPa))
  }

  # Do wood related against common state / flux variables
  if (any(check_list == "wood_gCm2")) {
      # Determine time varying difference since t=1
      dCwood = rowMeans(states_all$wood_gCm2 - states_all$wood_gCm2[,1]) # difference in wood from initial
      # Determine the time varying difference
      ensWood = rowMeans(states_all$wood_gCm2)
      # Correlation between wood stock and parameters
      states_all$wood_parameter_correlation = cor(tmp,ensWood)
      # Correlations with LAI
      states_all$lai_m2m2_to_wood_gCm2_correlation = cor(ensLAI,ensWood)      
      states_all$lai_m2m2_to_dCwood_gCm2_correlation = cor(ensLAI,dCwood)      
      # Correlations with NBP
      states_all$NBP_gCm2day_to_wood_gCm2_correlation = cor(ensNBP,ensWood)
      states_all$NBP_gCm2day_to_dCwood_gCm2_correlation = cor(ensNBP,dCwood)
  }
  if (any(check_list == "som_gCm2")) {
      # Determine time varying difference since t=1
      dCsom = rowMeans(states_all$som_gCm2 - states_all$som_gCm2[,1]) # difference in som from initial  
      # Determine time varying difference
      ensSOM = rowMeans(states_all$som_gCm2)
      # Correlation between som stock and parameters
      states_all$som_parameter_correlation = cor(tmp,ensSOM)
      # Correlations with LAI
      states_all$lai_m2m2_to_som_gCm2_correlation = cor(ensLAI,ensSOM)
      states_all$lai_m2m2_to_dCsom_gCm2_correlation = cor(ensLAI,dCsom)      
      # Correlation with NBP
      states_all$NBP_gCm2day_to_dCsom_gCm2_correlation = cor(ensNBP,dCsom)      
      states_all$NBP_gCm2day_to_som_gCm2_correlation = cor(ensNBP,ensSOM)
  }

  # Correlations between LAI and common gross and net fluxes
  states_all$lai_m2m2_to_GPP_gCm2day_correlation = cor(ensLAI,ensGPP)
  states_all$lai_m2m2_to_NEE_gCm2day_correlation = cor(ensLAI,ensNEE)
  states_all$lai_m2m2_to_NBP_gCm2day_correlation = cor(ensLAI,ensNBP)
  states_all$lai_m2m2_to_Rauto_gCm2day_correlation = cor(ensLAI,ensRauto)
  states_all$lai_m2m2_to_Rhet_gCm2day_correlation = cor(ensLAI,ensRhet)

  # If harvest is estimated
  if (any(check_list == "harvest_gCm2day")) {
      states_all$lai_m2m2_to_harvest_gCm2day_correlation = cor(ensLAI,rowMeans(states_all$harvest_gCm2day))
  }

  # Correlations between NBP and key gross and net fluxes
  states_all$NBP_gCm2day_to_GPP_gCm2day_correlation = cor(ensNBP,ensGPP)
  states_all$NBP_gCm2day_to_NEE_gCm2day_correlation = cor(ensNBP,ensNEE)
  states_all$NBP_gCm2day_to_lai_m2m2_correlation = cor(ensNBP,ensLAI)
  states_all$NBP_gCm2day_to_Rauto_gCm2day_correlation = cor(ensNBP,ensRauto)
  states_all$NBP_gCm2day_to_Rhet_gCm2day_correlation = cor(ensNBP,ensRhet)

  # If harvest is estimated
  if (any(check_list == "harvest_gCm2day")) {
      states_all$NBP_gCm2day_to_harvest_gCm2day_correlation = cor(ensNBP,rowMeans(states_all$harvest_gCm2day))
  }

  # Determine whether have have both mean transit time and allocation to wood
  if (any(check_list == "MTT_wood_years") && (any(check_list == "alloc_wood_gCm2day") || any(check_list == "labile_to_wood_gCm2day"))) {
      
      # Determine where the carbon inputs are coming from
      if (any(check_list == "alloc_wood_gCm2day") & any(check_list == "alloc_wolabile_to_wood_gCm2dayod_gCm2day")) {
          # Multi-use variable
          ensAwood = rowMeans(states_all$alloc_wood_gCm2day+states_all$labile_to_wood_gCm2day)
      } else if (any(check_list == "alloc_wood_gCm2day")) {
          # Multi-use variable
          ensAwood = rowMeans(states_all$alloc_wood_gCm2day)
      } else if (any(check_list == "labile_to_wood_gCm2day")){
          # Multi-use variable
          ensAwood = rowMeans(states_all$labile_to_wood_gCm2day)
      }

      # As both exist determine their correlations with parameters...
      states_all$MTT_wood_years_parameter_correlation = cor(tmp,states_all$MTT_wood_years)
      states_all$NPP_wood_gCm2day_parameter_correlation = cor(tmp,ensAwood)
      states_all$NPP_wood_fraction_parameter_correlation = cor(tmp,ensAwood)
      # ...against key other observables
      # NPP wood - Flux
      states_all$NPP_wood_gCm2day_to_GPP_gCm2day_correlation = cor(ensAwood,ensGPP)
      states_all$NPP_wood_gCm2day_to_NEE_gCm2day_correlation = cor(ensAwood,ensNEE)
      states_all$NPP_wood_gCm2day_to_Rauto_gCm2day_correlation = cor(ensAwood,ensRauto)
      states_all$NPP_wood_gCm2day_to_Rhet_gCm2day_correlation = cor(ensAwood,ensRhet)   
      states_all$NPP_wood_gCm2day_to_wood_gCm2_correlation = cor(ensAwood,ensWood)   
      states_all$NPP_wood_gCm2day_to_som_gCm2_correlation = cor(ensAwood,ensSOM)   
      states_all$NPP_wood_gCm2day_to_lai_m2m2_correlation = cor(ensAwood,ensLAI)   
      states_all$NPP_wood_gCm2day_to_dCwood_gCm2_correlation = cor(ensAwood,dCwood)         
      states_all$NPP_wood_gCm2day_to_dCsom_gCm2_correlation = cor(ensAwood,dCsom)         
      # NPP wood - fraction
      states_all$NPP_wood_fraction_to_GPP_gCm2day_correlation = cor(states_all$NPP_wood_fraction,ensGPP)
      states_all$NPP_wood_fraction_to_NEE_gCm2day_correlation = cor(states_all$NPP_wood_fraction,ensNEE)
      states_all$NPP_wood_fraction_to_Rauto_gCm2day_correlation = cor(states_all$NPP_wood_fraction,ensRauto)
      states_all$NPP_wood_fraction_to_Rhet_gCm2day_correlation = cor(states_all$NPP_wood_fraction,ensRhet)   
      states_all$NPP_wood_fraction_to_wood_gCm2_correlation = cor(states_all$NPP_wood_fraction,ensWood)   
      states_all$NPP_wood_fraction_to_som_gCm2_correlation = cor(states_all$NPP_wood_fraction,ensSOM)   
      states_all$NPP_wood_fraction_to_lai_m2m2_correlation = cor(states_all$NPP_wood_fraction,ensLAI)     
      states_all$NPP_wood_fraction_to_dCwood_gCm2_correlation = cor(states_all$NPP_wood_fraction,dCwood)         
      states_all$NPP_wood_fraction_to_dCsom_gCm2_correlation = cor(states_all$NPP_wood_fraction,dCsom)                   
      # MTT wood
      states_all$MTT_wood_years_to_GPP_gCm2day_correlation = cor(states_all$MTT_wood_years,ensGPP)
      states_all$MTT_wood_years_to_NEE_gCm2day_correlation = cor(states_all$MTT_wood_years,ensNEE)
      states_all$MTT_wood_years_to_Rauto_gCm2day_correlation = cor(states_all$MTT_wood_years,ensRauto)
      states_all$MTT_wood_years_to_Rhet_gCm2day_correlation = cor(states_all$MTT_wood_years,ensRhet)
      states_all$MTT_wood_years_to_wood_gCm2_correlation = cor(states_all$MTT_wood_years,ensWood)   
      states_all$MTT_wood_years_to_som_gCm2_correlation = cor(states_all$MTT_wood_years,ensSOM)   
      states_all$MTT_wood_years_to_lai_m2m2_correlation = cor(states_all$MTT_wood_years,ensLAI)      
      states_all$MTT_wood_years_to_dCwood_gCm2_correlation = cor(states_all$MTT_wood_years,dCwood)         
      states_all$MTT_wood_years_to_dCsom_gCm2_correlation = cor(states_all$MTT_wood_years,dCsom)                            
      # ...and with each other
      states_all$MTT_wood_years_to_NPP_wood_gCm2day_correlation = cor(states_all$MTT_wood_years,ensAwood)
      states_all$MTT_wood_years_to_NPP_wood_fraction_correlation = cor(states_all$MTT_wood_years,states_all$NPP_wood_fraction)
      states_all$MTT_wood_years_to_MTT_som_years_correlation = cor(states_all$MTT_wood_years,states_all$MTT_som_years)      
      # As MTT wood exists we must also be able to work out the correlations for change in wood over time
      states_all$dCwood_gCm2_to_gpp_gCm2day_correlation = cor(dCwood,ensGPP) 
      states_all$dCwood_gCm2_to_rauto_gCm2day_correlation = cor(dCwood,ensRauto) 
      states_all$dCwood_gCm2_to_nee_gCm2day_correlation = cor(dCwood,ensNEE)    
      states_all$dCwood_gCm2_to_rhet_gCm2day_correlation = cor(dCwood,ensRhet)    
      states_all$dCwood_gCm2_to_wood_gCm2_correlation = cor(dCwood,ensWood)          
      states_all$dCwood_gCm2_to_som_gCm2_correlation = cor(dCwood,ensSOM)      
      states_all$dCwood_gCm2_to_dCsom_gCm2_correlation = cor(dCwood,dCsom)          
      # Remove multi-use variable
      rm(ensAwood)
  } else {

      # Both are not present, so we will determine whether we can generate one of the correlation estimates

      # If Mean transit time for wood is provided generate a correlation estimate
      if (any(check_list == "MTT_wood_years")) {
          states_all$MTT_wood_years_parameter_correlation = cor(tmp,states_all$MTT_wood_years)
      }
      # If Mean mean allocation to wood is provided generate a correlation estimate
      if (any(check_list == "alloc_wood_gCm2day")) {
          states_all$NPP_wood_gCm2day_parameter_correlation = cor(tmp,rowMeans(states_all$alloc_wood_gCm2day))
      }
      # If Mean mean allocation to wood is provided generate a correlation estimate
      if (any(check_list == "labile_to_wood_gCm2day")) {
          states_all$NPP_wood_gCm2day_parameter_correlation = cor(tmp,rowMeans(states_all$labile_to_wood_gCm2day))
      }

  } # Both MTT wood and alloc_wood present?

  if (any(check_list == "MTT_som_years") == TRUE) {
      states_all$MTT_som_years_parameter_correlation = cor(tmp,states_all$MTT_som_years)  
      # Assess within pixel correlations with soil turnover
      states_all$MTT_som_years_to_GPP_gCm2day_correlation = cor(states_all$MTT_som_years,ensGPP)
      states_all$MTT_som_years_to_NEE_gCm2day_correlation = cor(states_all$MTT_som_years,ensNEE)
      states_all$MTT_som_years_to_Rauto_gCm2day_correlation = cor(states_all$MTT_som_years,ensRauto)
      states_all$MTT_som_years_to_Rhet_gCm2day_correlation = cor(states_all$MTT_som_years,ensRhet)
      states_all$MTT_som_years_to_som_gCm2_correlation = cor(states_all$MTT_som_years,ensSOM)   
      states_all$MTT_som_years_to_lai_m2m2_correlation = cor(states_all$MTT_som_years,ensLAI)      
      states_all$MTT_som_years_to_dCsom_gCm2_correlation = cor(states_all$MTT_som_years,dCsom)
      # As MTT som exists we must also be able to work out the correlations for change in som over time
      # Note reuse of dCbio
      states_all$dCsom_gCm2_to_gpp_gCm2day_correlation = cor(dCsom,ensGPP) 
      states_all$dCsom_gCm2_to_rauto_gCm2day_correlation = cor(dCsom,ensRauto) 
      states_all$dCsom_gCm2_to_nee_gCm2day_correlation = cor(dCsom,ensNEE)    
      states_all$dCsom_gCm2_to_rhet_gCm2day_correlation = cor(dCsom,ensRhet)    
      states_all$dCsom_gCm2_to_som_gCm2_correlation = cor(dCsom,ensSOM)
      # Special case for the managed grassland model
      if (exists("ensWood")) {
          states_all$MTT_som_years_to_dCwood_gCm2_correlation = cor(states_all$MTT_som_years,dCwood) 
          states_all$MTT_som_years_to_wood_gCm2_correlation = cor(states_all$MTT_som_years,ensWood)   
          states_all$dCsom_gCm2_to_wood_gCm2_correlation = cor(dCsom,ensWood)          
          tmp = rowMeans(states_all$wood_to_litter_gCm2day + states_all$litter_to_som_gCm2day)
          states_all$dCsom_gCm2_to_som_input_gCm2_correlation = cor(dCsom,tmp)
      } else {
          tmp = rowMeans(states_all$litter_to_som_gCm2day)
          states_all$dCsom_gCm2_to_som_input_gCm2_correlation = cor(dCsom,tmp)
      }

  }   

  # Tidy multi-use variables
  rm(ensLAI,ensNBP,ensNEE,ensGPP,ensRauto,ensRhet)

  # Return back to user
  return(states_all)
          
} # end function post_process_dalec
## Use byte compile
post_process_dalec<-cmpfun(post_process_dalec)

###
# Quantify the proportion of ensemble members within the uncertainty bounds of the calibration datasets
###

# This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE)

assess_ensemble_fit_to_calibration_data<-function(states_all,parameters,drivers,PROJECT) {

  ###
  ## Comparison with assimilated observation - to what extent does the ensemble overlap?

  # Check the list variables in states_all which we will be searching
  check_list = names(states_all)


  ## Parameter prior informatoin
  # Initise array
  states_all$priors_assim_data_overlap_fraction = rep(NA, max(PROJECT$model$nopars))
  for (p in seq(1, max(PROJECT$model$nopars))) {
       if (drivers$parpriors[p] != -9999) {
           # Loop through time to assess model overlap with observations
           states_all$priors_assim_data_overlap_fraction[p] = 0
           # Estimate the min / max values for the observations
           obs_max = drivers$parpriors[p] + drivers$parpriorunc[p]
           obs_min = drivers$parpriors[p] - drivers$parpriorunc[p]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = as.vector(parameters[p,,]))
           # Estimate average model ensemble within observated range
           states_all$priors_assim_data_overlap_fraction[p] = ensemble_within_range(hist_list$o,hist_list$m)
           # Average the overlap
       } # was the obs assimilated?
  } # parameter loop

  ## GPP (gC/m2/day)
  obs_id = 1 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$gpp_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$gpp_gCm2day[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$gpp_gCm2day[,tt:t],1,mean))
           }                                            
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction / nobs
      } else {
          states_all$gpp_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## LAI (m2/m2)
  obs_id = 4 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$lai_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$lai_m2m2[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$lai_m2m2[,tt:t],1,mean))
           }
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
           } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction / nobs
      } else {
          states_all$lai_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## NEE (gC/m2/day)
  obs_id = 7 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$nee_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$nee_gCm2day[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$nee_gCm2day[,tt:t],1,mean))
           }
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction / nobs
      } else {
          states_all$nee_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Reco (gC/m2/day)
  obs_id = 13 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$reco_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$reco_gCm2day[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$reco_gCm2day[,tt:t],1,mean))
           }           
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$reco_assim_data_overlap_fraction = states_all$reco_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$reco_assim_data_overlap_fraction = states_all$reco_assim_data_overlap_fraction / nobs
      } else {
          states_all$reco_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Wood (gC/m2)
  obs_id = 19 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  # If there is a prior assign it to the first timestep of the observation timeseries
  if (drivers$parpriors[21] > 0) { drivers$obs[1,obs_id] = drivers$parpriors[21] ; drivers$obs[1,unc_id] = drivers$parpriorunc[21] }
  if (any(drivers$obs[,obs_id] != -9999) && any(check_list == "wood_gCm2")) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$wood_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$wood_gCm2[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$wood_gCm2[,tt:t],1,mean))
           }           
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction / nobs
      } else {
          states_all$wood_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Soil (gC/m2)
  obs_id = 28 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  # If there is a prior assign it to the first timestep of the observation timeseries
  if (drivers$parpriors[23] > 0) { 
      drivers$obs[1,obs_id] = drivers$parpriors[23] 
      drivers$obs[1,unc_id] = drivers$parpriorunc[23] 
      drivers$obs[1,lag_id] = 0
  }
  if (any(drivers$obs[,obs_id] != -9999) && any(check_list == "som_gCm2")) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$soil_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$som_gCm2[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$som_gCm2[,tt:t],1,mean))
           }                      
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction / nobs
      } else {
          states_all$soil_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## fAPAR (0-1)
  obs_id = 34 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$fapar_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           # NOTE 0.5 is assumed fraction of shortwave radiation that is PAR
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$APAR_MJm2day[,t]/(drivers$met[t,4]*0.5))
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$APAR_MJm2day[,tt:t]/(drivers$met[tt:t,4]*0.5),1,mean))
           }                      
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$fapar_assim_data_overlap_fraction = states_all$fapar_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
           } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$fapar_assim_data_overlap_fraction = states_all$fapar_assim_data_overlap_fraction / nobs
      } else {
          states_all$fapar_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## ET (kgH2O/m2/day)
  obs_id = 40 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$et_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$ET_kgH2Om2day[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$ET_kgH2Om2day[,tt:t],1,mean))
           }                                 
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction / nobs
      } else {
          states_all$et_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## NBE (gC/m2/day)
  obs_id = 46 ; unc_id = obs_id+1  ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$nbe_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$nbe_gCm2day[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$nbe_gCm2day[,tt:t],1,mean))
           }                                 
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction / nobs
      } else {
          states_all$nbe_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Fire (gC/m2/day)
  obs_id = 10 ; unc_id = obs_id+1 ; lag_id = unc_id+1
  if (any(drivers$obs[,obs_id] != -9999)) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$fire_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in 1:length(to_do)) {
           # Assign correct time step
           t = to_do[a] ; tt = max(1,t - drivers$obs[t,lag_id])
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           if (t == tt) {
               hist_list = list(o = c(obs_min,obs_max), m = states_all$fire_gCm2day[,t])
           } else {
               hist_list = list(o = c(obs_min,obs_max), m = apply(states_all$fire_gCm2day[,tt:t],1,mean))
           }                                            
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction / nobs
      } else {
          states_all$fire_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  # Return back to user
  return(states_all)

} # end function assess_ensemble_fit_to_calibration_data
## Use byte compile
assess_ensemble_fit_to_calibration_data<-cmpfun(assess_ensemble_fit_to_calibration_data)

