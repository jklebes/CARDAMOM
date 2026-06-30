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
# Function to direct the creation of site level timeseries plots of CARDAMOM analyes
# 
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory). Translation to R and subsequent 
# modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

single_site_plotting_control<-function(n,PROJECT) {

   # Function deals with the control of site level plotting of parameters and
   # model stock / fluxes estimates

   # find relevant parameter information first
   # output is order dimensions(npar+1,iter,chain)
   parameters = read_parameter_chains(PROJECT,n)
   # If an analysis has been carried out for this location (parameters[1] != -9999)
   if (parameters[1] != -9999) {
       # Determine whether chains have converged (true/false)
       converged = have_chains_converged(parameters)
       plot_parameters(PROJECT,parameters,converged,n)
       # generate file name of the output file created in stage 3
       loadfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")
       if (file.exists(loadfile)) {
           # model state and flux plotting with uncertainty
           uncertainty_figures(n,PROJECT,loadfile)
       }
   } # parameters[1] != -9999

} # end function plot_each_site

generate_uncertainty_figures<-function(PROJECT) {

   # Function is responsible for using either serial or parallel
   # approaches to generating site specific plots of retrieved parameter
   # ensembles and time series information on ecosystem states and fluxes

   # Request the creation of the plots
   if (use_parallel & PROJECT$nosites > 1) {
       cl <- makeCluster(min(PROJECT$nosites,numWorkers), type = "PSOCK")
       # load R libraries in cluster
       clusterExport(cl,c("load_r_libraries","rmse","have_chains_converged",
                          "read_parameter_chains","plotconfidence","psrf",
                          "uncertainty_figures","plot_parameters",
                          "use_parallel"))
       clusterEvalQ(cl, load_r_libraries())
       dummy = parLapply(cl,1:PROJECT$nosites,fun=single_site_plotting_control,PROJECT=PROJECT)
       stopCluster(cl)
   } else {
       # or use serial
       dummy = lapply(1:PROJECT$nosites,FUN=single_site_plotting_control,PROJECT=PROJECT)
   } # parallel option

} # end function generate_uncertainty_figures
## Use byte compile
generate_uncertainty_figures<-cmpfun(generate_uncertainty_figures)
