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
# Function to determine convergence of chains
# 
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory). Translation to R and subsequent modifications 
# by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

have_chains_converged<-function (param_sets) {

      # input is order dimensions(npar+1,iter,chain)

      # lets not beat about things lets see if the models have converged
      # need to re-arrange the array to the same as the function needs
      var = array(NA,dim=c(dim(param_sets)[2],dim(param_sets)[1],dim(param_sets)[3]))
      for (i in seq(1, dim(param_sets)[1])) {
	         var[1:dim(param_sets)[2],i,1:dim(param_sets)[3]] = param_sets[i,1:dim(param_sets)[2],1:dim(param_sets)[3]]
      }

      # pass parameters to convergence function
      converged = psrf(var) #; print(converged$R[length(converged$R)])
      # assume critical value of 1.2 (default is 1.1)
      bob = which(converged$R > 1.2)
      # start with all pass assumption
      converged = array("PASS",length(converged$R))
      # changed passes to fails if needed
      converged[bob] = "FAIL"
      # return the result
      return(converged)

} # end function have_chains_converged

## Use byte compile
have_chains_converged<-cmpfun(have_chains_converged)
