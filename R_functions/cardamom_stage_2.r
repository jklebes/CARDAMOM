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
# Function to carry out stage 2 processes,
# i.e. submitting CARDAMOM jobs to the relevant computer
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

cardamom_stage_2 <-function(PROJECT) {
      print('Welcome to Stage 2 - running CARDAMOM')
      print('The code will be run on cluster or your local machine');

      if (PROJECT$ecdf) {
          # submit files to eddie
          submit_processes_to_cluster(PROJECT)
      } else {
          if (PROJECT$request_use_local_slurm) {
              submit_processes_to_local_slurm_machine(PROJECT)
          } else {
              # submit to local machine
              submit_processes_to_local_machine(PROJECT)
          }
      }
      # report to the user
      return(paste("CARDAMOM Report: 2 completed", sep=""))

} # end function cardamom_stage_2

## Use byte compile
cardamom_stage_2<-cmpfun(cardamom_stage_2)