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
# Function to handle running the run_each_site task via slurm job submission.
# The intension is that this script will be piped into an R process via the slurm scheduler.
# Unlike the normal R batch parallel approach the memory management must be done here to allow
# the task to run totally independently and still provide information back to the main job once completed.
# 
# Author: T. Luke Smallman (12/11/2024)
# Exceptions states below in specific functions
#
#########################################################################################

# read arguments taken from the command line
args = commandArgs(trailingOnly=TRUE)

# Check that the command line contains the right number of arguments
if (length(args) == 3) {
   
    # First argument should be the project infofile.RData for the current task
    load(args[1])
    # Second should be the repair status of the task
    repair = as.numeric(args[2])
    # Third should be the site number from the project
    n = as.numeric(args[3])

    # The looping structure of the slurm script currently tries to 
    # run the largest number of sites per task. SO in the final task
    # there will be more requests than there are sites.
    # This guards against running sites that do not exist.
    if (n <= PROJECT$nosites) {

        # Load needed libraries for the processes to occur
        library(compiler) # Just-in-time byte compiler for R
        library(zoo)      # Library of rolling average functions

        # Any default conditions we can assume for this process
        use_parallel = TRUE # In this instance TRUE supresses print statements to the screen 
        grid_override = FALSE # Assume that if we are using the cluster to do this, 
                              # we are not crazy enough to want keep full ensemble information

        # Load any user specified functions needed for the job
        source(paste(PROJECT$paths$cardamom,"R_functions/read_binary_file_format.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/ensemble_within_range.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/simulate_all.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/psrf_function.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/have_chains_converged.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/read_parameter_chains.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/read_parameter_covariance.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/run_mcmc_results.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/post_process_dalec.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/post_process_for_grid.r",sep=""))
        source(paste(PROJECT$paths$cardamom,"R_functions/run_each_site.r",sep=""))

        # Run analysis, output is a character object of the output file if it was successful
        # otherwise various negative numerical values are returned.
        site_output_all = run_each_site(n,PROJECT,repair,grid_override)

        # How to write the names to file somewhere for the R code to check against to continue the normal process?
        
    } # guard against trying to run a site that doesn't exist

} else {

    #stop("Incorrect number of command line arguments past into R. There should be two. See run_each_site_slurm.r for details")

}
