
###
## Function to handle running the run_each_site task via slurm job submission
###

# This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE)

# The intension is that this script will be piped into an R process via the slurm schedular.
# Unlike the normal R batch parallel approach the memory management must be done here to allow
# the task to run totally independently and still provide information back to the main job once completed.

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

} else {

    #stop("Incorrect number of command line arguments past into R. There should be two. See run_each_site_slurm.r for details")

}
