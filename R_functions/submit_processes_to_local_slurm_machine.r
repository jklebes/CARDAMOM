
###
## Function to submit processes to eddie
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (UoE).

submit_processes_to_local_slurm_machine<-function (PROJECT_in) {

    print('PREPARING TO SUBMIT MCMC TO LOCAL CLUSTER MACHINE (SLURM scheduler assumed)')
    print('This function should be valid for all CARDAMOM compatible DALEC MCMC functions')

    ## Some housekeeping

    # Remove any previous output files?
    delete_old = readline("Delete any previous output files for this project name?(y/n)")
    if (delete_old == "y") {
        system(paste("rm ",PROJECT_in$resultspath,"/*",sep=""))
    }

    # CARDAMOM typically uses a multi-phase MCMC process. First, an EDC searching phase to
    # to find parameters with an EDC compliant starting point. Second, a pre-mcmc during which
    # the likelihood scores for each data stream are normalised by their sample size. Third, the 
    # main analysis during which the likelihood scores are weighted based on the 'cost_function_scaling'
    # However, if this is an extended run, i.e. going beyond the parameter proposals originally requested
    # the pre-mcmc must be turned off to maintain consistency in the likelihood scores being assessed.
    if (request_extended_mcmc) {
        # In an extended run the number of proposals (samples) is added to by the requested extension
        nsamples = as.integer(PROJECT_in$nsamples) + request_nos_extended_samples
        # Assume pre-mcmc is not to be used if this is an extended run
        pre_mcmc = 0
    } else {
        # In a normal run the number of proposals is passed into the local variable unchanged
        nsamples = as.integer(PROJECT_in$nsamples)
        # Asse the pre-mcmc is used as default
        pre_mcmc = 1
    }

    # Check presence of PROJECT_in$cost_function_scaling
    if (exists(x = "cost_function_scaling", where = PROJECT_in) == FALSE) {
        # If not, assume default cost function
        PROJECT_in$cost_function_scaling = 0
    }

    ## Create the two files needed, one which contains the list of jobs to be ran
    ## and a second which write the correct slurm job submission script

    # create the new file name in the correct location
    outfile = paste(PROJECT_in$exepath,"CARDAMOM_ECDF_EXECUTABLES_LIST.txt",sep="")
    # begin writing out the file contents
    # construct the file now
    first_pass=TRUE
    for (c in seq(1, PROJECT_in$nochains)) {
         for (n in seq(1, PROJECT_in$nosites)) {
              infile = paste(PROJECT_in$datapath,PROJECT_in$name,"_",PROJECT_in$sites[n],".bin",sep="")
              output = paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_",c,"_",sep="")
              if (first_pass) {
                  write(paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                              infile," ",
                              output," ",
                              as.integer(nsamples),
                              " 0 ",
                              as.integer(PROJECT_in$samplerate)," ",
                              as.integer(pre_mcmc)," ",
                              as.integer(PROJECT_in$cost_function_scaling),sep=""),sep=" ", ncolumn=1,file=outfile,append="F")
                  first_pass=FALSE
              } else {
                  write(paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                              infile," ",
                              output," ",
                              as.integer(nsamples),
                              " 0 ",
                              as.integer(PROJECT_in$samplerate)," ",
                              as.integer(pre_mcmc)," ",
                              as.integer(PROJECT_in$cost_function_scaling),sep=""),sep=" ", ncolumn=1,file=outfile,append="T")
              } # first pass or not
         } # chain no
    } # nosite

    # Create the shell script for submitting the job to slurm on the local cluster
    slurm_file = paste(PROJECT_in$exepath,"/slurm_submission.sh",sep="")
    col_sep = "" ; nos_cols = 20
    write(    c("#!/bin/bash"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = FALSE)
    write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("# Slurm directives"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("#SBATCH --account=geos_extra"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("#SBATCH --job-name="), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("#SBATCH --ntasks=1"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("#SBATCH --cpus-per-task=1"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("#SBATCH --mem=1G "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c(paste("#SBATCH --output=",PROJECT_in$oestreampath,"/slurm-%A_%a.out ",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c(paste("#SBATCH --time=",as.numeric(PROJECT_in$chain_runtime),":00:00",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c(paste("#SBATCH --array=1-",PROJECT_in$nosites*PROJECT_in$nochains,sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("# THIS SCRIPT MUST BE ACCOMPANIED BY CARDAMOM_ECDF_EXECUTABLES_LIST.txt IN THE SAME DIRECTORY"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("# arguments are start and end lines!"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c(paste("task=$( cat $1CARDAMOM_ECDF_EXECUTABLES_LIST.txt | sed $SLURM_ARRAY_TASK_ID\\!d )",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
    write(    c("command ${task}"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)

    # Record directory to change back in a moment
    cwd = getwd()
    # Set working directory to the location of the executable we want to run
    setwd(PROJECT_in$exepath)
    # Submit jobs to the local slurm cluster
    system(paste("sbatch ",slurm_file,sep=""))
    # Return back to normal working directory
    setwd(cwd) ; rm(cwd)
    
    # Inform the user
    print("Command issued to local slurm machine")

} # end of function submit_processes_to_local_slurm_machine

## Use byte compile
submit_processes_to_local_slurm_machine<-cmpfun(submit_processes_to_local_slurm_machine)
