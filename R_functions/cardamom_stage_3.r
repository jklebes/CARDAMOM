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
# Function to carry out stage 3 processes,
# i.e. post processing CARDAMOM retrieved parameters
# back through DALEC and storing in R binary file formats
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

cardamom_stage_3 <-function(PROJECT,PROJECTfile) {

   print("Stage 3 will copy files back from cluster and begin postprocessing")
   print("NOTE: this will only be effective if cluster has completed its tasks")

   if (PROJECT$ecdf) {
       failed = TRUE
       while(failed) {
          # do we copy back the files?
          copy_back = readline("Copy results back from cluster? (y/n)")
          if (copy_back != "y" & copy_back != "n") { failed=TRUE } else { failed=FALSE }
       }
       # Are we copying back the files
       if (copy_back == "y") {
           #home_computer=Sys.info()["nodename"]
           # If yes, then we mist delete the existing files to ensure we do not mix analysis versions
           if (length(list.files(paste(PROJECT$resultspath,"/",sep=""))) > 0) {
               system(paste("find ",PROJECT$resultspath," -type f -name '*' -delete",sep=""))
           }
           # Ensure any existing zip directory has been deleted before creating a new one.
           command = paste("rm ",PROJECT$eresultspath,"cardamom_outputs*.zip",sep="")
           # Prepare copy back command
           # Compress all existing files into zip directory
           # There is a limit on how many files (based on the command length) that can be added at once using zip alone.
           # However, we can get around this by listing all files using find and then piping these into zip
           for (i in seq(1,PROJECT$nochains)) {
                command = c(command,paste("zip -j -r -q ",PROJECT$eresultspath,"cardamom_outputs_",i,".zip ",PROJECT$eresultspath," -i '*_",i,"_PARS'",sep=""))
           }
           command = c(command,paste("scp -r -q ",PROJECT$eresultspath,"cardamom_outputs*.zip ",username,"@",home_computer,":",PROJECT$resultspath,sep=""))
           # Execute on remote server
           ecdf_execute(command,PROJECT$paths$cardamom_cluster)
       } # copy back
       # Locally check that we have the cardamom_outputs.zip copied back from remote server.
       # Assuming they are present unzip and delete the zip directory
      for (i in seq(1,PROJECT$nochains)) {
           if (file.exists(paste(PROJECT$resultspath,"cardamom_outputs_",i,".zip",sep=""))) {
               # Unzip
               system(paste("unzip -qq -o ",PROJECT$resultspath,"cardamom_outputs_",i,".zip -d ",PROJECT$resultspath, sep=""))
               # Delete file now
               system(paste("rm ",PROJECT$resultspath,"cardamom_outputs_",i,".zip", sep=""))
           }
      }
      # Check that we have any output files
      tmp = list.files(PROJECT$resultspath, pattern="PARS")
      if (length(tmp) == 0) {return(paste("No files can be found in the RESULTS directory")) }
   } # ecdf condition
   # do we run the parameters yet for analysis
   # Changed to a hardcoded run of the analysis
   run_all = "y"#readline("Run all parameter vectors to generate confidence intervals? (y/n)")
   failed = TRUE
   while(failed) {
      if (run_all != "y" & run_all != "n") {run_all = readline("Run all parameter vectors to generate confidence intervals? (y/n)") ; failed=TRUE} else {failed = FALSE}
   }
   # If we are running
   if (run_all == "y") {
       # Specify how much of the storged parameter set to run
       # NOTE: this is deprecated in favour of hardcoded last 100 of each parameter file.
       # This is more resilient to running incomplete chains
       PROJECT$latter_sample_frac = 0.75 #0.5 # 0.75 #readline("What (latter) fraction of accepted parameters to use (e.g. 0.5)?")
       # Run the parameter back through DALEC
       run_mcmc_results(PROJECT,repair,grid_override)
   }

   # now save the project
   save(PROJECT,file=PROJECTfile)

   # report to the user
   return(paste("CARDAMOM Report: 3 completed", sep=""))

} # end function cardamom_stage_3

## Use byte compile
cardamom_stage_3<-cmpfun(cardamom_stage_3)