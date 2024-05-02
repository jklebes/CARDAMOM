
###
## Function to carry out stage 3 processes,
## i.e. post processing CARDAMOM retrieved parameters
## back through DALEC and storing in R binary file formats
###

# Author: T. Luke Smallman (02/05/2024)

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
           #command = c(command,paste("rm ",PROJECT$eresultspath,"cardamom_outputs.zip",sep=""))
           #command = paste("scp -r ",PROJECT$eresultspath,"* ",username,"@",home_computer,":",PROJECT$resultspath,sep="")
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
       run_mcmc_results(PROJECT,stage,repair,grid_override)
   }

   # now save the project
   save(PROJECT,file=PROJECTfile)

   # report to the user
   return(paste("CARDAMOM Report: 3 completed", sep=""))

} # end function cardamom_stage_3

## Use byte compile
cardamom_stage_3<-cmpfun(cardamom_stage_3)