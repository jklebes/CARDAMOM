
###
## Function to carry out stage 2 processes,
## i.e. submitting CARDAMOM jobs to the relevant computer
###

# Author: T. Luke Smallman (02/05/2024)

cardamom_stage_2 <-function(PROJECT) {
      print('Welcome to Stage 2 - running CARDAMOM')
      print('The code will be run on cluster or your local machine');

      if (PROJECT$ecdf) {
          # submit files to eddie
          submit_processes_to_cluster(PROJECT)
      } else {
          # submit to local machine
          submit_processes_to_local_machine(PROJECT)
      }
      # report to the user
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))

} # end function cardamom_stage_2

## Use byte compile
cardamom_stage_2<-cmpfun(cardamom_stage_2)