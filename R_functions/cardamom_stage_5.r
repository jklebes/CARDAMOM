
###
## Function to carry out stage 5 processes,
## i.e. dump of all outputs to netcdf files
###

# Author: T. Luke Smallman (02/05/2024)

cardamom_stage_5<-function(PROJECT) {

   # Inform the user
   print("Stage 5 write a netcdf dump of the CARDAMOM output")

   if (PROJECT$spatial_type == "grid") {
       # Create the netcdf file
       create_grid_output_nc(PROJECT)
   } else if (PROJECT$spatial_type == "site") {
       create_states_all_nc(PROJECT)
   } else {
       print("PROJECT$spatial_type does not have valid value")
   }

   # report to the user
   return(paste("CARDAMOM Report: 5 completed", sep=""))

} # end function cardamom_stage_5

## Use byte compile
cardamom_stage_5<-cmpfun(cardamom_stage_5)