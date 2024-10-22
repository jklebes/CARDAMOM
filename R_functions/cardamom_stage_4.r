
###
## Function to carry out stage 4 processes,
## i.e. creating generic figures of the analysis
###

# Author: T. Luke Smallman (02/05/2024)

cardamom_stage_4<-function(PROJECT) {

      print("Beginning stage 4: generating stardard outputs")

      # Generating site level plots or gridded
      if (PROJECT$spatial_type == "site" | grid_override) {
          # Generate figures of parameters and model values including
          # uncertainty information
          generate_uncertainty_figures(PROJECT)
      } else if (PROJECT$spatial_type == "grid") {
          # will generate spatial maps instead
          generate_parameter_maps(PROJECT)
          generate_simplified_stock_and_flux_maps(PROJECT)
      } else {
          stop('missing spatial_type definition (i.e. grid or site)')
      } # grid or site run

      # report to the user
      return(paste("CARDAMOM Report: 4 completed", sep=""))

} # end function cardamom_stage_4

## Use byte compile
cardamom_stage_4<-cmpfun(cardamom_stage_4)