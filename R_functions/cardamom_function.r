
###
## CARDAMOM function
## from here all other components are called
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom <-function (projname,model,method,stage) {

# Some useful hardcoding optinons
#stage <<- 4 ; repair <<- 1 ; use_parallel <<- FALSE

  ## load needed functions into R environment
  paths = load_paths()

  # Check that the control file has minimum default values and variables created 
  check_control_file_defaults()

  # Use this function to ensure that if the short model name has been provided that we translate
  # this into the full internal code version
  tmp = cardamom_model_details(model,"global",1)
  model = tmp$name

  # define file name for PROJECT file
  # this file will contain all information relating the the PROJECT
  PROJECTfile = paste(paths$cardamom_outputs,model,"_",method,"/",projname,"/infofile.RData",sep="")
  PROJECTtype = paste(model,"_",method,sep="")
  # information to the user
  print(paste("When this all began ",Sys.time(),sep=""))

  if (stage == -1 & model == "ACM") {
      # check the user understand what is about to happen
      failed=TRUE
      while (failed){
        understands=readline("Does the user understand that using the ACM model all datastreams except GPP MUST be set to '  ' and that 'path_to_site_obs' file has all needed information (yes/no)?")
        if (understands != "yes") {failed = TRUE} else {failed = FALSE}
        if (understands == "no") {stop('then you need to read the code to figure out what it assumes....')}
      } # while loop
  } # model ACM

  ###
  ## Begin Stage -1
  ## create or repair PROJECT info file, PROJECT initialiation

  if (file.exists(PROJECTfile) == FALSE | stage == -1){

      # Carry out stage 1 processes
      dummy = cardamom_stage_minus_1(PROJECTfile,PROJECTtype,paths,model,method,projname)
      # report to the user
      return(dummy)

  } else {

      # Inform the user
      print(paste("Loading the infofile = ",PROJECTfile,sep=""))

      # load PROJECT file from binary file
      load(PROJECTfile)
      #if (cardamom_type == "site") {print(PROJECT)}

  } # if (stage == -1)

  ###
  ## Stage 0 is a compile only option. 
  ## This assumes that the PROJECT has already been created

  if (stage == 0) {

      # call compile of the model
      # define PROJECT set up
      PROJECT = cardamom_project_setup(paths,PROJECT)

  } # if (stage == 0)

  ###
  ## Begin Stage 1 - creating CARDAMOM binary input files

  if (stage == 1) {

      # Carry out stage 1 processes
      dummy = cardamom_stage_1(PROJECT)
      # report to the user
      return(dummy)

  } # if stage == 1

  ###
  ## Begin Stage 2

  if (stage == 2) {

      # Carry out stage 2 processes
      dummy = cardamom_stage_2(PROJECT)
      # report to the user
      return(dummy)

  } # stage == 2

  ###
  ## Begin Stage 3 - post processing parameter sets

  if (stage == 3) {

      # Carry out stage 3 processes
      dummy = cardamom_stage_3(PROJECT,PROJECTfile)
      # report to the user
      return(dummy)

  } # stage == 3

  ###
  ## Begin Stage 4 - generating standard figures

  if (stage == 4) {

      # Carry out stage 4 processes
      dummy = cardamom_stage_4(PROJECT)
      # report to the user
      return(dummy)

  } # stage == 4

  ###
  ## Begin Stage 5

  # Currently empty
  if (stage == 5) {

      # Carry out stage 5 processes
      dummy = cardamom_stage_5(PROJECT)
      # report to the user
      return(dummy)

  } # stage == 5

} # end function cardamom

## Use byte compile
cardamom<-cmpfun(cardamom)
