
###
## Function to read the parameter prior distribution definition source code files.
## The function then outputs a model specific vector of maximum and minimum parameter values
## for use in determining the posterior:prior uncertainty reduction
##
## Created by T. L. Smallman (t.l.smallman@ed.ac.uk, UoE)
## First version: 03/09/2021
## Last modification: 06/09/2021
###

read_src_model_priors<- function(PROJECT) {

  # Note this code has been created to only work with the default DALEC versions,
  # i.e. it will not return the crop model parameters

  # Define the prior max / min bounds
  parmin = rep(NA, max(PROJECT$model$nopars))
  parmax = rep(NA, max(PROJECT$model$nopars))
  p_index = 1:max(PROJECT$model$nopars)

  # Determine file path for the model specific parameter file
  src_par_file = paste(PROJECT$paths$cardamom_library,"/CARDAMOM_F/model/",PROJECT$model$name,"/src/",PROJECT$model$name,"_PARS.f90",sep="")
  # Open connection to the file
  src_par = file(src_par_file, open="r")
  # Begin reading the source code for the specific model
  keep_going = TRUE ; mn_done = 1 ; mx_done = 1 ; n = 0
  while (keep_going) {
       # Read the current line
       cur_line = readLines(src_par, n=1) ; n = n + 1
       # Check the current line is valid
       if (length(cur_line) > 0) {
           # Check this is not a comment line, i.e. is the first character a '!'
           if (grepl("^!", cur_line) == FALSE) {
              # This isn't a comment line, ok so now check to see if this is line contains parameter bounds information
              # Check whether min parameter bound can be found
              if (grepl(pattern = paste("parmin(*)",sep=""), x = cur_line, fixed=FALSE)) {
                  # Determine what is the parmeter number we are dealing with
                  mn_search = FALSE ; mn = 0
                  while (mn_search == FALSE) {
                     mn = mn + 1
                     mn_search = grepl(pattern = paste("parmin(",p_index[mn],")",sep=""), x = cur_line, fixed=TRUE)
                  }
                  # Now extract the information
                  tmp = unlist(strsplit(cur_line,"="))
                  tmp = tmp[length(tmp)]
                  # Swap double precision definition (fortran) to R compatible
                  tmp = gsub("d","e",tmp)
                  # Check for further comments
                  if (grepl(pattern = "!", x = tmp)) {
                      tmp = unlist(strsplit(tmp,"!"))
                      tmp = tmp[1]
                  } # are their further comments to remove
                  parmin[mn] = eval(parse(text = tmp))
                  # Increment counter
                  mn_done = mn_done + 1
              } # is parmin?
              if (grepl(pattern = paste("parmax(*)",sep=""), x = cur_line, fixed=FALSE)) {
                  # Determine what is the parmeter number we are dealing with
                  mx_search = FALSE ; mx = 0
                  while (mx_search == FALSE) {
                     mx = mx + 1
                     mx_search = grepl(pattern = paste("parmax(",p_index[mx],")",sep=""), x = cur_line, fixed=TRUE)
                  }
                  # Now extract the information
                  tmp = unlist(strsplit(cur_line,"="))
                  tmp = tmp[length(tmp)]
                  # Swap double precision definition (fortran) to R compatible
                  tmp = gsub("d","e",tmp)
                  # Check for further comments
                  if (grepl(pattern = "!", x = tmp)) {
                      tmp = unlist(strsplit(tmp,"!"))
                      tmp = tmp[1]
                  } # are their further comments to remove
                  parmax[mx] = eval(parse(text = tmp))
                  # Increment counter
                  mx_done = mx_done + 1
             } # is parmin?

           } # is this a commented line?
       } else {
           # None valid line read, assume end of file, which is a problem
           close(src_par) ; stop(paste("read_src_model_priors: Error reading ",src_par_file,sep=""))
       } # Is the current line valid

       # Check whether we have got all our parameters
       if (mx_done > length(parmax) & mn_done > length(parmin)) {
           # Stop searching
           keep_going = FALSE
       } # Got all my parameters

  } # keep_going == TRUE

  # Close the parameter file
  close(src_par) ; gc()

  # Return to the user
  return(list(parmin = parmin, parmax = parmax))

} # end function read_src_model_priors
