
###
## Function to read final parameter covariance matrix from file
###

read_parameter_covariance<- function(PROJECT_in,n) {

  # search for all output files
  cfile = list.files(paste(PROJECT_in$resultspath,sep=""), full.names=TRUE)
  # select the correct project
  is_it = grepl(PROJECT_in$name,cfile) ; cfile = cfile[is_it]
  # select the PARS files only
  is_it = grepl("COV",cfile) ; cfile = cfile[is_it]
  # need to duplicate the list at this point to ensure that we can be certain we do not confuse the chain number and site numbers
  cfile_tmp = gsub(c("_COV"),"",cfile)
  # select the correct site
  is_it = grepl(paste("_",PROJECT_in$sites[n],"_",sep=""),cfile_tmp) ; cfile = cfile[is_it] ; rm(cfile_tmp)

  # just in case
  if (length(cfile) <= 1) {return(parameter_covariance = -9999, info = "file not present")}

  # calculate the number of chains
  chains = seq(1,length(cfile))
  # which site are we on now
  print("Beginning parameter covariance extraction")
  print(paste("Site = ",PROJECT_in$sites[n]," ",n," of ",PROJECT_in$nosites," ",Sys.time(),sep=""))
  # create error flag, initial value zero
  status=array(0,dim=c(length(chains)))

  # define output variable
  parameter_covariance = array(-9999,dim=c(PROJECT_in$model$nopars[n],PROJECT_in$model$nopars[n],length(chains)))
  # loop through each chain
  for (c in seq(1,length(chains))) {
    print(paste("...chain ",c," of ",length(chains),sep=""))
    # open this chains binary file into R, instructing 'r' to read and 'b' for binary
    bob = file(paste(cfile[c],sep=""),'rb') ; nos_var = 1e6
    set1 = readBin(bob, double(),nos_var) ; temp = 0
    if (length(set1) > 0) {
       # keep reading until we have read all that can be read
       while (length(temp) > 0) {
         temp = readBin(bob, double(),nos_var)
         set1 = append(set1,temp)
       }
       # re-arrange into correct structure
       parameter_covariance[,,c] = array(set1,dim=c((PROJECT_in$model$nopars[n]+1),(PROJECT_in$model$nopars[n]+1)))
    }
    # now close this chain
    close(bob) ; set1 = 0 ; rm(set1)

  } # end of chains loop

  # return the parameter solutions
  return(parameter_covariance = parameter_covariance, info = "dim = npar,npar,chain")

} # end of function
