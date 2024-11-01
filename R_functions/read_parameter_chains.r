###
## Function to read parameter chains and determine which will be ran, based on their convergence criterion
###

# This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

determine_parameter_chains_to_run<-function(PROJECT,n) {

  # load only the desired latter fraction of the parameter vectors
  # output is order dimensions(npar+1,iter,chain)
  parameters = read_parameter_chains(PROJECT,n)
  #parameter_covariance = read_parameter_covariance(PROJECT,n)

  # determine whether we have any actual completed chains and whether they include EDC consistent value only
  error_check = FALSE
  if (length(parameters) == 1 & parameters[1] == -9999) {
      error_check = TRUE
      #print("Site not available / parameters file empty")
      dummy = -1 ; return(dummy)
  } else {
      if (length(which(as.vector(is.na(parameters)))) > 0 ) {
          error_check = TRUE
          print("NA found in likelihood score")
          dummy = -2 ; return(dummy)
      } else if (min(as.vector(parameters)) == -Inf) {
          error_check = TRUE
          print("Inf found in likelihood score")
          #print("Site not available / parameters file empty")
          dummy = -3 ; return(dummy)
      } # NaN / Inf check
  } # error check

  # test for convergence and whether or not there is any single chain which can be removed in they do not converge
  notconv = TRUE ; converged = rep("TRUE", times = max(PROJECT$model$nopars))
  while (dim(parameters)[3] > 2 & notconv) {
     if (use_parallel == FALSE) {print("begin convergence checking")}
     converged = have_chains_converged(parameters)
     # if log-likelihood has passed then we are not interested
     if (converged[length(converged)] == "FAIL") {
         #if (use_parallel == FALSE) {print("...not converged begin removing potential parameter vectors")}
         i = 1 ; max_likelihood = rep(NA, length.out=dim(parameters)[3]) ; CI90 = rep(NA,length.out=c(2))
         while (notconv){
            #if (use_parallel == FALSE) {print(paste("......trying removal of parameter vector ",i,sep=""))}
            # Track the maximum likelihood across each chain.
            max_likelihood[i] = max(parameters[dim(parameters)[1],,i])
            converged = have_chains_converged(parameters[,,-i]) ; i = i + 1
            # if removing one of the chains get convergence then great
            if (converged[length(converged)] == "PASS") {
                #if (use_parallel == FALSE) {print(".........convergence found on chain removal")}
                # likelihoods converge now but we need to check for the possibility that the chain we have removed is actually better than the others
                CI90[1] = quantile(parameters[dim(parameters)[1],,(i-1)], prob=c(0.10))
                CI90[2] = quantile(parameters[dim(parameters)[1],,-(i-1)], prob=c(0.90))
                # if the rejected chain is significantly better (at 90 % CI) than the converged chains then we have a problem
                if (CI90[1] > CI90[2]) {
                    # rejected chain (while others converge) is actually better and the others have gotten stuck in a local minima.
                    # we will now assume that we use the single good chain instead...
                    parameters = array(parameters[,,(i-1)],dim=c(dim(parameters)[1:2],2))
                    notconv = FALSE ; i = (i-1) * -1
                    if (use_parallel == FALSE) {print(paste("............chain ",i*-1," only has been accepted",sep=""))}
                } else {
                    # if the non-converged chain is worse or just the same in likelihood terms as the others then we will ditch it
                    notconv = FALSE ; i = i-1 # converged now?
                    parameters = parameters[,,-i]
                    if (use_parallel == FALSE) {print(paste("............chain rejected = ",i,sep=""))}
                }
            }
            # If removing one chain does not lead to convergence then lowest average likelihood chain could be removed
            # NOTE: this should be made optional, such that non-converging locations are excluded from the analysis instead
            if (i > dim(parameters)[3] & notconv) {
                # Which is lowest likelihood
                i = which(max_likelihood == min(max_likelihood))
                # Remove from array
                parameters = parameters[,,-i]
                # Update the maximum likelihood vector also
                max_likelihood = max_likelihood[-i]
                # Update the user
                if (use_parallel == FALSE) {print(paste(".........single chain removal couldn't find convergence; removing lowest likelihood chain = ",i,sep=""))}
                # reset counter
                i = 1
                # If we have removed chains down to 2 (or kept just one) we need to abort
                if (dim(parameters)[3] < 3 & notconv) {
                    notconv = FALSE
                    if (use_parallel == FALSE) {print(".........have removed all low likelihood chains without successful convergence")}
                }
                 }
            } # while to removing chains
     } else {
         #if (use_parallel == FALSE) {print("All chains converge")}
         # we have conveged
         notconv = FALSE
     } # if likelihood not converged
  } # if more than 2 chains

  # Return parameters to user
  return(list(parameters = parameters,converged = converged))

} # end function determine_parameter_chains_to_run
## Use byte compile
determine_parameter_chains_to_run<-cmpfun(determine_parameter_chains_to_run)

###
## Function to read parameter chains info
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

dump_binary_files <-function(infile) {

  # Function will read in the contents of a raw fortran binary file format and
  # dump into a big vector for the user to sort out themselves

  # Check file exists
  if (file.exists(infile)) {

      # open connection to the file
      # 'r' to read and 'b' for binary
      bob = file(infile,'rb') ; nos_var = 1e6
      # Read from the file, double() is the data type expected
      set1 = readBin(bob, double(),nos_var) ; temp = 0
      # keep reading until we have read all that can be read
      while (length(temp) > 0) {
         temp = readBin(bob, double(),nos_var)
         set1 = append(set1,temp)
      }
      # now close this chain
      close(bob)

      # Tidy up
      rm(temp)

      # Back to user
      return(set1)

  } # file exists

} # end function

read_parameter_chains<- function(PROJECT_in,n) {

  # Determine the intended name for the parmeter files
  pfile=paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_",c(1:PROJECT_in$nochains),"_PARS",sep="")

  # Find and remove any files which have no data in them
  is_it = file.size(pfile) ; is_it = which(is_it > 0) ; pfile = pfile[is_it]
  # Return if no files
  if (length(pfile) == 0) {return(-9999)}

  # calculate the number of chains
  chains = seq(1, length(pfile))
  # How many parameter sets to take from the end of the available
  par_vector_length = 100 
  # which site are we on now
  if (use_parallel == FALSE) {
      print("Beginning parameter extraction and chain merge")
      #print(paste("Site = ",PROJECT_in$sites[n]," ",n," of ",PROJECT_in$nosites," ",Sys.time(),sep=""))
  }

  # create error flag, initial value zero
  status = array(0,dim=c(length(chains)))

  # define output variable
  param_sets_out = array(NA,dim=c((PROJECT_in$model$nopars[n]+1),par_vector_length,length(chains)))
  # loop through each chain
  for (c in seq(1,length(chains))) {
       if (use_parallel == FALSE) {print(paste("...chain ",c," of ",length(chains),sep=""))}
       # open this chains binary file into R, instructing 'r' to read and 'b' for binary
       bob = file(paste(pfile[c],sep=""),'rb') ; nos_var = 1e6
       set1 = readBin(bob, double(),nos_var) ; temp = 0
       # keep reading until we have read all that can be read
       while (length(temp) > 0) {
          temp = readBin(bob, double(),nos_var)
          set1 = append(set1,temp)
       }
       # now close this chain
       close(bob)

       # re-arrange into correct structure
       param_sets = array(set1,dim=c((PROJECT_in$model$nopars[n]+1),(length(set1)/(PROJECT_in$model$nopars[n]+1))))
       set1 = 0 ; rm(set1)

       # check for inconsistencies, is there a different number of parameter sets
       # stored than expected - regardless or of more or less
       if (abs(dim(param_sets)[2] - PROJECT_in$nsubsamples) > 0) {
           # Check if more parameter sets are more...
           if (dim(param_sets)[2] > PROJECT_in$nsubsamples) {
               print('*************************************************************')
               print(paste('Warning! Too many parameter vectors in ',pfile[c],sep=""))
               print('*************************************************************')
               print('Likely cause is appended solutions in previously existing file')
               print(paste('To solve this, only the last ',PROJECT_in$nsubsamples,' will be used',sep=""))
               # keep only the end of the parameter sets
               param_sets = param_sets[,((dim(param_sets)[2]-PROJECT_in$nsubsamples):dim(param_sets)[2])]
               status[c] = 2
               # ...or fewer parameter sets than expected
           } else if (dim(param_sets)[2] < PROJECT_in$nsubsamples) {
               print('*************************************************************')
               print(paste('Warning! Missing parameter vectors in ',pfile[c],sep=""))
               print('*************************************************************')
               print('Likely cause is incomplete chain: it may have been stopped by ECDF, or')
               print('MCMC chain got stuck (happens sometimes). This may also result as an')
               print('inconsistency between requested number of samples vs PROJECT$nsubsamples')
               print('CONSIDER DELETING (or re-running) THIS CHAIN!!')
               print('Likely error will occur next!')
               if (dim(param_sets)[2] < par_vector_length) {status[c] = 1}
               # just in case, below we will need to remove this chain from the output
           } else {
               status[c] = 0
           }
       } # mismatch between expected and actual parameter outputs > 1

       # Do not assign if too small
       if (status[c] != 1) {
           # keep a sample form the end of the chain
           param_sets = param_sets[,(((dim(param_sets)[2]-par_vector_length)+1):dim(param_sets)[2])]
           # add these output to the final output variable
           param_sets_out[1:(PROJECT_in$model$nopars[n]+1),,c] = param_sets
       }
       # clean up
       param_sets = 0 ; rm(param_sets)

  } # end of chains loop
    
  # Check status of each chain that has been read in
  filter = which(status == 1)
  if (length(filter) > 0) {
      if (length(filter) == length(chains)) {
          # All chains are being removed as too small, therefore we must return
          # empty vector
          return(-9999)
      }
      # Assuming we have spare chains to work with we can proceed
      param_sets_out = param_sets_out[,,-filter]
      status = status[-filter]
      chains = chains[-filter]
      # Ensure that the shape of the output array is still correct
      param_sets_out = array(param_sets_out, dim=c(dim(param_sets_out)[1:2],length(chains)))
  }

  # Potentially dangerous hack, take modulus of parameters which are nominally 1-365, 
  # but retrieved using broader parameter ranges to aid searching.
  if (PROJECT_in$model$name == "DALEC_1005" || PROJECT_in$model$name == "DALEC_1005a" ||
      PROJECT_in$model$name == "DALEC.C1.D1.F2.P1.#" || PROJECT_in$model$name == "DALEC.A1.C2.D2.F2.H2.P2.R3.#" ||
      PROJECT_in$model$name == "DALEC.A1.C1.D2.F2.H1.P1.#" || PROJECT_in$model$name == "DALEC.A1.C1.D2.F2.H2.P1.#" ||
      PROJECT_in$model$name == "DALEC.A1.C1.D2.F2.H3.P1" ||
      PROJECT_in$model$name == "DALEC.A1.C1.D2.F2.H2.P1.R1.#" || PROJECT_in$model$name == "DALEC.A1.C2.D2.F2.H2.P1.R1.#" ||
      PROJECT_in$model$name == "DALEC.A1.C2.D2.F2.H2.P2.R1.#" || PROJECT_in$model$name == "DALEC.A2.C1.D2.F2.H2.P1.#" ||
      PROJECT_in$model$name == "DALEC.A1.C1.D2.F2.H2.P2.#" || PROJECT_in$model$name == "DALEC.A1.C1.D2.F2.H2.P5.#" ||
      PROJECT_in$model$name == "DALEC.A4.C6.D2.F2.H2.P11.#") {
      param_sets_out[c(12,15),,] = ((param_sets_out[c(12,15),,]-1)%%365.25)+1
  }
  if (PROJECT_in$model$name == "DALEC.C5.D1.F2.P1.#") {
      param_sets_out[c(8,11),,] = ((param_sets_out[c(8,11),,]-1)%%365.25)+1
  }

  # return the parameter solutions
  return(param_sets_out)

} # end of function
## Use byte compile
read_parameter_chains<-cmpfun(read_parameter_chains)
