## Alternative cardamom main loop as R script !
#
# jklebes 2025
#
# The function model_loglikelihood, containing
# all the cardamom model integration and data comparison, 
# is wrapped to an R function vector -> loglikelihood
# which can be called by any standard R sampler .
# One can try MCMC, DEz samplers implemented in R on it (instead of cardamom_samplers)
# and take advantage of their logging, plotting tools 
# (instead of cardamom main loop file writing and output file writing),
# for comparison and isolation testing.
# 
# Disadvantages: 
# - The Roberts and Rosenthal adaptive method with beta is not implemented in R
# - We have standard R parallelism, i.e. probably threads, not OMP or MPI
# 
# I estimate performance should be ok as the model integration is still in fortran;
# top level sampling loop was not a big factor in performance. 

library(BayesianTools) # if not found install.packages("BayesianTools")

# Run the "cmake ..", "make" of cardamom to generate the shared library
dyn.load("/home/jklebes/CARDAMOM/build/LIBRARY/CARDAMOM_F/libcardamom-lib.so")
print(paste("function modellikelihood is loaded:", is.loaded("C_modellikelihood")))


# command line args : infile, outfile, solution_wanted_char, freq_print_char, &
#                   freq_write_char, do_inflate_char, cost_func_scaling_char


#derived variables:  5 of these as int, 
# do_inflate flag 

# some checks on input args 

# ======= INIT ===========  
# random seed - for R library sampler 


# read_pari_data (also get npars and allocates arrays based on npars - shouldnt 
#                 should instead just check npars in data matches npars in model)
out <- .C("C_TMP_initialize")
print(out)

# read_options 
# check files for restart
# ~~initialise_mcmc_output~~
# ~~open_output_files~~
# ~~ buffering to prepare output streams~~ ... all to be outsourced to R mcmc data collection

## ----- Pass model to R ----------

# TODO keep R wrappers in a different file
#model_name <- .C("getmodelname")
out <- as.integer(0)
model_npars <- .C("C_getmodelnpars", out)[[1]]
#print(paste("detected model name:", model_name))
print(paste("detected model npars:", model_npars))

out <- rep(as.numeric(0), model_npars)
model_parmin <- .C("C_getmodelparmin", n=model_npars, out)[[2]]
print("Fetched model parmin")
print(model_parmin)
model_parmax <- .C("C_getmodelparmax", n=model_npars, out)[[2]]
print("Fetched model parmax")
print(model_parmax)

print("Generating npars random values TODO with bounds in parinfo")
initial <- runif(model_npars)
print("initial loglikelihood:")
ll0 <- .C("C_modellikelihood", initial, model_npars)

print("Running R MCMC")

## Check - call modellikelihood function for initial logloikehood
ll <- modellikelihood_R(initial) #error unless everything is propoerly initialized on fortran side

bayesianSetup = createBayesianSetup(modellikelihood_R, parmin, parmax)

# ============= RUN MCMC ===========

iter = 1000
settings = list(iterations=iter)
out <- runMCMC(bayesianSetup, sampler="Metropolis", settings=settings)



