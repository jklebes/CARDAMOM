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
cardamom_dll = "/home/jklebes/CARDAMOM/build/LIBRARY/CARDAMOM_F/libCARDAMOM.so"
dyn.load(cardamom_dll)


# command line args : infile, outfile, solution_wanted_char, freq_print_char, &
#                   freq_write_char, do_inflate_char, cost_func_scaling_char


#derived variables:  5 of these as int, 
# do_inflate flag 

# some checks on input args 

# ======= INIT ===========  
# random seed - for R library sampler 


# read_pari_data (also get npars and allocates arrays based on npars - shouldnt 
#                 should instead just check npars in data matches npars in model)

# read_options 
# check files for restart
# ~~initialise_mcmc_output~~
# ~~open_output_files~~
# ~~ buffering to prepare output streams~~ ... all to be outsourced to R mcmc data collection

## ----- Pass model to R ----------
out_ <- .C("C_initialize_stresstest_circle")

# TODO keep R wrappers in a different file
#model_name <- .C("getmodelname")
model_npars <- as.integer(10)

out <- rep(as.numeric(0), model_npars)
model_parmin <- .C("C_getstresstestparmin", npars=model_npars, out)[[2]]
print("Fetched model parmin")
print(model_parmin)
model_parmax <- .C("C_getstresstestparmax", npars=model_npars, out)[[2]]
print("Fetched model parmax")
print(model_parmax)


get_initial <- function(){
    initial <- runif(model_npars)*(model_parmax-model_parmin) + model_parmin
}

# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)
#print(ll0)


#wrap that .C function to a more usual R function
cardamom_stresstestcirclelikelihood <- function(pars){
    #print("Calling")
    npars <- as.integer(11)
    out_ <- 0.0
    ll <- .C("C_stresstest_likelihood", pars, model_npars, out_)[[3]]
    #l <- -pars[3]**2
    #ll <- generateTestDensityMultiNormal(sigma = "no correlation")
}


print("initial loglikelihood:")
ll <- cardamom_stresstestcirclelikelihood(initial)
print(ll)

print("Best loglikelihood would be")
answers = c(3.14, 1.2, 3, 8, 10, 15, 200, 193, 88, 291, 1)
print(cardamom_stresstestcirclelikelihood(answers))

# ============= RUN MCMC ===========

print("Running R adaptive MCMC on Stresstest Circle")

ncores <- 4
bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=ncores, #auto parallelism - works with DEzs
                                     # make fortran/C function available on all cores
                                     parallelOptions = list(variables = "all", packages = "all", dlls = cardamom_dll),
                                     )
iter = 10000

settings = list(iterations = 20, startValue=4 , message = TRUE) # hacky: short run to have runMCMC start a cluster
out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings) 
# so that we can initialize cardamom stresstest (i.e. fill DATAin, PI) on each core
parallel::clusterEvalQ(out_parallel[["setup"]][["likelihood"]][["cl"]], out_ <- .C("C_initialize_stresstest_circle") )

settings = list(iterations = iter, startValue=4 , message = TRUE)
out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)


# compare 
bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=FALSE
)
out_serial <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

summary(out_parallel)
plot(out_parallel)