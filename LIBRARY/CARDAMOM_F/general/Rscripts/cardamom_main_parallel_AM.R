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
#   ; this is not exactly the same as cardamom-samplers adaptive MCMC
#
# This script runs parallel adaptive MCMC with parallelism implemented 
# at R script level rather than by fortran and ompenmp .  
# It demonstrates that the cardamom dalec model loglikelihood calculation 
# can be queried in parallel, i.e. is thread safe,
# and what the performance is like with a standard sampler.


library(BayesianTools) # if not found install.packages("BayesianTools")
library(assert)

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

nchains <- as.integer(4)
# TODO take from command line arg of Rscript
# changing filename here currently has no effect :(
filename <- "/home/jklebes/CARDAMOM/test/data/UK_baseline_sites_AliceHolt.bin"
print("File exists:")
print(file.exists(filename))

## ----- Pass model to R ----------
# this function triggers several initialization steps:
# load PI , DATAin info from model file, input file into the 
# type objects.  Plus initialize model_working_variables type objects x nchains 
# in the model module.
# WARNING filename ignored currently, hard-coded inside
 out_ <- .C("C_initialize_model", nchains)

# TODO keep R wrappers in a different file
#model_name <- .C("getmodelname")
out <- as.integer(0)
model_npars <- .C("C_getmodelnpars", out)[[1]]
assert(model_npars==32)

out <- rep(as.numeric(0), model_npars)
model_parmin <- .C("C_getmodelparmin", npars=model_npars, out)[[2]]
print("Fetched model parmin")
print(model_parmin)
model_parmax <- .C("C_getmodelparmax", npars=model_npars, out)[[2]]
print("Fetched model parmax")
print(model_parmax)


get_initial <- function(){
    initial <- runif(model_npars)*(model_parmax-model_parmin) + model_parmin
}

# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)


#wrap that .C function to a more usual R function
cardamom_edc_modellikelihood <- function(pars){
    out_ <- 0.0
    ll <- .C("C_edcmodellikelihood", pars, model_npars, out_, as.integer(1))[[3]]
}


print("initial loglikelihood:")
ll <- cardamom_edc_modellikelihood(initial)
print(ll)

cl <- parallel::makeCluster(nchains)
#p_likelihood <- function(param) parallel::parApply(cl = cl, X = param, MARGIN = 1, FUN = cardamom_stresstestcirclelikelihood)

# ============= RUN MCMC ===========

print("Running R adaptive MCMC on Stresstest Circle")

bayesianSetup <- createBayesianSetup(likelihood = cardamom_edc_modellikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=F
                                     )

# sampler AM works with this outer-level parallelization only.  We run N
# separate single-chain samplers on N cores.

parallel::clusterEvalQ(cl, library(BayesianTools))
parallel::clusterExport(cl, "cardamom_dll" )
parallel::clusterExport(cl, "model_npars" )
parallel::clusterExport(cl, "model_parmax" )
parallel::clusterExport(cl, "model_parmin" )
parallel::clusterExport(cl, "filename" )
parallel::clusterEvalQ(cl, dyn.load(cardamom_dll))
parallel::clusterExport(cl, "cardamom_edc_modellikelihood")
parallel::clusterEvalQ(cl, out_ <- .C("C_initialize_model") )
parallel::clusterExport(cl, "get_initial")
parallel::clusterEvalQ(cl,  intial <- get_initial())
iter = 1000000

settings = list(iterations = iter, nrChains=1, message = TRUE)
parallel::clusterExport(cl, "bayesianSetup")
parallel::clusterExport(cl, "settings")
# This will be useful for when you want to pass chainId X to function:
out <- parallel::parLapply(cl, 1:nchains, function(X, bayesianSetup, settings) runMCMC(
    bayesianSetup, settings, sampler = "DRAM") , bayesianSetup, settings)


out <- createMcmcSamplerList(out)
summary(out)
plot(out[[1]][["chain"]][999000:999999,'LL'])
