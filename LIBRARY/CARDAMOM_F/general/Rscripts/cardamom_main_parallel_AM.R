## Alternative cardamom main loop as R script !
#
# for `here()` to find CARDAMOM/ as project root, must be run from somewhere within CARDAMOM/
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
library(here)

nchains <- as.integer(4)

source(file.path(here(), "LIBRARY/CARDAMOM_F/general/Rscripts/load_cardamom_model.R"))

assert(model_npars==32)

# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)
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
parallel::clusterEvalQ(cl,  initial <- get_initial())
iter = 10000

settings = list(iterations = iter, nrChains=1, message = TRUE)
parallel::clusterExport(cl, "bayesianSetup")
parallel::clusterExport(cl, "settings")
# This will be useful for when you want to pass chainId X to function:
out <- parallel::parLapply(cl, 1:nchains, function(X, bayesianSetup, settings) runMCMC(
    bayesianSetup, settings, sampler = "AM") , bayesianSetup, settings)
out <- createMcmcSamplerList(out)
summary(out)
plot(out[[3]][["chain"]][,'LL'])


