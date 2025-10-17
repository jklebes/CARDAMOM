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
# This script runs samplers from BayesianTools such as adaptive MCMC
# - single thread.  Not samplers requiring concurrency such as DEcz cannot 
# be run from this script.
# see _parallel_AM script for MCMC with multiple chains in parallel
# currently not sure if differential evolution family DEzs samplers
# can be called from R at all.

library(BayesianTools) # if not found install.packages("BayesianTools")
library(assert)
library(here)

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

iter = 10000

settings = list(iterations = iter, nrChains=1, message = TRUE)
out <- runMCMC(bayesianSetup, settings, sampler = "DRAM") 
# DRAM seems to work best from Metropolis family; others have no movement in 
# some parameters


summary(out)

