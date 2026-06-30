## Alternative cardamom stresstest main as R script 
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
# - We have standard R parallelism, i.e. probably threads, not testing OMP or MPI
# 
# runs stresstest circle, roughly equivalent to running cardamom with args
# 

library(BayesianTools) # if not found install.packages("BayesianTools")
library(assert)
library(here)

source(file.path(here(), "LIBRARY/CARDAMOM_F/general/Rscripts/load_stresstest.R"))


# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)
print("initial loglikelihood:")
ll <- cardamom_stresstestcirclelikelihood(initial)
print(ll)

#print("Best loglikelihood would be")
#answers = c(3.14, 1.2, 3, 8, 10, 15, 200, 193, 88, 291, 1)
#print(cardamom_stresstestcirclelikelihood(answers))

print("Running R adaptive MCMC on Stresstest Circle")


bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax)

# ============= RUN MCMC ===========
# This runs a single-thread adaptive metropolis sampling
# can try different BayesianTools samplers (with no parallelism) here
iter = 100000
settings = list(iterations = iter, startValue=initial , message = TRUE)
out <- runMCMC(bayesianSetup, sampler="AM", settings=settings)

plot(out)
