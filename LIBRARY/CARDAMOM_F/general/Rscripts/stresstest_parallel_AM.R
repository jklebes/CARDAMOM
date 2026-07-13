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
library(assert)
library(here)

source(file.path(here(), "LIBRARY/CARDAMOM_F/general/Rscripts/load_stresstest.R"))

# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)
#print(ll0)
print("initial loglikelihood:")
ll <- cardamom_stresstestcirclelikelihood(initial)
print(ll)

nchains <- 4
cl <- parallel::makeCluster(nchains)

# ============= RUN MCMC ===========

print("Running R adaptive MCMC on Stresstest Circle")

bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=F
                                     )

# sampler AM works with this outer-level parallelization only.  We run N
# separate single-chain samplers on N cores of a cluser.
# This is N completely independednt AM runs.  Could be handled by taskarray,
# here parallelism is instead handled by R parallel cluster.  This runs on 
# multiple CPUs .

parallel::clusterEvalQ(cl, library(BayesianTools))
parallel::clusterExport(cl, "cardamom_dll" )
parallel::clusterExport(cl, "model_npars" )
parallel::clusterEvalQ(cl, dyn.load(cardamom_dll))
parallel::clusterExport(cl, "cardamom_stresstestcirclelikelihood")
parallel::clusterEvalQ(cl, out_ <- .C("C_initialize_stresstest_circle") )
iter = 100000

# one chain within each core
settings = list(iterations = iter, nrChains=1, message = TRUE)
#parallel::clusterExport(cl, "bayesianSetup")
#parallel::clusterExport(cl, "settings")

# This parLapply will be useful for when you want to pass chainId X to function:
out <- parallel::parLapply(cl, 1:nchains, function(X, bayesianSetup, settings) runMCMC(
    bayesianSetup, settings, sampler = "AM") , bayesianSetup, settings)


out <- createMcmcSamplerList(out)
plot(out)

