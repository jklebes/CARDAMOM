## Alternative cardamom stresstest as R script !
#
# jklebes 2025
#
#
# this version demonstrates simplest automatic R parallelism using R bayesianTools DEMcz and 
# CARDAMOM stresstest, which unlike cardamom main model is not stateful.  For
# Cardamom main model different parallelism will be required.

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


# ============= RUN MCMC ===========

print("Running R adaptive MCMC on Stresstest Circle")

ncores <- 4
bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=ncores, #auto parallelism - works with DEzs
                                     # make same fortran/C function available on all cores
                                     parallelOptions = list(variables = "all", packages = "all", dlls = cardamom_dll),
                                     )
# with default parallelism the same likelihood function will be used on each chain
# i.e. this won't work for main CARDAMOM DALEC simulation, where models are stateful
# and need to look up the correct arrays via thread_id 


iter = 5000

settings = list(iterations = 20, startValue=4 , message = TRUE) # hacky: short run to have runMCMC start a cluster
out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings) 

# so that we can initialize cardamom stresstest (i.e. fill DATAin, PI) on each core
parallel::clusterEvalQ(out_parallel[["setup"]][["likelihood"]][["cl"]], out_ <- .C("C_initialize_stresstest_circle") )

settings = list(iterations = iter, startValue=4, message = TRUE)
out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

# Result : 60x slower with 4 cores than with parallel=FALSE.  
# There is syncronizing and information sharing between chains
# Not faster because the circle stresstest likleihood function is trivial

# For comparison - serial DEzs 
bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=FALSE
)
out_serial <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

summary(out_parallel)
plot(out_parallel)
