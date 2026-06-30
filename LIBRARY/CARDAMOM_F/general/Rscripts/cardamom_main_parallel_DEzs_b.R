## Alternative cardamom main loop as R script !
#
# jklebes 2025
#
#
# with top level R parallelism.  
# This does not deomnstrate that the model is thread-safe against 
# shared-mmeory parallelism in the same fortran simulation/

library(BayesianTools) # if not found install.packages("BayesianTools")
library(assert)
library(here)

# Run the "cmake ..", "make" of cardamom to generate the shared library
cardamom_dll = file.path(here(),"build/LIBRARY/CARDAMOM_F/libCARDAMOM.so")
dyn.load(cardamom_dll)

source(file.path(here(), "LIBRARY/CARDAMOM_F/general/Rscripts/load_cardamom_model.R"))

assert(model_npars==32)

# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)
#print(ll0)
print("initial loglikelihood:")
ll <- cardamom_edc_modellikelihood(initial)
print(ll)


# ============= RUN MCMC ===========

print("Running R adaptive MCMC on Stresstest Circle")

nchains <- 4

cl <- parallel::makeCluster(nchains)
parallel::clusterEvalQ(cl, library(BayesianTools))
parallel::clusterExport(cl, "cardamom_dll" )
parallel::clusterExport(cl, "model_npars" )
parallel::clusterExport(cl, "model_parmax" )
parallel::clusterExport(cl, "model_parmin" )
parallel::clusterExport(cl, "filename" )
parallel::clusterEvalQ(cl, dyn.load(cardamom_dll))
parallel::clusterExport(cl, "cardamom_edc_modellikelihood")
parallel::clusterEvalQ(cl, out_ <- .C("C_initialize_model", filename) )
parallel::clusterExport(cl, "get_initial")
parallel::clusterEvalQ(cl,  initial <- get_initial())


# see also parLapply etc for different-shaped parallelization mappings
# here we run the exact same function over a list of different parameter vectors 
# , because unlike cardamom main model stresstest evaluation is not stateful
plikelihood <- function(paramArr) {
    if (! is.matrix(paramArr)){
        paramArr = t(as.matrix(paramArr))
    }
    parallel::parApply(cl= cl, X=paramArr, MARGIN = 1, FUN = cardamom_modellikelihood)
}



initialMatrix = rbind(get_initial(), get_initial(), get_initial(), get_initial())
bayesianSetup <- createBayesianSetup(likelihood = plikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel='external', # use the cluster
                                     )
iter = 10000


settings = list(iterations = iter, nrChains=1, message = TRUE)
parallel::clusterExport(cl, "bayesianSetup")
parallel::clusterExport(cl, "settings")
# This will be useful for when you want to pass chainId X to function:
out <- parallel::parLapply(cl, 1:nchains, function(X, bayesianSetup, settings) runMCMC(
  bayesianSetup, settings, sampler = "AM") , bayesianSetup, settings)
out <- createMcmcSamplerList(out)

# make sure all chains got to EDC_loglikelihood 0 !

settings = list(iterations = iter, nrChains=1 , message = TRUE , 
                startValue = t ( sapply(c(1:nchains), function(x) (out[[x]]$current)) ))
out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

plot(out_parallel$chain[,'LL'])

