
## Alternative cardamom stresstest as R script !
#
# jklebes 2025
#
#
# R parallelism using R bayesianTools DEMcz - Option b
# similar to parallel_AM demo, this versio uses an R-cluster
# which is given parallel function over a matrix containting N
# prpoposed paramter vectors as loglikelihood function.

# very slow version - lots of communication overhead

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

nchains <- 4

cl <- parallel::makeCluster(nchains)
parallel::clusterEvalQ(cl, library(BayesianTools))
# pass the variables on to cluster
parallel::clusterExport(cl, "cardamom_dll" )
parallel::clusterExport(cl, "model_npars" )
parallel::clusterExport(cl, "model_parmax" )
parallel::clusterExport(cl, "model_parmin" )
parallel::clusterExport(cl, "filename" )
# load cardamom library on each cluster node
parallel::clusterEvalQ(cl, dyn.load(cardamom_dll))
parallel::clusterExport(cl, "cardamom_stresstestcirclelikelihood")
# initialize a model on each cluster node
parallel::clusterEvalQ(cl, out_ <- .C("C_initialize_stresstest_circle", filename) )
parallel::clusterExport(cl, "get_initial")
parallel::clusterEvalQ(cl,  initial <- get_initial())


# see also parLapply etc for different-shaped parallelization mappings
# here we run the exact same function over a list of different parameter vectors 
# , because unlike cardamom main model stresstest evaluation is not stateful
plikelihood <- function(paramArr) {
    if (! is.matrix(paramArr)){
        paramArr = t(as.matrix(paramArr))
    }
    parallel::parApply(cl= cl, X=paramArr, MARGIN = 1, FUN = cardamom_stresstestcirclelikelihood)
}



initialMatrix = rbind(get_initial(), get_initial(), get_initial(), get_initial())
bayesianSetup <- createBayesianSetup(likelihood = plikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel='external', # use the cluster
                                     )

# note this will be interpreted as iter/nchains iterations per chain
# while in cardamom-samplers each chain does n_out number of steps
iter = 10000

# only one chain per cluster core
settings = list(iterations = iter, nrChains=1 , message = TRUE , startValue = bayesianSetup$prior$sampler(nchains))
#test
plikelihood(initialMatrix)
plikelihood(bayesianSetup$prior$sampler(nchains))

out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)
plot(out_parallel$chain[,'LL'])

#expect convergence on approx. pi :
plot(out_parallel$chain[,1])
