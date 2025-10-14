## Alternative cardamom stresstest as R script !
#
# jklebes 2025
#
#
# with top level R parallelism.  
# This does not deomnstrate that the model is thread-safe against 
# shared-mmeory parallelism in the same fortran simulation/

library(BayesianTools) # if not found install.packages("BayesianTools")

# Run the "cmake ..", "make" of cardamom to generate the shared library
cardamom_dll = "/home/jklebes/CARDAMOM/build/LIBRARY/CARDAMOM_F/libCARDAMOM.so"
dyn.load(cardamom_dll)


## ----- Pass model to R ----------
filename <- "/home/jklebes/CARDAMOM/test/data/UK_baseline_sites_AliceHolt.bin"
out_ <- .C("C_initialize_model", filename)

# TODO keep R wrappers in a different file
out <- as.integer(0)
model_npars <- .C("C_getmodelnpars", out)[[1]]

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



#wrap that .C function to a more usual R function
cardamom_edc_modellikelihood <- function(pars){
    out_ <- 0.0
    ll <- .C("C_edcmodellikelihood", pars, model_npars, out_, as.integer(1))[[3]]
    #l <- -pars[3]**2
    #ll <- generateTestDensityMultiNormal(sigma = "no correlation")
}


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
    parallel::parApply(cl= cl, X=paramArr, MARGIN = 1, FUN = cardamom_edc_modellikelihood)
}



initialMatrix = rbind(get_initial(), get_initial(), get_initial(), get_initial())
bayesianSetup <- createBayesianSetup(likelihood = plikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel='external', # use the cluster
                                     )
iter = 10000

settings = list(iterations = iter, nrChains=1 , message = TRUE , startValue = bayesianSetup$prior$sampler(nchains))
#test
plikelihood(initialMatrix)
plikelihood(bayesianSetup$prior$sampler(nchains))

out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

plot(out_parallel$chain[,'LL'])

