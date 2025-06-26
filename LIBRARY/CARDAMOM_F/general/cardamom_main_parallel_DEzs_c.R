## Alternative cardamom stresstest as R script !
#
# jklebes 2025
#
#
# this version demonstrates simplest automatic R parallelism using R bayesianTools DEMcz and 
# CARDAMOM stresstest, which unlike cardamom main model is not stateful.  For
# Cardamom main model different parallelism will be required.

library(BayesianTools) # if not found install.packages("BayesianTools")

# Run the "cmake ..", "make" of cardamom to generate the shared library
cardamom_dll = "/home/jklebes/CARDAMOM/build/LIBRARY/CARDAMOM_F/libCARDAMOM.so"
dyn.load(cardamom_dll)


## ----- Pass model to R ----------
nchains <- 4
out_ <- .C("C_initialize_model", nchains)

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
cardamom_edc_modellikelihood <- function(pars, chain_id){
    #print("Calling")
    out_ <- 0.0
    chain_id = 1
    ll <- .C("C_edcmodellikelihood", pars, model_npars, out_, as.integer(chain_id))[[3]]
}


# call modellikelihood once as a test
print("random initial:")
initial <- get_initial()
print(initial)
#print(ll0)
print("initial loglikelihood:")
ll <- cardamom_edc_modellikelihood(initial, 1)
print(ll)
ll <- cardamom_edc_modellikelihood(initial, 4)
print(ll)


# ============= RUN MCMC ===========

print("Running R adaptive MCMC on Stresstest Circle")


# see also parLapply etc for different-shaped parallelization mappings
# here we run the exact same function over a list of different parameter vectors 
# , because unlike cardamom main model stresstest evaluation is not stateful
plikelihood <- function(paramArr) {
    if (! is.matrix(paramArr)){
        paramArr = t(as.matrix(paramArr))
    }
    print(dim(paramArr))
    apply(X=paramArr, MARGIN = 1, FUN = cardamom_edc_modellikelihood)
}



initialMatrix = rbind(get_initial(), get_initial(), get_initial(), get_initial())
bayesianSetup <- createBayesianSetup(likelihood = plikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel='external', # use the cluster
                                     )
iter = 1000

settings = list(iterations = iter, nrChains=1 , message = TRUE , startValue = bayesianSetup$prior$sampler(nchains))
#test
plikelihood(initialMatrix)
plikelihood(bayesianSetup$prior$sampler(nchains))

out_parallel <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

summary(out_parallel)
plot(out_parallel)

