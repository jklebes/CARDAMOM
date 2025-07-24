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
out_ <- .C("C_initialize_stresstest_circle")

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
cardamom_stresstestcirclelikelihood <- function(pars){
    #print("Calling")
    out_ <- 0.0
    ll <- .C("C_stresstest_likelihood", pars, model_npars, out_)[[3]]
    #l <- -pars[3]**2
    #ll <- generateTestDensityMultiNormal(sigma = "no correlation")
}


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

# compare 
bayesianSetup <- createBayesianSetup(likelihood = cardamom_stresstestcirclelikelihood, 
                                     lower = model_parmin,
                                     upper = model_parmax, 
                                     parallel=FALSE
)
out_serial <- runMCMC(bayesianSetup, sampler="DEzs", settings=settings)

summary(out_parallel)
plot(out_parallel)
