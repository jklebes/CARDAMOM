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

# if not present, run the "cmake ..", "make" of cardamom to generate the shared library
dyn.load("/home/jklebes/CARDAMOM/build/LIBRARY/CARDAMOM_F/libCARDAMOM.so")


## ----- Pass model to R ----------

# TODO wrap these functions from R side once 
# Prompt initialization, i.e. allocation of arrays for stresstest
out_ <- .C("C_initialize_stresstest_circle")

# Fetch variables from fortran internals to R : npars, PI%parmin, PI%parmax
out <- as.integer(0)
model_npars <- .C("C_getmodelnpars", out)[[1]]
assert(model_npars == 10)
out <- rep(as.numeric(0), model_npars)
model_parmin <- .C("C_getmodelparmin", npars=model_npars, out)[[2]]
print("Fetched model parmin")
print(model_parmin)
model_parmax <- .C("C_getmodelparmax", npars=model_npars, out)[[2]]
print("Fetched model parmax")
print(model_parmax)

# a function to generate initial values in uniform distribution in range
get_initial <- function(){
    initial <- runif(model_npars)*(model_parmax-model_parmin) + model_parmin
}

#wrap that .C loglikelihood function to a more usual-shaped R function
cardamom_stresstestcirclelikelihood <- function(pars){
    #print("Calling")
    npars <- as.integer(11)
    out_ <- 0.0
    ll <- .C("C_stresstest_likelihood", pars, model_npars, out_)[[3]]
}

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
# can try different BayesianTools samplers (with no parallelism) here
iter = 100000
settings = list(iterations = iter, startValue=initial , message = TRUE)
out <- runMCMC(bayesianSetup, sampler="AM", settings=settings)

plot(out)
