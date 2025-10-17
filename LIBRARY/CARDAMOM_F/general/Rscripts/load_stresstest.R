# jklebes 2025
#
# wrappers from the R side, loading functions from compiled CARDAMOM shared libary
# and wrapping them into more standard-shaped R functions

library(here)

# If run from within CARDAMOM/, here() is CARDAMOM/ , detected based on .git directory
# DLL: if missing run the "cmake ..", "make" of cardamom to generate the shared library
cardamom_dll = file.path(here(),"build/LIBRARY/CARDAMOM_F/libCARDAMOM.so")
dyn.load(cardamom_dll)

initialize_stresstest <- function(){
    out_ <- .C("C_initialize_stresstest_circle")
    return()
}

initialize_stresstest()

get_model_npars <-function(){
    out <- as.integer(0)
    model_npars <- .C("C_getmodelnpars", out)[[1]]
    return(model_npars)
}

model_npars <- get_model_npars()

get_model_parmin <- function(){
    #model_npars <- get_model_npars()
    out <- rep(as.numeric(0), model_npars)
    model_parmin <- .C("C_getmodelparmin", npars=model_npars, out)[[2]]
}

get_model_parmax <- function(){
    #model_npars <- get_model_npars()
    out <- rep(as.numeric(0), model_npars)
    model_parmax <- .C("C_getmodelparmax", npars=model_npars, out)[[2]]
}

model_parmin <- get_model_parmin()
model_parmax <- get_model_parmax()

get_initial <- function(){
    initial <- runif(model_npars)*(model_parmax-model_parmin) + model_parmin
}

## Likelihood functions

cardamom_stresstestcirclelikelihood <- function(pars){
    out_ <- 0.0
    ll <- .C("C_stresstest_likelihood", pars, model_npars, out_, as.integer(1))[[3]]
}
