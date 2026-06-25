# jklebes 2025
#
# wrappers from the R side, loading functions from compiled CARDAMOM shared libary
# and wrapping them into more standard-shaped R functions

library(here)
library(inline)
#library(Rcpp)

# If run from within CARDAMOM/, here() is CARDAMOM/ , detected based on .git directory
# DLL: if missing run the "cmake ..", "make" of cardamom to generate the shared library
cardamom_dll = file.path(here(),"build/LIBRARY/CARDAMOM_F/libCARDAMOM.so")
dyn.load(cardamom_dll)

# TODO this filename actually has no effect, there's a hardcoded one in fortran initialize_model
# because of difficulty passing strings R to C to fortran
filename <- file.path(here(),"test/data/UK_baseline_sites_AliceHolt.bin")


call_initialize_model <- cfunction(c(bin="SEXP", nbins="SEXP"), '
  C_initialize_model("/home/jklebes/CARDAMOM/test/data/UK_baseline_sites_AliceHolt.bin");
  return nbins;

', includes = '#include "/home/jklebes/CARDAMOM/LIBRARY/CARDAMOM_F/general/Rscripts/libCARDAMOM.h" ',
libargs = c("-L/home/jklebes/CARDAMOM/build/LIBRARY/CARDAMOM_F/libCARDAMOM"),
verbose=TRUE )

print(call_initialize_model(3,9))

#sourceCpp(file.path(here(),"LIBRARY/CARDAMOM_F/general/Rscripts/R_string_fcts.cpp"))
#x <- runif(10)
#print(meanC(x))

initialize_model <- function(filename){
    out_ <- .C("C_initialize_model", filename)
    return()
}

initialize_model(filename)

get_model_npars <-function(){
    out <- as.integer(0)
    model_npars <- .C("C_getmodelnpars", out)[[1]]
    return(model_npars)
}

model_npars <- get_model_npars()

get_model_parmin <- function(){
    #model_npars <- get_model_npars()
    print(model_npars)
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

cardamom_edc_modellikelihood <- function(pars){
    out_ <- 0.0
    ll <- .C("C_edcmodellikelihood", pars, model_npars, out_, as.integer(1))[[3]]
}

cardamom_modellikelihood <- function(pars){
    out_ <- 0.0
    ll <- .C("C_modellikelihood", pars, model_npars, out_, as.integer(1))[[3]]
}

# TODO scaled sqrt variant
