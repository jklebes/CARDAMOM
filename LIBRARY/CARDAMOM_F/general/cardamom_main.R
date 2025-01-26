dyn.load("../../build/LIBRARY/CARDAMOM_F/libcardamom-lib.so")
print(paste("function modellikelihood is loaded:", is.loaded("modellikelihood")))


model_name <- .C("getmodelname")
model_npars <- .C("getmodelnpars")
model_parmin <- .C("getparmin")
model_parmax <- .C("getparmax")

print(paste("detected model name:", model_name))
print(paste("detected model npars:", model_npars))

print("Generating npars random values TODO with bounds in parinfo")
model_npars = 36
initial <- runif(model_npars)

print("initial loglikelihood:")
ll0 <- .C("modellikelihood", initial, model_npars)

print("Running R MCMC")
