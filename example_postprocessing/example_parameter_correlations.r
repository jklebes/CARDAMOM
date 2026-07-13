
###
## Example creating parameter correlations
###

# Load the CARDAMOM analysis you want to work with
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A3.C1.D2.F2.H2.P1.#_MHMCMC/FI-Hyy_example/RESULTS_PROCESSED/FI-Hyy.RData")

# Extract SoilGrids v2 input values
drivers$top_sand ; drivers$bot_sand
drivers$top_clay ; drivers$bot_clay

# Determine the total number of timesteps
finish = dim(states_all$gpp_gCm2day)[2]
start = finish-12+1

# Estimate the ensemble of annual GPPs
annualGPP = apply(states_all$gpp_gCm2day[,start:finish],1,mean)

# Correlate with parameters
# p = parameter counter
for (p in seq(1, dim(parameters)[1])) {
     print(paste("p = ",p," cor = ",round(cor(annualGPP,as.vector(parameters[p,,])), digits = 4), sep=""))
}
