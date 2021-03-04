
###
## Script to demonstrate running of DALEC-CDEA (aka DALEC2, C1; Bloom & Williams 2015)
###

# To run this example you should only need to modify the working directory below.
# All other commands work based on having the correct starting point.
# NOTE: that this code has been written to work on Linux Systems only

# Example information:
# The example draws data from an eddy covariance site in the central Amazon
# FLUXNET validation site BR-Sa1 (2002-2011); old growth site
# lat (-90/90) = -2.85667 ; long (-180/180) = -54.9588900
# doi: 10.18140/FLX/1440032, (http://sites.fluxdata.org/BR-Sa1/)
#
# CARDAMOM has been used to retrieve location specific parameters for this site.
# Assimilated data are:
# 1 km resolution Copernicus LAI (2001-2017)
# Prior on the initial state of the soil C pool drawn from SoilGrids
# Eddy covariance derived NEE (2002-2011)
# CRUJRP meteorological reanalysis (2001-2017)
#
# Time period: 2001-2017 (inclusive) at monthly time step

# Set working directory to the the location of the example files
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/example_files/DALEC_C1_example_pack")

###
## No need to modify below this line

# Load needed libraries and functions
library(compiler)
source("./src/simulate_all.r")

# Load site specific information needed for run
PROJECT = list(model = list(name = "DALEC_CDEA_LU_FIRES",
                            nopars = 23, nofluxes = 28,
                            nopools = 8),
               start_year = 2001, end_year = 2017,
               parameter_type = "global",
               ctessel_pft = 0,
               exepath = "./")
# Counter is used in mainframework to loop through many sites in a grid but here we just have one example, so hardcode
n = 1

# Load driving data
drivers_input = read.csv("./inputs/drivers.csv", header=TRUE)
# Load example parameter set
parameters_input = read.csv("./inputs/parameters.csv", header=TRUE)
# Load observational constraints used
observations_input = read.csv("./inputs/observations.csv", header=TRUE)

# Arrange drivers into array structure needed by DALEC
met = array(NA, dim=dim(drivers_input))
for (i in seq(1,dim(met)[2])) {
     met[,i] = drivers_input[,i]
}
drivers = list(met = met,
               lat = -2.85667,
               top_sand = 45.9, bot_sand = 50.6,
               top_clay = 45.1, bot_clay = 40.0)
# Arrange parameters into array structure needed by DALEC
parameters = array(NA, dim=c(PROJECT$model$nopars[n],dim(parameters_input)[1]))
for (i in seq(1, PROJECT$model$nopars[n])) {
     parameters[i,] = parameters_input[,i]
}

# Compile DALEC shared object for use
system(paste("gfortran -O2 -shared ./src/model/",PROJECT$model$name,"/src/",PROJECT$model$name,".f90 ",
             "./src/model/",PROJECT$model$name,"/src/",PROJECT$model$name,"_CROP.f90 ",
             "./src/model/",PROJECT$model$name,"/src/",PROJECT$model$name,"_R_interface.f90 ","-o dalec.so -fPIC",sep=""))

# Run example parameter files
print("running model ensemble")
# NOTE: soil_info is not actually used in this example version of DALEC but values are needed for code consistency
soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
states_all = simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,parameters[1:PROJECT$model$nopars[n],],
                          drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                          PROJECT$exepath,soil_info)

# Remoe the '#' from the below command to see what else is stored in the output
#names(states_all)

# Do something with the output...?
par(mfrow=c(2,3))
plot(apply(states_all$lai_m2m2,2,median), ylab="LAI (m2/m2)", xlab = "Month of analysis", col="green", lwd=3, type="l") ; points(observations_input$LAI_m2m2, pch=16)
plot(apply(states_all$nee_gCm2day,2,median), ylab="NEE (gC/m2/d)", xlab = "Month of analysis", col="red", lwd=3, type="l") ; points(observations_input$NEE_gCm2day, pch=16)
plot(apply(states_all$wood_gCm2,2,median), ylab="Wood (gC/m2)", xlab = "Month of analysis", col="brown", lwd=3, type="l")
plot(apply(states_all$gpp_gCm2day,2,median), ylab="GPP (gC/m2/d)", xlab = "Month of analysis", col="blue", lwd=3, type="l")
plot(apply(states_all$fire_gCm2day,2,median), ylab="Fire (gC/m2/d)", xlab = "Month of analysis", col="yellow", lwd=3, type="l")

# Write quantiles to file for above fluxes
write.table(t(apply(states_all$lai_m2m2,2,quantile, prob=c(0.025,0.25,0.50,0.75,0.975))),
            file = "./outputs/LAI_m2m2_TimeQuantile.csv", append = FALSE, row.name=FALSE, col.name=c("2.5%","25%","50%","75%","97.5%"), sep=",")
write.table(t(apply(states_all$nee_gCm2day,2,quantile, prob=c(0.025,0.25,0.50,0.75,0.975))),
            file = "./outputs/NEE_gCm2day_TimeQuantile.csv", append = FALSE, row.name=FALSE, col.name=c("2.5%","25%","50%","75%","97.5%"), sep=",")
write.table(t(apply(states_all$gpp_gCm2day,2,quantile, prob=c(0.025,0.25,0.50,0.75,0.975))),
            file = "./outputs/GPP_gCm2day_TimeQuantile.csv", append = FALSE, row.name=FALSE, col.name=c("2.5%","25%","50%","75%","97.5%"), sep=",")
write.table(t(apply(states_all$fire_gCm2day,2,quantile, prob=c(0.025,0.25,0.50,0.75,0.975))),
            file = "./outputs/FIRE_gCm2day_TimeQuantile.csv", append = FALSE, row.name=FALSE, col.name=c("2.5%","25%","50%","75%","97.5%"), sep=",")
write.table(t(apply(states_all$wood_gCm2,2,quantile, prob=c(0.025,0.25,0.50,0.75,0.975))),
            file = "./outputs/WOOD_gCm2_TimeQuantile.csv", append = FALSE, row.name=FALSE, col.name=c("2.5%","25%","50%","75%","97.5%"), sep=",")
