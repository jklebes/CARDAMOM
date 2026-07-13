
# Read and compare the likelihood score plus other traits(?) between the default and alternate model structure

load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/AU-How_default/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> default_states_all ; parameters->default_parameters
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/AU-How_alt3/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> alt1_states_all ; parameters->alt1_parameters

load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/FI-Hyy_default/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> default_states_all ; parameters->default_parameters
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/FI-Hyy_alt3/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> alt1_states_all ; parameters->alt1_parameters

load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/FR-LBr_default/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> default_states_all ; parameters->default_parameters
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/FR-LBr_alt3/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> alt1_states_all ; parameters->alt1_parameters

load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/FR-Pue_default/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> default_states_all ; parameters->default_parameters
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/FR-Pue_alt3/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> alt1_states_all ; parameters->alt1_parameters

load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/GF-Guy_default/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> default_states_all ; parameters->default_parameters
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/GF-Guy_alt3/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> alt1_states_all ; parameters->alt1_parameters

load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/US-Ha1_default/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> default_states_all ; parameters->default_parameters
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_GSI_BUCKET_MHMCMC/US-Ha1_alt3/infofile.RData")
load(paste(PROJECT$results_processed,"/",PROJECT$site[1],".RData", sep=""))
states_all-> alt1_states_all ; parameters->alt1_parameters

# Tidy
rm(states_all,parameters)

summary(as.vector(default_parameters[dim(default_parameters)[1],,]))
summary(as.vector(alt1_parameters[dim(alt1_parameters)[1],,]))

# Then create new colours with high 'alpha', i.e. transparency
c_blue =  rgb(0,0,255, max = 255, alpha = 80, names = "lt.blue")
c_green = rgb(0,255,0, max = 255, alpha = 80, names = "lt.pink")
c_black = rgb(0,0,0,   max = 255, alpha = 80, names = "lt.black")
c_yellow = rgb(255,255,0, max = 255, alpha = 80, names = "lt.yellow")

# Overlay histograms

par(mfrow=c(7,7), mar = c(3,3,3,1))
for (p in seq(1, dim(alt1_parameters)[1])) {
     # Set to local variables
     default = as.vector(default_parameters[p,,])
     alternate = as.vector(alt1_parameters[p,,])
     # Determine the x axis range and breakpoints
     b <- min(c(default,alternate), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(default,alternate), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.1) ; e = e + abs(mean(b,e)*0.1) # add a buffer
     ax <- pretty(c(b,e), n = 12) # Make a neat vector for the breakpoints
     # Plot the seperate histograms and store them in an object, do not save them yet
     hgA <- hist(default, breaks = ax, plot = FALSE) # Save first histogram data
     hgB <- hist(alternate, breaks = ax, plot = FALSE) # Save 2nd histogram data
     # Now plot them together
     ymax = max(c(hgA$counts,hgB$counts))
     plot(hgA, col = c_blue, main=paste("Parameter = ",p,sep=""), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with FCP
     plot(hgB, col = c_green, add = TRUE) # Add Alternate
     if (p == 1) {legend("topright",legend = c("Default", "Alternate"), col = c(c_blue,c_green), pch=16, bty = "n")}
}

# RaGPP
default = as.vector(apply(default_states_all$rauto_gCm2day/default_states_all$gpp_gCm2day,1,mean))
alternate = as.vector(apply(alt1_states_all$rauto_gCm2day/alt1_states_all$gpp_gCm2day,1,mean))
# Determine the x axis range and breakpoints
b <- min(c(default,alternate), na.rm=TRUE) # Set the minimum for the breakpoints
e <- max(c(default,alternate), na.rm=TRUE) # Set the maximum for the breakpoints
b = b - abs(mean(b,e)*0.1) ; e = e + abs(mean(b,e)*0.1) # add a buffer
ax <- pretty(c(b,e), n = 12) # Make a neat vector for the breakpoints
# Plot the seperate histograms and store them in an object, do not save them yet
hgA <- hist(default, breaks = ax, plot = FALSE) # Save first histogram data
hgB <- hist(alternate, breaks = ax, plot = FALSE) # Save 2nd histogram data
# Now plot them together
ymax = max(c(hgA$counts,hgB$counts))
plot(hgA, col = c_blue, main="RaGPP", xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with FCP
plot(hgB, col = c_green, add = TRUE) # Add Alternate
if (p == 1) {legend("topright",legend = c("Default", "Alternate"), col = c(c_blue,c_green), pch=16, bty = "n")}

# Foliar MRT
default = as.vector(default_states_all$MTT[,1])
alternate = as.vector(alt1_states_all$MTT[,1])
# Determine the x axis range and breakpoints
b <- min(c(default,alternate), na.rm=TRUE) # Set the minimum for the breakpoints
e <- max(c(default,alternate), na.rm=TRUE) # Set the maximum for the breakpoints
b = b - abs(mean(b,e)*0.1) ; e = e + abs(mean(b,e)*0.1) # add a buffer
ax <- pretty(c(b,e), n = 12) # Make a neat vector for the breakpoints
# Plot the seperate histograms and store them in an object, do not save them yet
hgA <- hist(default, breaks = ax, plot = FALSE) # Save first histogram data
hgB <- hist(alternate, breaks = ax, plot = FALSE) # Save 2nd histogram data
# Now plot them together
ymax = max(c(hgA$counts,hgB$counts))
plot(hgA, col = c_blue, main="Foliar MRT", xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with FCP
plot(hgB, col = c_green, add = TRUE) # Add Alternate
if (p == 1) {legend("topright",legend = c("Default", "Alternate"), col = c(c_blue,c_green), pch=16, bty = "n")}

