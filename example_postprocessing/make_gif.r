
# Set path to project file
project_file = "~/gcel_ceph/cardamom_analyses/lsmallma/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.004_MHMCMC/global_0.5deg_dalec4_trendyv14_LCA_TWB_GPP_fAPAR/infofile.RData"
# Select variable name to process
variable_name = "nbp_gCm2day" # Note the variable selected is assumed to have dimensions of site,quantiles,time

# Load needed functions 
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/load_all_cardamom_functions.r")

# Load project file
load(project_file)
# load the CARDAMOM files
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

###
## Begin creating information for processing and subsequent saving to files

# Time information
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
steps_per_year = dim(grid_output$lai_m2m2)[3] / nos_years

# create lat / long axes, assumes regular WGS-84 grid
grid_lat = grid_output$lat ; grid_long = grid_output$long

# load colour palette
colour_choices = colorRampPalette(brewer.pal(11,"Spectral"))
colour_choices_CI = colorRampPalette(brewer.pal(9,"Purples"))
# array sizes are always the same so
colour_choices = colour_choices(100)
colour_choices_CI = colour_choices_CI(100)   

# Extract the available quantiles 
quantiles_wanted = grid_output$num_quantiles
nos_quantiles = length(quantiles_wanted)
# Check that the quantiles we want to use are available
# Minimum quantiles
if (length(which(quantiles_wanted == 0.025)) == 1) {
    q1_quant = which(quantiles_wanted == 0.025)
    q1_quant_lab = "2.5pc"
    q1_quant_longlab = "2.5 % quantile"
} else {
    stop("Desired low quantile cannot be found")
}
# A lower quantile
if (length(which(quantiles_wanted == 0.05)) == 1) {
    q2_quant = which(quantiles_wanted == 0.05)
    q2_quant_lab = "5pc"
    q2_quant_longlab = "5 % quantile"
} else {
    stop("Desired low quantile cannot be found")
}
# A lower quartile
if (length(which(quantiles_wanted == 0.16)) == 1) {
    q3_quant = which(quantiles_wanted == 0.16)
    q3_quant_lab = "16pc"
    q3_quant_longlab = "16 % quantile"
} else {
    stop("Desired lower quartile cannot be found")
}
# The median estimate
if (length(which(quantiles_wanted == 0.5)) == 1) {
    mid_quant = which(quantiles_wanted == 0.5)
} else {
    stop("Median quantile cannot be found")
}
# A upper quartile
if (length(which(quantiles_wanted == 0.84)) == 1) {
    q4_quant = which(quantiles_wanted == 0.84)
    q4_quant_lab = "84pc"
    q4_quant_longlab = "84 % quantile"
} else {
    stop("Desired upper quartile cannot be found")
}
# A upper quantile
if (length(which(quantiles_wanted == 0.95)) == 1) {
    q5_quant = which(quantiles_wanted == 0.95)
    q5_quant_lab = "95pc"
    q5_quant_longlab = "95 % quantile"
} else {
    stop("Desired high quantile cannot be found")
}
# Maximum quantile
if (length(which(quantiles_wanted == 0.975)) == 1) {
    q6_quant = which(quantiles_wanted == 0.975)
    q6_quant_lab = "97.5pc"
    q6_quant_longlab = "97.5 % quantile"
} else {
    stop("Desired max quantile cannot be found")
}

# Create variable to populate
variable = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))

# Fill the output arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n] ; j = grid_output$j_location[n]
         # load into the variable
         variable[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output[[variable_name]][n,,] 

     } # Does the file exist / has it been processed

} # site loop

###

fig_height = 3000*0.65 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height
if (grepl("27700",PROJECT$grid_type)) { fig_height = 8000*0.65 ; fig_width = 7200*0.65 }
fig_height = fig_height * 0.5 ; fig_width = fig_width * 0.5

# create_pngs for converting for the actual variable
png(file=paste(PROJECT$figpath,variable_name,"_%03d.png",sep=""), width=fig_width, height=fig_height)
   z_axis = c(-1,1)*max(abs(quantile(as.vector(variable[,,mid_quant,]), prob=c(0.01,0.99),na.rm=TRUE)))
   variable[,,mid_quant,][which(variable[,,mid_quant,] < z_axis[1])] = z_axis[1]
   variable[,,mid_quant,][which(variable[,,mid_quant,] > z_axis[2])] = z_axis[2]
   for (t in seq(1, length(PROJECT$model$timestep_days))){
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.6), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(x = grid_long, y = grid_lat, z = variable[,,mid_quant,t],col = colour_choices
                  ,main = variable_name, zlim = z_axis, axes=FALSE
                  ,cex.main = 0.9, legend.width = 3.0, cex = 1.5, axis.args = list(cex.axis = 1.8,hadj = 0.1))
        map(add=TRUE, lwd = 2)     
   }
dev.off()
# convert pngs to one gif using ImageMagick
system(paste("convert -delay 45 ",PROJECT$figpath,variable_name,"*.png ",PROJECT$figpath,variable_name,".gif", sep=""))
# cleaning up
#remove_list = list.files(PROJECT$figpath,pattern=".png", full.names=TRUE) ; remove_list = remove_list[grepl(paste(PROJECT$figpath,variable_name,sep=""), remove_list)]
#file.remove(remove_list)

# create_pngs for converting for the uncertainty layer (95 % CI)
uncertainty = variable[,,q6_quant,]-variable[,,q1_quant,]
png(file=paste(PROJECT$figpath,variable_name,"_uncertainty_%03d.png",sep=""), width=fig_width, height=fig_height)
   z_axis = c(0,quantile(as.vector(uncertainty), prob=c(0.99),na.rm=TRUE))
   uncertainty[,,][which(uncertainty[,,] < z_axis[1])] = z_axis[1]
   uncertainty[,,][which(uncertainty[,,] > z_axis[2])] = z_axis[2]
   for (t in seq(1, length(PROJECT$model$timestep_days))){
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.6), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(x = grid_long, y = grid_lat, z = uncertainty[,,t],col = colour_choices_CI
                  ,main = variable_name, zlim = z_axis, axes=FALSE
                  ,cex.main = 0.9, legend.width = 3.0, cex = 1.5, axis.args = list(cex.axis = 1.8,hadj = 0.1))
        map(add=TRUE, lwd = 2)     
   }
dev.off()
# convert pngs to one gif using ImageMagick
system(paste("convert -delay 45 ",PROJECT$figpath,variable_name,"_uncertainty_*.png ",PROJECT$figpath,variable_name,"_uncertainty.gif", sep=""))
# cleaning up
#remove_list = list.files(PROJECT$figpath,pattern=".png", full.names=TRUE) ; remove_list = remove_list[grepl(paste(PROJECT$figpath,variable_name,"_uncertainty",sep=""), remove_list)]
#file.remove(remove_list)

