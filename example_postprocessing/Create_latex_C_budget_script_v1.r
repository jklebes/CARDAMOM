
###
## Script for the creation of a latex script to generate C-budgets for DALEC
## Author: T L Smallman (t.l.smallman@ed.ac.uk)
## Created: 12/07/2022
## Last updated: 14/07/2022
## Contributing authors: 
## NOTES: 
## 1) This version presents a combined labile and foliage pool. 
## 2) Assumes a CDEA style phenology with both direct GPP to foliage and via labile pathways
###

###
## Prepare the work space

# Set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# Load R libraries

# Load user defined functions
#source("./R_functions/load_all_cardamom_functions.r")

# Define the output file name for the created script
outfilename = "latex_C_budget.tex"

###
## Load files from which C-budget is extracted, prepare C-budget values

# Load information file
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Forests2020_Mexico_Kiuic_chronosequence/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/ODA_extension_Africa_agb/infofile.RData")

# If this is a site analysis 
# which site are we using?
site_nos = 1

###
## Define figure labels

# The caption must be written in correct latex syntex
# NOTE: if the latex language requires use of "\" ensure it is a "\\".
figure_caption = "Steady state C-budget for semi-evergreen FCP. Numbers show median estimate of fluxes (alongside arrows) and of stocks (in boxes). Units are gC m^\\({-2}\\) for stocks and gC m\\(^{-2}\\) yr^\\({-1}\\) for fluxes. 95\\% confidence intervals are shown in a fractional form with 2.5 and 97.5 percentiles as numerator and denominator. Black fluxes are biogenic, including net primary production (\\(NPP\\)), mortality (\\(Mort\\)), autotrophic respiration (\\(Ra\\)) and heterotrophic respiration (\\(Rh\\)). \\(NEE = Ra + Rh - GPP\\). \\(NBE = NEE + E_{total}\\). Red fluxes are fire-driven emissions (\\(E\\))."
# The label will be used for referencing the figure in the latex document
figure_label = "SIFig:fcp_budget"
# Desired precision, i.e. decimal places
dp = 1

if (PROJECT$spatial_type == "grid") {
    # A gridded analysis
    
    # Specifiy quantiles we want to extract, should only be 3 (low/median/high).
    # From gridded analysis this must be from thos available in the file. 
    # Check grid_output$num_quantiles for available
    quantiles_wanted = c(1,4,7)

    # Load the processed site file
    load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

    # Extract or calculate required derived values
    # NOTE: unit conversion from gC/m2/day to gC/m2/yr

    # NATURAL FLUXES
    gpp_gCm2yr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    rauto_gCm2yr = format(round(apply(grid_output$mean_rauto_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    rhet_litter_gCm2yr = format(round(apply(grid_output$mean_rhet_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    rhet_som_gCm2yr = format(round(apply(grid_output$mean_rhet_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    rhet_gCm2yr = format(round(apply(grid_output$mean_rhet_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    npp_gCm2yr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    npp_labilefoliage_gCm2yr = format(round(apply(grid_output$mean_combined_alloc_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    npp_roots_gCm2yr = format(round(apply(grid_output$mean_alloc_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    npp_wood_gCm2yr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    foliage_to_litter_gCm2yr = format(round(apply(grid_output$mean_foliage_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    roots_to_litter_gCm2yr = format(round(apply(grid_output$mean_roots_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    wood_to_litter_gCm2yr = format(round(apply(grid_output$mean_wood_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    litter_to_som_gCm2yr = format(round(apply(grid_output$mean_litter_to_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    # FIRE FLUXES
    FIRElitter_labile_foliage_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_labile_gCm2day[,,quantiles_wanted]+grid_output$mean_FIRElitter_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIRElitter_roots_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIRElitter_wood_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIRElitter_litter_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_labile_foliage_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_labile_gCm2day[,,quantiles_wanted]+grid_output$mean_FIREemiss_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_roots_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_wood_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_litter_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_som_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    fire_gCm2yr = format(round(apply(grid_output$mean_fire_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    # Net fluxes
    nbe_gCm2yr = format(round(apply(grid_output$mean_nbe_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    nee_gCm2yr = format(round(apply(grid_output$mean_nee_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
    # STOCKS
    labile_foliage_gCm2 = format(round(apply(grid_output$mean_labile_gCm2[,,quantiles_wanted]+grid_output$mean_foliage_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    roots_gCm2 = format(round(apply(grid_output$mean_roots_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    wood_gCm2 = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    litter_gCm2 = format(round(apply(grid_output$mean_litter_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    som_gCm2 = format(round(apply(grid_output$mean_som_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
  
} else if (PROJECT$spatial_type == "site") {
    # A site analysis
    
    # Specifiy quantiles we want to extract, should only be 3 (low/median/high)
    quantiles_wanted = c(0.025,0.50,0.975)    
    
    # Load the processed site file
    load(paste(PROJECT$results_processedpath,PROJECT$sites[site_nos],".RData",sep=""))

    # Extract or calculate required derived values
    # NOTE: unit conversion from gC/m2/day to gC/m2/yr

    # NATURAL FLUXES
    gpp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    rauto_gCm2yr = format(round(quantile(apply(states_all$rauto_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    rhet_litter_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    rhet_som_gCm2yr = format(round(quantile(apply(states_all$rhet_som_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    rhet_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    npp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day-states_all$rauto_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    npp_labilefoliage_gCm2yr = format(round(quantile(apply(states_all$alloc_foliage_gCm2day+states_all$labile_to_foliage_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    npp_roots_gCm2yr = format(round(quantile(apply(states_all$alloc_roots_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    npp_wood_gCm2yr = format(round(quantile(apply(states_all$alloc_wood_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    foliage_to_litter_gCm2yr = format(round(quantile(apply(states_all$foliage_to_litter_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    roots_to_litter_gCm2yr = format(round(quantile(apply(states_all$roots_to_litter_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    wood_to_litter_gCm2yr = format(round(quantile(apply(states_all$wood_to_litter_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    litter_to_som_gCm2yr = format(round(quantile(apply(states_all$litter_to_som_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    # FIRE FLUXES
    FIRElitter_labile_foliage_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_labile_gCm2day + states_all$FIRElitter_foliage_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIRElitter_roots_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_roots_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIRElitter_wood_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_wood_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIRElitter_litter_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_litter_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_labile_foliage_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_labile_gCm2day + states_all$FIREemiss_foliage_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_roots_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_roots_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_wood_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_wood_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_litter_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_litter_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    FIREemiss_som_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_som_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    fire_gCm2yr = format(round(quantile(apply(states_all$fire_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    # Net fluxes
    nbe_gCm2yr = format(round(quantile(apply((states_all$fire_gCm2day+states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)-states_all$gpp_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    nee_gCm2yr = format(round(quantile(apply((states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)-states_all$gpp_gCm2day,1,mean), prob=quantiles_wanted) * 365.25, digits = dp), nsmall = dp)
    # STOCKS
    labile_foliage_gCm2 = format(round(quantile(apply(states_all$labile_gCm2 + states_all$foliage_gCm2,1,mean), prob=quantiles_wanted), digits = dp), nsmall = dp)
    roots_gCm2 = format(round(quantile(apply(states_all$roots_gCm2,1,mean), prob=quantiles_wanted), digits = dp), nsmall = dp)
    wood_gCm2 = format(round(quantile(apply(states_all$wood_gCm2,1,mean), prob=quantiles_wanted), digits = dp), nsmall = dp)
    litter_gCm2 = format(round(quantile(apply(states_all$litter_gCm2,1,mean), prob=quantiles_wanted), digits = dp), nsmall = dp)
    som_gCm2 = format(round(quantile(apply(states_all$som_gCm2,1,mean), prob=quantiles_wanted), digits = dp), nsmall = dp)
    
} else {
    # We have a compatibility problem
    stop("PROJECT$spatial_type does not have a compatible value, i.e. grid or site")
} # end if is_grid

###
## Begin writting the latex code

col_sep = ""
nos_cols = 20

write(    c("% Biogenic flux and emissions figure with budgets"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = FALSE)
write(    c("\\begin{figure}[]"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("% NOTE: \\put(x coord,y coord){ ... } where to put ..."), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\vector(x slope,y slope){length} used for arrows"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\line(x slope,y slope){length} used for lines"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\framebox(x,y){...} puts ... at box centre"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("   \\centering"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("   \\fbox{"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\setlength{\\unitlength}{0.60cm}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\begin{picture}(28,12)"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % GPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,4.50){\\vector(1,0){3.25}} % arrow for GPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(1.0,4.7){$GPP$}               % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(0.25,3.9){\\small ",gpp_gCm2yr[2],"}                 % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(1.90,3.925){\\scriptsize $\\frac{",gpp_gCm2yr[1],"}{",gpp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Ra"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.0,0.78){$Ra$}              % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.2,0.1){\\small ",rauto_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(26.9,0.1){\\scriptsize $\\frac{",rauto_gCm2yr[1],"}{",rauto_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % CUE and associated arrows         "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.67,4.35){$CUE$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.50,3.70){\\dashbox{0.2}(1.75,1.65)} % Box around CUE"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.375,3.70){\\line(0,-1){3.05}} % Vertical line down from box"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.375,0.65){\\vector(1,0){21.625}} % horizontal line to Ra"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)         
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add labels and box for the internal C-cycle dynamics"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.30,7.70){\\textit{Internal carbon rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.25,0){\\dashbox{0.2}(20.25,8.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add labels for input / output"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,11.5){\\textit{Input}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,11.0){\\textit{carbon}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,10.5){\\textit{rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26,11.5){\\textit{Output}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26,11.0){\\textit{carbon}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26,10.5){\\textit{rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Total NPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(5.25,4.5){\\vector(1,0){3.0}} % arrow for NPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.0,4.7){$NPP$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(5.35,3.9){\\small ",npp_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.9,3.925){\\scriptsize $\\frac{",npp_gCm2yr[1],"}{",npp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add partitioning point in graphic"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(8.70,4.5){\\small \\rotatebox[origin=c]{90}{NPP allocation}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(8.25,1.5){\\dashbox{0.2}(1.2,6)} % Add box around label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to labile + foliage - change by model?"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.45,7){\\vector(1,0){3.3}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.65,7.2){\\small $NPP_{fol+lab}$}    % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.88,6.4){\\small ",npp_labilefoliage_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.02,6.45){\\scriptsize $\\frac{",npp_labilefoliage_gCm2yr[1],"}{",npp_labilefoliage_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to fine root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.45,4.5){\\vector(1,0){3.3}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.75,4.7){\\small $NPP_{root}$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.85,3.9){\\small ",npp_roots_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.07,3.95){\\scriptsize $\\frac{",npp_roots_gCm2yr[1],"}{",npp_roots_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.45,2){\\vector(1,0){3.3}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.75,2.2){\\small $NPP_{wood}$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.85,1.45){\\small ",npp_wood_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.07,1.50){\\scriptsize $\\frac{",npp_wood_gCm2yr[1],"}{",npp_wood_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Foliage + labile C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,6.0){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.95,7.45){$C_{fol+lab}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,6.80){\\small ",labile_foliage_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,6.20){\\scriptsize $\\frac{",labile_foliage_gCm2[1],"}{",labile_foliage_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine root C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,3.5){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.20,4.95){$C_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,4.35){\\small ",roots_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,3.75){\\scriptsize $\\frac{",roots_gCm2[1],"}{",roots_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,1){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.1,2.45){$C_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.10,1.85){\\small ",wood_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.20,1.25){\\scriptsize $\\frac{",wood_gCm2[1],"}{",wood_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % litter/mortality, natural and fire driven (red)"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile + foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,7.2){\\small $Mort_{fol+lab}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,6.4){\\small ",foliage_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,6.45){\\scriptsize $\\frac{",foliage_to_litter_gCm2yr[1],"}{",foliage_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,6.4){\\color{red}{\\small ",FIRElitter_labile_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,6.45){\\color{red}{\\scriptsize $\\frac{",FIRElitter_labile_foliage_gCm2yr[1],"}{",FIRElitter_labile_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,7.0){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,7.0){\\vector(1,-1){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,4.65){\\small $Mort_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,3.9){\\small ",roots_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,3.95){\\scriptsize $\\frac{",roots_to_litter_gCm2yr[1],"}{",roots_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,3.9){\\color{red}{\\small ",FIRElitter_roots_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,3.95){\\color{red}{\\scriptsize $\\frac{",FIRElitter_roots_gCm2yr[1],"}{",FIRElitter_roots_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,4.5){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,4.5){\\vector(1,2){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.70,2.2){\\small $Mort_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,1.45){\\small ",wood_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,1.50){\\scriptsize $\\frac{",wood_to_litter_gCm2yr[1],"}{",wood_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,1.45){\\color{red}{\\small ",FIRElitter_wood_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,1.50){\\color{red}{\\scriptsize $\\frac{",FIRElitter_wood_gCm2yr[1],"}{",FIRElitter_wood_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Csom"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,2.0){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,2.0){\\vector(1,1){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fire emission fluxes"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile + foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.37,9.4){\\small \\rotatebox[origin=c]{90}{$E_{fol+lab}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(14.0,8.1){\\color{red}{\\line(0,1){2.9}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(14.1,9.7){\\color{red}{\\small ",FIREemiss_labile_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(14.1,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_labile_foliage_gCm2yr[1],"}{",FIREemiss_labile_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine roots"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(11.85,9.4){\\small \\rotatebox[origin=c]{90}{$E_{root}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.45,5.0){\\color{red}{\\line(1,0){0.25}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.45,5.0){\\color{red}{\\line(0,1){6.0}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(12.50,9.7){\\color{red}{\\small ",FIREemiss_roots_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(12.45,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_roots_gCm2yr[1],"}{",FIREemiss_roots_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.05,9.4){\\small \\rotatebox[origin=c]{90}{$E_{wood}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.2,2.7){\\color{red}{\\line(1,0){0.35}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.55,2.7){\\color{red}{\\line(0,1){8.3}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,9.7){\\color{red}{\\small ",FIREemiss_wood_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.70,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_wood_gCm2yr[1],"}{",FIREemiss_wood_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.25,9.4){\\small \\rotatebox[origin=c]{90}{$E_{litter}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.75,7.2){\\color{red}{\\line(0,1){3.8}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.85,9.7){\\color{red}{\\small ",FIREemiss_litter_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.90,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_litter_gCm2yr[1],"}{",FIREemiss_litter_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % SOM"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.5,9.4){\\small \\rotatebox[origin=c]{90}{$E_{som}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,3.4){\\color{red}{\\line(1,0){0.3}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.05,3.4){\\color{red}{\\line(0,1){7.6}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.1,9.7){\\color{red}{\\small ",FIREemiss_som_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.1,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_som_gCm2yr[1],"}{",FIREemiss_som_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fire emissions total"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.4,11.2){$E_{total}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.45,11.0){\\color{red}{\\vector(1,0){10.95}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.4,10.5){\\color{red}{\\small ",fire_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(24.2,10.5){\\color{red}{\\scriptsize $\\frac{",fire_gCm2yr[1],"}{",fire_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Litter C pool"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.25,5.2){\\framebox(2.5,2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.75,6.45){$C_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.4,5.6){\\small ",litter_gCm2[2],"}               % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.3,5.6){\\scriptsize $\\frac{",litter_gCm2[1],"}{",litter_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Som C pool"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.25,1.0){\\framebox(2.5,2.75)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.65,2.9){$C_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.55,2.1){\\small ",som_gCm2[2],"}                 % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.75,1.45){\\scriptsize $\\frac{",som_gCm2[1],"}{",som_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - natural and fire"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(21.5,5.2){\\vector(0,-1){1.40}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - natural"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.7,4.65){\\small ",litter_to_som_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.65,4.1){\\scriptsize $\\frac{",litter_to_som_gCm2yr[1],"}{",litter_to_som_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - combusted litter to som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.4,4.65){\\color{red}{\\small ",FIRElitter_litter_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.5,4.1){\\color{red}{\\scriptsize $\\frac{",FIRElitter_litter_gCm2yr[1],"}{",FIRElitter_litter_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Heterotrophic respiration of litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,6.0){\\line(1,0){2.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.65,6.3){\\small $Rh_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.65,5.45){\\small ",rhet_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.25,5.50){\\scriptsize $\\frac{",rhet_litter_gCm2yr[1],"}{",rhet_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.0,6.0){\\vector(1,-2){0.70}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Heterotrophic respiration of som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,2.75){\\line(1,0){2.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.65,3.0){\\small $Rh_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.65,2.2){\\small ",rhet_som_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.25,2.25){\\scriptsize $\\frac{",rhet_som_gCm2yr[1],"}{",rhet_som_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.0,2.75){\\vector(1,1){0.85}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Total heterotrophic respiration"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.0,4.6){$Rh$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.2,4.0){\\small ",rhet_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(26.9,4.0){\\scriptsize $\\frac{",rhet_gCm2yr[1],"}{",rhet_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Net Biome Exchange "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.25,8.5){\\dashbox{0.2}(2.2,2.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.75,10.05){NBE}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.75,9.4){\\small ",nbe_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.75,8.8){\\scriptsize $\\frac{",nbe_gCm2yr[1],"}{",nbe_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Net Ecosystem Exchange "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.25,8.5){\\dashbox{0.2}(2.2,2.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.75,10.05){NEE}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.75,9.4){\\small ",nee_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.75,8.8){\\scriptsize $\\frac{",nee_gCm2yr[1],"}{",nee_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\end{picture}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("        }"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(paste("   \\caption{",figure_caption,"}", sep="")), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(paste("   \\label{",figure_label,"}"), sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\end{figure}   "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)

# DONE!
