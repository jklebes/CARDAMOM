###
## This readme briefly describes the CARDAMOM methodology,
## the creation of Global Carbon Project and TRENDY v14 compatible 
## ascii and netcdf files.
##
## Created: 16/03/2026
## Author: T. Luke Smallman (t.l.smallman@ed.ac.uk)
###

###
## Analysis scope

Spatial domain: longitude = -180/180, latitude = -90/90 on a regular WGS84 grid.
Spatial resolution: 0.5 x 0.5 degree
Temporal domain: 2003-2024 (23 years)
Temporal resolution: monthly

###
## CARDAMOM methodology

The CARbon DAta MOdel fraMework (CARDAMOM, Bloom & Williams 2015; Bloom et al., 2016) uses a Bayesian approach within an 
Adpative Proposal - Markov Chain Monte Carlo (AP-MCMC, Haario et al., 2001) algorithm to retrieve ensembles of parameters 
for the Data Assimilation Linked ECosystem (DALEC) model. The retrieval of parameter ensembles is a function of the information 
contained in observations and their uncertainties, ecological theory contained within DALEC's model structure 
(i.e  process representation) and ecological theory contained within the ecological and dynamical constraints (EDCs).
CARDAMOM's AP-MCMC proposes 100 million parameter sets for each analysis. Three analyses, known as chains, are carried out at each 
pixel/location across the terrestrial land surface. Creating multiple chains allows assessment of how robust the analyses are by 
ensuring that the chains at each pixel are consistent with each other (Gelman-Rubin convergence criterion, Gelman and Ruben 1992). 
From each chain a sub-sample of 1000 parameter sets are drawn from across the AP-MCMC analysis are stored for future analysis / 
diagnosis. In practice, the parameter posteriors, ecosystem state and flux estimates are derived from the final 100 subsampled values
per chain, i.e. at each pixel 300 parameter sets (ensemble members) will be processed from which the uncertainty information can 
directly estimated.

CARDAMOM estimates parameters for each pixel independently of all others, i.e. as a function of local / location specific information.
As a result CARDAMOM does not impose plant functional type (PFT) specific parameter values or parameter priors. Therefore, the 
resulting location specific parameters and any coherent spatial variation found is a function of local observations and broad
ecological theory imposed by the EDCs.

The EDCs reduce the likelihood of CARDAMOM retrieving ecologically unrealistic parameters (an inherent challenge due to equifinality). 
EDCs have two broad groupings. The first grouping rejects proposed parameter combinations which are unrealistic, e.g. the prior range
for root and wood (i.e. structural) carbon (C) residence times overlap but for any given combination the root residence time should 
always be shorter. The second grouping restricts unrealistic C dynamics such as exponential change in C pool change in the absence of 
disturbance (e.g. fire or harvest). For further details see Bloom et al., (2016).

DALEC is a model of intermediate complexity which represents C stored and exchanges between both live (labile, foliage, roots and wood), 
dead organic matter (litter and soil) and the atmosphere (Bloom and Williams 2015; Smallman and Williams 2019). Plant photosynthesis 
and water cycling are simulated by the ACM-GPP-ET model as a function of air temperature, absorbed photosynthetically active radiation 
(a function of incoming radiation and DALEC simulated leaf area) and atmospheric CO2 concentration (Smallman and Williams 2019). Leaf
internal CO2 concentration is determined by stomatal conductance which is itself estimated based on the C return on stomatal opening 
within plant hydraulic constraints of water supply from the roots and atmospheric demand for water. Allocation of photosynthate to 
plant tissues is determined by fixed fractions. Canopy phenology is simulated by a day of year model controlling leaf growth 
(i.e. allocating labile C to foliage) and leaf senescence. Turnover of the root and wood pools follows first order kinetics, while 
turnover of the litter and soil pools follows first order kinetics modified by an exponential temperature response function. Each C 
pool and fluxes of carbon between pools and with the ecosystems environment is influences by a minimum of one parameter which are
estimated by CARDAMOM at pixel scale.

###
## Forcing information

Forcing information is used by the DALEC model as a driver, i.e. there is no uncertainty assumed, and thus are inputs to DALEC.

Meteorology: CRU_TS 4.08 - min / max temperature, vapour pressure deficit, precipitation
             CRU-JRA v2.5.5 - short wave radiation, wind speed
Atmospheric CO2: Trendy v14 provided
Fire: GFEDv5.1 Burned area
Deforestation: Global Forest Watch, with fire induced losses removed.
Soil sand / clay: SoilGrids version 2

###
## Observational constraints

Observational constrains are information with uncertainty estimates, used to estimate the likelihood of each proposed parameter set.

Leaf area index: MODIS MCD15A2H (2003-2024) product
Fraction absorbed photosynthetically active radiation: MODIS MCD15A2H (2003-2024) product
Wood stock: A combined map of total above and below ground woody biomass (2007, 2010, 2015-2019) from 
            Xu et al., (2021), Xu, L. et al. Sci. Adv. 7, DOI: 10.1126/sciadv.abe9829
            Santoro, M.; Cartus, O. (2025): ESA Biomass Climate Change Initiative (Biomass_cci): Global datasets of forest above-ground biomass for the years 2007, 2010, 2015, 2016, 2017,             
                                            2018, 2019, 2020, 2021 and 2022, v6.0. NERC EDS Centre for Environmental Data Analysis, 17 April 2025. 
                                            doi:10.5285/95913ffb6467447ca72c4e9d8cf30501.
            The CCI biomass map was upscaled to total wood following Saatchi et al., 2011 allometry.
            The CCI maps were the baseline map, however, non-forest areas are filled with zero values which is unrealistic. The filled areas are replaced with Xu et al.
Soil C prior: SoilGrids version 2
LCA prior: Butler et al., (2017), PNAS, https://doi.org/10.1073/pnas.1708984114
Globally applied prior on the ratio of autotrophic respiration to gross primary productivity (0.54 +/- 0.12, Collalti and Prentice 2019)
 
###
## Creation of TRENDY compatible output files

CARDAMOM's output are converted to netcdf files with variable names and units converted to match those outlined in the TRENDYv13
variables description Excel spreadsheet. Each file contains a variable defining the pixel area and land cover fraction. All CARDAMOM
estimates are provided in units per pixel area. Thus when scaling to global totals each pixel should be multiplied by the grid area
and land cover fraction. 

In contrast to other models contributing to TRENDYv13, CARDAMOM has an explicit estimate of uncertainty which comes from the 
underlying ensembles of parameters estimated for each location in the analysis grid. Storing full ensembles of all variables is 
impractical due to the very large data storage requirements. Instead, for each pixel and time step the ensemble median, 2.5 % and 
97.5 % quantiles are estimates and stored in the appropriate variable file. Note that we do no include information on temporal 
correlations in these files, which means that uncertainty propagation from monthly to annual time steps will lead to an overestimate 
of uncertainty typically on the order of >30 %. Annual estimates with explicit uncertainty propagation can be provided on request. 

###
## Creation of Global Carbon Project NBP files

Annual totals for net biome productivity (NBP, PgC/yr) are provided in ascii files with white space delimited, as requested in the 
TRENDYv13 protocol. Estimates are generated by first determining the annual NBP per pixel, i.e. including explicit propagation of 
ensemble uncertainty from monthly to annual time scales. These pixel estimates are scaled by pixel area and land cover fraction (provided in file). 
We currently lack a robust understanding of how errors / uncertainty are correlated in space, therefore we conservatively assume 
uncertainties are fully correlated between pixels when aggregating across pixel. The resulting uncertainty estimates are provided in 
the ascii file along with the median estimate. We note, however, that this approach means our reported uncertainty is an overestimate.

###
## Creation of TRENDYv14 simulations

The TRENDYv14 protocol targets the creation of 4 simulations to partition the natural C-balance versus the impacts of changing 
climate, CO2 concentrations and land use change / management, i.e. S0-S3. CARDAMOM is inherently dependent on available observations 
to constrain its analyses, as a result does not simulate the pre-industrial period or initialise via a 'spin-up' to steady state. 

As a result we provide simulations S2 and S3 but not S0 and S1. 
S3, i.e. all forcings is CARDAMOM's default analysis. We provide S2, exclusion of land use change, by re-running the already retrieved
parameter sets, but ignoring any deforestation information drawn from the global forest watch dataset. Burned area is assumed to be
imposed as normal.

###
## References
###

Bloom, A. A. and M. Williams (2015). Constraining ecosystem carbon dynamics in a data-limited world: integrating ecological 
"common sense" in a model data fusion framework., Biogeosciences, 12, 1299-1315, doi:10.5194/bg-12-1299-2015.

Bloom, A. A., Exbrayat, J.-F., van~der Velde, I. R., Feng, L., and Williams, M. (2016). The decadal state of the terrestrial carbon 
cycle: Global retrievals of terrestrial carbon allocation, pools, and residence times, P. Natl. Acad. Sci. USA, 113, 1285-1290, 
doi: 10.1073/pnas.1515160113.

Collalti, A. and Prentice, I. C. (2019). Is NPP proportional to GPP? Waring’s hypothesis 20 years on, Tree Physiol., 39, 1473-1483, 
doi10.1093/treephys/tpz034

Haario, H., Saksman, E. and Tamminen, J. (2001). An adaptive Metropolis algorithm. Bernoulli, 7 (2), 223-242, 2001.

Gelman, A., and Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. Stat. Sci., 7, 457–472. 

Smallman, T. L., and Williams, M. (2019). Description and validation of an intermediate complexity model for ecosystem photosynthesis 
and evapotranspiration: ACM-GPP-ETv1, Geosci. Model Dev., 12 (6), doi:10.5194/gmd-12-2227-2019.

