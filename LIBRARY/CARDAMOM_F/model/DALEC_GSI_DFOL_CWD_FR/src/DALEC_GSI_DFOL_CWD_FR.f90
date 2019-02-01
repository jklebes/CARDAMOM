
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL        &
         ,dble_one,dble_zero  &
         ,vsmall              &
         ,arrhenious          &
         ,acm_gpp             &
         ,acm_albedo_gc       &
         ,acm_meteorological_constants  &
         ,calculate_shortwave_balance   &
         ,calculate_longwave_isothermal &
         ,calculate_daylength           &
         ,freeze              &
         ,co2comp_saturation  &
         ,co2comp_half_sat_conc &
         ,co2_half_saturation &
         ,co2_compensation_point &
         ,kc_saturation       &
         ,kc_half_sat_conc    &
         ,calculate_Rtot      &
         ,calculate_aerodynamic_conductance &
         ,linear_model_gradient &
         ,seconds_per_day     &
         ,seconds_per_hour    &
         ,seconds_per_step    &
         ,root_biomass        &
         ,root_reach          &
         ,min_root            &
         ,max_depth           &
         ,root_k              &
         ,top_soil_depth      &
         ,mid_soil_depth      &
         ,soil_depth          &
         ,previous_depth      &
         ,nos_root_layers     &
         ,deltat_1            &
         ,water_flux          &
         ,layer_thickness     &
         ,min_layer        &
         ,nos_soil_layers  &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,meant            &
         ,meant_K          &
         ,stomatal_conductance &
         ,avN              &
         ,iWUE             &
         ,NUE                           &
         ,pn_max_temp                   &
         ,pn_opt_temp                   &
         ,pn_kurtosis                   &
         ,e0                            &
         ,co2_half_sat                  &
         ,co2_comp_point                &
         ,minlwp                        &
         ,max_lai_lwrad_transmitted     &
         ,lai_half_lwrad_transmitted    &
         ,max_lai_nir_reflection        &
         ,lai_half_nir_reflection       &
         ,max_lai_par_reflection        &
         ,lai_half_par_reflection       &
         ,max_lai_par_transmitted       &
         ,lai_half_par_transmitted      &
         ,max_lai_nir_transmitted       &
         ,lai_half_nir_transmitted      &
         ,max_lai_lwrad_reflected       &
         ,lai_half_lwrad_reflected      &
         ,soil_swrad_absorption         &
         ,max_lai_lwrad_release         &
         ,lai_half_lwrad_release        &
         ,leafT                         &
         ,mint             &
         ,maxt             &
         ,swrad            &
         ,co2              &
         ,doy              &
         ,wind_spd         &
         ,vpd_pa           &
         ,lai              &
         ,days_per_step    &
         ,days_per_step_1  &
         ,dayl_seconds_1   &
         ,dayl_seconds     &
         ,dayl_hours       &
         ,disturbance_residue_to_litter &
         ,disturbance_residue_to_cwd    &
         ,disturbance_residue_to_som    &
         ,disturbance_loss_from_litter  &
         ,disturbance_loss_from_cwd     &
         ,disturbance_loss_from_som     &
         ,itemp,ivpd,iphoto&
         ,extracted_C      &
         ,dim_1,dim_2      &
         ,nos_trees        &
         ,nos_inputs       &
         ,leftDaughter     &
         ,rightDaughter    &
         ,nodestatus       &
         ,xbestsplit       &
         ,nodepred         &
         ,bestvar

!!!!!!!!!!
! Random Forest GPP emulator
!!!!!!!!!!

! arrays for the emulator, just so we load them once and that is it cos they be
! massive
integer ::    dim_1, & ! dimension 1 of response surface
              dim_2, & ! dimension 2 of response surface
          nos_trees, & ! number of trees in randomForest
         nos_inputs    ! number of driver inputs
double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                    rightDaughter, & ! right daughter for forets
                                                       nodestatus, & ! nodestatus for forests
                                                       xbestsplit, & ! for forest
                                                         nodepred, & ! prediction value for each tree
                                                          bestvar    ! for randomForests
!!!!!!!!!
! Parameters
!!!!!!!!!

! useful technical parameters
double precision, parameter :: xacc = 1d-4        & ! accuracy parameter for zbrent bisection proceedure ! 0.0001
                              ,dble_zero = 0d0    &
                              ,dble_one = 1d0     &
                              ,vsmall = tiny(0d0)

integer, parameter :: nos_root_layers = 3, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: pi = 3.1415927d0,    &
                             pi_1 = pi**(-1d0),     &
                              pi2 = pi**2d0,        &
                           two_pi = pi*2d0,         &
                       deg_to_rad = pi/180d0,       &
              sin_dayl_deg_to_rad = sin( 23.45d0 * deg_to_rad ), & ! repeated function in acm
                          gravity = 9.8067d0,       & ! acceleration due to gravity, ms-1
                            boltz = 5.670400d-8,    & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96d0,         &
                      emiss_boltz = emissivity * boltz, &
                  sw_par_fraction = 0.5d0,          & ! fraction of short-wave radiation which is PAR
                           freeze = 273.15d0,       &
                       gs_H2O_CO2 = 1.646259d0,     & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                     gs_H2O_CO2_1 = gs_H2O_CO2 ** (-dble_one), &
                       gb_H2O_CO2 = 1.37d0,         & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
          partial_molar_vol_water = 18.05d-6,       & ! partial molar volume of water, m3 mol-1 at 20C
                       umol_to_gC = 1d-6*12d0,      & ! conversion of umolC -> gC
                       gC_to_umol = umol_to_gC**(-dble_one), & ! conversion of gC -> umolC
                 mmol_to_kg_water = 1.8d-5,         & ! milli mole conversion to kg
                   mol_to_g_water = 18d0,           & ! molecular mass of water
                     mol_to_g_co2 = 12d0,           & ! molecular mass of CO2 (g)
                     g_to_mol_co2 = 1d0/12d0,       &
!snowscheme       density_of_water = 998.9d0,         & ! density of !water kg.m-3
                   gas_constant_d = 287.04d0,       & ! gas constant for dry air (J.K-1.mol-1)
                             Rcon = 8.3144d0,       & ! Universal gas constant (J.K-1.mol-1)
                        vonkarman = 0.41d0,         & ! von Karman's constant
                      vonkarman_2 = vonkarman**2d0, & ! von Karman's constant^2
                            cpair = 1004.6d0          ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

! photosynthesis / respiration parameters
double precision, parameter :: &
                    kc_saturation = 310d0,        & ! CO2 half saturation, saturation value
                 kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, half sat
               co2comp_saturation = 36.5d0,       & ! CO2 compensation point, saturation
            co2comp_half_sat_conc = 9.46d0,       & ! CO2 comp point, half sat
                                                    ! Each of these are temperature
                                                    ! sensitivty
              leaf_life_weighting = 1d0/2d0,      & ! inverse of averaging period of lagged effects
                                                    ! probably should be an actual parmeter
                      Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                    ! lost as growth respiration
                                                    ! (i.e. 0.28 .eq. xNPP)
                  one_Rg_fraction = dble_one - Rg_fraction
! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5d0,        & ! tortuosity
                           gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                        max_depth = 2d0,          & ! max root depth (m)
                           root_k = 100d0,        & ! root biomass needed to reach 50% depth (gbiomass/m2)
                      root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                    !                 Williams et al 1996 = 0.0001
                    root_radius_1 = root_radius**(-dble_one), &
              root_cross_sec_area = 3.141593d-08, & ! root cross sectional area (m2)
                                                    ! = pi * root_radius * root_radius
                     root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                    ! 0.5e6 Williams et al 1996
                                                    ! 0.31e6 Bonan et al 2014
          root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1d0), &
               const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                             head = 0.009807d0,   & ! head of pressure (MPa/m)
                           head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

! structural parameters
double precision, parameter :: &
                    canopy_height = 9d0,          & ! canopy height assumed to be 9 m
                     tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                         min_wind = 0.1d0,        & ! minimum wind speed at canopy top
                        min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                      soil_roughl = 0.05d0,       & ! soil roughness length (m)
                   top_soil_depth = 0.1d0,        & ! thickness of the top soil layer (m)
                   mid_soil_depth = 0.2d0,        & ! thickness of the second soil layer (m)
                         min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                          min_lai = 1.5d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                  min_throughfall = 0.2d0           ! minimum fraction of precipitation which
                                                    ! is through fall

! timing parameters
double precision, parameter :: &
                 seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                  seconds_per_day = 86400d0,        & ! Number of seconds per day
                seconds_per_day_1 = 1d0/seconds_per_day       ! Inverse of seconds per day

!!!!!!!!!
! Module level variables
!!!!!!!!!

! management and gsi related values
integer :: gsi_lag_remembered
! local variables for GSI phenology model
double precision :: Tfac,Photofac,VPDfac & ! oC, seconds, Pa
                   ,Rg_from_labile       &
                   ,delta_gsi,tmp,tmp_lai,gradient         &
                   ,fol_turn_crit,lab_turn_crit    &
                   ,gsi_history(22),just_grown

double precision, allocatable, dimension(:) :: extracted_C,itemp,ivpd,iphoto, &
                                               disturbance_residue_to_litter, &
                                               disturbance_residue_to_som,    &
                                               disturbance_residue_to_cwd,    &
                                               disturbance_loss_from_litter,  &
                                               disturbance_loss_from_cwd,     &
                                               disturbance_loss_from_som,     &
                                               tmp_x, tmp_m

! Phenological choices
! See source below for details of these variables
double precision :: root_cost,root_life,       &
                    leaf_cost,leaf_life

! hydraulic model variables
integer :: water_retention_pass, soil_layer, sunrise, sunset
double precision, dimension(nos_soil_layers+1) :: layer_thickness    ! thickness of soil layers (m)
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil
double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                         demand, & ! maximum potential canopy hydraulic demand
                                                     water_flux    ! potential transpiration flux (mmol.m-2.s-1)

double precision :: root_reach, root_biomass,soil_depth, &
                    soilRT, &
  new_depth,previous_depth, & ! depth of bottom of soil profile
               canopy_wind, & ! wind speed (m.s-1) at canopy top
                     ustar, & ! friction velocity (m.s-1)
                  ustar_Uh, &
            air_density_kg, & ! air density kg/m3
                    roughl, & ! roughness length (m)
              displacement, & ! zero plane displacement (m)
                max_supply, & ! maximum water supply (mmolH2O/m2/day)
                     meant, & ! mean air temperature (oC)
                   meant_K, & ! mean air temperature (K)
        canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
          soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
            soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
             sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
      stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.day-1)
   aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
         convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                    lambda, & ! latent heat of vapourisation (J.kg-1)
                     psych, & ! psychrometric constant (kPa K-1)
                     slope    ! Rate of change of saturation vapour pressure
                              ! with temperature (kPa.K-1)

! Module level variables for ACM_GPP_ET parameters
double precision :: delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                       avN, & ! average foliar N (gN/m2)
                      iWUE, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
                  NUE_mean, & !
               NUE_optimum, & ! NUE but without age correction applied
                       NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                              ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
               pn_max_temp, & ! Maximum temperature for photosynthesis (oC)
               pn_opt_temp, & ! Optimum temperature fpr photosynthesis (oC)
               pn_kurtosis, & ! Kurtosis of photosynthesis temperature response
                        e0, & ! Quantum yield gC/MJ/m2/day PAR
              co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
            co2_comp_point, & ! CO2 at which photosynthesis > 0 (ppm)
                    minlwp, & ! min leaf water potential (MPa)
 max_lai_lwrad_transmitted, & ! Max fraction of LW from sky transmitted by canopy
lai_half_lwrad_transmitted, & ! LAI at which canopy LW transmittance = 50 %
    max_lai_nir_reflection, & ! Max fraction of NIR reflected by canopy
   lai_half_nir_reflection, & ! LAI at which canopy NIR refection = 50 %
    max_lai_par_reflection, & ! Max fraction of PAR refected by canopy
   lai_half_par_reflection, & ! LAI at which canopy PAR reflection = 50 %
   max_lai_par_transmitted, & ! minimum transmittance = 1-par
  lai_half_par_transmitted, & ! LAI at which 50 %
   max_lai_nir_transmitted, & ! minimum transmittance = 1-par
  lai_half_nir_transmitted, & ! LAI at which 50 %
   max_lai_lwrad_reflected, & !
  lai_half_lwrad_reflected, & ! LAI at which 50 % LW is reflected back to sky
     soil_swrad_absorption, & ! Fraction of SW rad absorbed by soil
     max_lai_lwrad_release, & ! 1-Max fraction of LW emitted from canopy to be
    lai_half_lwrad_release    ! LAI at which LW emitted from canopy to be released at 50 %

! Module level variables for step specific met drivers
double precision :: mint, & ! minimum temperature (oC)
                    maxt, & ! maximum temperature (oC)
                   leafT, & ! leaf temperature (oC)
                   swrad, & ! incoming short wave radiation (MJ/m2/day)
                     co2, & ! CO2 (ppm)
                     doy, & ! Day of year
                wind_spd, & ! wind speed (m/s)
                  vpd_pa, & ! Vapour pressure deficit (Pa)
                     lai    ! leaf area index (m2/m2)

! Module level varoables for step specific timing information
double precision :: cos_solar_zenith_angle, &
                    seconds_per_step, & !
                       days_per_step, & !
                     days_per_step_1, & !
                        dayl_seconds, & ! day length in seconds
                      dayl_seconds_1, &
                          dayl_hours    ! day length in hours

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                                meant_time, &
                                       co2_half_saturation, & ! (ppm)
                                    co2_compensation_point    ! (ppm)

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out,NEE_out,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP_out)

    ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
    ! Index using DeltaGSI to determine FOLiar phenology- Forest Rotation (DALEC_GSI_DFOL_FR) model.
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2

    implicit none

    ! declare input variables
    integer, intent(in) :: start &
                          ,finish &
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout) :: lai_out & ! leaf area index
                                               ,GPP_out & ! Gross primary productivity
                                               ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    double precision :: Rtot    ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)
    integer :: p,f,nxp,n,test,m

    ! local fire related variables
    double precision :: burnt_area &
                       ,CFF(7) = dble_zero, CFF_res(4) = dble_zero    & ! combusted and non-combustion fluxes
                       ,NCFF(7) = dble_zero, NCFF_res(4) = dble_zero  & ! with residue and non-residue seperates
                       ,combust_eff(5)                                & ! combustion efficiency
                       ,rfac                                            ! resilience factor

    ! local deforestation related variables
    double precision, dimension(5) :: post_harvest_burn & ! how much burning to occur after
                                     ,foliage_frac_res  &
                                     ,roots_frac_res    &
                                     ,rootcr_frac_res   &
                                     ,stem_frac_res     &
                                     ,branch_frac_res   &
                                     ,Cbranch_part      &
                                     ,Crootcr_part      &
                                     ,soil_loss_frac

    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   &
                       ,wood_pellets,C_total         &
                       ,labile_frac_res              &
                       ,Cstem,Cbranch,Crootcr        &
                       ,stem_residue,branch_residue  &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    ! variables for phenology model update / adjustments
    double precision :: lai_save, &
                  canopy_lw_save, &
                    soil_lw_save, &
                  canopy_sw_save, &
                 canopy_par_save, &
                    soil_sw_save, &
                         gs_save, &
                         ga_save

    integer :: reforest_day, harvest_management, restocking_lag, gsi_lag

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kgH2O.m-2.s-1)
    ! 8th deforestation fraction
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature (oC)
    ! 11th 21 day average photoperiod (seconds)
    ! 12th 21 day average VPD (kPa)
    ! 13th Forest management practice to accompany any clearing
    ! 14th avg daily temperature (oC)
    ! 15th avg daily wind speed (m.s-1)
    ! 16th vapour pressure deficit (Pa)

    ! POOLS are:
    ! 1 = labile (p18)
    ! 2 = foliar (p19)
    ! 3 = root   (p20)
    ! 4 = wood   (p21)
    ! 5 = litter (p22)
    ! 6 = som    (p23)
    ! 7 = cwd    (p37)
    ! 8 = soil water content (currently assumed to field capacity)

    ! p(30) = labile replanting
    ! p(31) = foliar replanting
    ! p(32) = fine root replanting
    ! p(33) = wood replanting

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile production
    ! 9 = leaffall factor
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som
    ! 16 = labrelease factor
    ! 17 = carbon flux due to fire
    ! 18 = growing season index
    ! 19 = CWD turnover to litter
    ! 20 = C extracted as harvest
    ! 21 = labile loss due to disturbance
    ! 22 = foliage loss due to disturbance
    ! 23 = root loss due to disturbance
    ! 24 = wood loss due to disturbance

    ! PARAMETERS
    ! 31 process parameters; 7 C pool initial conditions

    ! p(1) = Litter to SOM conversion rate (fraction)
    ! p(2) = RmGPP fraction
    ! p(3) = GSI sensitivity for leaf growth
    ! p(4) = Max labile turnover to roots (fraction)
    ! p(5) = Max leaf turnover (GSI; fraction)
    ! p(6) = Turnover rate of wood (fraction)
    ! p(7) = Turnover rate of roots (fraction)
    ! p(8) = Litter turnover rate to heterotrophic respiration (fraction)
    ! p(9) = SOM turnover rate to heterotrophic respiration (fraction)
    ! p(10) = Exponential coefficient for temperature response for heterotrophic respiration
    ! p(11) = Average foliar nitrogen content (log10(gN/m2))
    ! p(12) = Max labile turnover to leaves (GSI; fraction)
    ! p(13) = Max labile turnover to wood (fraction)
    ! p(14) = Min temp threshold (GSI; Kelvin)
    ! p(15) = Max temp threshold (GSI; Kelvin)
    ! p(16) = Min photoperiod threshold (GSI; seconds)
    ! p(17) = Leaf Mass per unit Area (gC.m-2)
    ! p(18) = Initial labile pool (gC/m2)
    ! p(19) = Initial foliage pool (gC/m2)
    ! p(20) = Initial root pool (gC/m2)
    ! p(21) = Initial wood pool (gC/m2)
    ! p(22) = Initial litter pool (gC/m2)
    ! p(23) = Initial som pool (gC/m2)
    ! p(24) = Max photoperiod threshold (GSI; seconds)
    ! p(25) = Min vapour pressure deficit threshold (GSI; Pa)
    ! p(26) = Max vapour pressure deficit (GSI; Pa)
    ! p(27) = Fraction return threshold on GPP for LAI growth
    ! p(28) = Fraction of Cwood which is Cbranch
    ! p(29) = Fraction of Cwood which is Ccoarseroot
    ! p(37) = Initial CWD pool (gC/m2)
    ! p(38) = CWD turnover fraction (fraction)
    ! p(39) = Nitrogen use efficiency (gC/gN/day)

    ! variables related to deforestation
    ! labile_loss = total loss from labile pool from deforestation
    ! foliar_loss = total loss form foliar pool from deforestation
    ! roots_loss = total loss from root pool from deforestation
    ! wood_loss = total loss from wood pool from deforestation
    ! labile_residue = harvested labile remaining in system as residue
    ! foliar_residue = harested foliar remaining in system as residue
    ! roots_residue = harvested roots remaining in system as residue
    ! wood_residue = harvested wood remaining in system as residue
    ! coarse_root_residue = expected coarse woody root left in system as residue

    ! parameters related to deforestation
    ! labile_frac_res = fraction of labile harvest left as residue
    ! foliage_frac_res = fraction of foliage harvest left as residue
    ! roots_frac_res = fraction of roots harvest left as residue
    ! wood_frac_res = fraction of wood harvest left as residue
    ! Crootcr_part = fraction of wood pool expected to be coarse root
    ! Crootcr_frac_res = fraction of coarse root left as residue
    ! soil_loss_frac = fraction determining Csom expected to be physically
    ! removed along with coarse roots

! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0
!real :: Rtot_times=0, aero_time=0 , soilwater_time=0 , acm_et_time = 0
!call cpu_time(begin)
!call cpu_time(finish)

    ! load ACM-GPP-ET parameters
    NUE                        = pars(39)      ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                               ! ,unlimited by CO2, light and
                                               ! photoperiod (gC/gN/m2leaf/day)
    pn_max_temp                = 5.357174d+01  ! Maximum temperature for photosynthesis (oC)
    pn_opt_temp                = 3.137242d+01  ! Optimum temperature for photosynthesis (oC)
    pn_kurtosis                = 1.927458d-01  ! Kurtosis of photosynthesis temperature response
    e0                         = 5.875662d+00  ! Quantum yield gC/MJ/m2/day PAR
    max_lai_lwrad_transmitted  = 7.626683d-01  ! Max fractional reduction of LW from sky transmitted through canopy
    lai_half_lwrad_transmitted = 7.160363d-01  ! LAI at which canopy LW transmittance reduction = 50 %
    max_lai_nir_reflection     = 4.634860d-01  ! Max fraction of NIR reflected by canopy
    lai_half_nir_reflection    = 1.559148d+00  ! LAI at which canopy NIR reflected = 50 %
    minlwp                     =-1.996830d+00  ! minimum leaf water potential (MPa)
    max_lai_par_reflection     = 1.623013d-01  ! Max fraction of PAR reflected by canopy
    lai_half_par_reflection    = 1.114360d+00  ! LAI at which canopy PAR reflected = 50 %
    lai_half_lwrad_reflected   = 1.126214d+00  ! LAI at which 50 % LW is reflected back to sky
    iWUE                       = 1.602503d-06  ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
    soil_swrad_absorption      = 6.643079d-01  ! Fraction of SW rad absorbed by soil
    max_lai_par_transmitted    = 8.079519d-01  ! Max fractional reduction in PAR transmittance by canopy
    lai_half_par_transmitted   = 9.178784d-01  ! LAI at which PAR transmittance reduction = 50 %
    max_lai_nir_transmitted    = 8.289803d-01  ! Max fractional reduction in NIR transmittance by canopy
    lai_half_nir_transmitted   = 1.961831d+00  ! LAI at which NIR transmittance reduction = 50 %
    max_lai_lwrad_release      = 9.852855d-01  ! Max fraction of LW emitted (1-par) from canopy to be released
    lai_half_lwrad_release     = 7.535450d-01  ! LAI at which LW emitted from canopy to be released at 50 %
    max_lai_lwrad_reflected    = 1.955832d-02  ! LAI at which 50 % LW is reflected back to sky

    ! load some values
    avN = 10d0**pars(11)  ! foliar N gN/m2

    ! initial values for deforestation variables
    labile_loss = dble_zero    ; foliar_loss = dble_zero
    roots_loss = dble_zero     ; wood_loss = dble_zero
    labile_residue = dble_zero ; foliar_residue = dble_zero
    roots_residue = dble_zero  ; wood_residue = dble_zero
    stem_residue = dble_zero   ; branch_residue = dble_zero
    reforest_day = 0
    soil_loss_with_roots = dble_zero
    coarse_root_residue = dble_zero
    post_harvest_burn = dble_zero

    ! now load the hardcoded forest management parameters into their locations

    ! Parameter values for deforestation variables
    ! scenario 1
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = dble_one
    roots_frac_res(1)   = dble_one
    rootcr_frac_res(1) = dble_one
    branch_frac_res(1) = dble_one
    stem_frac_res(1)   = dble_zero !
    ! wood partitioning (fraction)
    Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(1) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(1) = dble_one

    !## scen 2
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = dble_one
    roots_frac_res(2)   = dble_one
    rootcr_frac_res(2) = dble_one
    branch_frac_res(2) = dble_one
    stem_frac_res(2)   = dble_zero !
    ! wood partitioning (fraction)
    Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(2) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(2) = dble_zero

    !## scen 3
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5d0
    roots_frac_res(3)   = dble_one
    rootcr_frac_res(3) = dble_one
    branch_frac_res(3) = dble_zero
    stem_frac_res(3)   = dble_zero !
    ! wood partitioning (fraction)
    Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(3) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(3) = dble_zero

    !## scen 4
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5d0
    roots_frac_res(4)   = dble_one
    rootcr_frac_res(4) = dble_zero
    branch_frac_res(4) = dble_zero
    stem_frac_res(4)   = dble_zero
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(4) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(4) = dble_zero

    !## scen 5 (grassland grazing / cutting)
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(5) = 0.1d0
    roots_frac_res(5)   = dble_zero
    rootcr_frac_res(5)  = dble_zero
    branch_frac_res(5)  = 0.1d0
    stem_frac_res(5)    = 0.1d0
    ! wood partitioning (fraction)
    Crootcr_part(5) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(5) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(5) = dble_zero ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(5) = dble_zero

    ! for the moment override all paritioning parameters with those coming from
    ! CARDAMOM
    Cbranch_part = pars(28)
    Crootcr_part = pars(29)

    ! declare fire constants (labile, foliar, roots, wood, litter)
    combust_eff(1) = 0.1d0 ; combust_eff(2) = 0.9d0
    combust_eff(3) = 0.1d0 ; combust_eff(4) = 0.5d0
    combust_eff(5) = 0.3d0 ; rfac = 0.5d0

    ! assigning initial conditions
    POOLS(1,1)=pars(18)
    POOLS(1,2)=pars(19)
    POOLS(1,3)=pars(20)
    POOLS(1,4)=pars(21)
    POOLS(1,5)=pars(22)
    POOLS(1,6)=pars(23)
    POOLS(1,7)=pars(37)

    if (.not.allocated(disturbance_residue_to_som)) then
        allocate(disturbance_residue_to_litter(nodays), &
                 disturbance_residue_to_som(nodays),    &
                 disturbance_residue_to_cwd(nodays),    &
                 disturbance_loss_from_litter(nodays),  &
                 disturbance_loss_from_som(nodays),     &
                 disturbance_loss_from_cwd(nodays))
    endif
    disturbance_residue_to_litter = dble_zero ; disturbance_loss_from_litter = dble_zero
    disturbance_residue_to_som = dble_zero ; disturbance_loss_from_som = dble_zero
    disturbance_residue_to_cwd = dble_zero ; disturbance_loss_from_cwd = dble_zero

    ! calculate some values once as these are invarient between DALEC runs
    if (.not.allocated(tmp_x)) then
        ! 21 days is the maximum potential so we will fill the maximum potential
        ! + 1 for safety
        allocate(tmp_x(22),tmp_m(nodays))
        do f = 1, 22
           tmp_x(f) = f
        end do
        do n = 1, nodays
          ! calculate the gradient / trend of GSI
          if (sum(deltat(1:n)) < 21d0) then
              tmp_m(n) = n-1
          else
             ! else we will try and work out the gradient to see what is
             ! happening
             ! to the system over all. The default assumption will be to
             ! consider
             ! the averaging period of GSI model (i.e. 21 days). If this is not
             ! possible either the time step of the system is used (if step
             ! greater
             ! than 21 days) or all available steps (if n < 21).
             m = 0 ; test = 0
             do while (test < 21d0)
                m = m + 1 ; test = nint(sum(deltat((n-m):n)))
                if (m > (n-1)) test = 21
             end do
             tmp_m(n) = m
          endif ! for calculating gradient
        end do ! calc daily values once
        ! allocate GSI history dimension
        gsi_lag_remembered=nint(max(2d0,maxval(tmp_m)))
    end if ! .not.allocated(tmp_x)
    ! assign our starting value
    gsi_history = pars(36)-dble_one
     just_grown = pars(35)

    ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
    if (.not.allocated(deltat_1)) then
       allocate(deltat_1(nodays),meant_time(nodays), &
                co2_half_saturation(nodays),co2_compensation_point(nodays))
       ! inverse of time step (days-1) to avoid divisions
       deltat_1 = deltat**(-dble_one)
       ! meant time step temperature
       meant_time = (met(2,1:nodays)+met(3,1:nodays)) * 0.5d0
       do n = 1, nodays
          ! Temperature adjustments for Michaelis-Menten coefficients
          ! for CO2 (kc) and O2 (ko) and CO2 compensation point.
          co2_compensation_point(n) = arrhenious(co2comp_saturation,co2comp_half_sat_conc,met(3,n))
          co2_half_saturation(n) = arrhenious(kc_saturation,kc_half_sat_conc,met(3,n))
       end do

    else

    endif ! has SWP already been determined?

    ! initialise root reach based on initial conditions
    root_biomass = max(min_root,POOLS(1,3)*2d0)
    ! needed to initialise soils
    call calculate_Rtot(Rtot)

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit=pars(34)-1d0
    lab_turn_crit=pars(3)-1d0

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      leafT = maxt
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = met(6,n)   ! Day of year
      meant = meant_time(n)  ! mean air temperature (oC)
      meant_K = meant + freeze
      wind_spd = met(15,n) ! wind speed (m/s)
      vpd_pa = met(16,n)  ! Vapour pressure deficit (Pa)

      ! states needed for module variables
      lai_out(n) = POOLS(n,2)/pars(17)
      lai = lai_out(n) ! leaf area index (m2/m2)

      ! Temperature adjustments for Michaelis-Menten coefficients
      ! for CO2 (kc) and O2 (ko) and CO2 compensation point
      ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
      co2_half_sat   = co2_half_saturation(n)
      co2_comp_point = co2_compensation_point(n)

      ! extract timing related values
      call calculate_daylength((doy-(deltat(n)*0.5d0)),lat)
      dayl_seconds_1 = dayl_seconds ** (-dble_one)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      !!!!!!!!!!
      ! calculate soil water potential and total hydraulic resistance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2d0)
      call calculate_Rtot(Rtot)

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

      ! calculate some temperature dependent meteorologial properties
      call acm_meteorological_constants(maxt)
      ! calculate aerodynamic using consistent approach with SPA
      call calculate_aerodynamic_conductance

      !!!!!!!!!!
      ! Determine net shortwave and isothermal longwave energy balance
      !!!!!!!!!!

      call calculate_shortwave_balance
      call calculate_longwave_isothermal(leafT,maxt)

      ! calculate variables used commonly between ACM_GPP and ACM_ET
      call acm_albedo_gc(abs(minlwp),Rtot)

      !!!!!!!!!!
      ! GPP (gC.m-2.day-1)
      !!!!!!!!!!

      if (stomatal_conductance > dble_zero .and. lai > vsmall) then
          FLUXES(n,1) = max(dble_zero,acm_gpp(stomatal_conductance))
      else
          FLUXES(n,1) = dble_zero
      endif

      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*meant)
      ! (maintenance) autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = dble_zero !(FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(13)
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))*pars(4)
      ! wood production
      FLUXES(n,7) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5)-FLUXES(n,6)

      ! GSI added to fortran version by TLS 24/11/2014
      ! /* 25/09/14 - JFE
      ! Here we calculate the Growing Season Index based on
      ! Jolly et al. A generalized, bioclimatic index to predict foliar
      ! phenology in response to climate Global Change Biology, Volume 11, page 619-632,
      ! 2005 (doi: 10.1111/j.1365-2486.2005.00930.x)
      ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
      ! Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res.,
      ! 116, G03020, doi:10.1029/2010JG001545.

      ! It is the product of 3 limiting factors for temperature, photoperiod and
      ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated
      ! min and max value. Photoperiod, VPD and avgTmin are direct input

      ! temperature limitation, then restrict to 0-1; correction for k-> oC
      Tfac = (met(10,n)-(pars(14)-freeze)) / (pars(15)-pars(14))
      Tfac = min(dble_one,max(dble_zero,Tfac))
      ! photoperiod limitation
      Photofac = (met(11,n)-pars(16)) / (pars(24)-pars(16))
      Photofac = min(dble_one,max(dble_zero,Photofac))
      ! VPD limitation
      VPDfac = 1d0 - ( (met(12,n)-pars(25)) / (pars(26)-pars(25)) )
      VPDfac = min(dble_one,max(dble_zero,VPDfac))

      ! calculate and store the GSI index
      FLUXES(n,18) = Tfac*Photofac*VPDfac

      ! we will load up some needed variables
      m = tmp_m(n)
      ! update gsi_history for the calculation
      if (n == 1) then
          ! in first step only we want to take the initial GSI value only
          gsi_history(gsi_lag) = FLUXES(n,18)
      else
          gsi_history((gsi_lag-m):gsi_lag) = FLUXES((n-m):n,18)
      endif
      ! calculate gradient
      gradient = linear_model_gradient(tmp_x(1:(gsi_lag)),gsi_history(1:gsi_lag),gsi_lag)
      ! adjust gradient to daily rate
      gradient = gradient /  nint((sum(deltat((n-m+1):n))) / (gsi_lag-1))
      gsi_lag_remembered = gsi_lag

      ! first assume that nothing is happening
      FLUXES(n,9) = dble_zero  ! leaf turnover
      FLUXES(n,16) = dble_zero ! leaf growth

      ! save original values for re-allocation later
      canopy_lw_save = canopy_lwrad_Wm2 ; soil_lw_save  = soil_lwrad_Wm2
      canopy_sw_save = canopy_swrad_MJday ; canopy_par_save  = canopy_par_MJday
      soil_sw_save = soil_swrad_MJday ; gs_save = stomatal_conductance
      ga_save = aerodynamic_conductance ; lai_save = lai

      ! now update foliage and labile conditions based on gradient calculations
      if (gradient < fol_turn_crit .or. FLUXES(n,18) == dble_zero) then
         ! we are in a decending condition so foliar turnover
         FLUXES(n,9) = pars(5)*(dble_one-FLUXES(n,18))
         just_grown = 0.5d0
      else if (gradient > lab_turn_crit) then
         ! we are in a assending condition so labile turnover
         FLUXES(n,16) = pars(12)*FLUXES(n,18)
         just_grown = 1.5d0
         ! check carbon return
         tmp = POOLS(n,1)*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))/deltat(n)
         lai = (POOLS(n,2)+tmp)/pars(17)
         tmp = lai / lai_save
         aerodynamic_conductance = aerodynamic_conductance * tmp
         stomatal_conductance = stomatal_conductance * tmp
         call calculate_shortwave_balance
         if (lai_save < vsmall) then
             call calculate_aerodynamic_conductance
             call acm_albedo_gc(abs(minlwp),Rtot)
         endif
         tmp = max(dble_zero,acm_gpp(stomatal_conductance))
         ! determine if increase in LAI leads to an improvement in GPP greater
         ! than critical value, if not then no labile turnover allowed
         if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
             FLUXES(n,16) = dble_zero
         endif
      else
         ! probably we want nothing to happen, however if we are at the seasonal
         ! maximum we will consider further growth still
         if (just_grown >= dble_one) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            FLUXES(n,16) = pars(12)*FLUXES(n,18)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = POOLS(n,1)*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))/deltat(n)
            lai = (POOLS(n,2)+tmp)/pars(17)
            tmp = lai / lai_save
            aerodynamic_conductance = aerodynamic_conductance * tmp
            stomatal_conductance = stomatal_conductance * tmp
            call calculate_shortwave_balance
            if (lai_save < vsmall) then
                call calculate_aerodynamic_conductance
                call acm_albedo_gc(abs(minlwp),Rtot)
            endif
            tmp = max(dble_zero,acm_gpp(stomatal_conductance))
            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than critical value, if not then no labile turnover allowed
            if ( (tmp - FLUXES(n,1)) < (pars(27)*FLUXES(n,1)) ) then
                FLUXES(n,16) = dble_zero
            endif

         end if ! Just grown?
      endif ! gradient choice

      ! restore original value back from memory
      lai = lai_save
      !canopy_lwrad_Wm2 = canopy_lw_save ; soil_lwrad_Wm2 = soil_lw_save
      canopy_swrad_MJday = canopy_sw_save ; canopy_par_MJday = canopy_par_save
      soil_swrad_MJday = soil_sw_save ; stomatal_conductance = gs_save
      aerodynamic_conductance = ga_save

      ! these allocated if post-processing
      if (allocated(itemp)) then
         itemp(n) = Tfac
         ivpd(n) = VPDfac
         iphoto(n) = Photofac
      endif

      !
      ! those with time dependancies
      !

      ! total labile release
      FLUXES(n,8)  = POOLS(n,1)*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*min(dble_one,dble_one-(dble_one-FLUXES(n,9))**deltat(n))/deltat(n)
      ! total wood litter production
      FLUXES(n,11) = POOLS(n,4)*min(dble_one,dble_one-(dble_one-pars(6))**deltat(n))/deltat(n)
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*min(dble_one,dble_one-(dble_one-pars(7))**deltat(n))/deltat(n)

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(8)))**deltat(n))/deltat(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(9)))**deltat(n))/deltat(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(1)))**deltat(n))/deltat(n)
      ! CWD to litter
      FLUXES(n,20) = POOLS(n,7)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(38)))**deltat(n))/deltat(n)

      ! mass balance check for decomposition of litter as occurs via two
      ! pathways
      if ( (FLUXES(n,13) + FLUXES(n,15))*deltat(n) >= POOLS(n,5) ) then
           tmp = FLUXES(n,13) / (FLUXES(n,13) + FLUXES(n,15))
           FLUXES(n,13) = tmp * (POOLS(n,5)/deltat(n))
           FLUXES(n,15) = (dble_one - tmp) * (POOLS(n,5)/deltat(n))
      end if

      ! calculate growth respiration and adjust allocation to pools assuming
      ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
      ! foliage
      Rg_from_labile = FLUXES(n,8)*0.21875d0 ; FLUXES(n,8) = FLUXES(n,8) * 0.78125d0
      ! now update the Ra flux
      FLUXES(n,3) = FLUXES(n,3) + Rg_from_labile
      ! roots
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,6)*0.21875d0) ; FLUXES(n,6) = FLUXES(n,6) * 0.78125d0
      ! wood
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,7)*0.21875d0) ; FLUXES(n,7) = FLUXES(n,7) * 0.78125d0

      ! calculate the NEE
      NEE_out(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP_out(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8)-Rg_from_labile)*deltat(n)
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10)+FLUXES(n,8))*deltat(n)
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)+FLUXES(n,20)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14))*deltat(n)
      ! cwd pool
      POOLS(n+1,7) = POOLS(n,7) + (FLUXES(n,11)-FLUXES(n,20))*deltat(n)

      !!!!!!!!!!
      ! deal first with deforestation
      !!!!!!!!!!

      if (n == reforest_day) then
          POOLS(n+1,1) = pars(30) ! labile
          POOLS(n+1,2) = pars(31) ! foliar
          POOLS(n+1,3) = pars(32) ! roots
          POOLS(n+1,4) = pars(33) ! wood
      endif

      ! reset values
      FLUXES(n,17) = dble_zero ; FLUXES(n,21:25) = dble_zero
      harvest_management = 0 ; burnt_area = dble_zero

      if (met(8,n) > dble_zero) then

          ! pass harvest management to local integer
          harvest_management = int(met(13,n))

          ! assume that labile is proportionally distributed through the plant
          ! root and wood and therefore so is the residual fraction
          C_total = POOLS(n+1,3) + POOLS(n+1,4)
          ! partition wood into its components
          Cbranch = POOLS(n+1,4)*Cbranch_part(harvest_management)
          Crootcr = POOLS(n+1,4)*Crootcr_part(harvest_management)
          Cstem   = POOLS(n+1,4)-(Cbranch + Crootcr)
          ! now calculate the labile fraction of residue
          if (C_total > dble_zero) then
              labile_frac_res = ((POOLS(n+1,3)/C_total) * roots_frac_res(harvest_management)  ) &
                              + ((Cbranch/C_total)      * branch_frac_res(harvest_management) ) &
                              + ((Cstem/C_total)        * stem_frac_res(harvest_management)   ) &
                              + ((Crootcr/C_total)      * rootcr_frac_res(harvest_management) )
          else
              labile_frac_res = dble_zero
          endif

          ! you can't remove any biomass if there is none left...
          if (C_total > vsmall) then

              ! Loss of carbon from each pools
              labile_loss = POOLS(n+1,1)*met(8,n)
              foliar_loss = POOLS(n+1,2)*met(8,n)
              ! roots are not removed under grazing
              if (harvest_management /= 5) then
                 roots_loss = POOLS(n+1,3)*met(8,n)
              else
                 roots_loss = dble_zero
              endif
              wood_loss   = (Cbranch+Crootcr+Cstem)*met(8,n)
              ! estimate labile loss explicitly from the loss of their storage
              ! tissues
              labile_loss = POOLS(n+1,1) * ((roots_loss+wood_loss) / (POOLS(n+1,3)+POOLS(n+1,4)))

              ! For output / EDC updates
              if (met(8,n) <= 0.99d0) then
                  FLUXES(n,22) = labile_loss * deltat_1(n)
                  FLUXES(n,23) = foliar_loss * deltat_1(n)
                  FLUXES(n,24) = roots_loss * deltat_1(n)
                  FLUXES(n,25) = wood_loss * deltat_1(n)
              endif
              ! Transfer fraction of harvest waste to litter or som pools
              ! easy pools first
              labile_residue = labile_loss*labile_frac_res
              foliar_residue = foliar_loss*foliage_frac_res(harvest_management)
              roots_residue  = roots_loss*roots_frac_res(harvest_management)
              ! Explicit calculation of the residues from each fraction
              coarse_root_residue  = Crootcr*met(8,n)*rootcr_frac_res(harvest_management)
              branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
              stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
              ! Now finally calculate the final wood residue
              wood_residue = stem_residue + branch_residue + coarse_root_residue
              ! Mechanical loss of Csom due to coarse root extraction
              soil_loss_with_roots = Crootcr*met(8,n)*(dble_one-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

              ! Update pools
              POOLS(n+1,1) = max(dble_zero,POOLS(n+1,1)-labile_loss)
              POOLS(n+1,2) = max(dble_zero,POOLS(n+1,2)-foliar_loss)
              POOLS(n+1,3) = max(dble_zero,POOLS(n+1,3)-roots_loss)
              POOLS(n+1,4) = max(dble_zero,POOLS(n+1,4)-wood_loss)
              POOLS(n+1,5) = max(dble_zero, POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue))
              POOLS(n+1,6) = max(dble_zero, POOLS(n+1,6) - soil_loss_with_roots)
              POOLS(n+1,7) = max(dble_zero, POOLS(n+1,7) + wood_residue)

              ! Some variable needed for the EDCs
              ! reallocation fluxes for the residues
              disturbance_residue_to_litter(n) = labile_residue+foliar_residue+roots_residue
              disturbance_loss_from_litter(n)  = dble_zero
              disturbance_residue_to_cwd(n)    = wood_residue
              disturbance_loss_from_cwd(n)     = dble_zero
              disturbance_residue_to_som(n)    = dble_zero
              disturbance_loss_from_som(n)     = soil_loss_with_roots
              ! Convert all to rates to be consistent with the FLUXES in EDCs
              disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) * deltat_1(n)
              disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) * deltat_1(n)
              disturbance_residue_to_cwd(n)    = disturbance_residue_to_cwd(n) * deltat_1(n)
              disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) * deltat_1(n)
              disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) * deltat_1(n)
              disturbance_loss_from_som(n)     = disturbance_loss_from_som(n) * deltat_1(n)
              ! estimate total C extraction
              ! NOTE: this calculation format is to prevent precision error in calculation
              FLUXES(n,21) = wood_loss + labile_loss + foliar_loss + roots_loss
              FLUXES(n,21) = FLUXES(n,21) - (wood_residue + labile_residue + foliar_residue + roots_residue)
              ! Convert to daily rate
              FLUXES(n,21) = FLUXES(n,21) * deltat_1(n)

          end if ! C_total > vsmall

          ! Total carbon loss from the system
          C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                  - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

          ! If total clearance occured then we need to ensure some minimum
          ! values and reforestation is assumed one year forward
          if (met(8,n) > 0.99d0) then
              m = 0 ; test = nint(sum(deltat(n:(n+m))))
              ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
              restocking_lag = 365*2
              do while (test < restocking_lag)
                 m = m + 1 ; test = nint(sum(deltat(n:(n+m))))
                 !  get out clause for hitting the end of the simulation
                 if (m+n >= nodays) test = restocking_lag
              enddo
              reforest_day = min((n+m), nodays)
          endif ! if total clearance

      endif ! end deforestation info

      !!!!!!!!!!
      ! then deal with fire
      !!!!!!!!!!

      if (met(9,n) > dble_zero .or.(met(8,n) > dble_zero .and. harvest_management > 0)) then

          burnt_area = met(9,n)
          if (met(8,n) > dble_zero .and. burnt_area > dble_zero) then
              ! pass harvest management to local integer
              burnt_area = min(dble_one,burnt_area + post_harvest_burn(harvest_management))
          else if (met(8,n) > dble_zero .and. burnt_area <= dble_zero) then
              burnt_area = post_harvest_burn(harvest_management)
          endif

          if (burnt_area > dble_zero) then

              !/*first fluxes*/
              !/*LABILE*/
              CFF(1) = POOLS(n+1,1)*burnt_area*combust_eff(1)
              NCFF(1) = POOLS(n+1,1)*burnt_area*(dble_one-combust_eff(1))*(dble_one-rfac)
              !/*foliar*/
              CFF(2) = POOLS(n+1,2)*burnt_area*combust_eff(2)
              NCFF(2) = POOLS(n+1,2)*burnt_area*(dble_one-combust_eff(2))*(dble_one-rfac)
              !/*root*/
              CFF(3) = dble_zero !POOLS(n+1,3)*burnt_area*combust_eff(3)
              NCFF(3) = dble_zero !POOLS(n+1,3)*burnt_area*(dble_one-combust_eff(3))*(dble_one-rfac)
              !/*wood*/
              CFF(4) = POOLS(n+1,4)*burnt_area*combust_eff(4)
              NCFF(4) = POOLS(n+1,4)*burnt_area*(dble_one-combust_eff(4))*(dble_one-rfac)
              !/*litter*/
              CFF(5) = POOLS(n+1,5)*burnt_area*combust_eff(5)
              NCFF(5) = POOLS(n+1,5)*burnt_area*(dble_one-combust_eff(5))*(dble_one-rfac)
              ! CWD; assume same as live wood (should be improved later)
              CFF(7) = POOLS(n+1,7)*burnt_area*combust_eff(4)
              NCFF(7) = POOLS(n+1,7)*burnt_area*(dble_one-combust_eff(4))*(dble_one-rfac)
              !/*fires as daily averages to comply with units*/
              FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)) * deltat_1(n)
              !/*update net exchangep*/
              NEE_out(n) = NEE_out(n) + FLUXES(n,17)
              ! determine the as daily rate impact on live tissues for use in EDC and
              ! MTT calculations
              FLUXES(n,22) = FLUXES(n,22) + ((CFF(1) + NCFF(1)) * deltat_1(n)) ! labile
              FLUXES(n,23) = FLUXES(n,23) + ((CFF(2) + NCFF(2)) * deltat_1(n)) ! foliar
              FLUXES(n,24) = FLUXES(n,24) + ((CFF(3) + NCFF(3)) * deltat_1(n)) ! root
              FLUXES(n,25) = FLUXES(n,25) + ((CFF(4) + NCFF(4)) * deltat_1(n)) ! wood

              !// update pools
              !/*Adding all fire pool transfers here*/
              POOLS(n+1,1) = max(dble_zero,POOLS(n+1,1)-CFF(1)-NCFF(1))
              POOLS(n+1,2) = max(dble_zero,POOLS(n+1,2)-CFF(2)-NCFF(2))
              POOLS(n+1,3) = max(dble_zero,POOLS(n+1,3)-CFF(3)-NCFF(3))
              POOLS(n+1,4) = max(dble_zero,POOLS(n+1,4)-CFF(4)-NCFF(4))
              POOLS(n+1,5) = max(dble_zero,POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3))
              POOLS(n+1,6) = max(dble_zero,POOLS(n+1,6)+NCFF(4)+NCFF(5)+NCFF(7))
              POOLS(n+1,7) = max(dble_zero,POOLS(n+1,7)-CFF(7)-NCFF(7))

              ! some variable needed for the EDCs
              ! reallocation fluxes for the residues
              disturbance_residue_to_litter(n) = (NCFF(1)+NCFF(2)+NCFF(3))
              disturbance_residue_to_som(n)    = (NCFF(4)+NCFF(5)+NCFF(7))
              disturbance_loss_from_litter(n)  = CFF(5)+NCFF(5)
              disturbance_loss_from_cwd(n)     = CFF(7) - NCFF(7)
              ! convert to daily rate for consistency with the EDCs
              disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n)  * deltat_1(n)
              disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) * deltat_1(n)
              disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) * deltat_1(n)
              disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) * deltat_1(n)

          endif ! burn area > 0

      endif ! fire activity

      do nxp = 1, nopools
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < dble_zero) then
            print*,"step",n,"POOL",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"POOLS+1",POOLS(n+1,:)
            stop
         endif
      enddo

    end do ! nodays loop
!call cpu_time(done)
!print*,done-begin
  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp(gs)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pn, pd, pp, qq, ci, mult, pl &
                       ,gc ,gs_mol, gb_mol

    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
!    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,leafT)
!    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,leafT)

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    pn = lai*avN*NUE*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,leafT)

    !
    ! Diffusion limited photosynthesis
    !

    ! daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2
    gs_mol = gs * 1d-3 * seconds_per_day * gs_H2O_CO2
    ! canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance * seconds_per_day * convert_ms1_mol_1 * gb_H2O_CO2
    ! Combining in series the stomatal and boundary layer conductances
    gc = (gs_mol ** (-dble_one) + gb_mol ** (-dble_one)) ** (-dble_one)

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = (pn*gC_to_umol)/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))
    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    pd = (gc * (co2-ci)) * umol_to_gC
    ! scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    pd = pd * dayl_hours * 0.04166667d0

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    pl = e0 * canopy_par_MJday

    !
    ! CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis
    acm_gpp = pl*pd/(pl+pd)
    ! sanity check
    if (acm_gpp /= acm_gpp) acm_gpp = dble_zero
    ! don't forget to return
    return

  end function acm_gpp
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE (iWUE).

    ! Arguments
    double precision, intent(in) :: gs_in

    ! Local variables
    double precision :: gs_high, gpp_high, gpp_low

    ! Estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! Estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE - ((gpp_high - gpp_low)/lai)

  end function find_gs_iWUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance

    !
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    !

    implicit none

    ! local variables
    double precision :: mixing_length_momentum, & ! mixing length parameter for momentum (m)
                         length_scale_momentum    ! length scale parameter for momentum (m)

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
!    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh
    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    if (lai > min_lai) then
        length_scale_momentum = (4d0*canopy_height) / lai
        mixing_length_momentum = max(canopy_height*0.02d0, 2d0*(ustar_Uh**3)*length_scale_momentum)
    else
        length_scale_momentum = vonkarman * tower_height
        mixing_length_momentum = canopy_height * vonkarman
    endif

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! now we are interested in the within canopy wind speed,
    ! here we assume that the wind speed just inside of the canopy is most important.
    canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*0.75d0)-canopy_height))/mixing_length_momentum)

    ! calculate leaf level conductance (m/s) for water vapour under forced convective conditions
    call average_leaf_conductance(aerodynamic_conductance)

  end subroutine calculate_aerodynamic_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine average_leaf_conductance(gv_forced)

    !
    ! Subroutine calculates the forced conductance of water vapour for non-cylinder within canopy leaves (i.e. broadleaf)
    ! Free convection (i.e. that driven by energy balance) is negelected here due to the lack of an energy balance
    ! calculation in DALEC. Should a energy balance be added then this code could be expanded include free conductance
    ! Follows a simplified approach to that used in SPA (Smallman et al 2013).
    !

    implicit none

    ! arguments
    double precision, intent(out) :: gv_forced ! canopy conductance (m/s) for water vapour under forced convection

    ! local parameters
    double precision, parameter :: leaf_width = 0.08d0   & ! leaf width (m)
                                          ,Pr = 0.72d0     ! Prandtl number
    ! local variables
    double precision :: Dwv & ! Diffusion coefficient of water in air (m2.s-1); air temperature and pressure dependant
                              ! variables for the more advanced
                              ! boundary conditions
            ,nusselt_forced & ! Nusselt value under forced convection
         ,dynamic_viscosity & ! dynamic viscosity (kg.m-2.s-1)
       ,kinematic_viscosity & ! kinematic viscosity (m2.s-1)
                 ,Sh_forced & ! Sherwood number under forced convection
                        ,Re   ! Reynolds number

    ! Determine diffusion coefficient (m2s-1), temperature dependant (pressure dependence neglected). Jones p51;
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    Dwv = 0.0000242d0*(((maxt+freeze)/293.15d0)**1.75d0)
    ! Calculate the dynamic viscosity of air
    dynamic_viscosity = (((maxt+freeze)**1.5d0)/((maxt+freeze)+120d0))*1.4963d-6
    kinematic_viscosity = dynamic_viscosity/air_density_kg
    Re = (leaf_width*canopy_wind)/kinematic_viscosity
    ! calculate nusselt value under forced convection conditions
    nusselt_forced = (1.18d0*(Pr**(0.33d0))*(sqrt(Re)))
    ! update specific Sherwood numbers
    Sh_forced = 0.962d0*nusselt_forced
    ! This is the forced conductance of water vapour for the current leaf
    gv_forced = ((Dwv*Sh_forced)/leaf_width)*0.5d0
    ! apply lai correction
    gv_forced = gv_forced * lai

  end subroutine average_leaf_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine acm_albedo_gc(deltaWP,Rtot)

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! arguments
    double precision, intent(in) :: deltaWP, & ! minlwp-wSWP (MPa)
                                       Rtot    ! total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    ! local variables
    double precision :: denom, isothermal, deltaTemp, deltaR
    double precision, parameter :: max_gs = 2500d0, & ! mmolH2O.m-2.s-1
                                   min_gs = 0.0001d0, & !
                                   tol_gs = 4d0!4d0       !

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (deltaWP > vsmall) then
       ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
       max_supply = (deltaWP/Rtot) * seconds_per_day
    else
       ! set minimum (computer) precision level flow
       max_supply = vsmall
    end if

    if (lai > vsmall .and. aerodynamic_conductance > vsmall) then

        ! there is lai therefore we have have stomatal conductance

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * ((canopy_swrad_MJday * 1d6 * dayl_seconds_1) + canopy_lwrad_Wm2) &
              + (air_density_kg * cpair * vpd_pa * 1d-3 * aerodynamic_conductance)
        denom = (denom / (lambda * max_supply * mmol_to_kg_water * dayl_seconds_1)) - slope
        denom = denom / psych
        stomatal_conductance = aerodynamic_conductance / denom

        ! convert m.s-1 to mmolH2O.m-2.s-1
        stomatal_conductance = stomatal_conductance * 1d3 * convert_ms1_mol_1
        ! if conditions are dew forming then set conductance to maximum as we are not going to be limited by water demand
        if (stomatal_conductance <= dble_zero .or. stomatal_conductance > max_gs) stomatal_conductance = max_gs

        ! if we are potentially limited by stomatal conductance or we are using instrinsic water use efficiency (rather than WUE)
        ! then iterate to find optimum gs otherwise just go with the max...
        if (stomatal_conductance /= max_gs) then
            ! If there is a positive demand for water then we will solve for photosynthesis limits on gs through iterative solution
            delta_gs = 1d-3*lai ! mmolH2O/m2leaf/day
            ! intrinsic WUE optimisation
            stomatal_conductance = zbrent('acm_albedo_gc:find_gs_iWUE',find_gs_iWUE,min_gs,stomatal_conductance,tol_gs)
        end if

    else

        ! if no LAI then there can be no stomatal conductance
        stomatal_conductance = dble_zero

    endif ! if LAI > vsmall

  end subroutine acm_albedo_gc
  !
  !------------------------------------------------------------------
  !
  subroutine acm_meteorological_constants(input_temperature)

    ! Determine some multiple use constants

    implicit none

    ! arguments
    double precision, intent(in) :: input_temperature

    ! local variables
    double precision :: s, mult

    ! Density of air (kg/m3)
    air_density_kg = 353d0/(input_temperature+freeze)
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / ((input_temperature+freeze)*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
!    if (input_temperature < dble_one) then
!        lambda = 2.835d6
!    else
        lambda = 2501000d0-2364d0*input_temperature
!    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*input_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = input_temperature+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    s = 2502.935945d0*exp(17.269d0*input_temperature/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)

  end subroutine acm_meteorological_constants
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal(canopy_temperature,soil_temperature)

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! arguments
    double precision, intent(in) :: canopy_temperature, soil_temperature ! oC

    ! local variables
    double precision :: lwrad, & ! downward long wave radiation from sky (W.m-2)
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy, & ! assuming isothermal condition (W.m-2)
            trans_lw_fraction, &
        reflected_lw_fraction, &
         absorbed_lw_fraction, &
      canopy_release_fraction, & ! fraction of longwave emitted from within the canopy to ultimately be released
   canopy_absorption_from_sky, & ! canopy absorbed radiation from downward LW (W.m-2)
  canopy_absorption_from_soil, & ! canopy absorbed radiation from soil surface (W.m-2)
                  canopy_loss, & ! longwave radiation released from canopy surface (W.m-2).
                                 ! i.e. this value is released from the top and the bottom
     soil_absorption_from_sky, & ! soil absorbed radiation from sky (W.m-2)
  soil_absorption_from_canopy    ! soil absorbed radiation emitted from canopy (W.m-2)

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! calculate fraction of longwave radiation coming from the sky to pentrate to the soil surface
    trans_lw_fraction = dble_one - (max_lai_lwrad_transmitted*lai)/(lai+lai_half_lwrad_transmitted)
    ! calculate the fraction of longwave radiation from sky which is reflected back into the sky
    reflected_lw_fraction = (max_lai_lwrad_reflected*lai) / (lai+lai_half_lwrad_reflected)
    ! calculate absorbed longwave radiation coming from the sky
    absorbed_lw_fraction = dble_one - trans_lw_fraction - reflected_lw_fraction
    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky
    canopy_release_fraction = dble_one - (max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release)

    !!!!!!!!!!
    ! Distribute longwave from sky
    !!!!!!!!!!

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * absorbed_lw_fraction
    ! Long wave absorbed by soil from the sky, soil absorption assumed to be equal to emissivity
    soil_absorption_from_sky = trans_lw_fraction * lwrad * emissivity
    ! Long wave reflected directly back into sky
    sky_lwrad_Wm2 = lwrad * reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! First, calculate longwave radiation coming up from the soil plus the radiation which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (trans_lw_fraction * lwrad * (dble_one-emissivity))
    ! Second, use this total to estimate the longwave returning to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * trans_lw_fraction)
    ! Third, now calculate the longwave from the soil surface absorbed by the canopy
    canopy_absorption_from_soil = canopy_absorption_from_soil * absorbed_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave originating from the canopy itself
    !!!!!!!!!!

    ! calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = longwave_release_canopy * lai * canopy_release_fraction
    ! Calculate longwave absorbed by soil which is released by the canopy itself
    soil_absorption_from_canopy = canopy_loss * emissivity
    ! Canopy released longwave returned to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + canopy_loss

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    canopy_lwrad_Wm2 = (canopy_absorption_from_sky + canopy_absorption_from_soil) - (canopy_loss + canopy_loss)
    ! determine isothermal net soil
    soil_lwrad_Wm2 = (soil_absorption_from_sky + soil_absorption_from_canopy) - longwave_release_soil

  end subroutine calculate_longwave_isothermal
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

    ! local variables
    double precision :: absorbed_nir_fraction_soil &
                       ,absorbed_par_fraction_soil &
                       ,soil_par_MJday             &
                       ,soil_nir_MJday             &
                       ,trans_nir_MJday            &
                       ,trans_par_MJday            &
                       ,canopy_nir_MJday           &
                       ,refl_par_MJday             &
                       ,refl_nir_MJday             &
                       ,reflected_nir_fraction     & !
                       ,reflected_par_fraction     & !
                       ,absorbed_nir_fraction      & !
                       ,absorbed_par_fraction      & !
                       ,trans_nir_fraction         & !
                       ,trans_par_fraction

    !!!!!!!!!!
    ! Determine canopy absorption / reflectance as function of LAI
    !!!!!!!!!!

    ! Canopy transmitted of PAR & NIR radiation towards the soil
    trans_par_fraction = dble_one - (lai*max_lai_par_transmitted) &
                       / (lai+lai_half_par_transmitted)
    trans_nir_fraction = dble_one - (lai*max_lai_nir_transmitted) &
                       / (lai+lai_half_nir_transmitted)
    ! Canopy reflected of near infrared and photosynthetically active radiation
    reflected_nir_fraction = (lai*max_lai_nir_reflection) &
                          / (lai+lai_half_nir_reflection)
    reflected_par_fraction = (lai*max_lai_par_reflection) &
                          / (lai+lai_half_par_reflection)
    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = dble_one - reflected_nir_fraction - trans_nir_fraction
    absorbed_par_fraction = dble_one - reflected_par_fraction - trans_par_fraction

    !!!!!!!!!!
    ! Estimate canopy absorption of incoming shortwave radiation
    !!!!!!!!!!

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = (sw_par_fraction * swrad * absorbed_par_fraction)
    canopy_nir_MJday = ((dble_one - sw_par_fraction) * swrad * absorbed_nir_fraction)
    trans_par_MJday = (sw_par_fraction * swrad * trans_par_fraction)
    trans_nir_MJday = ((dble_one - sw_par_fraction) * swrad * trans_nir_fraction)
    refl_par_MJday = (sw_par_fraction * swrad * reflected_par_fraction)
    refl_nir_MJday = ((dble_one - sw_par_fraction) * swrad * reflected_nir_fraction)

    !!!!!!!!!
    ! Estimate soil absorption of shortwave passing through the canopy
    !!!!!!!!!

    ! Update soil reflectance based on snow cover
    absorbed_par_fraction_soil = soil_swrad_absorption
    absorbed_nir_fraction_soil = soil_swrad_absorption

    ! Then the radiation incident and ultimately absorbed by the soil surface itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday * absorbed_par_fraction_soil
    soil_nir_MJday = trans_nir_MJday * absorbed_nir_fraction_soil
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = soil_nir_MJday + soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection / absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! Update the canopy radiation absorption based on the reflected radiation (MJ.m-2.day-1)
    canopy_par_MJday = canopy_par_MJday + ((trans_par_MJday-soil_par_MJday) * absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + ((trans_nir_MJday-soil_nir_MJday) * absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + ((trans_par_MJday-soil_par_MJday) * trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + ((trans_nir_MJday-soil_nir_MJday) * trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = canopy_par_MJday + canopy_nir_MJday

  end subroutine calculate_shortwave_balance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! log law decay
    canopy_wind = (ustar / vonkarman) * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot(Rtot)

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! declare inputs
    double precision,intent(inout) :: Rtot ! MPa.s-1.m-2.mmol-1

    ! local variables
    integer :: i
    double precision :: bonus, sum_water_flux, &
                        transpiration_resistance,root_reach_local, &
                        root_depth_50
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,root_length  &
                                                   ,ratio
    double precision, parameter :: root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located
    ! reset water flux
    water_flux = dble_zero
    ratio = dble_zero ; ratio(1) = dble_one ; root_mass = dble_zero
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    transpiration_resistance = canopy_height / (gplant * max(min_lai,lai))

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! initialise root reach based on initial conditions
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! Determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
    layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
    layer_thickness(4) = max_depth - sum(layer_thickness(1:3))
    layer_thickness(5) = top_soil_depth

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

    ! top 25 % of root profile
    root_depth_50 = root_reach * root_depth_frac_50
    if (root_depth_50 <= layer_thickness(1)) then
        ! Greater than 50 % of the fine root biomass can be found in the top
        ! soil layer

        ! Start by assigning all 50 % of root biomass to the top soil layer
        root_mass(1) = root_biomass * 0.5d0
        ! Then quantify how much additional root is found in the top soil layer
        ! assuming that the top 25 % depth is found somewhere within the top
        ! layer
        bonus = (root_biomass-root_mass(1)) &
              * (layer_thickness(1)-root_depth_50) / (root_reach - root_depth_50)
        root_mass(1) = root_mass(1) + bonus
        ! partition the remaining root biomass between the seconds and third
        ! soil layers
        if (root_reach > sum(layer_thickness(1:2))) then
            root_mass(2) = (root_biomass - root_mass(1)) &
                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
            root_mass(3) = root_biomass - sum(root_mass(1:2))
        else
            root_mass(2) = root_biomass - root_mass(1)
        endif
    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then
        ! Greater than 50 % of fine root biomass found in the top two soil
        ! layers. We will divide the root biomass uniformly based on volume,
        ! plus bonus for the second layer (as done above)
        root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/root_depth_50)
        root_mass(2) = root_biomass * 0.5d0 * ((root_depth_50-layer_thickness(1))/root_depth_50)
        ! determine bonus for the seconds layer
        bonus = (root_biomass-sum(root_mass(1:2))) &
              * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach-root_depth_50))
        root_mass(2) = root_mass(2) + bonus
        root_mass(3) = root_biomass - sum(root_mass(1:2))
    else
        ! Greater than 50 % of fine root biomass stock spans across all three
        ! layers
        root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/root_depth_50)
        root_mass(2) = root_biomass * 0.5d0 * (layer_thickness(2)/root_depth_50)
        root_mass(3) = root_biomass - sum(root_mass(1:2))
    endif
    ! now convert root mass into lengths
    root_length = root_mass * root_mass_length_coef_1

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    demand = abs(minlwp)+head*canopy_height
    ! now loop through soil layers, where root is present
    do i = 1, nos_root_layers
      if (root_mass(i) > dble_zero) then
          ! if there is root then there is a water flux potential...
          root_reach_local = min(root_reach,layer_thickness(i))
          ! calculate and accumulate steady state water flux in mmol.m-2.s-1
          water_flux(i) = plant_soil_flow(i,root_length(i),root_mass(i) &
                         ,demand(i),root_reach_local,transpiration_resistance)
      else
          ! ...if there is not then we wont have any below...
          exit
      end if ! root present in current layer?
    end do ! nos_root_layers

    ! if freezing then assume soil surface is frozen
    if (meant < dble_one) then
        water_flux(1) = dble_zero
        ratio(1) = dble_zero
        ratio(2:nos_root_layers) = layer_thickness(2:nos_root_layers) / sum(layer_thickness(2:nos_root_layers))
    else
        ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif

    ! calculate sum value
    sum_water_flux = sum(water_flux)
    ! calculate uptake fraction
    uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux

    ! sanity check in case of zero flux
    if (sum_water_flux <= vsmall) then
        uptake_fraction = dble_zero ; uptake_fraction(1) = dble_one
    endif

    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum_water_flux

    ! finally convert transpiration flux (mmol.m-2.s-1)
    ! into kg.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength(doy,lat)

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate
    ! 1) day length in hours and seconds
    ! 2) hour of sunrise and sunset
    ! 3) cosine of solar zenith angle to allow scaling of evaporation over course of day

    implicit none

    ! arguments
    double precision, intent(in) :: doy, lat

    ! local variables
    double precision :: dec, mult, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! Declination
    ! NOTE: 0.002739726d0 = 1/365
!    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
!    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) * 0.002739726d0 ) )

    ! latitude in radians
    mult = lat * deg_to_rad
    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-dble_one,min(dble_one,sinld / cosld))

    ! estimate day length in hours and seconds and upload to module variables
    dayl_hours = 12d0 * ( dble_one + 2d0 * asin( aob ) * pi_1 )
    dayl_seconds = dayl_hours * seconds_per_hour

    ! estimate sun rise and run set hours
!    sunrise = 12 - nint(dayl_hours*0.5d0) ; sunset = sunrise + nint(dayl_hours)

    ! estimate the solar cosine zenith angle for 12 noon
!    cos_solar_zenith_angle = sinld + cosld

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    ! local variables
    double precision  sqrt_cd1_lai &
                     ,local_lai &
                     ,phi_h       ! roughness sublayer influence function
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
!                          ustar_Uh_max = 0.3,   & ! Maximum observed ratio of
                                                   ! (friction velocity / canopy top wind speed) (m.s-1)
                          ustar_Uh_max = 1d0, ustar_Uh_min = 0.2d0, &
                                    Cw = 2d0      ! Characterises roughness sublayer depth (m)

    ! assign new value to min_lai to avoid max min calls
    local_lai = max(min_lai,lai)
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (dble_one-((dble_one-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    phi_h = 0.19314718056d0
!    phi_h = log(Cw)-dble_one+Cw**(-dble_one) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl = ((dble_one-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

    ! sanity check
!    if (roughl /= roughl) then
!        write(*,*)"TLS:  ERROR roughness length calculations"
!        write(*,*)"Roughness lenght", roughl, "Displacement", displacement
!        write(*,*)"canopy height", canopy_height, "lai", lai
!    endif

  end subroutine z0_displacement
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split & d0 has been used. !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    ! local variables..
    double precision :: answer, denominator, numerator

    numerator   = t - 25d0
    denominator = t + freeze
    answer      = a * exp( b * dble_one * numerator / denominator )
    arrhenious  = answer

  end function arrhenious
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
    double precision :: sum_x, sum_y, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    sum_product_xy = sum(x*y)
    ! calculate the gradient
    linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val , optimum , kurtosis , current )

    ! estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis

    implicit none

    ! arguments..
    double precision,intent(in) :: max_val, optimum, kurtosis, current

    ! local variables..
    double precision :: dummy

    if ( current >= max_val ) then
       opt_max_scaling = dble_zero
    else
       dummy     = ( max_val - current ) / ( max_val - optimum )
       dummy     = exp( log( dummy ) * kurtosis * ( max_val - optimum ) )
       opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
    end if

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function root_resistance (root_mass,thickness)

   !
   ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
   !

   implicit none

   ! arguments
   double precision :: root_mass, & ! root biomass in layer (gbiomass)
                       thickness    ! thickness of soil zone roots are in

   ! calculate root hydraulic resistance
   root_resistance = root_resist / (root_mass*thickness)
   ! return
   return

  end function root_resistance
  !
  !------------------------------------------------------------------
  !
  double precision function plant_soil_flow(root_layer,root_length,root_mass &
                                           ,demand,root_reach_in,transpiration_resistance)

   !
   ! Calculate soil layer specific water flow form the soil to canopy (mmolH2O.m-2.s-1)
   ! Accounting for soil, root and plant resistance, and canopy demand
   !

   ! calculate and accumulate steady state water flux in mmol.m-2.s-1
   ! From the current soil layer given an amount of root within the soil layer.
   ! NOTE: that soil moisture impacts are not included in this case

   implicit none

   ! arguments
   integer, intent(in) :: root_layer
   double precision, intent(in) :: root_length, &
                                     root_mass, &
                                        demand, &
                                 root_reach_in, &
                      transpiration_resistance

   ! local arguments
   double precision :: soilR2

   ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
   soilR2 = root_resistance(root_mass,root_reach_in)
   plant_soil_flow = demand/(transpiration_resistance + soilR2)

   ! return
   return

  end function plant_soil_flow
  !
  !------------------------------------------------------------------
  !
  double precision function zbrent( called_from , func , x1 , x2 , tol )

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    !  reference to a particular function and also two values which bound !
    !  the arguments for the function of interest. ZBRENT finds a root of !
    !  the function (i.e. the point where the function equals zero), that !
    !  lies between the two bounds.                                       !
    ! For a full description see Press et al. (1986).                     !

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    double precision,intent(in) :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision ,intent(in) :: xval
       end function func
    end interface

    ! local variables..
    integer            :: iter
    integer, parameter :: ITMAX = 30
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision, parameter :: EPS = 3d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0d0 ) then
      zbrent = a
      return
    elseif ( fb .eq. 0d0 ) then
      zbrent = b
      return
    end if
    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
    if ( sign(1d0,fa) .eq. sign(1d0,fb) ) then
       fa = func( a )
       fb = func( b )
       ! tell me otherwise what is going on
!       print*,"Supplied values must bracket the root of the function.",new_line('x'),  &
!         "     ","You supplied x1:",x1,new_line('x'),                     &
!         "     "," and x2:",x2,new_line('x'),                             &
!         "     "," which give function values of fa :",fa,new_line('x'),  &
!         "     "," and fb:",fb," .",new_line('x'),                        &
!         " zbrent was called by: ",trim(called_from)
!       fa = func( a )
!       fb = func( b )
    end if
    c = b
    fc = fb

    do iter = 1 , ITMAX

       ! If the new value (f(c)) doesn't bracket
       ! the root with f(b) then adjust it..
       if ( sign(1d0,fb) .eq. sign(1d0,fc) ) then
          c  = a
          fc = fa
          d  = b - a
          e  = d
       end if
       if ( abs(fc) .lt. abs(fb) ) then
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2d0 * EPS * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
          zbrent = b
          return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2d0 * xm * s
             q = 1d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2d0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1d0 ) )
             q = ( q - 1d0 ) * ( r - 1d0 ) * ( s - 1d0 )
          end if
          if ( p .gt. 0d0 ) q = -q
          p = abs( p )
          if ( (2d0*p) .lt. min( 3d0*xm*q-abs(tol1*q) , abs(e*q) ) ) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a  = b
       fa = fb
       if ( abs(d) .gt. tol1 ) then
          b = b + d
       else
          b = b + sign( tol1 , xm )
       end if
       fb = func(b)
    enddo

!    print*,"zbrent has exceeded maximum iterations",new_line('x'),&
!           "zbrent was called by: ",trim(called_from)

    zbrent = b

  end function zbrent
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
