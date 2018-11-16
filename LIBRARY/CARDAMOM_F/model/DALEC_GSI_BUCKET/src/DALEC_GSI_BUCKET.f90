
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL                  &
         ,dble_one,dble_zero            &
         ,vsmall                        &
         ,arrhenious                    &
         ,acm_gpp                       &
         ,calculate_transpiration       &
         ,calculate_wetcanopy_evaporation &
         ,calculate_soil_evaporation    &
         ,calculate_shortwave_balance   &
         ,calculate_longwave_isothermal &
         ,acm_albedo_gc                 &
         ,acm_meteorological_constants  &
         ,update_soil_initial_conditions&
         ,calculate_daylength           &
         ,daylength_hours               &
         ,daylength_seconds             &
         ,opt_max_scaling               &
         ,freeze                        &
         ,co2comp_saturation            &
         ,co2comp_half_sat_conc         &
         ,co2_half_saturation           &
         ,co2_compensation_point        &
         ,kc_saturation                 &
         ,kc_half_sat_conc              &
         ,calculate_update_soil_water   &
         ,calculate_Rtot                &
         ,calculate_aerodynamic_conductance &
         ,saxton_parameters             &
         ,initialise_soils              &
         ,linear_model_gradient         &
         ,seconds_per_day               &
         ,seconds_per_step              &
         ,root_biomass                  &
         ,root_reach                    &
         ,min_root                      &
         ,max_depth                     &
         ,root_k                        &
         ,top_soil_depth                &
         ,mid_soil_depth                &
         ,previous_depth                &
         ,nos_root_layers               &
         ,wSWP                          &
         ,SWP                           &
         ,SWP_initial                   &
         ,deltat_1                      &
         ,water_flux                    &
         ,layer_thickness               &
         ,waterloss,watergain           &
         ,potA,potB                     &
         ,cond1,cond2,cond3             &
         ,soil_conductivity             &
         ,soil_waterfrac                &
         ,soil_waterfrac_initial        &
         ,porosity                      &
         ,porosity_initial              &
         ,field_capacity                &
         ,field_capacity_initial        &
         ,drythick                      &
         ,min_drythick                  &
         ,min_layer                     &
         ,soilwatermm                   &
         ,wSWP_time                     &
         ,soil_frac_clay                &
         ,soil_frac_sand                &
         ,nos_soil_layers               &
         ,meant                         &
         ,meant_K                       &
         ,stomatal_conductance          &
         ,iWUE                          &
         ,avN                           &
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
         ,mint                          &
         ,maxt                          &
         ,leafT                         &
         ,swrad                         &
         ,co2                           &
         ,doy                           &
         ,rainfall                      &
         ,airt_zero_fraction            &
         ,snowfall                      &
         ,snow_melt                     &
         ,wind_spd                      &
         ,vpd_pa                        &
         ,lai                           &
         ,days_per_step                 &
         ,days_per_step_1               &
         ,dayl_seconds                  &
         ,dayl_seconds_1                &
         ,dayl_hours                    &
         ,snow_storage                  &
         ,canopy_storage                &
         ,intercepted_rainfall          &
         ,disturbance_residue_to_litter &
         ,disturbance_residue_to_cwd    &
         ,disturbance_residue_to_som    &
         ,disturbance_loss_from_litter  &
         ,disturbance_loss_from_cwd     &
         ,disturbance_loss_from_som     &
         ,itemp,ivpd,iphoto             &
         ,dim_1,dim_2                   &
         ,nos_trees                     &
         ,nos_inputs                    &
         ,leftDaughter                  &
         ,rightDaughter                 &
         ,nodestatus                    &
         ,xbestsplit                    &
         ,nodepred                      &
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
logical :: do_iWUE = .true.  ! Use iWUE or WUE for stomatal optimisation
double precision, parameter :: dble_zero = 0d0    &
                              ,dble_one = 1d0     &
                              ,vsmall = tiny(0d0)*1d3 ! *1d3 to add a little breathing room

integer, parameter :: nos_root_layers = 3, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: pi = 3.1415927d0,  &
                             pi_1 = pi**(-dble_one),   &
                              pi2 = pi**2,        &
                           two_pi = pi*2d0,       &
                       deg_to_rad = pi/180d0,     &
              sin_dayl_deg_to_rad = sin( 23.45d0 * deg_to_rad ), & ! repeated function in acm
                          gravity = 9.8067d0,     & ! acceleration due to gravity, ms-1
                            boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96d0,       &
                      emiss_boltz = emissivity * boltz, &
                  sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                           freeze = 273.15d0,     &
                       gs_H2O_CO2 = 1.646259d0,   & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                     gs_H2O_CO2_1 = gs_H2O_CO2 ** (-dble_one), &
                       gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
          partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                       umol_to_gC = 1d-6*12d0,    & ! conversion of umolC -> gC
                 mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                   mol_to_g_water = 18d0,         & ! molecular mass of water (g)
                     mol_to_g_co2 = 12d0,         & ! molecular mass of CO2 (g)
                     g_to_mol_co2 = 1d0/12d0,     &
!snowscheme       density_of_water = 998.9d0,         & ! density of !water kg.m-3
                   gas_constant_d = 287.04d0,     & ! gas constant for dry air (J.K-1.mol-1)
                             Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                        vonkarman = 0.41d0,       & ! von Karman's constant
                      vonkarman_2 = vonkarman**2, & ! von Karman's constant^2
                            cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

! photosynthesis / respiration parameters
double precision, parameter :: &
                    kc_saturation = 310d0,        & ! CO2 half saturation, saturation value
                 kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, half sat
               co2comp_saturation = 36.5d0,       & ! CO2 compensation point, saturation
            co2comp_half_sat_conc = 9.46d0,       & ! CO2 comp point, half sat
                                                    ! Each of these are temperature
                                                    ! sensitivty
                      Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                    ! lost as growth respiration
                                                    ! (i.e. 0.28 .eq. xNPP)
                  one_Rg_fraction = dble_one - Rg_fraction

! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5d0,        & ! tortuosity
                           gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                      root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                    ! Williams et al 1996 = 0.0001
                    root_radius_1 = root_radius**(-dble_one), &
              root_cross_sec_area = pi * root_radius**2, & ! root cross sectional area (m2)
                                                           ! = pi * root_radius * root_radius
                     root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                    ! 0.5e6 Williams et al 1996
                                                    ! 0.31e6 Bonan et al 2014
          root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-dble_one), &
               const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                             head = 0.009807d0,   & ! head of pressure (MPa/m)
                           head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

! structural parameters
double precision, parameter :: &
                    canopy_height = 9d0,          & ! canopy height assumed to be 9 m
                     tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                         min_wind = 0.1d0,        & ! minimum wind speed at canopy top
                     min_drythick = 0.01d0,       & ! minimum dry thickness depth (m)
                        min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                      soil_roughl = 0.05d0,       & ! soil roughness length (m)
                   top_soil_depth = 0.1d0,        & ! thickness of the top soil layer (m)
                   mid_soil_depth = 0.2d0,        & ! thickness of the second soil layer (m)
                         min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                          min_lai = 1.5d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                  min_throughfall = 0.2d0,        & ! minimum fraction of precipitation which
                                                    ! is through fall
                      min_storage = 0.2d0           ! minimum canopy water (surface) storage (mm)

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
                   ,SLA & ! Specific leaf area
                   ,avail_labile,Rg_from_labile    &
                   ,Cwood_labile_release_gradient  &
                   ,Cwood_labile_half_saturation   &
                   ,Croot_labile_release_gradient  &
                   ,Croot_labile_half_saturation   &
                   ,Cwood_hydraulic_gradient       &
                   ,Cwood_hydraulic_half_saturation&
                   ,Cwood_hydraulic_limit          &
                   ,delta_gsi,tmp,gradient         &
                   ,fol_turn_crit,lab_turn_crit    &
                   ,gsi_history(22),just_grown

double precision, allocatable, dimension(:) :: itemp,ivpd,iphoto, &
                                               disturbance_residue_to_litter, &
                                               disturbance_residue_to_som,    &
                                               disturbance_residue_to_cwd,    &
                                               disturbance_loss_from_litter,  &
                                               disturbance_loss_from_cwd,     &
                                               disturbance_loss_from_som,     &
                                               tmp_x, tmp_m

! hydraulic model variables
integer :: water_retention_pass, soil_layer, sunrise, sunset
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil
double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                         demand, & ! maximum potential canopy hydraulic demand
                                                     water_flux    ! potential transpiration flux (mmol.m-2.s-1)
double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa)
                                          SWP_initial, &
                                    soil_conductivity, & ! soil conductivity
                                            waterloss, & ! water loss from specific soil layers (m)
                                            watergain, & ! water gained by specfic soil layers (m)
                                       field_capacity, & ! soil field capacity (m3.m-3)
                               field_capacity_initial, &
                                       soil_waterfrac, & ! soil water content (m3.m-3)
                               soil_waterfrac_initial, &
                                             porosity, & ! soil layer porosity, (fraction)
                                     porosity_initial, &
                                      layer_thickness, & ! thickness of soil layers (m)
                      cond1, cond2, cond3, potA, potB    ! Saxton equation values

double precision :: root_reach, root_biomass, &
                  drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                      wSWP, & ! weighted soil water potential (MPa) used in GSI calculate.
                              ! Removes / limits the fact that very low root density in young plants
                              ! give values too large for GSI to handle.
                 max_depth, & ! maximum possible root depth (m)
                    root_k, & ! biomass to reach half max_depth
                    runoff, & ! runoff (kgH2O.m-2.day-1)
                 underflow, & ! drainage from the bottom of soil column (kgH2O.m-2.day-1)
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
                 maxt_lag1, &
                     leafT, & ! canopy temperature (oC)
          mean_annual_temp, &
        canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
          soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
            soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
             sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                 ci_global, & ! internal CO2 concentration (ppm or umol/mol)
      stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.s-1)
   aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
          soil_conductance, & ! soil surface conductance (m.s-1)
         convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
       air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                    lambda, & ! latent heat of vapourisation (J.kg-1)
                     psych, & ! psychrometric constant (kPa K-1)
                     slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
              snow_storage, & ! snow storage (kgH2O/m2)
            canopy_storage, & ! water storage on canopy (kgH2O.m-2)
      intercepted_rainfall    ! intercepted rainfall rate equivalent (kgH2O.m-2.s-1)

! Module level variables for ACM_GPP_ET parameters
double precision :: delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                       avN, & ! average foliar N (gN/m2)
                      iWUE, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
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
      airt_zero_fraction, & ! fraction of air temperature above freezing
                   swrad, & ! incoming short wave radiation (MJ/m2/day)
                     co2, & ! CO2 (ppm)
                     doy, & ! Day of year
                rainfall, & ! rainfall (kgH2O/m2/s)
                snowfall, &
               snow_melt, & ! snow melt (kgH2O/m2/s)
                wind_spd, & ! wind speed (m/s)
                  vpd_pa, & ! Vapour pressure deficit (Pa)
                     lai    ! leaf area index (m2/m2)

! Module level varoables for step specific timing information
double precision :: cos_solar_zenith_angle, &
                  mean_days_per_step, & !
                    seconds_per_step, & !
                       days_per_step, & !
                     days_per_step_1, & !
                        dayl_seconds, & ! day length in seconds
                      dayl_seconds_1, &
                          dayl_hours    ! day length in hours

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                           daylength_hours, &
                                         daylength_seconds, &
                                       co2_half_saturation, & ! (ppm)
                                    co2_compensation_point, & ! (ppm)
                                 Cwood_labile_release_coef, & ! time series of labile release to wood
                                 Croot_labile_release_coef, & ! time series of labile release to root
                                               soilwatermm, &
                                                 wSWP_time

save

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out,NEE_out,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP_out)

    !
    ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
    ! Index - BUCKET (DALEC_GSI_BUCKET) model.
    !
    ! The Aggregated Canopy Model for Gross Primary Productivity and Evapotranspiration (ACM-GPP-ET)
    ! simulates coupled photosynthesis-transpiration (via stomata), soil and intercepted canopy evaporation and
    ! soil water balance (4 layers).
    !
    ! Carbon allocation based on fixed fraction and turnover follows first order kinetics with the following exceptions.
    ! 1) Foliar allocation and turnover is based on the Growing Season Index framework.
    ! 2) Turnover of litter and soil includes an exponential temperature dependency.
    !
    ! This version was coded by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! Version 1: 15/11/18

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   &
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

    ! declare general local variables
    double precision ::  infi &
                ,Tfac_range_1 &
            ,Photofac_range_1 &
              ,VPDfac_range_1 &
                     ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
                        ,Rtot & ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)
               ,act_pot_ratio &
               ,transpiration &
             ,soilevaporation &
              ,wetcanopy_evap &
            ,snow_sublimation

    integer :: f,nxp,n,test,m

    ! local fire related variables
    double precision :: burnt_area &
                       ,CFF(7) = dble_zero, CFF_res(4) = dble_zero    & ! combusted and non-combustion fluxes
                       ,NCFF(7) = dble_zero, NCFF_res(4) = dble_zero  & ! with residue and non-residue seperates
                       ,combust_eff(5)                                & ! combustion efficiency
                       ,rfac                                            ! resilience factor

    ! local deforestation related variables
    double precision, dimension(5) :: post_harvest_burn   & ! how much burning to occur after
                                     ,foliage_frac_res    &
                                     ,roots_frac_res      &
                                     ,rootcr_frac_res     &
                                     ,stem_frac_res       &
                                     ,branch_frac_res     &
                                     ,Cbranch_part        &
                                     ,Crootcr_part        &
                                     ,soil_loss_frac

    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   &
                       ,C_total,labile_frac_res      &
                       ,Cstem,Cbranch,Crootcr        &
                       ,stem_residue,branch_residue  &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    integer :: reforest_day, harvest_management,restocking_lag, gsi_lag

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
    ! 12th 21 day average VPD (Pa)
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
    ! 19 = Evapotranspiration (kgH2O.m-2.day-1)

    ! PARAMETERS
    ! 23 process parameters; 7 C pool initial conditions; 1 soil water initial condition

    ! p(1) Litter to SOM conversion rate  - m_r
    ! p(2) Fraction of GPP respired - f_a
    ! p(3) Fraction of NPP allocated to foliage - f_f
    ! p(4) Fraction of NPP allocated to roots - f_r
    ! p(5) max leaf turnover (GSI)
    ! p(6) Turnover rate of wood - t_w
    ! p(7) Turnover rate of roots - t_r
    ! p(8) Litter turnover rate - t_l
    ! p(9) SOM turnover rate  - t_S
    ! p(10) Parameter in exponential term of temperature - \theta
    ! p(11) mean foliar nitrogen content (gN/m2)
    ! p(12) = max labile turnover(GSI)
    ! p(13) = Fraction allocated to Clab - f_l
    ! p(14) = min temp threshold (GSI)
    ! p(15) = max temp threshold (GSI)
    ! p(16) = min photoperiod threshold (GIS)
    ! p(17) = LMA
    ! p(24) = max photoperiod threshold (GSI)
    ! p(25) = min VPD threshold (GSI)
    ! p(26) = max VPD threshold (GSI)
    ! p(27) = minimum GPP benefit of increased LAI for labile allocation to be allowed
    ! p(28) = fraction of Cwood which is Cbranch
    ! p(29) = fraction of Cwood which is Ccoarseroot
    ! p(37) = Initial CWD pool
    ! p(38) = CWD turnover fraction
    ! p(39) = Fine root (gbiomass.m-2) needed to reach 50% of max depth
    ! p(40) = Maximum rooting depth (m)
    ! p(41) = Initial root profile water content (m3/m3)

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
!call cpu_time(start)
!call cpu_time(finish)

    ! load ACM-GPP-ET parameters
    ! load ACM-GPP-ET parameters
    NUE                        = 1.182549d+01  ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
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

    ! Other initial values for ACM_GPP_ET
    avN = 10d0**pars(11)  ! foliar N gN/m2
    deltaWP = minlwp ! leafWP-soilWP (i.e. -2-0)
    Rtot = dble_one

    ! plus ones being calibrated
    root_k = pars(39) ; max_depth = pars(40)

    ! reset values
    intercepted_rainfall = dble_zero ; canopy_storage = dble_zero ; snow_storage = dble_zero

    ! initial values for deforestation variables
    labile_loss = 0d0    ; foliar_loss = 0d0
    roots_loss = 0d0     ; wood_loss = 0d0
    labile_residue = 0d0 ; foliar_residue = 0d0
    roots_residue = 0d0  ; wood_residue = 0d0
    stem_residue = 0d0   ; branch_residue = 0d0
    reforest_day = 0
    soil_loss_with_roots = 0d0
    coarse_root_residue = 0d0
    post_harvest_burn = 0d0

    ! now load the hardcoded forest management parameters into their locations

    ! Parameter values for deforestation variables
    ! scenario 1
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = 1.0d0
    roots_frac_res(1)   = 1.0d0
    rootcr_frac_res(1) = 1.0d0
    branch_frac_res(1) = 1.0d0
    stem_frac_res(1)   = 0d0 !
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
    post_harvest_burn(1) = 1d0

    !## scen 2
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = 1.0d0
    roots_frac_res(2)   = 1.0d0
    rootcr_frac_res(2) = 1.0d0
    branch_frac_res(2) = 1.0d0
    stem_frac_res(2)   = 0d0 !
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
    post_harvest_burn(2) = 0d0

    !## scen 3
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5d0
    roots_frac_res(3)   = 1.0d0
    rootcr_frac_res(3) = 1.0d0
    branch_frac_res(3) = 0d0
    stem_frac_res(3)   = 0d0 !
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
    post_harvest_burn(3) = 0d0

    !## scen 4
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5d0
    roots_frac_res(4)   = 1.0d0
    rootcr_frac_res(4) = 0d0
    branch_frac_res(4) = 0d0
    stem_frac_res(4)   = 0d0 !
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(4) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(4) = 0d0

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
    ! POOL(1,8) assigned later

    if (.not.allocated(deltat_1)) then
        ! allocate variables dimension which are fixed per site only the once
        allocate(disturbance_residue_to_litter(nodays), &
                 disturbance_residue_to_cwd(nodays),    &
                 disturbance_residue_to_som(nodays),    &
                 disturbance_loss_from_litter(nodays),  &
                 disturbance_loss_from_cwd(nodays),     &
                 disturbance_loss_from_som(nodays),     &
                 Cwood_labile_release_coef(nodays),     &
                 Croot_labile_release_coef(nodays),     &
                 co2_compensation_point(nodays),        &
                 co2_half_saturation(nodays),           &
                 deltat_1(nodays),soilwatermm(nodays),  &
                 wSWP_time(nodays))
        ! calculate inverse of number of days per time step, as this will be useful to avoid divisions later
        deltat_1 = deltat**(-1d0)

        ! first those linked to the time period of the analysis
        do n = 1, nodays
          ! Temperature adjustments for Michaelis-Menten coefficients
          ! for CO2 (kc) and O2 (ko) and CO2 compensation point.
          co2_compensation_point(n) = arrhenious(co2comp_saturation,co2comp_half_sat_conc,met(3,n))
          co2_half_saturation(n) = arrhenious(kc_saturation,kc_half_sat_conc,met(3,n))
        end do

        ! Temperature limitiation parameters on wood and fine root growth.
        ! Parmeters generated on the assumption of 5 % / 95 % activation at key
        ! temperature values. Roots 1oC/30oC, wood 5oC/30oC.
        ! NOTE: Foliage and root potential turnovers use the same temperature curve
        Croot_labile_release_gradient = 0.1962d0 ; Croot_labile_half_saturation = 15.0d0
        Cwood_labile_release_gradient = 0.2355d0 ; Cwood_labile_half_saturation  = 17.5d0
        ! calculate temperature limitation on potential wood/root growth
        Cwood_labile_release_coef = (dble_one+exp(-Cwood_labile_release_gradient* &
                                  (((met(3,:)+met(2,:))*0.5d0)-Cwood_labile_half_saturation)))**(-dble_one)
        Croot_labile_release_coef = (dble_one+exp(-Croot_labile_release_gradient* &
                                  (((met(3,:)+met(2,:))*0.5d0)-Croot_labile_half_saturation)))**(-dble_one)
        ! hydraulic limitation parameters for wood cell expansion, i.e. growth
        Cwood_hydraulic_gradient = 5d0 ; Cwood_hydraulic_half_saturation = -1.5d0

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
              if (sum(deltat(1:n)) < 21) then
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
                 do while (test < 21)
                    m = m+1 ; test = nint(sum(deltat((n-m):n)))
                    if (m > (n-1)) test = 21
                 end do
                 tmp_m(n) = m
              endif ! for calculating gradient
            end do ! calc daily values once
            ! allocate GSI history dimension
            gsi_lag_remembered = max(2,maxval(nint(tmp_m)))
        end if ! .not.allocated(tmp_x)
        ! assign our starting value
        gsi_history = pars(36)-1d0
         just_grown = pars(35)

        ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS

       !
       ! Initialise the water model
       !

       ! zero variables not done elsewhere
       water_flux = dble_zero
       ! initialise some time invarient parameters
       call saxton_parameters(soil_frac_clay,soil_frac_sand)
       call initialise_soils(soil_frac_clay,soil_frac_sand)
       call update_soil_initial_conditions(pars(41))
       ! save the initial conditions for later
       soil_waterfrac_initial = soil_waterfrac
       SWP_initial = SWP
       field_capacity_initial = field_capacity
       porosity_initial = porosity

    else
       !
       ! Load initial soil water conditions from memory
       !

       water_flux = dble_zero
       field_capacity = field_capacity_initial
       porosity = porosity_initial

       ! input initial soil water fraction then
       ! update SWP and soil conductivity accordingly
       call update_soil_initial_conditions(pars(41))

    endif

    ! Initialise root reach based on initial conditions
    root_biomass = max(min_root,POOLS(1,3)*2d0)
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! Determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
    layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
    layer_thickness(4) = max_depth - sum(layer_thickness(1:3))
    layer_thickness(5) = top_soil_depth
    previous_depth = max(top_soil_depth,root_reach)
    ! Needed to initialise soils
    call calculate_Rtot(Rtot)
    ! Used to initialise soils
    call calculate_update_soil_water(dble_zero,dble_zero,dble_zero,FLUXES(1,19)) ! assume no evap or rainfall
    ! Store soil water content of the surface zone (mm)
    POOLS(1,8) = 1d3*soil_waterfrac(1)*layer_thickness(1)

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit=pars(34)-dble_one
    lab_turn_crit=pars(3)-dble_one
    Tfac_range_1 = (pars(15)-pars(14))**(-1d0)
    Photofac_range_1 = (pars(24)-pars(16))**(-1d0)
    VPDfac_range_1 = abs(pars(26)-pars(25))**(-1d0)
    SLA = pars(17)**(-1d0)

    infi = 0d0

    ! clear disturbance arrays to avoid throw overs
    disturbance_residue_to_litter = 0d0 ; disturbance_loss_from_litter = 0d0
    disturbance_residue_to_som = 0d0 ; disturbance_loss_from_som = 0d0
    disturbance_residue_to_cwd = 0d0 ; disturbance_loss_from_cwd = 0d0


    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! set lag information using previous time step value for temperature
      maxt_lag1 = maxt
      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      leafT = maxt     ! initial canopy temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = met(6,n)   ! Day of year
      rainfall = met(7,n) ! rainfall (kgH2O/m2/s)
      meant = (maxt + mint) * 0.5d0  ! mean air temperature (oC)
      meant_K = meant + freeze
      airt_zero_fraction = (maxt-dble_zero) / (maxt-mint) ! fraction of temperture period above freezing
      wind_spd = met(15,n) ! wind speed (m/s)
      vpd_pa = met(16,n)  ! Vapour pressure deficit (Pa)

      ! states needed for module variables
      lai_out(n) = POOLS(n,2)*SLA
      lai = lai_out(n) ! leaf area index (m2/m2)

      ! Temperature adjustments for Michaelis-Menten coefficients
      ! for CO2 (kc) and O2 (ko) and CO2 compensation point
      ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
      co2_half_sat   = co2_half_saturation(n)
      co2_comp_point = co2_compensation_point(n)

      ! calculate daylength in hours and seconds
      call calculate_daylength((doy-(deltat(n)*0.5d0)),lat)
      ! extract timing related values
      dayl_seconds_1 = dayl_seconds ** (-dble_one)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      ! snowing or not...?
      snow_melt = dble_zero ; snowfall = dble_zero
      if (mint < dble_zero .and. maxt > dble_zero) then
          ! if minimum temperature is below freezing point then we weight the
          ! rainfall into snow or rain based on proportion of temperature below
          ! freezing
          snowfall = dble_one - airt_zero_fraction
          snowfall = rainfall * snowfall ; rainfall = rainfall - snowfall
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)

          ! Also melt some of the snow
          snow_melt = airt_zero_fraction
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_melt * snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      elseif (maxt < dble_zero) then
          ! if whole day is below freezing then we should assume that all
          ! precipitation is snowfall
          snowfall = rainfall ; rainfall = dble_zero
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)
      else if (mint > dble_zero) then
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      end if

      !!!!!!!!!!
      ! calculate soil water potential and total hydraulic resistance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2d0)
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(dble_one,dble_one - (soil_waterfrac(1) / porosity(1))))
      call calculate_Rtot(Rtot)
      ! Pass wSWP to output variable and update deltaWP between minlwp and
      ! current weighted soil WP
      wSWP_time(n) = wSWP ; deltaWP = min(dble_zero,minlwp-wSWP)

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

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

      !!!!!!!!!!
      ! Estimate evaporative and photosynthetic fluxes
      !!!!!!!!!!

      ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
      call calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,canopy_storage,dble_zero)
      ! calculate radiation absorption and estimate stomatal conductance
      call acm_albedo_gc(abs(deltaWP),Rtot)

      ! if snow present assume that soilevaporation is sublimation of soil first
      snow_sublimation = dble_zero
      if (snow_storage > dble_zero) then
          snow_sublimation = soilevaporation
          if (snow_sublimation*deltat(n) > snow_storage) snow_sublimation = snow_storage * deltat_1(n)
          soilevaporation = soilevaporation - snow_sublimation
          snow_storage = snow_storage - (snow_sublimation * deltat(n))
      end if

      ! Note that soil mass balance will be calculated after phenology
      ! adjustments

      ! reset output variable
      ci_global = dble_zero
      if (lai > vsmall .and. stomatal_conductance > vsmall) then
          ! Gross primary productivity (gC/m2/day)
          FLUXES(n,1) = max(dble_zero,acm_gpp(stomatal_conductance))
          ! Canopy transpiration (kgH2O/m2/day)
          call calculate_transpiration(transpiration)
          ! restrict transpiration to positive only
          transpiration = max(dble_zero,transpiration)
      else
          ! assume zero fluxes
          FLUXES(n,1) = dble_zero
          transpiration = dble_zero
      endif

      ! Soil surface (kgH2O.m-2.day-1)
      call calculate_soil_evaporation(soilevaporation)

      !!!!!!!!!!
      ! GPP allocation
      !!!!!!!!!!

      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*meant)
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = 0d0 !(FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))

      !!!!!!!!!!
      ! calculate canopy phenology
      !!!!!!!!!!

      ! Determine leaf growth and turnover based on GSI model + some economics
      ! NOTE: that turnovers will be bypassed in favour of mortality turnover
      ! should available labile be exhausted
      call calculate_leaf_dynamics(n,deltat,nodays        &
                                  ,pars(14),pars(16),pars(25)       &
                                  ,Tfac_range_1,Photofac_range_1    &
                                  ,VPDfac_range_1,pars(5),pars(12)  &
                                  ,met(10,n),met(11,n),met(12,n),deltaWP,Rtot &
                                  ,FLUXES(n,1),POOLS(n,2),pars(27)  &
                                  ,FLUXES(:,18),FLUXES(n,9),FLUXES(n,16))

      ! Total labile release to foliage
      FLUXES(n,8) = avail_labile*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))*deltat_1(n)
      ! Retrict based on available labile stores
      FLUXES(n,8) = min(avail_labile*deltat_1(n),FLUXES(n,8))
      ! Update available labile supply for fine roots and wood
      avail_labile = avail_labile - (FLUXES(n,8)*deltat(n))

      ! calculate allocation of labile to roots and wood including, where appropriate, marginal return calculations
      call calculate_wood_root_growth(n,pars(4),pars(13),deltaWP,Rtot,FLUXES(n,1) &
                                     ,POOLS(n,3),POOLS(n,4),FLUXES(n,6),FLUXES(n,7))

      !
      ! litter creation with time dependancies
      !

      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**deltat(n))*deltat_1(n)
      ! total wood litter production
      FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-pars(6))**deltat(n))*deltat_1(n)
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*(1d0-(1d0-pars(7))**deltat(n))*deltat_1(n)

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))*deltat_1(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1d0-(1d0-FLUXES(n,2)*pars(9))**deltat(n))*deltat_1(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(1))**deltat(n))*deltat_1(n)
      ! CWD to litter
      FLUXES(n,20) = POOLS(n,7)*(1d0-(1d0-FLUXES(n,2)*pars(38))**deltat(n))*deltat_1(n)

      !!!!!!!!!!
      ! calculate growth respiration and adjust allocation to pools assuming
      ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
      !!!!!!!!!!

      ! foliage
      Rg_from_labile =                   FLUXES(n,8)*Rg_fraction  ; FLUXES(n,8) = FLUXES(n,8) * (one_Rg_fraction)
      ! roots
      Rg_from_labile = Rg_from_labile + (FLUXES(n,6)*Rg_fraction) ; FLUXES(n,6) = FLUXES(n,6) * (one_Rg_fraction)
      ! wood
      Rg_from_labile = Rg_from_labile + (FLUXES(n,7)*Rg_fraction) ; FLUXES(n,7) = FLUXES(n,7) * (one_Rg_fraction)
      ! now update the Ra flux with Rg
      FLUXES(n,3) = FLUXES(n,3) + Rg_from_labile

      ! calculate the NEE
      NEE_out(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP_out(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8)-FLUXES(n,6)-FLUXES(n,7)-Rg_from_labile)*deltat(n)
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
      ! Update soil water balance
      !!!!!!!!!!

      ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
      rainfall = rainfall + (snow_melt / seconds_per_step)
      ! do mass balance (i.e. is there enough water to support ET)
      call calculate_update_soil_water(transpiration,soilevaporation,((rainfall-intercepted_rainfall)*seconds_per_day) &
                                      ,FLUXES(n,19))
      ! now that soil mass balance has been updated we can add the wet canopy
      ! evaporation (kg.m-2.day-1)
      FLUXES(n,19) = FLUXES(n,19) + wetcanopy_evap
      ! store soil water content of the surface zone (mm)
      POOLS(n,8) = 1d3*soil_waterfrac(1)*layer_thickness(1)

      !
      ! deal first with deforestation
      !

      if (n == reforest_day) then
          POOLS(n+1,1) = pars(30)
          POOLS(n+1,2) = pars(31)
          POOLS(n+1,3) = pars(32)
          POOLS(n+1,4) = pars(33)
      endif

      ! reset values
      FLUXES(n,17) = dble_zero ; FLUXES(n,22:25) = dble_zero
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
              wood_loss   = POOLS(n+1,4)*met(8,n)
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
              labile_residue = POOLS(n+1,1)*met(8,n)*labile_frac_res
              foliar_residue = POOLS(n+1,2)*met(8,n)*foliage_frac_res(harvest_management)
              roots_residue  = POOLS(n+1,3)*met(8,n)*roots_frac_res(harvest_management)
              ! Explicit calculation of the residues from each fraction
              coarse_root_residue  = Crootcr*met(8,n)*rootcr_frac_res(harvest_management)
              branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
              stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
              ! Now finally calculate the final wood residue
              wood_residue = stem_residue + branch_residue + coarse_root_residue
              ! Mechanical loss of Csom due to coarse root extraction
              soil_loss_with_roots = Crootcr*met(8,n)*(dble_one-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

              ! Update living pools directly
              POOLS(n+1,1) = max(dble_zero,POOLS(n+1,1)-labile_loss)
              POOLS(n+1,2) = max(dble_zero,POOLS(n+1,2)-foliar_loss)
              POOLS(n+1,3) = max(dble_zero,POOLS(n+1,3)-roots_loss)
              POOLS(n+1,4) = max(dble_zero,POOLS(n+1,4)-wood_loss)

              ! Set burn related values
              FLUXES(n,17) = dble_zero
              CFF = dble_zero ; NCFF = dble_zero
              CFF_res = dble_zero ; NCFF_res = dble_zero

              ! Update all pools this time
              POOLS(n+1,1) = max(dble_zero, POOLS(n+1,1) - CFF(1) - NCFF(1) )
              POOLS(n+1,2) = max(dble_zero, POOLS(n+1,2) - CFF(2) - NCFF(2) )
              POOLS(n+1,3) = max(dble_zero, POOLS(n+1,3) - CFF(3) - NCFF(3) )
              POOLS(n+1,4) = max(dble_zero, POOLS(n+1,4) - CFF(4) - NCFF(4) )
              POOLS(n+1,5) = max(dble_zero, POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue) &
                                                  + (NCFF(1)+NCFF(2)+NCFF(3))-CFF(5)-NCFF(5) )
              POOLS(n+1,6) = max(dble_zero, POOLS(n+1,6) - soil_loss_with_roots + (NCFF(4)+NCFF(5)+NCFF(7)))
              POOLS(n+1,7) = max(dble_zero, POOLS(n+1,7) + wood_residue - CFF(7) - NCFF(7) )
              ! Some variable needed for the EDCs
              ! reallocation fluxes for the residues
              disturbance_residue_to_litter(n) = (labile_residue+foliar_residue+roots_residue) &
                                               + (NCFF(1)+NCFF(2)+NCFF(3))
              disturbance_loss_from_litter(n)  = CFF(5)+NCFF(5)
              disturbance_residue_to_cwd(n)    = wood_residue
              disturbance_loss_from_cwd(n)     = CFF(7) - NCFF(7)
              disturbance_residue_to_som(n)    = NCFF(4)+NCFF(5)+NCFF(7)
              disturbance_loss_from_som(n)     = soil_loss_with_roots
              ! Convert all to rates to be consistent with the FLUXES in EDCs
              disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) * deltat_1(n)
              disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) * deltat_1(n)
              disturbance_residue_to_cwd(n)    = disturbance_residue_to_cwd(n) * deltat_1(n)
              disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) * deltat_1(n)
              disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) * deltat_1(n)
              disturbance_loss_from_som(n)     = disturbance_loss_from_som(n) * deltat_1(n)
              ! This is intended for use with the R interface for subsequent post
              ! processing
              FLUXES(n,21) =  (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                               + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                               + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                               + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))
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
!              !/*update net exchangep*/
!              NEE(n)=NEE(n)+FLUXES(n,17)
              ! determine the as daily rate impact on live tissues for use in EDC and
              ! MTT calculations
              FLUXES(n,22) = FLUXES(n,22) + ((CFF(1) + NCFF(1)) * deltat_1(n))
              FLUXES(n,23) = FLUXES(n,23) + ((CFF(2) + NCFF(2)) * deltat_1(n))
              FLUXES(n,24) = FLUXES(n,24) + ((CFF(3) + NCFF(3)) * deltat_1(n))
              FLUXES(n,25) = FLUXES(n,25) + ((CFF(4) + NCFF(4)) * deltat_1(n))

              !// update pools
              !/*Adding all fire pool transfers here*/
              POOLS(n+1,1)=max(dble_zero,POOLS(n+1,1)-CFF(1)-NCFF(1))
              POOLS(n+1,2)=max(dble_zero,POOLS(n+1,2)-CFF(2)-NCFF(2))
              POOLS(n+1,3)=max(dble_zero,POOLS(n+1,3)-CFF(3)-NCFF(3))
              POOLS(n+1,4)=max(dble_zero,POOLS(n+1,4)-CFF(4)-NCFF(4))
              POOLS(n+1,5)=max(dble_zero,POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3))
              POOLS(n+1,6)=max(dble_zero,POOLS(n+1,6)+NCFF(4)+NCFF(5)+NCFF(7))
              POOLS(n+1,7)=max(dble_zero,POOLS(n+1,7)-CFF(7)-NCFF(7))
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
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0) then
            print*,"step",n, nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"wSWP",wSWP
            print*,"waterfrac",soil_waterfrac
            stop
         endif
      enddo

    end do ! nodays loop

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
    pp = (pn/umol_to_gC)/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))
    ci = min(ci,co2) ! C3 can't have more CO2 than is in the atmosphere
    ci_global = ci
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
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in) :: gs_in

    ! local variables
    double precision :: tmp
    double precision :: gs_high, gs_store, &
                        gpp_high, gpp_low

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)

    ! determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE - ((gpp_high - gpp_low)/lai)

    ! remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_WUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to WUE.

    ! arguments
    double precision, intent(in) :: gs_in

    ! local variables
    double precision :: tmp
    double precision :: gs_high, gs_store, &
                        gpp_high, gpp_low, &
                        evap_high, evap_low

    !!!!!!!!!!
    ! Optimise water use efficiency
    !!!!!!!!!!

    ! Globally stored upper stomatal conductance estimate in memory
    gs_store = stomatal_conductance
    ! now assign the current estimate
    stomatal_conductance = gs_in
    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)
    call calculate_transpiration(evap_low)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! now assign the incremented estimate
    stomatal_conductance = gs_high
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)
    call calculate_transpiration(evap_high)

    ! estimate marginal return on GPP for water loss, less water use efficiency criterion (gC.kgH2O-1.m-2.s-1)
    find_gs_WUE = ((gpp_high - gpp_low)/(evap_high - evap_low)) / lai
    find_gs_WUE = find_gs_WUE - iWUE

    ! return original stomatal value back into memory
    stomatal_conductance = gs_store

    ! remember to return back to the user
    return

  end function find_gs_WUE  !
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
    double precision :: denom
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
        if (stomatal_conductance /= max_gs .or. do_iWUE ) then
            ! If there is a positive demand for water then we will solve for photosynthesis limits on gs through iterative solution
            delta_gs = 1d-3*lai ! mmolH2O/m2leaf/day
            if (do_iWUE) then
                ! intrinsic WUE optimisation
                stomatal_conductance = zbrent('acm_albedo_gc:find_gs_iWUE',find_gs_iWUE,min_gs,stomatal_conductance,tol_gs)
            else
                ! WUE optimisation
                stomatal_conductance = zbrent('acm_albedo_gc:find_gs_WUE',find_gs_WUE,min_gs,stomatal_conductance,tol_gs)
            endif
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
    if (input_temperature < dble_one) then
        lambda = 2.835d6
    else
        lambda = 2501000d0-2364d0*input_temperature
    endif
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
  subroutine calculate_transpiration(transpiration)

    ! Models leaf cnaopy transpiration based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out) :: transpiration ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation & ! isothermal net radiation (W/m2)
                           ,water_supply & ! Potential water supply to canopy from soil (kgH2O.m-2.day-1)
                                  ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! calculate potential water supply (kgH2O.m-2.day-1)
    ! provided potential upper bound on evaporation
    water_supply = max_supply * mmol_to_kg_water

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.s-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
    gs = stomatal_conductance / (convert_ms1_mol_1 * 1d3)
    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kg.m-2.day-1)
    transpiration = (slope*canopy_radiation) + (air_density_kg*cpair*vpd_pa*1d-3*gb)
    ! Calculate the transpiration flux and restrict by potential water supply
    ! over the day
    transpiration = min(water_supply,(transpiration / (lambda*(slope+(psych*(dble_one+gb/gs)))))*dayl_seconds)

  end subroutine calculate_transpiration
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,storage,transpiration)

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(in) :: transpiration      ! kgH2O/m2/day
    double precision, intent(inout) :: storage         ! canopy water storage kgH2O/m2
    double precision, intent(out) :: wetcanopy_evap, & ! kgH2O.m-2.day-1
                                      act_pot_ratio    ! Ratio of potential evaporation to actual

    ! local variables
    double precision :: canopy_radiation, & ! isothermal net radiation (W/m2)
                                      gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * seconds_per_day_1)

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    wetcanopy_evap = (slope*canopy_radiation) + (air_density_kg*cpair*vpd_pa*1d-3*gb)
    ! Calculate the potential wet canopy evaporation, limited by energy used for
    ! transpiration
    wetcanopy_evap = (wetcanopy_evap / (lambda*(slope+psych))) * seconds_per_day !dayl_seconds
    ! substract transpiration from potential surface evaporation
    wetcanopy_evap = wetcanopy_evap - transpiration

    ! dew is unlikely to occur (if we had energy balance) if mint > 0
    if (wetcanopy_evap < dble_zero .and. mint > dble_zero) wetcanopy_evap = dble_zero
    ! Sublimation is unlikely to occur (if we had energy balance) if maxt < 0
    if (wetcanopy_evap > dble_zero .and. maxt < dble_zero) wetcanopy_evap = dble_zero

    ! Remember potential evaporation to later calculation of the potential
    ! actual ratio
    act_pot_ratio = wetcanopy_evap

    ! assuming there is any rainfall, currently water on the canopy or dew formation
    if (rainfall > dble_zero .or. storage > dble_zero .or. wetcanopy_evap < dble_zero) then
        ! Update based on canopy water storage
        call canopy_interception_and_storage(wetcanopy_evap,storage)
    else
        ! there is no water movement possible
        intercepted_rainfall = dble_zero ; wetcanopy_evap = dble_zero
    endif

    ! now calculate the ratio of potential to actual evaporation
    if (act_pot_ratio == dble_zero) then
        act_pot_ratio = dble_zero
    else
        act_pot_ratio = abs(wetcanopy_evap / act_pot_ratio)
    endif

  end subroutine calculate_wetcanopy_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_evaporation(soilevap)

    ! Estimate soil surface evaporation based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out) :: soilevap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: soil_radiation & ! isothermal net radiation (W/m2)
                      ,water_diffusion & ! Diffusion of water through soil matrix (m.s-1)
                                ,esurf & ! see code below
                                 ,esat & ! soil air space saturation vapour pressure
                                  ,gws & ! water vapour conductance through soil air space (m.s-1)
                                   ,Qc

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday * 1d6 * dayl_seconds_1)
    ! estimate ground heat flux from statistical approximation, positive if energy moving up profile
    ! NOTE: linear coefficient estimates from SPA simulations
    Qc = -0.4108826d0 * (maxt - maxt_lag1)
    soil_radiation = soil_radiation + Qc

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * (maxt+freeze) - 4717.306081d0 ) / ( maxt+freeze - 35.86d0 ) )
    air_vapour_pressure = esat - (vpd_pa * 1d-3)

    ! Estimate water diffusion rate (m2.s-1) Jones (2014) appendix 2
    water_diffusion = 24.2d-6 * ( (maxt+freeze) / 293.2d0 )**1.75d0
    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_diffusion / (tortuosity*drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water
    esurf = esat * exp( 1d6 * SWP(1) * partial_molar_vol_water / ( Rcon * (maxt+freeze) ) )
    ! now difference in vapour pressure between soil and canopy air spaces
    esurf = esurf - air_vapour_pressure

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
    soilevap = (slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)
    soilevap = soilevap / (lambda*(slope+(psych*(dble_one+soil_conductance/gws))))
    soilevap = soilevap * dayl_seconds

    ! dew is unlikely to occur (if we had energy balance) if mint > 0
!    if (soilevap < dble_zero .and. mint > dble_one) soilevap = dble_zero
    ! Sublimation is unlikely to occur (if we had energy balance) if maxt < 0
!    if (soilevap > dble_zero .and. maxt < dble_one) soilevap = dble_zero

  end subroutine calculate_soil_evaporation
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

    ! calculate soil surface conductance
    call calculate_soil_conductance(mixing_length_momentum)

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
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer        :: i
    double precision :: x1, x2

    x1 = 0.1d0 ; x2 = 0.7d0 ! low/high guess
    do i = 1 , nos_soil_layers+1
       water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns , x1 , x2 , 0.001d0 )
    enddo

  end subroutine calculate_field_capacity
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

    ! declination
!    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
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
    sunrise = 12 - nint(dayl_hours*0.5d0) ; sunset = sunrise + nint(dayl_hours)

    ! estimate the solar cosine zenith angle for 12 noon
    cos_solar_zenith_angle = sinld + cosld

    ! return to user
    return

  end subroutine calculate_daylength
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
    double precision :: balance                    &
                       ,absorbed_nir_fraction_soil &
                       ,absorbed_par_fraction_soil &
                       ,fsnow                      &
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

    ! local parameters
    double precision, parameter :: newsnow_nir_abs = 0.27d0 & ! NIR absorption fraction
                                  ,newsnow_par_abs = 0.05d0   ! PAR absorption fraction

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
    if (snow_storage > dble_zero) then
        fsnow = dble_one - exp( - snow_storage * 1d-2 )  ! fraction of snow cover on the ground
        absorbed_par_fraction_soil = ((dble_one - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_par_abs)
        absorbed_nir_fraction_soil = ((dble_one - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_nir_abs)
    else
        absorbed_par_fraction_soil = soil_swrad_absorption
        absorbed_nir_fraction_soil = soil_swrad_absorption
    endif

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

    ! check energy balance
    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
!    if ((balance - swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual = ",balance,"SW in = ",swrad
!    endif

  end subroutine calculate_shortwave_balance
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
    water_flux = dble_zero ; wSWP = dble_zero
    ratio = dble_zero ; ratio(1) = dble_one ; root_mass = dble_zero
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
!    transpiration_resistance = (gplant * lai)**(-dble_one)
    transpiration_resistance = canopy_height / (gplant * max(min_lai,lai))

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

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
!    root_length = root_mass / (root_density * root_cross_sec_area)

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
    demand = abs(minlwp-SWP(1:nos_root_layers))+head*canopy_height
    ! now loop through soil layers, where root is present
    do i = 1, nos_root_layers
      if (root_mass(i) > dble_zero) then
          ! if there is root then there is a water flux potential...
          root_reach_local = min(root_reach,layer_thickness(i))
          ! calculate and accumulate steady state water flux in mmol.m-2.s-1
          water_flux(1) = plant_soil_flow(i,root_length(i),root_mass(i) &
                         ,demand(i),root_reach_local,transpiration_resistance)
      else
          ! ...if there is not then we wont have any below...
          exit
      end if ! root present in current layer?
    end do ! nos_root_layers
    ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))

    ! if freezing then assume soil surface is frozen
    if (meant < dble_one) then
        water_flux(1) = dble_zero
        ratio(1) = dble_zero
        ratio(2:nos_root_layers) = layer_thickness(2:nos_root_layers) / sum(layer_thickness(2:nos_root_layers))
    endif
    ! calculate sum value
    sum_water_flux = sum(water_flux)

    ! calculate weighted SWP and uptake fraction
    wSWP = sum(SWP(1:nos_root_layers) * water_flux(1:nos_root_layers))
    uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
    wSWP = wSWP / sum_water_flux

    ! sanity check in case of zero flux
    if (sum_water_flux <= vsmall) then
        wSWP = -20d0
        uptake_fraction = dble_zero ; uptake_fraction(1) = dble_one
    endif

    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum(water_flux)

    ! finally convert transpiration flux (mmol.m-2.s-1)
    ! into kg.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(potential_evaporation,storage)

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

    ! arguments
    double precision, intent(inout) :: storage, & ! canopy water storage (kgH2O/m2)
                         potential_evaporation    ! wet canopy evaporation (kgH2O.m-2.day-1),
                                                  ! enters as potential but leaves as water balance adjusted
    ! local variables
    integer :: i, hr
    double precision :: a, through_fall, max_storage, max_storage_1, daily_addition, wetcanopy_evaporation &
                       ,potential_drainage_rate ,drain_rate, evap_rate, initial_canopy, co_mass_balance, dx, tmp1, tmp2, tmp3
    ! local parameters
    double precision, parameter :: CanIntFrac = -0.5d0,  & ! Coefficient scaling rainfall interception fraction with LAI
                                  CanStorFrac = 0.1d0,   & ! Coefficient scaling canopy water storage with LAI
                                 RefDrainRate = 0.002d0, & ! Reference drainage rate (mm/min; Rutter et al 1975)
                                  RefDrainLAI = 1.05d0,  & ! Reference drainage LAI (m2/m2; Rutter et al 1975)
                                 RefDrainCoef = 3.7d0,   & ! Reference drainage Coefficient (Rutter et al 1975)
                               RefDrainCoef_1 = RefDrainCoef ** (-dble_one)

    ! hold initial canopy storage in memory
    initial_canopy = storage
    ! determine maximum canopy storage & through fall fraction
    through_fall = max(min_throughfall,exp(CanIntFrac*lai))
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood
    max_storage = max(min_storage,CanStorFrac*lai) ; max_storage_1 = max_storage**(-dble_one)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    intercepted_rainfall = rainfall * (dble_one - through_fall)

    ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
    ! 0.002 is canopy specific coefficient modified by 0.002*(max_storage/1.05)
    ! where max_storage is the canopy maximum capacity (mm) (LAI based) and
    ! 1.05 is the original canopy capacitance
    a = log( RefDrainRate * ( max_storage / RefDrainLAI ) ) - RefDrainCoef * max_storage

    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = intercepted_rainfall * seconds_per_day

    ! reset cumulative variables
    through_fall = dble_zero ; wetcanopy_evaporation = dble_zero
    drain_rate = dble_zero ; evap_rate = dble_zero

    ! deal with rainfall additions first
    do i = 1, int(days_per_step)

       ! add rain to the canopy and overflow as needed
       storage = storage + daily_addition

       if (storage > max_storage) then

           if (potential_evaporation > dble_zero) then

               ! assume co-access to available water above max_storage by both drainage and
               ! evaporation. Water below max_storage is accessable by evaporation only.

               ! Trapezium rule for approximating integral of drainage rate.
               ! Allows estimation of the mean drainage rate between starting point and max_storage,
               ! thus the time period appropriate for co-access can be quantified
               dx = storage - ((storage + max_storage)*0.5d0)
               tmp1 = exp(a + (RefDrainCoef * storage))
               tmp2 = exp(a + (RefDrainCoef * max_storage))
               tmp3 = exp(a + (RefDrainCoef * (storage+dx)))
               potential_drainage_rate = 0.5d0 * dx * ((tmp1 + tmp2) + 2d0 * tmp3)
               potential_drainage_rate = potential_drainage_rate * 1440d0

               ! restrict evaporation and drainage to the quantity above max_storage
               evap_rate = potential_evaporation ; drain_rate = min(potential_drainage_rate,storage-max_storage)

               ! limit based on available water if total demand is greater than excess
               co_mass_balance = ((storage-max_storage) / (evap_rate + drain_rate))
               evap_rate = evap_rate * co_mass_balance ; drain_rate = drain_rate * co_mass_balance

               ! estimate evaporation from remaining water, less that already removed from storage and evaporation energy used
               evap_rate = evap_rate + min(potential_evaporation - evap_rate, storage - evap_rate - drain_rate)

           else

               ! load dew formation to the current local evap_rate variable
               evap_rate = potential_evaporation
               ! restrict drainage the quantity above max_storage, adding dew formation too
               drain_rate = (storage - evap_rate) - max_storage

           endif

       else
           ! no drainage just apply evaporation / dew formation fluxes directly
           evap_rate = potential_evaporation
           drain_rate = dble_zero
           if (evap_rate > dble_zero) then
               ! evaporation restricted by fraction of surface actually covered
               ! in water
               evap_rate = evap_rate * min(dble_one,storage * max_storage_1)
               ! and the total amount of water
               evap_rate = min(evap_rate,storage)
           else
               ! then dew formation has occurred, if this pushes storage > max_storage add it to drainage
               drain_rate = max(dble_zero,(storage - evap_rate) - max_storage)
           endif ! evap_rate > 0
       endif ! storage > max_storage

       ! update canopy storage with water flux
       storage = max(dble_zero,storage - evap_rate - drain_rate)
       wetcanopy_evaporation = wetcanopy_evaporation + evap_rate
       through_fall = through_fall + drain_rate

    end do ! days

    ! correct intercepted rainfall rate to kgH2O.m-2.s-1
    intercepted_rainfall = intercepted_rainfall - ((through_fall * days_per_step_1) * seconds_per_day_1)

    ! sanity checks; note 1e-8 prevents precision errors causing flags
    if (intercepted_rainfall > rainfall .or. storage < dble_zero &
   .or. (wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day)) ) then
       print*,"Condition 1",intercepted_rainfall > rainfall
       print*,"Condition 2",storage < dble_zero
       print*,"Condition 3",(wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day))
       print*,"storage (kgH2O/m2)",storage,"max_storage (kgH2O/m2)",max_storage,"initial storage (kgH2O/m2)", initial_canopy
       print*,"rainfall (kgH2O/m2/day)", rainfall*seconds_per_day, "through_fall (kgH2O/m2/day)", (through_fall * days_per_step_1)
       print*,"through_fall_total (kgH2O/m2/step)",through_fall
       print*,"potential_evaporation (kgH2O/m2/day)",potential_evaporation
       print*,"actual evaporation    (kgH2O/m2/day)",wetcanopy_evaporation * days_per_step_1
       stop
    endif

    ! average evaporative flux to daily rate (kgH2O/m2/day)
    potential_evaporation = wetcanopy_evaporation * days_per_step_1

    ! final clearance of canopy storage of version small values at the level of system precision
    if (storage < 10d0*vsmall) storage = dble_zero

  end subroutine canopy_interception_and_storage
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_update_soil_water(ET_leaf,ET_soil,rainfall_in,corrected_ET)

   !
   ! Function updates the soil water status and layer thickness
   ! Soil water profile is updated in turn with evaporative losses,
   ! rainfall infiltration and gravitational drainage
   ! Root layer thickness is updated based on changes in the rooting depth from
   ! the previous step
   !

   implicit none

   ! arguments
   double precision, intent(in) :: ET_leaf,ET_soil & ! evapotranspiration estimate (kgH2O.m-2.day-1)
                                      ,rainfall_in   ! rainfall (kgH2O.m-2.day-1)
   double precision, intent(out) :: corrected_ET     ! water balance corrected evapotranspiration (kgH2O/m2/day)

   ! local variables
   integer :: day
   double precision ::  depth_change, water_change
   double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses

   ! reset soil water exchanges
   underflow = dble_zero ; runoff = dble_zero ; corrected_ET = dble_zero

   ! to allow for smooth water balance integration carry this out at daily time step
   do day = 1, nint(days_per_step)

      !!!!!!!!!!
      ! Evaporative losses
      !!!!!!!!!!

      ! Assume leaf transpiration is drawn from the soil based on the
      ! update_fraction estimated in calculate_Rtot
      evaporation_losses = ET_leaf * uptake_fraction
      ! Assume all soil evaporation comes from the soil surface only
      evaporation_losses(1) = evaporation_losses(1) + ET_soil
      ! can not evaporate from soil more than is available (m -> mm)
      avail_flux = soil_waterfrac(1:nos_root_layers) * layer_thickness(1:nos_root_layers) * 1d3
      where (evaporation_losses > avail_flux) evaporation_losses = avail_flux * 0.999d0

      ! this will update the ET estimate outside of the function
      ! days_per_step corrections happens outside of the loop below
      corrected_ET = corrected_ET + sum(evaporation_losses)

      ! pass information to waterloss variable and zero watergain
      ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
      waterloss = dble_zero ; watergain = dble_zero
      waterloss(1:nos_root_layers) = evaporation_losses(1:nos_root_layers)*1d-3
      ! update soil water status with evaporative losses
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                           + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = dble_zero ; watergain = dble_zero

      !!!!!!!!!!
      ! Gravitational drainage
      !!!!!!!!!!

      ! determine drainage flux between surface -> sub surface and sub surface
      call gravitational_drainage

      ! update soil water status with drainage
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                           + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = dble_zero ; watergain = dble_zero

      !!!!!!!!!!
      ! Rainfall infiltration drainage
      !!!!!!!!!!

      ! determine infiltration from rainfall (kgH2O/m2/step),
      ! if rainfall is probably liquid / soil surface is probably not frozen
      if (rainfall_in > dble_zero) then
          call infiltrate(rainfall_in)
      else
          runoff = runoff + (rainfall_in * days_per_step_1)
      endif ! is there any rain to infiltrate?
      ! update soil profiles. Convert fraction into depth specific values (rather than m3/m3) then update fluxes
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                           + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = dble_zero ; watergain = dble_zero

      ! mass balance check, at this point do not try and adjust evaporation to
      ! correct for lack of supply. Simply allow for drought in next time step
      ! instead...
      where (soil_waterfrac <= dble_zero)
             soil_waterfrac = vsmall
      end where

   end do ! days_per_step

   ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
   corrected_ET = corrected_ET * days_per_step_1
   underflow = underflow * days_per_step_1
   runoff = runoff * days_per_step_1

   !!!!!!!!!!
   ! Update soil layer thickness
   !!!!!!!!!!

   depth_change = dble_zero ; water_change = dble_zero
   ! if roots extent down into the bucket
   if (root_reach > (top_soil_depth+mid_soil_depth) .or. previous_depth > (top_soil_depth+mid_soil_depth)) then
      ! how much has root depth extended since last step?
      depth_change = root_reach - previous_depth

      ! if there has been an increase
      if (depth_change > dble_zero .and. root_reach > sum(layer_thickness(1:2))+min_layer) then

         ! determine how much water is within the new volume of soil
         water_change = soil_waterfrac(nos_soil_layers) * depth_change
         ! now assign that new volume of water to the deep rooting layer
         soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers) * layer_thickness(nos_root_layers)) &
                                            + water_change) / (layer_thickness(nos_root_layers)+depth_change)
         ! explicitly update the soil profile if there has been rooting depth
         ! changes
         layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
         layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
         layer_thickness(4) = max_depth - sum(layer_thickness(1:3))


      elseif (depth_change < dble_zero .and. root_reach > layer_thickness(1)+min_layer) then

         ! determine how much water is lost from the old volume of soil
         water_change = soil_waterfrac(nos_root_layers) * abs(depth_change)
         ! now assign that new volume of water to the deep rooting layer
         soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers) * layer_thickness(nos_soil_layers)) &
                                            + water_change) / (layer_thickness(nos_soil_layers)+abs(depth_change))

         ! explicitly update the soil profile if there has been rooting depth
         ! changes
         layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
         layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
         layer_thickness(4) = max_depth - sum(layer_thickness(1:3))

      else

         ! we don't want to do anything, just recycle the previous depth

      end if ! depth change

   end if ! root reach beyond top layer

   ! in all cases keep track of the previous rooted depth
   previous_depth = root_reach

   ! finally update soil water potential
   call soil_water_potential

!   ! sanity check for catastrophic failure
!   do soil_layer = 1, nos_soil_layers
!      if (soil_waterfrac(soil_layer) < 0d0 .and. soil_waterfrac(soil_layer) > -0.01d0) then
!          soil_waterfrac(soil_layer) = 0d0
!      endif
!      if (soil_waterfrac(soil_layer) < 0d0 .or. soil_waterfrac(soil_layer) /= soil_waterfrac(soil_layer)) then
!         print*,'ET',ET,"rainfall",rainfall_in
!         print*,'evaporation_losses',evaporation_losses
!         print*,"watergain",watergain
!         print*,"waterloss",waterloss
!         print*,'depth_change',depth_change
!         print*,"soil_waterfrac",soil_waterfrac
!         print*,"porosity",porosity
!         print*,"layer_thicknes",layer_thickness
!         print*,"Uptake fraction",uptake_fraction
!         print*,"max_depth",max_depth,"root_k",root_k,"root_reach",root_reach
!         print*,"fail" ; stop
!      endif
!   end do

   ! explicit return needed to ensure that function runs all needed code
   return

  end subroutine calculate_update_soil_water
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kg.m-2.day-1)

    ! local argumemts
    integer :: i
    double precision    :: add   & ! surface water available for infiltration (m)
                          ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in * 1d-3

    do i = 1 , nos_soil_layers
       ! determine the available pore space in current soil layer
       wdiff = max(dble_zero,(porosity(i)-soil_waterfrac(i))*layer_thickness(i)-watergain(i)+waterloss(i))
       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer
       if (add > wdiff) then
          ! if so fill and subtract from input and move on to the next layer
          watergain(i) = watergain(i) + wdiff
          add = add - wdiff
       else
          ! otherwise infiltate all in the current layer
          watergain(i) = watergain(i) + add
          add = dble_zero
       end if
       ! if we have added all available water we are done
       if (add <= dble_zero) then
           add = dble_zero
           exit
       end if

    end do ! nos_soil_layers

    ! if after all of this we have some water left assume it is runoff (kgH2O.m-2.day-1)
    ! NOTE that runoff is reset outside of the daily soil loop
    runoff = runoff + (add * 1d3)

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage

    ! integrator for soil gravitational drainage !

    implicit none

    ! local variables..
    integer :: d, nos_integrate
    double precision  :: liquid & ! liquid water in local soil layer (m3/m3)
                    ,drainlayer & ! field capacity of local soil layer (m3/m3)
                         ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3)
                        ,change & ! absolute volume of water drainage in current layer (m3)
                      ,drainage & ! drainage rate of current layer (m/day)
                   ,local_drain & ! drainage of current layer (m/nos_minutes)
      ,iceprop(nos_soil_layers)

    ! local parameters
    integer, parameter :: nos_hours_per_day = 1440, nos_minutes = 360

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = dble_zero
    ! except the surface layer in the mean daily temperature is < 0oC
    if (meant < dble_one) iceprop(1) = dble_one

    do soil_layer = 1, nos_soil_layers

       ! soil water capacity of the current layer
       drainlayer = field_capacity( soil_layer )
       ! liquid content of the soil layer
       liquid     = soil_waterfrac( soil_layer ) &
                  * ( dble_one - iceprop( soil_layer ) )

       ! initial conditions; i.e. is there liquid water and more water than
       ! layer can hold
       if ( liquid > drainlayer ) then

          ! unsaturated volume of layer below (m3 m-2)..
          unsat = max( dble_zero , ( porosity( soil_layer+1 ) - soil_waterfrac( soil_layer+1 ) ) &
                             * layer_thickness( soil_layer+1 ) / layer_thickness( soil_layer ) )

          d = 1 ; nos_integrate = nos_hours_per_day / nos_minutes
          drainage = dble_zero ; local_drain = dble_zero
          do while (d <= nos_integrate .and. liquid > drainlayer)
              ! estimate drainage rate (m/s)
              call calculate_soil_conductivity(soil_layer,liquid,local_drain)
              ! scale to total number of seconds in increment
              local_drain = local_drain * dble(nos_minutes * 60)
              local_drain = min(liquid-drainlayer,local_drain)
              liquid = liquid - local_drain
              drainage = drainage + local_drain
              d = d + 1
          end do ! integrate over time

          ! layer below cannot accept more water than unsat
          if ( drainage > unsat ) drainage = unsat
          ! water loss from this layer (m3)
          change = drainage * layer_thickness(soil_layer)
          ! update soil layer below with drained liquid
          watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
          waterloss( soil_layer     ) = waterloss( soil_layer     ) + change

       end if ! some liquid water and drainage possible

    end do ! soil layers

    ! estimate drainage from bottom of soil column (kgH2O/m2/day)
    ! NOTES: that underflow is reset outside of the daily soil loop
    underflow = underflow + (waterloss(nos_soil_layers) * 1d3)

  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay,soil_frac_sand)

   ! Porosity is estimated from Saxton equations. !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    double precision, parameter :: H = 0.332d0, &
                                   J = -7.251d-4, &
                                   K = 0.1276d0

    ! loop over soil layers..
    porosity(1:nos_soil_layers) = H + J * soil_frac_sand(1:nos_soil_layers) + &
                                  K * log10(soil_frac_clay(1:nos_soil_layers))
    ! then assign same to core layer
    porosity(nos_soil_layers+1) = porosity(nos_soil_layers)

  end subroutine soil_porosity
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay,soil_frac_sand)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i

    ! include some hardcoded boundaries for the Saxton equations
    where (soil_frac_sand < 5d0) soil_frac_sand = 5d0
    where (soil_frac_clay < 5d0) soil_frac_clay = 5d0
    where (soil_frac_clay > 60d0) soil_frac_clay = 60d0
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity

    ! final sanity check for porosity
    where (porosity <= field_capacity) porosity = field_capacity + 0.01d0

  end subroutine initialise_soils
  !
  !---------------------------------------------------------------------
  !
  subroutine update_soil_initial_conditions(input_soilwater_frac)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision :: input_soilwater_frac

    ! local variables
    integer :: i

    ! Default assumption to be field capacity
    soil_waterfrac = field_capacity

    ! if prior value has been given
    if (input_soilwater_frac > -9998d0) then
        ! calculate initial soil water fraction
        soil_waterfrac(1:nos_soil_layers) = input_soilwater_frac
    endif

    ! calculate initial soil water potential
    SWP = dble_zero
    call soil_water_potential
    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers
    ! but apply the lowest soil layer to the core as well in initial conditions
    soil_conductivity(nos_soil_layers+1) = soil_conductivity(nos_soil_layers)

  end subroutine update_soil_initial_conditions
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity(soil_layer,waterfrac,conductivity)

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

    ! arguments
    integer, intent(in) :: soil_layer
    double precision, intent(in) :: waterfrac
    double precision, intent(out) :: conductivity

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    conductivity = cond1(soil_layer) * exp(cond2(soil_layer)+cond3(soil_layer)/waterfrac)

    ! protection against floating point error
    if (waterfrac < 0.05d0) conductivity = 1d-30

  end subroutine calculate_soil_conductivity
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay,soil_frac_sand)

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    double precision, parameter :: A = -4.396d0,  B = -0.0715d0,   CC = -4.880d-4, D = -4.285d-5, &
                                   E = -3.140d0,  F = -2.22d-3,     G = -3.484d-5, H = 0.332d0,   &
                                   J = -7.251d-4, K = 0.1276d0,     P = 12.012d0,  Q = -7.551d-2, &
                                   R = -3.895d0,  T = 3.671d-2,     U = -0.1103d0, V = 8.7546d-4, &
                                   mult1 = 100d0, mult2 = 2.778d-6

    ! layed out in this manor to avoid memory management issues in module
    ! variables
    potA(1:nos_soil_layers) = A + (B * soil_frac_clay) + &
                            (CC * soil_frac_sand * soil_frac_sand) + &
                            (D * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    potA(1:nos_soil_layers) = exp(potA(1:nos_soil_layers))
    potA(1:nos_soil_layers) = potA(1:nos_soil_layers) * mult1

    potB(1:nos_soil_layers) = E + (F * soil_frac_clay * soil_frac_clay) + &
                             (G * soil_frac_sand * soil_frac_sand * soil_frac_clay)

    cond1(1:nos_soil_layers) = mult2
    cond2(1:nos_soil_layers) = P + (Q * soil_frac_sand)
    cond3(1:nos_soil_layers) = R + (T * soil_frac_sand) + (U * soil_frac_clay) + &
                              (V * soil_frac_clay * soil_frac_clay)

    ! assign bottom of soil column value to core
    potA(nos_soil_layers+1)  = potA(nos_soil_layers)
    potB(nos_soil_layers+1)  = potB(nos_soil_layers)
    cond1(nos_soil_layers+1) = mult2
    cond2(nos_soil_layers+1) = cond2(nos_soil_layers)
    cond3(nos_soil_layers+1) = cond3(nos_soil_layers)

  end subroutine saxton_parameters
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_conductance(lm)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: lm

    ! local variables
    double precision :: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0 ! foliage drag coefficient

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht=vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = (((foliage_drag*canopy_height*max(min_lai,lai))/lm)**0.5d0)

    ! approximation of integral for soil resistance
    soil_conductance = canopy_height/(canopy_decay*Kh_canht) &
                     * (exp(canopy_decay*(dble_one-(soil_roughl/canopy_height)))- &
                        exp(canopy_decay*(dble_one-((roughl+displacement)/canopy_height))))

    ! convert resistance (s.m-1) to conductance (m.s-1)
    soil_conductance = soil_conductance ** (-dble_one)

  end subroutine calculate_soil_conductance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_potential

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    where (SWP(1:nos_soil_layers) < -20d0) SWP(1:nos_soil_layers) = -20d0

  end subroutine soil_water_potential
  !
  !------------------------------------------------------------------
  !
!  subroutine update_net_radiation(isothermal,tempC,area_scaling,act_pot_ratio &
!                                 ,sfc_exchange,aero_exchange,vapour_gradient,deltaTemp,deltaR)
!
!    ! Use steady state solution of evaporation, convective (sensible) and radiative heat loss
!    ! to update isothermal net radiation to net.
!    ! Area scaling (e.g. lai) is an input to allow for common useage for soil (neglecting ground heat) and canopy
!
!    ! arguments
!    double precision, intent(in) ::      tempC, & ! input surface / air temperature (oC)
!                                    isothermal, & ! isothermal net radiation (SW+LW; W/m2)
!                                  area_scaling, & ! area scaling to apply (m2/m2)
!                                 act_pot_ratio, & ! ratio of potential to actual evaporation, i.e. (avail / potenial)
!                                  sfc_exchange, & ! surface exchange conductance (m/s; e.g. stomatal conductance)
!                                 aero_exchange, & ! aerodynamic exchange conductance (m/s; e.g. aerodynamic conductance)
!                               vapour_gradient    ! vapour pressure gradient (Pa; either VPD or between air and soil)
!    double precision, intent(out) :: deltaTemp, & ! surface temperature difference (K)
!                                     deltaR    ! surface longwave radiation difference (W/m2); subtract from isothermal longwave
!
!    ! local variables
!    double precision ::  tempK, & ! ambient temperature as K
!          heat_loss_resistance, & ! resistance to heat loss from radiative and convection (s/m)
!        aerodynamic_resistance, & ! aerodynamic resistance to water or heat exchangce (s/m)
!           stomatal_resistance, & ! stomatal resistance to water exchangce (s/m)
!              water_resistance, & ! serial combination of resistances to water evaporation
! thermal_gains, thermal_losses
!
!    ! ambient temperature C -> K
!    tempK = tempC + freeze
!
!    !
!    ! Calculate resistance to heat loss (s/m)
!    !
!
!    ! First estimate radiative loss term, initially calculated as conductance)
!    heat_loss_resistance = 4d0 * emissivity * boltz * tempK ** 3 / (air_density_kg * cpair)
!    ! Combine in parallel with convective conductances with area correction
!    heat_loss_resistance = heat_loss_resistance + (2d0 * aero_exchange / area_scaling)
!    ! convert from conductance m/s to s/m
!    heat_loss_resistance = heat_loss_resistance ** (-1d0)
!
!    !
!    ! Convert aerodynamic and stomatal conductances to reisistance of water flux
!    !
!
!    aerodynamic_resistance = (aero_exchange/area_scaling) ** (-1d0)
!    if (sfc_exchange == dble_zero) then
!        ! if being used for surface water flux
!        stomatal_resistance = dble_zer!o
!    else
!        ! if used for transpiration
!        stomatal_resistance = (sfc_exchange/area_scaling) ** (-1d0)
!    endif
!
!    !
!    ! Estimate thermal gains and losses (K) to calculate temperature difference
!    !
!
!    water_resistance = (aerodynamic_resistance + stomatal_resistance)
!    thermal_gains = (heat_loss_resistance * water_resistance * psych * (isothermal/area_scaling)) &
!                  / (air_density_kg * cpair * ((psych*water_resistance) + (slope*heat_loss_resistance)))
!    thermal_losses = (heat_loss_resistance * vapour_gradient * 1d-3) &
!                   / ((psych*water_resistance) + (slope*heat_loss_resistance))
!    ! determine surface temperature difference (K); should be added to the canopy temperature
!    deltaTemp = thermal_gains - thermal_losses
!    ! apply actual potential ratio to scale wet surface evaporation when the
!    ! supply of water is limited
!    deltaTemp = deltaTemp * act_pot_ratio
!
!    ! estimate update between isothermal to net radiation (W/m2), including area correction
!    ! note that this MUST be added from the longwave component outside of this function
!    deltaR = -4d0 * emissivity * boltz * tempK ** 3 * ( deltaTemp )
!    deltaR = deltaR * area_scaling
!
!    ! return to user
!    return
!
!  end subroutine update_net_radiation
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
  subroutine calculate_leaf_dynamics(current_step,deltat,nodays           &
                                    ,Tfac_min,Photofac_min,VPDfac_min     &
                                    ,Tfac_range_1,Photofac_range_1        &
                                    ,VPDfac_range_1,pot_leaf_fall         &
                                    ,pot_leaf_growth,mean_min_airt        &
                                    ,mean_daylength,mean_vpd,deltaWP,Rtot &
                                    ,GPP_current,foliage,gpp_crit_frac    &
                                    ,GSI,leaf_fall,leaf_growth)

      ! Subroutine determines whether leaves are growing or dying.
      ! 1) Calculate the Growing Season Index (GSI)
      ! 2) Determines whether conditions are improving or declining
      ! 3) Performes marginal return calculation

      ! GSI added by JFE and TLS.
      ! Refs Jolly et al., 2005, doi: 10.1111/j.1365-2486.2005.00930.x)
      !      Stoeckli et al., 2010, doi:10.1029/2010JG001545.

      implicit none

      ! declare arguments
      integer, intent(in) :: nodays, current_step
      double precision, intent(in) :: deltat(nodays) & !
                                            ,foliage & !
                                        ,GPP_current & !
                                      ,gpp_crit_frac & !
                                      ,mean_min_airt & !
                                     ,mean_daylength & !
                                           ,mean_vpd & !
                                            ,deltaWP & !
                                               ,Rtot &
                                           ,Tfac_min & !
                                       ,Photofac_min & !
                                         ,VPDfac_min & !
                                       ,Tfac_range_1 & !
                                   ,Photofac_range_1 & !
                                     ,VPDfac_range_1 & !
                                      ,pot_leaf_fall & !
                                    ,pot_leaf_growth

      double precision, intent(inout) :: GSI(nodays) &
                                        ,leaf_fall,leaf_growth

      ! declare local variables
      integer :: gsi_lag, m
      double precision :: infi      &
                         ,tmp       &
                         ,deltaGPP  &
                         ,lai_save  &
                         ,canopy_lw_save &
                         ,canopy_sw_save &
                         ,canopy_par_save &
                         ,soil_lw_save &
                         ,soil_sw_save, gs_save

      ! save original values for re-allocation later
      canopy_lw_save = canopy_lwrad_Wm2 ; soil_lw_save = soil_lwrad_Wm2
      canopy_sw_save = canopy_swrad_MJday ; canopy_par_save  = canopy_par_MJday
      soil_sw_save = soil_swrad_MJday ; gs_save = stomatal_conductance
      gsi_lag = gsi_lag_remembered
      lai_save = lai

      ! for infinity checks
      infi = 0d0

      ! It is the product of 3 limiting factors for temperature, photoperiod and
      ! vapour pressure deficit that grow linearly from 0 to 1 between a
      ! calibrated min and max value.
      ! Photoperiod, VPD and avgTmin are direct input

      ! temperature limitation, then restrict to 0-1; correction for k-> oC
      Tfac = (mean_min_airt-(Tfac_min-freeze)) * Tfac_range_1
      Tfac = min(dble_one,max(dble_zero,Tfac))
      ! photoperiod limitation
      Photofac = (mean_daylength-Photofac_min) * Photofac_range_1
      Photofac = min(dble_one,max(dble_zero,Photofac))
      ! VPD limitation (kPa)
      VPDfac = dble_one - ((mean_VPD-VPDfac_min) * VPDfac_range_1)
      VPDfac = min(dble_one,max(dble_zero,VPDfac))

      ! calculate and store the GSI index
      GSI(current_step) = Tfac*Photofac*VPDfac

      ! we will load up some needed variables
      m = nint(tmp_m(current_step))
      ! update gsi_history for the calculation
      if (current_step == 1) then
          ! in first step only we want to take the initial GSI value only
          gsi_history(gsi_lag) = GSI(current_step)
      else
          gsi_history((gsi_lag-m):gsi_lag) = GSI((current_step-m):current_step)
      endif
      ! calculate gradient
      gradient = linear_model_gradient(tmp_x(1:(gsi_lag)),gsi_history(1:gsi_lag),gsi_lag)
      ! adjust gradient to daily rate
      gradient = gradient / dble(nint((sum(deltat((current_step-m+1):current_step))) / dble(gsi_lag-1)))
      gsi_lag_remembered = gsi_lag

      ! first assume that nothing is happening
      leaf_fall = dble_zero   ! leaf turnover
      leaf_growth = dble_zero ! leaf growth

      ! everything else in here was needed to keep track of GSI values but
      ! ultimately if there is not labile available no growth can occur and loss
      ! should have been managed else where as mortality
      if (avail_labile > dble_zero) then

          ! now update foliage and labile conditions based on gradient calculations
          if (gradient <= fol_turn_crit .or. GSI(current_step) == dble_zero) then

             ! we are in a decending condition so foliar turnover
             leaf_fall = pot_leaf_fall*(dble_one-GSI(current_step))
             just_grown = 0.5d0

          else if (gradient >= lab_turn_crit .and. deltaWP < dble_zero) then

             ! we are in an assending condition so labile turnover
             leaf_growth = pot_leaf_growth*GSI(current_step)
             just_grown = 1.5d0

             ! calculate potential C allocation to leaves
             tmp = avail_labile * &
                   min(dble_one,dble_one-(dble_one-leaf_growth)**deltat(current_step))*deltat_1(current_step)
             ! calculate new leaf area, GPP return
             lai = (foliage+tmp) * SLA
             call calculate_shortwave_balance 
             ! calculate stomatal conductance of water
             if (lai > vsmall .and. stomatal_conductance > vsmall) then
                 tmp = acm_gpp(stomatal_conductance)
             else
                 tmp = dble_zero
             endif
             deltaGPP = tmp - GPP_current
             ! is the marginal return for GPP (over the mean life of leaves)
             ! less than increase in maintenance respiration and C required to
             ! growth?

             if (deltaGPP < gpp_crit_frac*GPP_current) leaf_growth = dble_zero

          else if (gradient < lab_turn_crit .and. gradient > fol_turn_crit .and. &
                   deltaWP < dble_zero ) then

             ! probaly we want nothing to happen,

             ! However if we are at the seasonal
             ! maximum we will consider further growth still
             if (just_grown >= dble_one) then

                ! we have recently grown so we will not be losing leaves, but we
                ! might want to grow some more depending on the marginal return

                leaf_growth = pot_leaf_growth*GSI(current_step)

                ! calculate potential C allocation to leaves
                tmp = avail_labile * &
                      min(dble_one,dble_one-(dble_one-leaf_growth)**deltat(current_step))*deltat_1(current_step)
                ! calculate new leaf area, GPP return
                lai = (foliage+tmp) * SLA
                call calculate_shortwave_balance
                ! calculate stomatal conductance of water
                if (lai > vsmall .and. stomatal_conductance > vsmall) then
                    tmp = acm_gpp(stomatal_conductance)
                else
                    tmp = dble_zero
                endif
                deltaGPP = tmp - GPP_current
                ! is the marginal return for GPP (over the mean life of leaves)
                ! less than increase in maintenance respiration and C required to
                ! growth?

                if (deltaGPP < gpp_crit_frac*GPP_current) leaf_growth = dble_zero

             end if ! Just grown?

          endif ! gradient choice

      endif ! avail_labile > 0

      ! restore original value back from memory
      lai = lai_save
      canopy_lwrad_Wm2 = canopy_lw_save ; soil_lwrad_Wm2 = soil_lw_save
      canopy_swrad_MJday = canopy_sw_save ; canopy_par_MJday = canopy_par_save
      soil_swrad_MJday = soil_sw_save ; stomatal_conductance = gs_save

  end subroutine calculate_leaf_dynamics
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wood_root_growth(n,lab_to_roots,lab_to_wood &
                                       ,deltaWP,Rtot,current_gpp,Croot,Cwood  &
                                       ,root_growth,wood_growth)
    implicit none

    ! Premise of wood and root phenological controls

    ! Assumption 1:
    ! Based on plant physiology all cell expansion can only occur if there is
    ! sufficient water pressure available to drive the desired expansion.
    ! Moreover, as there is substantial evidence that shows wood and root growth
    ! do not follow the same phenology as leaves or GPP availability.
    ! Therefore, their phenological constrols should be separate from both that of
    ! the GSI model driving canopy phenology or GPP. Wood growth is limited by a
    ! logisitic temperature response assuming <5 % growth potential at 5oC and
    ! >95 % growth potential at 30 oC. Wood growth is also limited by a
    ! logistic response to water availability. When deltaWP (i.e. minleaf-wSWP)
    ! is less than -1 MPa wood growth is restricted to <5 % of potential.
    ! See review Fatichi et al (2013). Moving beyond phtosynthesis from carbon
    ! source to sink driven vegetation modelling. New Phytologist,
    ! https://doi.org/10.1111/nph.12614 for further details.

    ! As with wood, root phenology biologically speaking is independent of
    ! observed foliar phenological dynamics and GPP availabilty and thus has
    ! a separate phenology model. Similar to wood, a logistic temperature
    ! response is applied such that root growth is <5 % of potential at 0oC and
    ! >95 % of potential at 30oC. The different temperature minimua between
    ! wood and root growth is due to observed root growth when ever the soil
    ! is not frozen. We also assume that root growth is less sensitive to
    ! available hydraulic pressure, see assumption 3.

    ! Assumption 2:
    ! Actual biological theory suggests that roots support demands for resources
    ! made by the rest of the plant in this current model this is water only.
    ! Therefore there is an implicit assumption that roots should grow so long as
    ! growth is environmentally possible as growth leads to an improvement in
    ! C balance over their life time greater than their construction cost.

    ! Assumption 3:
    ! Determining when root growth should stop is poorly constrained.
    ! Similar to wood growth, here we assume root expansion is also dependent on water availability,
    ! but is less sensitive than wood. Root growth is assumed to stop when deltaWP approaches 0,
    ! determined by marginal return on root growth and temperature limits.

    ! arguments
    integer, intent(in) :: n
    double precision, intent(in) :: lab_to_roots,lab_to_wood &
                                   ,deltaWP,Rtot,current_gpp,Croot,Cwood
    double precision, intent(out) :: root_growth,wood_growth

    ! reset allocation to roots and wood
    root_growth = dble_zero ; wood_growth = dble_zero

    ! Is it currently hydraulically possible for cell expansion (i.e. is soil
    ! water potential more negative than min leaf water potential).
    if ( avail_labile > dble_zero .and. deltaWP < dble_zero ) then

        ! Assume potential root growth is dependent on hydraulic and temperature conditions.
        ! Actual allocation is only allowed if the marginal return on GPP,
        ! averaged across the life span of the root is greater than the rNPP and Rg_root.

        ! Temperature limited turnover rate of labile -> roots
        root_growth = lab_to_roots*Croot_labile_release_coef(n)

        ! calculate hydraulic limits on wood growth.
        ! NOTE: PARAMETERS NEED TO BE CALIBRATRED
        Cwood_hydraulic_limit = (dble_one+exp(Cwood_hydraulic_gradient*(deltaWP-Cwood_hydraulic_half_saturation)))**(-dble_one)
        ! determine wood growth based on temperature and hydraulic limits
        wood_growth = lab_to_wood*Cwood_labile_release_coef(n)*Cwood_hydraulic_limit

        ! cost of wood construction and maintenance not accounted for here due
        ! to no benefit being determined

    endif ! grow root and wood?

    ! track labile reserves to ensure that fractional losses are applied
    ! sequencially in assumed order of importance (leaf->root->wood)

    ! root production (gC.m-2.day-1)
    root_growth = avail_labile*(dble_one-(dble_one-root_growth)**days_per_step)*days_per_step_1
    root_growth = min(avail_labile*days_per_step_1,root_growth)
    avail_labile = avail_labile - (root_growth*days_per_step)
    ! wood production (gC.m-2.day-1)
    wood_growth = avail_labile*(dble_one-(dble_one-wood_growth)**days_per_step)*days_per_step_1
    wood_growth = min(avail_labile*days_per_step_1,wood_growth)
    avail_labile = avail_labile - (wood_growth*days_per_step)

    return

  end subroutine calculate_wood_root_growth
  !
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  ! Functions other than the primary ACM and ACM ET are stored
  ! below this line.
  !------------------------------------------------------------------
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
  double precision function root_resistance(root_mass,thickness)

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
  !-----------------------------------------------------------------
  !
  double precision function soil_resistance(root_length,thickness,soilC)

    !
    ! Calculates the soil hydraulic resistance (MPa m2 s mmol-1) for a given
    ! soil-root zone
    !

    implicit none

    ! arguments
    double precision :: root_length, & ! root length in soil layer (m)
                          thickness, & ! thickness of soil layer (m)
                              soilC    ! soil conductivity m2.s-1.MPa-1

    ! local variables
    double precision :: rs, rs2

    ! calculate
    rs  = (root_length*pi)**(-0.5d0)
    rs2 = log( rs * root_radius_1 ) / (two_pi*root_length*thickness*soilC)
    ! soil water resistance
    soil_resistance = rs2*1d-9*mol_to_g_water

    ! return
    return

  end function soil_resistance
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

   implicit none

   ! arguments
   integer, intent(in) :: root_layer
   double precision, intent(in) :: root_length, &
                                     root_mass, &
                                        demand, &
                                 root_reach_in, &
                      transpiration_resistance

   ! local arguments
   double precision :: soilR1, &
                       soilR2

   ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
   soilR1 = soil_resistance(root_length,root_reach_in,soil_conductivity(root_layer)*head_1)
   soilR2 = root_resistance(root_mass,root_reach_in)
   plant_soil_flow = demand/(transpiration_resistance + soilR1 + soilR2)

   ! return
   return

  end function plant_soil_flow
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    implicit none

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision ::soil_wp

    ! calculate the soil water potential (MPa)..
    ! note that some modifications to scaling values have been made compared to
    ! SPA src to reduce computational cost
!    soil_wp = -0.001 * potA( water_retention_pass ) * xin**potB( water_retention_pass )
!    water_retention_saxton_eqns = -1000.0 * soil_wp + 10.0    ! 10 kPa represents air-entry swp
    soil_wp = potA( water_retention_pass ) * xin**potB( water_retention_pass )
    water_retention_saxton_eqns = -1d0 * soil_wp + 10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
  !
  !--------------------------------------------------------------------------
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
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  ! Generic mathematical functions such as bisection and intergrator proceedures
  ! are stored below here
  !------------------------------------------------------------------
  !
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
    double precision,intent(in)             :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision ,intent(in) :: xval
       end function func
    end interface

    ! local variables..
    integer            :: iter
    integer,parameter  :: ITMAX = 30
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision,parameter :: EPS = 3e-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0. ) then
      zbrent = a
      return
    elseif ( fb .eq. 0. ) then
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
       tol1 = 2. * EPS * abs(b) + 0.5 * tol
       xm   = .5 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0. ) ) then
          zbrent = b
          return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2. * xm * s
             q = 1. - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2. * xm * q * ( q - r ) - ( b - a ) * ( r - 1. ) )
             q = ( q - 1. ) * ( r - 1. ) * ( s - 1. )
          end if
          if ( p .gt. 0. ) q = -q
          p = abs( p )
          if ( (2.*p) .lt. min( 3.*xm*q-abs(tol1*q) , abs(e*q) ) ) then
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
!  subroutine ode_int( called_from , ystart , nvar , x1 , x2 , eps , h1 , hmin , nok , nbad , derivs )
!
!    ! This is an integrator for ordinary differential equations. We use !
!    ! the RUNGE_KUTTA to track the dynamic behaviour of various model   !
!    ! state variables.. uch a leaf water potential, and water stored on !
!    ! leaf surfaces. The RUNGE_KUTTA finds the time step that ensures   !
!    ! dynamics are smooth and that the relevant feedbacks are properly  !
!    ! incorporated. For a full description see Press et al. (1986).     !
!
!    implicit none
!
!    ! arguments..
!    character(len=*),intent(in)    :: called_from    ! name of procedure calling (used to pass through for errors)
!    integer,intent(in)             :: nvar
!    double precision,intent(in)    :: h1, hmin, x1, x2
!    double precision,intent(inout) :: ystart(nvar), eps
!    integer,intent(out)            :: nbad, nok
!
!    ! Interfaces are the correct way to pass procedures as arguments.
!    ! (note 'derivs' is actually one of
!    ! canopy_water_store,soil_water_store,lwp_diff_eqn)
!!    external :: derivs
!    interface
!       subroutine derivs(time,y,dydt, max_iter)
!         integer, intent(in) :: max_iter
!         double precision,intent(in)    :: time
!         double precision,intent(in)    :: y(max_iter)
!         double precision,intent(out)   :: dydt(max_iter)
!       end subroutine derivs
!    end interface
!
!    ! local variables..
!    integer :: i, kount, nstp
!    double precision :: h, hdid, hnext, x, xsav ,         &
!                        dydx(max_nos_iterations), y(max_nos_iterations), &
!                        yscal(max_nos_iterations), xp(kmaxx), yp(max_nos_iterations,kmaxx)
!
!    ! calculations..
!    x = x1
!    h = sign( h1 , x2-x1 )
!    nok = 0
!    nbad = 0
!    kount = 0
!    do i = 1 , nvar
!       y(i) = ystart(i)
!    enddo
!    if ( kmax .gt. 0 ) xsav = x - 2. * dxsav
!    do nstp = 1 , MAXSTP
!       call derivs( x , y , dydx , max_nos_iterations)
!       do i = 1 , nvar
!          yscal(i) = abs(y(i))+abs(h*dydx(i))+TINY
!       enddo
!       if ( kmax .gt. 0 ) then
!          if ( abs( x - xsav ) .gt. abs( dxsav ) ) then
!             if ( kount .lt. kmax - 1 ) then
!                kount = kount + 1
!                xp(kount) = x
!                do i = 1 , nvar
!                   yp(i,kount) = y(i)
!                enddo
!                xsav = x
!             end if
!          end if
!       end if
!       if ( (x+h-x2) * (x+h-x1) .gt. 0. ) h = x2 - x
!
!       call runge_kutta_q_step( trim(called_from)//":ode_int" , y , dydx , nvar , x , h , eps , yscal , hdid , hnext , derivs )
!       if (hdid.eq.h) then
!          nok = nok+1
!       else
!          nbad = nbad+1
!       end if
!       if ( ( x-x2 ) * ( x2-x1 ) .ge. 0. ) then
!          do i = 1 , nvar
!            ystart(i) = y(i)
!          enddo
!          if ( kmax .ne. 0 ) then
!            kount = kount + 1
!            xp(kount) = x
!            do i = 1 , nvar
!              yp(i,kount) = y(i)
!            enddo
!          end if
!          return
!       end if
!       if ( abs(hnext) .lt. hmin ) then
!         print*,"stepsize smaller than permitted minimum in ode_int",new_line('x'),&
!                          "ode_int was called by: ",trim(called_from)
!       end if
!       h  =  hnext
!    enddo
!
!  end subroutine ode_int
  !
  !------------------------------------------------------------------
  !
!  subroutine runge_kutta_check( y , dydx , n , x , h , yout , yerr , derivs )
!
!    ! > subroutine summary? < !
!
!    implicit none
!
!    ! arguments..
!    integer,intent(in) :: n
!    double precision,intent(in)    :: h,x,dydx(n),y(n)
!    double precision,intent(out)   :: yerr(n),yout(n)
!
!    ! Interfaces are the correct way to pass procedures as arguments.
!    ! (note 'derivs' is actually one of
!    ! canopy_water_store,soil_water_store,lwp_diff_eqn)
!    external :: derivs
!    interface
!       subroutine derivs( time , y , dydt , max_iter)
!         integer, intent(in) :: max_iter
!         double precision,intent(in)   :: y(max_iter)
!         double precision,intent(in)   :: time
!         double precision,intent(out)  :: dydt(max_iter)
!       end subroutine derivs
!    end interface
!
!    ! local variables..
!    integer        :: i
!    double precision  :: ak2(max_nos_iterations),ak3(max_nos_iterations) &
!                        ,ak4(max_nos_iterations),ak5(max_nos_iterations) &
!                        ,ak6(max_nos_iterations),ytemp(max_nos_iterations)
!    double precision,parameter :: &
!         A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,B32=9./40., &
!         B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,B53=-70./27.,B54=35./27., &
!         B61=1631./55296.,B62=175./512.,B63=575./13824.,B64=44275./110592., &
!         B65=253./4096.,C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771., &
!         DC1=C1-2825./27648.,DC3=C3-18575./48384.,DC4=C4-13525./55296., &
!         DC5=-277./14336.,DC6=C6-.25
!
!    ! calculations...
!    do i = 1 , n
!       ytemp(i) = y(i) + B21 * h * dydx(i)
!    enddo
!    call derivs( x + A2 * h , ytemp , ak2 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B31 * dydx(i) + B32 * ak2(i) )
!    enddo
!    call derivs( x + A3 * h , ytemp , ak3 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i) )
!    enddo
!    call derivs( x + A4 * h , ytemp , ak4 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + B54 * ak4(i) )
!    enddo
!    call derivs( x + A5 * h , ytemp , ak5 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + B64 * ak4(i) + B65 * ak5(i) )
!    enddo
!    call derivs( x + A6 * h , ytemp , ak6 , max_nos_iterations )
!    do i = 1 , n
!       yout(i) = y(i) + h * ( C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + C6 * ak6(i) )
!       yerr(i) = h * ( DC1 * dydx(i) + DC3 * ak3(i) + DC4 * ak4(i) + DC5 * ak5(i) + DC6 * ak6(i) )
!    enddo
!
!  end subroutine runge_kutta_check
  !
  !----------------------------------------------------------------------
  !
!  subroutine runge_kutta_q_step( called_from , y , dydx , n , x , htry , eps , yscal , hdid , hnext , derivs )
!
!    ! > subroutine summary? < !
!
!    implicit none
!
!    ! arguments..
!    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
!    integer,intent(in)          :: n
!    double precision,intent(in)             :: eps, htry, dydx(n), yscal(n)
!    double precision,intent(inout)          :: x
!    double precision,intent(out)            :: hdid, hnext, y(n)
!
!    ! Interfaces are the correct way to pass procedures as arguments.
!    ! canopy_water_store,soil_water_store,lwp_diff_eqn)
!    interface
!       subroutine derivs( time , y , dydt , max_iter)
!         integer, intent(in) :: max_iter
!         double precision,intent(in)  :: y(max_iter)
!         double precision,intent(in)  :: time
!         double precision,intent(out) :: dydt(max_iter)
!       end subroutine derivs
!    end interface
!
!    ! local variables..
!    integer        :: i
!    double precision           :: errmax, h, htemp, xnew, yerr(max_nos_iterations), ytemp(max_nos_iterations)
!    double precision,parameter :: ERRCON = 1.89e-4, & !
!                                  PGROW  = -0.2,    & !
!                                  PSHRNK = -0.25,   & !
!                                  SAFETY = 0.9        !
!
!    ! calculations...
!    h = htry
!1   call runge_kutta_check( y , dydx , n , x , h , ytemp , yerr , derivs )  ! !YES THERE SHOULD BE A ONE AT THE BEGINNING OF THIS LINE
!    errmax = 0.
!    do  i = 1 , n
!       errmax = max( errmax , abs( yerr(i) / yscal(i) ) )
!    enddo
!    errmax = errmax / eps
!    if ( errmax .gt. 1. ) then
!       htemp = SAFETY * h * ( errmax**PSHRNK )
!       h = sign( max( abs(htemp) , 0.1*abs(h) ) , h )
!       xnew = x + h
!       if ( xnew .eq. x ) then
!         print*,"stepsize underflow in runge_kutta_q_step",new_line('x'),&
!                          "runge_kutta_q_step called from: ",trim(called_from)
!       end if
!       goto 1
!    else
!       if ( errmax .gt. ERRCON ) then
!          hnext = SAFETY * h * ( errmax**PGROW )
!       else
!          hnext = 5. * h
!       end if
!       hdid = h
!       x = x + h
!       do i = 1 , n
!          y(i) = ytemp(i)
!       enddo
!
!    end if
!
!  end subroutine runge_kutta_q_step
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
