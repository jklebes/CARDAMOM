
module CARBON_MODEL_MOD

  implicit none

  ! make all private
  private

  ! explicit publics
  public :: CARBON_MODEL                  &
           ,vsmall                        &
           ,arrhenious                    &
           ,acm_gpp_stage_1               &
           ,acm_gpp_stage_2               &
           ,calculate_radiation_balance   &
           ,calculate_stomatal_conductance&
           ,meteorological_constants      &
           ,calculate_daylength           &
           ,opt_max_scaling               &
           ,freeze                        &
           ,co2comp_saturation            &
           ,co2comp_half_sat_conc         &
           ,kc_saturation                 &
           ,kc_half_sat_conc              &
           ,calculate_Rtot                &
           ,calculate_aerodynamic_conductance &
           ,linear_model_gradient         &
           ,seconds_per_day               &
           ,seconds_per_step              &
           ,fine_root_biomass             &
           ,root_biomass                  &
           ,root_reach                    &
           ,min_root                      &
           ,max_depth                     &
           ,root_k                        &
           ,top_soil_depth                &
           ,mid_soil_depth                &
           ,previous_depth                &
           ,nos_root_layers               &
           ,deltat_1                      &
           ,water_flux                    &
           ,layer_thickness               &
           ,min_layer                     &
           ,soil_frac_clay                &
           ,soil_frac_sand                &
           ,nos_soil_layers               &
           ,meant                         &
           ,stomatal_conductance          &
           ,iWUE                          &
           ,avN                           &
           ,NUE                           &
           ,ceff                          &
           ,pn_max_temp                   &
           ,pn_opt_temp                   &
           ,pn_kurtosis                   &
           ,e0                            &
           ,co2_half_sat                  &
           ,co2_comp_point                &
           ,minlwp                        &
           ,mint                          &
           ,maxt                          &
           ,leafT                         &
           ,swrad                         &
           ,co2                           &
           ,doy                           &
           ,wind_spd                      &
           ,vpd_kPa                       &
           ,lai                           &
           ,dayl_seconds                  &
           ,dayl_seconds_1                &
           ,dayl_hours                    &
           ,dayl_hours_fraction           &
           ,disturbance_residue_to_litter &
           ,disturbance_residue_to_cwd    &
           ,disturbance_residue_to_som    &
           ,disturbance_loss_from_litter  &
           ,disturbance_loss_from_cwd     &
           ,disturbance_loss_from_som     &
           ,Rg_from_labile                &
           ,Rm_from_labile                &
           ,Resp_leaf, Resp_wood_root     &
           ,Rm_leaf, Rm_wood_root         &
           ,Rg_leaf, Rg_wood_root         &
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
  logical :: do_iWUE = .true. ! Use iWUE or WUE for stomatal optimisation
  double precision, parameter :: vsmall = tiny(0d0)*1d3 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  integer, parameter :: nos_root_layers = 3, nos_soil_layers = nos_root_layers + 1
  double precision, parameter :: pi = 3.1415927d0,  &
                               pi_1 = 0.3183099d0,  & ! pi**(-1d0)
                             two_pi = 6.283185d0,   & ! pi*2d0
                         deg_to_rad = 0.01745329d0, & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0 * deg_to_rad )
                              boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                         emissivity = 0.96d0,       &
                        emiss_boltz = 5.443584d-08, & ! emissivity * boltz
                    sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                             freeze = 273.15d0,     &
                  gs_H2Ommol_CO2mol = 0.001646259d0,& ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
              gs_H2Ommol_CO2mol_day = 142.2368d0,   & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
                         gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
            partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                     mol_to_g_water = 18d0,         & ! molecular mass of water
                   mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                         umol_to_gC = 1.2d-5,       & ! conversion of umolC -> gC
                         gC_to_umol = 83333.33d0,   & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                               Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                          vonkarman = 0.41d0,       & ! von Karman's constant
                        vonkarman_1 = 2.439024d0,   & ! 1 / von Karman's constant
                              cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

  ! photosynthesis / respiration parameters
  double precision, parameter :: &
                      kc_saturation = 310d0,        & ! CO2 half saturation, saturation value
                   kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, half sat
                 co2comp_saturation = 36.5d0,       & ! CO2 compensation point, saturation
              co2comp_half_sat_conc = 9.46d0,       & ! CO2 comp point, half sat
                                                      ! Each of these are temperature sensitivty
                        Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                      ! lost as growth respiration
                                                      ! (i.e. 0.28 .eq. xNPP)
                    one_Rg_fraction = 1d0 - Rg_fraction

  ! hydraulic parameters
  double precision, parameter :: &
                         tortuosity = 2.5d0,        & ! tortuosity
                             gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                        root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                        root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                      ! Williams et al 1996 = 0.0001
                      root_radius_1 = root_radius**(-1d0), &
                root_cross_sec_area = pi * root_radius**2, & ! root cross sectional area (m2)
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
                           min_wind = 0.2d0,        & ! minimum wind speed at canopy top
                       min_drythick = 0.01d0,       & ! minimum dry thickness depth (m)
                          min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                        soil_roughl = 0.05d0,       & ! soil roughness length (m)
                     top_soil_depth = 0.1d0,        & ! thickness of the top soil layer (m)
                     mid_soil_depth = 0.2d0,        & ! thickness of the second soil layer (m)
                           min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                            min_lai = 1.0d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                    min_throughfall = 0.2d0,        & ! minimum fraction of precipitation which
                                                      ! is through fall
                        min_storage = 0.2d0           ! minimum canopy water (surface) storage (mm)

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                    seconds_per_day = 86400d0,        & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05      ! Inverse of seconds per day

  ! ACM-GPP-ET parameters
  double precision, parameter :: &
                        pn_max_temp = 5.978665d+01, & ! Maximum temperature for photosynthesis (oC)
                        pn_opt_temp = 3.150166d+01, & ! Optimum temperature fpr photosynthesis (oC)
                        pn_kurtosis = 1.730812d-01, & ! Kurtosis of photosynthesis temperature response
                                 e0 = 5.897208d+00, & ! Quantum yield gC/MJ/m2/day PAR
         lai_half_lwrad_transmitted = 1.639527d-01, & ! LAI at which LW transmittance to soil at 50 %
           lai_half_lwrad_reflected = 6.746538d-03, & ! LAI at which LW reflectance to sky at 50 %
             max_lai_nir_reflection = 1.457203d-01, & ! Coefficient relating NIR reflectance and LAI
            lai_half_nir_reflection = 1.720843d-02, & ! LAI at which canopy NIR reflection = 50 %
                     minlwp_default =-1.936111d+00, & ! minimum leaf water potential (MPa)
             max_lai_par_reflection = 1.774887d-01, & ! Max fraction of PAR reflected by canopy
            lai_half_par_reflection = 1.094594d-02, & ! LAI at which canopy PAR reflected = 50 %
                      max_lw_escape = 5.813763d-01, & ! Max LW which is released from canopy that escapes in one direction
                               iWUE = 3.312080d-06, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
              soil_swrad_absorption = 9.641337d-01, & ! Fraction of SW rad absorbed by soil
            max_lai_par_transmitted =-1.096639d-02, & ! Max fraction reduction in PAR transmittance by canopy
            max_lai_nir_transmitted =-1.255590d-02, & ! Max fraction reduction in NIR transmittance by canopy
              max_lai_lwrad_release = 9.989129d-01, & ! 1-Max fraction of LW emitted from canopy to be released
             lai_half_lwrad_release = 2.009765d+00, & ! LAI at which LW emitted from canopy to be released at 50 %
               soil_iso_to_net_coef =-9.475060d-01, & ! Coefficient relating soil isothermal net radiation to net.
              soil_iso_to_net_const =-2.443723d+00, & ! Constant relating soil isothermal net radiation to net
                max_par_transmitted = 1.392685d-01, & ! Max fraction of canopy incident PAR transmitted to soil
                max_nir_transmitted = 2.352473d-01, & ! Max fraction of canopy incident NIR transmitted to soil
                  max_par_reflected = 4.700132d-01, & ! Max fraction of canopy incident PAR reflected to sky
                  max_nir_reflected = 5.338683d-01    ! Max fraction of canopy incident NIR reflected to sky

  double precision :: minlwp = minlwp_default

  !!!!!!!!!
  ! Module level variables
  !!!!!!!!!

  ! management and gsi related values
  integer :: gsi_lag_remembered
  ! local variables for GSI phenology model
  double precision :: Tfac,Photofac,VPDfac & ! oC, seconds, Pa
                     ,Rm_leaf_per_gC, Rm_leaf_baseline     & 
                     ,mean_Q10_adjustment  &
                     ,leaf_life, SLA, NCE_smoother &
                     ,avail_labile                   &
                     ,Cwood_labile_release_gradient  &
                     ,Cwood_labile_half_saturation   &
                     ,Croot_labile_release_gradient  &
                     ,Croot_labile_half_saturation   &
                     ,Cwood_hydraulic_gradient       &
                     ,Cwood_hydraulic_half_saturation&
                     ,Cwood_hydraulic_limit          &
                     ,tmp,gradient,fol_turn_crit

  double precision, allocatable, dimension(:) :: &
                                              tmp_x,gsi_history, &
                                                 Q10_adjustment, &
                                                 Rg_from_labile, &
                                                 Rm_from_labile, & 
                                      Resp_leaf, Resp_wood_root, & ! Total respiration (gC/m2/day)
                                          Rm_leaf, Rm_wood_root, & ! Maintenance respiration (gC/m2/day)
                                          Rg_leaf, Rg_wood_root, &
                                              itemp,ivpd,iphoto, &
                                  disturbance_residue_to_litter, &
                                     disturbance_residue_to_som, &
                                     disturbance_residue_to_cwd, &
                                   disturbance_loss_from_litter, &
                                      disturbance_loss_from_cwd, &
                                      disturbance_loss_from_som 

  ! hydraulic model variables
  integer :: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil
  double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                                       water_flux    ! potential transpiration flux (mmol.m-2.s-1)
  double precision, dimension(nos_soil_layers+1) :: layer_thickness ! thickness of soil layers (m)

  double precision :: root_reach, root_biomass, fine_root_biomass, & ! root depth, coarse+fine, and fine root biomass
                                     max_depth, & ! maximum possible root depth (m)
                                        root_k, & ! biomass to reach half max_depth
                      new_depth,previous_depth, & ! depth of bottom of soil profile
                                   canopy_wind, & ! wind speed (m.s-1) at canopy top
                                         ustar, & ! friction velocity (m.s-1)
                                      ustar_Uh, &
                                air_density_kg, & ! air density kg/m3
                                ET_demand_coef, & ! air_density_kg * vpd_kPa * cpair
                                        roughl, & ! roughness length (m)
                                  displacement, & ! zero plane displacement (m)
                                    max_supply, & ! maximum water supply (mmolH2O/m2/day)
                                         meant, & ! mean air temperature (oC)
                                         leafT, & ! canopy temperature (oC)
                              mean_annual_temp, &
                            canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
                              soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
                                soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
                                 sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                          stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.s-1)
                       aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
                             convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                            convert_ms1_mmol_1, & ! Conversion ratio for m/s -> mmol/m2/s
                           air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                                        lambda, & ! latent heat of vapourisation (J.kg-1)
                                         psych, & ! psychrometric constant (kPa K-1)
                                         slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
                        water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
                           kinematic_viscosity    ! kinematic viscosity (m2.s-1)

  ! Module level variables for ACM_GPP_ET parameters
  double precision ::   delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                            ceff, & ! canopy efficency, ceff = avN*NUE
                             avN, & ! average foliar N (gN/m2)
                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                    ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
metabolic_limited_photosynthesis, & ! temperature, leaf area and foliar N limiterd photosynthesis (gC/m2/day)
    light_limited_photosynthesis, & ! light limited photosynthesis (gC/m2/day)
                          gb_mol, & ! Canopy boundary layer conductance (molCO2/m2/day)
                    co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
                  co2_comp_point    ! CO2 at which photosynthesis > 0 (ppm)

  ! Module level variables for step specific met drivers
  double precision :: mint, & ! minimum temperature (oC)
                      maxt, & ! maximum temperature (oC)
                     swrad, & ! incoming short wave radiation (MJ/m2/day)
                       co2, & ! CO2 (ppm)
                       doy, & ! Day of year
                  wind_spd, & ! wind speed (m/s)
                   vpd_kPa, & ! Vapour pressure deficit (kPa)
                     lai_1, & ! inverse of LAI
                       lai    ! leaf area index (m2/m2)

  ! Module level varoables for step specific timing information
  integer :: steps_per_year
  double precision :: cos_solar_zenith_angle, &
                            seconds_per_step, & !
                          mean_days_per_step, &
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                         dayl_hours_fraction, &
                                  dayl_hours    ! day length in hours

  double precision, dimension(:), allocatable :: deltat_1, & ! inverse of decimal days
                                               meant_time, &
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1, &
                                Cwood_labile_release_coef, & ! time series of labile release to wood
                                Croot_labile_release_coef    ! time series of labile release to root

  save

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out,NEE_out,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP_out)

    !
    ! The Data Assimilation Linked Ecosystem Carbon (DALEC) model.
    !
    ! The Aggregated Canopy Model for Gross Primary Productivity and Evapotranspiration (ACM-GPP-ET)
    ! simulates coupled photosynthesis-transpiration (via stomata), soil and intercepted canopy evaporation and
    ! soil water balance (4 layers). Note that in this version the soil is
    ! assumed to be saturated. As a result much of the water cycle code has been
    ! removed for efficiency.
    !
    ! Carbon allocation based on fixed fraction and turnover follows first order kinetics with the following exceptions.
    ! 1) Foliar allocation is based on marginal return calculation estimating the net canopy export of carbon
    ! 2) Foliar turnover is based on the Growing Season Index framework.
    ! 3) Turnover of litter and soil includes an exponential temperature dependency.
    !
    ! This version was coded by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! Version 1: 29/01/2020

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
                        ,Rtot   ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    integer :: f,nxp,n,test,m

    ! local fire related variables
    double precision :: burnt_area           &
                       ,CFF(7) = 0d0   & ! combusted and non-combustion fluxes
                       ,NCFF(7) = 0d0  & ! with residue and non-residue seperates
                       ,combust_eff(5) & ! combustion efficiency
                       ,rfac             ! resilience factor

    ! local deforestation related variables
    double precision, dimension(5) :: post_harvest_burn   & ! how much burning to occur after
                                     ,foliage_frac_res    &
                                     ,roots_frac_res      &
                                     ,rootcr_frac_res     &
                                     ,stem_frac_res       &
                                     ,Crootcr_part        &
                                     ,soil_loss_frac

    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   &
                       ,C_total,labile_frac_res      &
                       ,Cstem,Crootcr,stem_residue   &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    integer :: reforest_day, harvest_management,restocking_lag

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
    ! 4 = respiration het cwd
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
    ! 19 = NOT IN USE
    ! 20 = CWD turnover to litter
    ! 21 = C extracted as harvest
    ! 22 = labile loss due to disturbance
    ! 23 = foliage loss due to disturbance
    ! 24 = root loss due to disturbance
    ! 25 = wood loss due to disturbance
    ! 26 = NOT IN USE
    ! 27 = NOT IN USE
    ! 28 = NOT IN USE

    ! PARAMETERS
    ! 35 process parameters; 7 C pool initial conditions

    ! p(1) Decomposition efficiency (fraction to som)
    ! p(2) Fraction of GPP respired as maintenance for wood and root
    ! p(3) Baseline foliar turnover (frac/day)
    ! p(4) Fraction of NPP allocated to roots 
    ! p(5) Max leaf turnover (GSI)
    ! p(6) Turnover rate of wood 
    ! p(7) Turnover rate of roots 
    ! p(8) Litter turnover rate 
    ! p(9) SOM mineralisation rate 
    ! p(10) Parameter in exponential term of temperature 
    ! p(11) Mean foliar nitrogen content (gN/m2)
    ! p(12) = Max labile turnover(NCE,GSI)
    ! p(13) = Fraction allocated to Clab 
    ! p(14) = Min temp threshold (GSI)
    ! p(15) = Max temp threshold (GSI)
    ! p(16) = Min photoperiod threshold (GSI)
    ! p(17) = LMA
    ! p(24) = Max photoperiod threshold (GSI)
    ! p(25) = Min VPD threshold (GSI)
    ! p(26) = Max VPD threshold (GSI)
    ! p(27) = Net Canopy Export return on new Cfol investment (gCperNCE per gCnewfol)
    ! p(28) = Initial GSI
    ! p(29) = fraction of Cwood which is Ccoarseroot
    ! p(34) = GSI sensitivity for leaf scenescence
    ! p(35) = Reich Maintenance respiration N leaf exponent
    ! p(36) = Reich Maintenance respiration N leaf baseline
    ! p(37) = Initial CWD pool
    ! p(38) = CWD turnover fraction
    ! p(39) = Fine root (gbiomass.m-2) needed to reach 50% of max depth
    ! p(40) = Maximum rooting depth (m)
    ! p(41) = NOT IN USE
    ! p(42) = Nitrogen use efficiency (gC/gN/day)
    ! p(43) = Initial leaf lifespan (days)

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

    ! infinity check requirement
    infi = 0d0
    ! reset basic input / output variables
    FLUXES = 0d0 ; POOLS = 0d0

    ! load ACM-GPP-ET parameters
    NUE = pars(42)       ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                         ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
                         ! Other initial values for ACM_GPP_ET
    avN = 10d0**pars(11) ! foliar N gN/m2
    ceff = avN*NUE       ! canopy efficiency, used to avoid what in most cases is a reductance multiplication
                         ! NOTE: must be updated any time NUE or avN changes
    Rtot = 1d0

    ! Estimate time invarient N response for maintenance respiration
    ! Include scalings from nmolC/g/s -> gC/m2/day
    ! Note that the mean temperature Q10 will be estimates in loop below where
    ! meant_time calculated
    Rm_leaf_baseline = Rm_reich_N(pars(17)/avN,pars(35),pars(36)) * umol_to_gC * seconds_per_day * 2d-3
    ! set initial leaf lifespan
    leaf_life = pars(43)

    ! plus ones being calibrated
    root_k = pars(39) ; max_depth = pars(40)

    ! if either of our disturbance drivers indicate disturbance will occur then
    ! set up these components
    if (maxval(met(8,:)) > 0d0 .or. maxval(met(9,:)) > 0d0) then

        ! initial values for deforestation variables
        labile_loss = 0d0    ; foliar_loss = 0d0
        roots_loss = 0d0     ; wood_loss = 0d0
        labile_residue = 0d0 ; foliar_residue = 0d0
        roots_residue = 0d0  ; wood_residue = 0d0
        stem_residue = 0d0
        reforest_day = 0
        soil_loss_with_roots = 0d0
        coarse_root_residue = 0d0
        post_harvest_burn = 0d0

        ! now load the hardcoded forest management parameters into their locations

        ! Parameter values for deforestation variables
        ! scenario 1
        ! harvest residue (fraction); 1 = all remains, 0 = all removed
        foliage_frac_res(1) = 1d0
        roots_frac_res(1)   = 1d0
        rootcr_frac_res(1) = 1d0
        stem_frac_res(1)   = 0.20d0 !
        ! wood partitioning (fraction)
        Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
        ! was the forest burned after deforestation
        post_harvest_burn(1) = 1d0

        !## scen 2
        ! harvest residue (fraction); 1 = all remains, 0 = all removed
        foliage_frac_res(2) = 1d0
        roots_frac_res(2)   = 1d0
        rootcr_frac_res(2) = 1d0

        stem_frac_res(2)   = 0.20d0 !
        ! wood partitioning (fraction)
        Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
        ! was the forest burned after deforestation
        post_harvest_burn(2) = 0d0

        !## scen 3
        ! harvest residue (fraction); 1 = all remains, 0 = all removed
        foliage_frac_res(3) = 0.5d0
        roots_frac_res(3)   = 1d0
        rootcr_frac_res(3) = 1d0
        stem_frac_res(3)   = 0d0 !
        ! wood partitioning (fraction)
        Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
        ! was the forest burned after deforestation
        post_harvest_burn(3) = 0d0

        !## scen 4
        ! harvest residue (fraction); 1 = all remains, 0 = all removed
        foliage_frac_res(4) = 0.5d0
        roots_frac_res(4)   = 1d0
        rootcr_frac_res(4) = 0d0
        stem_frac_res(4)   = 0d0
        ! wood partitioning (fraction)
        Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
        ! was the forest burned after deforestation
        post_harvest_burn(4) = 0d0

        !## scen 5 (grassland grazing / cutting)
        ! harvest residue (fraction); 1 = all remains, 0 = all removed
        foliage_frac_res(5) = 0.1d0
        roots_frac_res(5)   = 0d0
        rootcr_frac_res(5)  = 0d0
        stem_frac_res(5)    = 0.12d0
        ! wood partitioning (fraction)
        Crootcr_part(5) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(5) = 0d0 ! actually between 1-3 %
        ! was the forest burned after deforestation
        post_harvest_burn(5) = 0d0

        ! for the moment override all paritioning parameters with those coming from
        ! CARDAMOM
        Crootcr_part = pars(29)

        ! declare fire constants (labile, foliar, roots, wood, litter)
        combust_eff(1) = 0.1d0 ; combust_eff(2) = 0.9d0
        combust_eff(3) = 0.1d0 ; combust_eff(4) = 0.5d0
        combust_eff(5) = 0.3d0 ; rfac = 0.5d0

    end if ! disturbance ?

    ! assigning initial conditions for the current iteration
    POOLS(1,1) = pars(18)
    POOLS(1,2) = pars(19)
    POOLS(1,3) = pars(20)
    POOLS(1,4) = pars(21)
    POOLS(1,5) = pars(22)
    POOLS(1,6) = pars(23)
    POOLS(1,7) = pars(37)
    ! POOL(1,8) assigned later

    if (.not.allocated(deltat_1)) then
        ! allocate variables dimension which are fixed per site only the once
        allocate(disturbance_residue_to_litter(nodays),disturbance_residue_to_cwd(nodays), &
                 disturbance_residue_to_som(nodays),disturbance_loss_from_litter(nodays), & 
                 disturbance_loss_from_cwd(nodays),disturbance_loss_from_som(nodays),&
                 Cwood_labile_release_coef(nodays),Croot_labile_release_coef(nodays), &
                 deltat_1(nodays),daylength_hours(nodays),daylength_seconds(nodays), & 
                 daylength_seconds_1(nodays),meant_time(nodays), &
                 Rg_from_labile(nodays),Rm_from_labile(nodays),Resp_leaf(nodays), &
                 Resp_wood_root(nodays),Rm_leaf(nodays),Rm_wood_root(nodays),Q10_adjustment(nodays), &
                 Rg_leaf(nodays),Rg_wood_root(nodays),itemp(nodays),ivpd(nodays),iphoto(nodays))

        !
        ! Timing variables which are needed first
        !

        deltat_1 = deltat**(-1d0)

        !
        ! Iteration independent variables using functions and thus need to be in a loop
        !

        ! first those linked to the time period of the analysis
        do n = 1, nodays
           ! calculate daylength in hours and seconds
           call calculate_daylength((met(6,n)-(deltat(n)*0.5d0)),lat)
           daylength_hours(n) = dayl_hours ; daylength_seconds(n) = dayl_seconds
        end do

        ! calculate inverse for each time step in seconds
        daylength_seconds_1 = daylength_seconds ** (-1d0)
        ! meant time step temperature
        meant_time = (met(2,:)+met(3,:)) * 0.5d0

        ! number of time steps per year
        steps_per_year = nint(dble(nodays)/(sum(deltat)*0.002737851d0))
        ! mean days per step
        mean_days_per_step = sum(deltat) / dble(nodays)

        ! Calculate timing components needed for GSI / NCE gradient calculations
        gsi_lag_remembered = max(2,nint(21d0/mean_days_per_step))
        allocate(tmp_x(gsi_lag_remembered),gsi_history(gsi_lag_remembered))
        do f = 1, gsi_lag_remembered
           tmp_x(f) = dble(f) * mean_days_per_step
        end do

    endif ! deltat_1 allocated

    ! assign our starting value
    gsi_history = pars(28) 

    ! specific leaf area (m2/gC)
    SLA = pars(17)**(-1d0)

    ! load some needed module level values
    lai = POOLS(1,2)*SLA
    mint = met(2,1)  ! minimum temperature (oC)
    maxt = met(3,1)  ! maximum temperature (oC)
    swrad = met(4,1) ! incoming short wave radiation (MJ/m2/day)
    co2 = met(5,1)   ! CO2 (ppm)
    doy = met(6,1)   ! Day of year
    wind_spd = met(15,1) ! wind speed (m/s)
    vpd_kPa = met(16,1)*1d-3 ! vapour pressure deficit (Pa->kPa)
    meant = meant_time(1)
    leafT = maxt     ! initial canopy temperature (oC)
    seconds_per_step = deltat(1) * seconds_per_day

    ! reset disturbance at the beginning of iteration
    disturbance_residue_to_litter = 0d0 ; disturbance_loss_from_litter = 0d0
    disturbance_residue_to_som = 0d0 ; disturbance_loss_from_som = 0d0
    disturbance_residue_to_cwd = 0d0 ; disturbance_loss_from_cwd = 0d0

    ! Initialise root reach based on initial coarse root biomass
    fine_root_biomass = max(min_root,POOLS(1,3)*2d0)
    root_biomass = fine_root_biomass + max(min_root,POOLS(1,4)*pars(29)*2d0)
    ! Needed to initialise soils
    call calculate_Rtot(Rtot)

    ! assign climate sensitivities
    fol_turn_crit = pars(34) 

    ! Calculate GSI ranges
    Tfac_range_1 = (pars(15)-pars(14))**(-1d0)
    Photofac_range_1 = (pars(24)-pars(16))**(-1d0)
    VPDfac_range_1 = (pars(26)-pars(25))**(-1d0)

    !!!!!!!!!!!!
    ! assign climate sensitivities
    !!!!!!!!!!!!

    FLUXES(:,2) = exp(pars(10)*meant_time)

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
       leafT = maxt     ! initial canopy temperature (oC)
       swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
       co2 = met(5,n)   ! CO2 (ppm)
       doy = met(6,n)   ! Day of year
       meant = meant_time(n) ! mean air temperature (oC)
       wind_spd = met(15,n) ! wind speed (m/s)
       vpd_kPa = met(16,n)*1d-3 ! vapour pressure deficit (Pa->kPa)

       ! states needed for module variables
       lai_out(n) = POOLS(n,2)*SLA
       lai = lai_out(n) ! leaf area index (m2/m2)

       ! extract timing related values
       dayl_hours = daylength_hours(n)
       dayl_hours_fraction = dayl_hours * 0.04166667d0 ! 1/24 = 0.04166667
       dayl_seconds = daylength_seconds(n) ; dayl_seconds_1 = daylength_seconds_1(n)
       seconds_per_step = seconds_per_day * deltat(n)

       !!!!!!!!!!
       ! Calculate surface exchange coefficients
       !!!!!!!!!!

       ! calculate some temperature dependent meteorologial properties
       call meteorological_constants(maxt,maxt+freeze,vpd_kPa)
       ! pass variables from memory objects
       convert_ms1_mmol_1 = convert_ms1_mol_1 * 1d3
       ! calculate aerodynamic using consistent approach with SPA
       call calculate_aerodynamic_conductance

       !!!!!!!!!!
       ! Determine net shortwave and isothermal longwave energy balance
       !!!!!!!!!!

       call calculate_radiation_balance

       !!!!!!!!!!
       ! Calculate physically total hydraulic resistance
       !!!!!!!!!!

       ! calculate the minimum soil & root hydraulic resistance based on total
       ! fine root mass ! *2*2 => *RS*C->Bio
       fine_root_biomass = max(min_root,POOLS(n,3)*2d0)
       root_biomass = fine_root_biomass + max(min_root,POOLS(n,4)*pars(29)*2d0)
       call calculate_Rtot(Rtot)

       ! calculate radiation absorption and estimate stomatal conductance
       call calculate_stomatal_conductance(abs(minlwp),Rtot)

       ! Note that soil mass balance will be calculated after phenology
       ! adjustments

       ! Reset output variable
       if (stomatal_conductance > vsmall) then
           ! Gross primary productivity (gC/m2/day)
           call acm_gpp_stage_1 ; FLUXES(n,1) = acm_gpp_stage_2(stomatal_conductance)
       else
           ! assume zero fluxes
           FLUXES(n,1) = 0d0 
       endif

       !!!!!!!!!!
       ! GPP allocation
       !!!!!!!!!!

       ! autotrophic respiration (gC.m-2.day-1)
       Rm_wood_root(n) = pars(2)*FLUXES(n,1)
       ! labile production (gC.m-2.day-1)
       FLUXES(n,5) = (FLUXES(n,1)-Rm_wood_root(n))*pars(13)
       ! root production (gC.m-2.day-1)
       FLUXES(n,6) = (FLUXES(n,1)-Rm_wood_root(n)-FLUXES(n,5))*pars(4)
       ! wood production
       FLUXES(n,7) = FLUXES(n,1)-Rm_wood_root(n)-FLUXES(n,5)-FLUXES(n,6)

       !!!!!!!!!!
       ! calculate maintenance respiration demands and mass balance
       !!!!!!!!!!

       ! Scale baseline Rm (gC/gCleaf/day at 20oC) for leaf with current temperature Q10 response
       ! giving maintenance respiration at maximum of current or mean annual temperature per gC leaf
       ! (gC/gCleaf/day).
       Q10_adjustment(n) = Rm_reich_Q10(meant)
       Rm_leaf_per_gC = Rm_leaf_baseline * Q10_adjustment(n)
       ! Estimate time step demand on labile pool to support canopy maintenance
       ! NOTE: current Rm costs must be related to current temperature
       Rm_leaf(n) = min(POOLS(n,1),Q10_adjustment(n) * Rm_leaf_baseline * POOLS(n,2) * deltat(n)) * deltat_1(n)

       !!!!!!!!!!
       ! Calculate canopy phenology
       !!!!!!!!!!

       ! assign labile C available in current time step, less that needed for
       ! maintenance respiration
       avail_labile = POOLS(n,1) - (Rm_leaf(n)*deltat(n))

       ! Determine leaf growth and turnover based on GSI model + some economics
       ! NOTE: that turnovers will be bypassed in favour of mortality turnover
       ! should available labile be exhausted
       call calculate_leaf_dynamics(n,deltat,nodays        &
                                   ,pars(14),pars(16),pars(25)       &
                                   ,Tfac_range_1,Photofac_range_1    &
                                   ,VPDfac_range_1,pars(3),pars(5),pars(12)  &
                                   ,met(10,n),met(11,n),met(12,n),minlwp,Rtot &
                                   ,FLUXES(n,1),POOLS(n,2),pars(27)  &
                                   ,FLUXES(:,18),FLUXES(n,9),FLUXES(n,16))

       ! Total labile release to foliage
       FLUXES(n,8) = avail_labile*(1d0-(1d0-FLUXES(n,16))**deltat(n))*deltat_1(n)

       ! if 12 months has gone by, update the leaf lifespan variable
       if (n > steps_per_year .and. met(6,n) < met(6,n-1)) then
           ! determine the mean life span (days)
           tmp = sum(POOLS((n-steps_per_year):(n-1),2)) &
               / sum(FLUXES((n-steps_per_year):(n-1),10) + FLUXES((n-steps_per_year):(n-1),23))
           ! i.e. we cannot / should not update the leaf lifespan if there has
           ! been no turnover and / or there is no foliar pool.
           ! 2933 = 365.25 * 8 years
           if (tmp > 0d0 .and. tmp < 2933d0) then
               ! We assume that leaf life span is weighted 50:50 between the
               ! previous year and its history
               leaf_life = (tmp + leaf_life) * 0.5d0
           end if
       endif ! n /= 1 and new calendar year

       !
       ! litter creation with time dependancies
       !

       ! total leaf litter production
       FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**deltat(n))*deltat_1(n)
       ! total wood litter production
       ! NOTE: we assume that the coarse root fraction (pars(29)) is passed to the soil directly
       FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-pars(6))**deltat(n))*deltat_1(n)
       ! total root litter production
       FLUXES(n,12) = POOLS(n,3)*(1d0-(1d0-pars(7))**deltat(n))*deltat_1(n)

       !
       ! those with temperature AND time dependancies
       !

       ! turnover of litter
       tmp = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))*deltat_1(n)
       ! respiration heterotrophic litter ; decomposition of litter to som
       FLUXES(n,13) = tmp * (1d0-pars(1)) ; FLUXES(n,15) = tmp * pars(1)

       ! respiration heterotrophic som
       FLUXES(n,14) = POOLS(n,6)*(1d0-(1d0-FLUXES(n,2)*pars(9))**deltat(n))*deltat_1(n)

       ! respiration heterotrophic cwd ; decomposition of CWD to som
       tmp = POOLS(n,7)*(1d0-(1d0-FLUXES(n,2)*pars(38))**deltat(n))*deltat_1(n)
       FLUXES(n,4) = tmp * (1d0-pars(1)) ; FLUXES(n,20) = tmp * pars(1)

       !!!!!!!!!!
       ! calculate growth respiration and adjust allocation to pools assuming
       ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
       !!!!!!!!!!

       ! calculate growth respiration and adjust allocation to pools assuming
       ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
       ! foliage. NOTE: that in the current version only Rg_leaf and Rm_leaf comes from
       ! Clabile
    
       ! Determine growth respiration...
       Rg_leaf(n) = FLUXES(n,8)*Rg_fraction
       Rg_wood_root(n) = (FLUXES(n,7)*Rg_fraction) + (FLUXES(n,6)*Rg_fraction) 
       ! ...and remove from the bulk flux to tissue
       FLUXES(n,8) = FLUXES(n,8) * one_Rg_fraction ! Leaf
       FLUXES(n,6) = FLUXES(n,6) * one_Rg_fraction ! Wood
       FLUXES(n,7) = FLUXES(n,7) * one_Rg_fraction ! Fine root
       ! Combine to give total tissue respiration fluxes
       Resp_leaf(n) = Rg_leaf(n) + Rm_leaf(n)
       Resp_wood_root(n) = Rg_wood_root(n) + Rm_wood_root(n)
       ! Total Autotrophic respiration
       FLUXES(n,3) = Resp_leaf(n) + Resp_wood_root(n)
       ! Total fluxes from labile pool
       Rg_from_labile(n) = Rg_leaf(n) ; Rm_from_labile(n) = Rm_leaf(n)

       !
       ! update pools for next timestep
       !

       ! labile pool
       POOLS(n+1,1) = POOLS(n,1) + ((FLUXES(n,5)-FLUXES(n,8)-Rg_from_labile(n)-Rm_from_labile(n))*deltat(n))
       ! foliar pool
       POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,8)-FLUXES(n,10))*deltat(n)
       ! wood pool
       POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
       ! root pool
       POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
       ! litter pool
       POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
       ! som pool
       POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)+FLUXES(n,20)-FLUXES(n,14))*deltat(n)
!       POOLS(n+1,6) = POOLS(n,6) + ((FLUXES(n,11)*pars(29))+FLUXES(n,15)+FLUXES(n,20)-FLUXES(n,14))*deltat(n)
       ! litter wood pool
       POOLS(n+1,7) = POOLS(n,7) + (FLUXES(n,11)-FLUXES(n,20)-FLUXES(n,4))*deltat(n)
!       POOLS(n+1,7) = POOLS(n,7) +  ((FLUXES(n,11)*(1d0-pars(29)))-FLUXES(n,20)-FLUXES(n,4))*deltat(n)

       !!!!!!!!!!
       ! deal first with deforestation
       !!!!!!!!!!

       if (n == reforest_day) then
           POOLS(n+1,1) = pars(30)
           POOLS(n+1,2) = pars(31)
           POOLS(n+1,3) = pars(32)
           POOLS(n+1,4) = pars(33)
       endif

       ! reset values
       harvest_management = 0 ; burnt_area = 0d0

       if (met(8,n) > 0d0) then

           ! pass harvest management to local integer
           harvest_management = int(met(13,n))

           ! assume that labile is proportionally distributed through the plant
           ! root and wood and therefore so is the residual fraction
           C_total = POOLS(n+1,3) + POOLS(n+1,4)
           ! partition wood into its components
           Crootcr = POOLS(n+1,4)*Crootcr_part(harvest_management)
           Cstem   = POOLS(n+1,4)-Crootcr !(Cbranch + Crootcr)
           ! now calculate the labile fraction of residue
           if (C_total > 0d0) then
               labile_frac_res = ((POOLS(n+1,3)/C_total) * roots_frac_res(harvest_management)  ) &
                               + ((Cstem/C_total)        * stem_frac_res(harvest_management)   ) &
                               + ((Crootcr/C_total)      * rootcr_frac_res(harvest_management) )
           else
               labile_frac_res = 0d0
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
                   roots_loss = 0d0
               endif
               wood_loss   = (Crootcr+Cstem)*met(8,n)
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
!               branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
               stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
               ! Now finally calculate the final wood residue
               wood_residue = stem_residue + coarse_root_residue
               ! Mechanical loss of Csom due to coarse root extraction
               soil_loss_with_roots = Crootcr*met(8,n)*(1d0-rootcr_frac_res(harvest_management)) &
                                    * soil_loss_frac(harvest_management)

               ! Update pools
               POOLS(n+1,1) = POOLS(n+1,1) - labile_loss
               POOLS(n+1,2) = POOLS(n+1,2) - foliar_loss
               POOLS(n+1,3) = POOLS(n+1,3) - roots_loss
               POOLS(n+1,4) = POOLS(n+1,4) - wood_loss
               POOLS(n+1,5) = POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue)
               POOLS(n+1,6) = POOLS(n+1,6) - soil_loss_with_roots
               POOLS(n+1,7) = POOLS(n+1,7) + wood_residue
               ! mass balance check
               where (POOLS(n+1,1:7) < 0d0) POOLS(n+1,1:7) = 0d0

               ! Some variable needed for the EDCs
               ! reallocation fluxes for the residues
               disturbance_residue_to_litter(n) = labile_residue+foliar_residue+roots_residue
               disturbance_loss_from_litter(n)  = 0d0
               disturbance_residue_to_cwd(n)    = wood_residue
               disturbance_loss_from_cwd(n)     = 0d0
               disturbance_residue_to_som(n)    = 0d0
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

       if (met(9,n) > 0d0 .or.(met(8,n) > 0d0 .and. harvest_management > 0)) then

           burnt_area = met(9,n)
           if (met(8,n) > 0d0 .and. burnt_area > 0d0) then
               ! pass harvest management to local integer
               burnt_area = min(1d0,burnt_area + post_harvest_burn(harvest_management))
           else if (met(8,n) > 0d0 .and. burnt_area <= 0d0) then
               burnt_area = post_harvest_burn(harvest_management)
           endif

           if (burnt_area > 0d0) then

               !/*first fluxes*/
               !/*LABILE*/
               CFF(1) = POOLS(n+1,1)*burnt_area*combust_eff(1)
               NCFF(1) = POOLS(n+1,1)*burnt_area*(1d0-combust_eff(1))*(1d0-rfac)
               !/*foliar*/
               CFF(2) = POOLS(n+1,2)*burnt_area*combust_eff(2)
               NCFF(2) = POOLS(n+1,2)*burnt_area*(1d0-combust_eff(2))*(1d0-rfac)
               !/*root*/
               CFF(3) = 0d0 !POOLS(n+1,3)*burnt_area*combust_eff(3)
               NCFF(3) = 0d0 !POOLS(n+1,3)*burnt_area*(1d0-combust_eff(3))*(1d0-rfac)
               !/*wood*/
               CFF(4) = POOLS(n+1,4)*burnt_area*combust_eff(4)
               NCFF(4) = POOLS(n+1,4)*burnt_area*(1d0-combust_eff(4))*(1d0-rfac)
               !/*litter*/
               CFF(5) = POOLS(n+1,5)*burnt_area*combust_eff(5)
               NCFF(5) = POOLS(n+1,5)*burnt_area*(1d0-combust_eff(5))*(1d0-rfac)
               ! CWD; assume same as live wood (should be improved later)
               CFF(7) = POOLS(n+1,7)*burnt_area*combust_eff(4)
               NCFF(7) = POOLS(n+1,7)*burnt_area*(1d0-combust_eff(4))*(1d0-rfac)
               !/*fires as daily averages to comply with units*/
               FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)) * deltat_1(n)
               !              !/*update net exchangep*/
               !              NEE_out(n) = NEE_out(n)+FLUXES(n,17)
               ! determine the as daily rate impact on live tissues for use in EDC and
               ! MTT calculations
               FLUXES(n,22) = FLUXES(n,22) + ((CFF(1) + NCFF(1)) * deltat_1(n)) ! labile
               FLUXES(n,23) = FLUXES(n,23) + ((CFF(2) + NCFF(2)) * deltat_1(n)) ! foliar
               FLUXES(n,24) = FLUXES(n,24) + ((CFF(3) + NCFF(3)) * deltat_1(n)) ! root
               FLUXES(n,25) = FLUXES(n,25) + ((CFF(4) + NCFF(4)) * deltat_1(n)) ! wood

               !// update pools
               !/*Adding all fire pool transfers here*/
               POOLS(n+1,1) = POOLS(n+1,1)-CFF(1)-NCFF(1)
               POOLS(n+1,2) = POOLS(n+1,2)-CFF(2)-NCFF(2)
               POOLS(n+1,3) = POOLS(n+1,3)-CFF(3)-NCFF(3)
               POOLS(n+1,4) = POOLS(n+1,4)-CFF(4)-NCFF(4)
               POOLS(n+1,5) = POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3)
               POOLS(n+1,6) = POOLS(n+1,6)+NCFF(4)+NCFF(5)+NCFF(7)
               POOLS(n+1,7) = POOLS(n+1,7)-CFF(7)-NCFF(7)
               ! mass balance check
               where (POOLS(n+1,1:7) < 0d0) POOLS(n+1,1:7) = 0d0

               ! Some variable needed for the EDCs
               ! Reallocation fluxes for the residues, remember to
               ! convert to daily rate for consistency with the EDCs
               ! NOTE: accumulation because fire and removal may occur concurrently...
               disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) &
                                                + ((NCFF(1)+NCFF(2)+NCFF(3)) * deltat_1(n))
               disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) &
                                                + ((NCFF(4)+NCFF(5)+NCFF(7)) * deltat_1(n))
               disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) &
                                                + ((CFF(5) + NCFF(5)) * deltat_1(n))
               disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) &
                                                + ((CFF(7) - NCFF(7)) * deltat_1(n))
           endif ! burn area > 0

       endif ! fire activity

!       do nxp = 1, nopools
!          if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0) then
!              print*,"step",n, nxp
!              print*,"met",met(:,n)
!              print*,"POOLS",POOLS(n,:)
!              print*,"FLUXES",FLUXES(n,:)
!              print*,"POOLS+1",POOLS(n+1,:)
!              stop
!          endif
!       enddo
!
!       if (FLUXES(n,1) < 0d0 .or. FLUXES(n,1) /= FLUXES(n,1)) then
!           print*,"step",n, nxp
!           print*,"met",met(:,n)
!           print*,"POOLS",POOLS(n,:)
!           print*,"FLUXES",FLUXES(n,:)
!           print*,"POOLS+1",POOLS(n+1,:)
!           stop
!       endif

    end do ! nodays loop

    !!!!!!!!!!
    ! Calculate Ecosystem diagnostics
    !!!!!!!!!!

    ! calculate NEE
    NEE_out = -FLUXES(:,1) & ! GPP
            + FLUXES(:,3)+FLUXES(:,13)+FLUXES(:,14)+FLUXES(:,4) ! Respiration
    ! load GPP
    GPP_out = FLUXES(:,1)

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
!    pn = lai*avN*NUE*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,leafT)
    metabolic_limited_photosynthesis = lai*ceff*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,leafT)

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    light_limited_photosynthesis = e0 * canopy_par_MJday

    !
    ! Stomatal conductance independent variables for diffusion limited
    ! photosynthesis
    !

    ! Canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance * seconds_per_day * convert_ms1_mol_1 * gb_H2O_CO2

    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,leafT)
    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,leafT)

    ! don't forget to return
    return

  end subroutine acm_gpp_stage_1
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp_stage_2(gs)

    ! Combine the temperature (pn) and light (pl) limited gross primary productivity
    ! estimates with CO2 supply limited via stomatal conductance (gs).
    ! See acm_gpp_stage_1() for additional details on pn and pl calculation.

    implicit none

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pp, qq, ci, mult, gc, pl, pn, pd

    ! load into local variables (for code readability)
    pn = metabolic_limited_photosynthesis
    pl = light_limited_photosynthesis

    !
    ! Diffusion limited photosynthesis
    !

    ! Daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2 then all multiplied by 86400 seconds
    ! 
    ! Combining in series the stomatal and boundary layer conductances
    gc = ((gs*gs_H2Ommol_CO2mol_day) ** (-1d0) + gb_mol ** (-1d0)) ** (-1d0)

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
    pd = pd * dayl_hours_fraction

    !
    ! Estimate CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis
    acm_gpp_stage_2 = pl*pd/(pl+pd)

    ! sanity check
    if (acm_gpp_stage_2 /= acm_gpp_stage_2 .or. acm_gpp_stage_2 < 0d0) acm_gpp_stage_2 = 0d0

    ! don't forget to return
    return

  end function acm_gpp_stage_2
  !
  !------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in) :: gs_in

    ! local variables
    double precision :: gs_high, gs_store, &
                        gpp_high, gpp_low

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp_stage_2(gs_in)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp_stage_2(gs_high)

    ! determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE - ((gpp_high - gpp_low)*lai_1)

    ! remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance(deltaWP,Rtot)

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! arguments
    double precision, intent(in) :: deltaWP, & ! minlwp-wSWP (MPa)
                                       Rtot    ! total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    ! local variables
    double precision :: denom
    double precision, parameter :: max_gs = 3500d0, &  ! mmolH2O.m-2.s-1
                                   min_gs = 0.001d0, & !
                                   tol_gs = 4d0

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (aerodynamic_conductance > vsmall .and. deltaWP > vsmall) then

        ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
        max_supply = (deltaWP/Rtot) * seconds_per_day

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * ((canopy_swrad_MJday * 1d6 * dayl_seconds_1) + canopy_lwrad_Wm2) &
              + (ET_demand_coef * aerodynamic_conductance)
        denom = (denom / (lambda * max_supply * mmol_to_kg_water * dayl_seconds_1)) - slope
        stomatal_conductance = aerodynamic_conductance / (denom / psych)

        ! Convert m.s-1 to mmolH2O.m-2.s-1
        stomatal_conductance = stomatal_conductance * convert_ms1_mmol_1
        ! If conditions are dew forming then set conductance to maximum as we
        ! are not going to be limited by water demand
        if (stomatal_conductance < vsmall .or. stomatal_conductance > max_gs) stomatal_conductance = max_gs

        ! If we are potentially limited by stomatal conductance or we are using
        ! instrinsic water use efficiency (rather than WUE),
        ! then iterate to find optimum gs otherwise just go with the max...
        if (stomatal_conductance /= max_gs .or. do_iWUE ) then
            ! If there is a positive demand for water then we will solve for
            ! photosynthesis limits on gs through iterative solution
            delta_gs = 1d-3*lai ! mmolH2O/m2leaf/day
            ! Estimate inverse of LAI to avoid division in optimisation
            lai_1 = lai**(-1d0)
            ! Calculate stage one acm, temperature and light limitation which
            ! are independent of stomatal conductance effects
            call acm_gpp_stage_1
            ! Intrinsic WUE optimisation
            stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE',find_gs_iWUE,min_gs,stomatal_conductance,tol_gs)
        end if

    else

        ! if no LAI then there can be no stomatal conductance
        stomatal_conductance = 0d0
        ! set minimum (computer) precision level flow
        max_supply = vsmall

    endif ! if aerodynamic conductance > vsmall

  end subroutine calculate_stomatal_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature,input_temperature_K,input_vpd_kPa)

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory...

    implicit none

    ! arguments
    double precision, intent(in) :: input_temperature, input_temperature_K, &
                                    input_vpd_kPa

    ! local variables
    double precision :: s, mult, &
              dynamic_viscosity    ! dynamic viscosity (kg.m-2.s-1)

    !
    ! Used for soil, canopy evaporation and transpiration
    !

    ! Density of air (kg/m3)
    air_density_kg = 353d0/input_temperature_K
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    lambda = 2501000d0-2364d0*input_temperature

    ! psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*input_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = input_temperature+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = (2502.935945d0*exp(17.269d0*input_temperature/mult)) / (mult*mult)

    ! estimate frequently used atmmospheric demand component
    ET_demand_coef = air_density_kg*cpair*input_vpd_kPa

    !
    ! Used for soil evaporation and leaf level conductance
    !

    ! Determine diffusion coefficient (m2.s-1), temperature dependant (pressure dependence neglected). Jones p51; appendix 2
    ! Temperature adjusted from standard 20oC (293.15 K), NOTE that 1/293.15 = 0.003411223
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    water_vapour_diffusion = 0.0000242d0*((input_temperature_K/293.15d0)**1.75d0)

    !
    ! Used for calculation of leaf level conductance
    !

    ! Calculate the dynamic viscosity of air (kg.m-2.s-1)
    dynamic_viscosity = ((input_temperature_K**1.5d0)/(input_temperature_K+120d0))*1.4963d-6
    ! and kinematic viscosity (m2.s-1)
    kinematic_viscosity = dynamic_viscosity/air_density_kg

  end subroutine meteorological_constants
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
    double precision :: local_lai, &
           mixing_length_momentum, & ! mixing length parameter for momentum (m)
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
    local_lai = max(min_lai,lai)
    length_scale_momentum = (4d0*canopy_height) / local_lai
    mixing_length_momentum = 2d0*(ustar_Uh**3)*length_scale_momentum

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
    double precision, parameter :: leaf_width = 0.04d0, & ! leaf width (m) (original 0.08)
                                 leaf_width_1 = leaf_width ** (-1d0), &
                                           Pr = 0.72d0, & ! Prandtl number
                                      Pr_coef = 1.05877d0 !1.18d0*(Pr**(0.33d0))
    ! local variables
    double precision :: &
         nusselt_forced & ! Nusselt value under forced convection
             ,Sh_forced & ! Sherwood number under forced convection
                    ,Re   ! Reynolds number

    ! Reynold number
!    Re = (leaf_width*canopy_wind)/kinematic_viscosity
    ! calculate nusselt value under forced convection conditions
!    nusselt_forced = (1.18d0*(Pr**(0.33d0))*(sqrt(Re)))
!    nusselt_forced = Pr_coef*(sqrt(Re))
    ! update specific Sherwood numbers
!    Sh_forced = 0.962d0*nusselt_forced
!    Sh_forced = (0.962d0*Pr_coef*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity)))
    ! Estimate the the forced conductance of water vapour
    gv_forced = 0.962d0*Pr_coef*sqrt((leaf_width*canopy_wind)/kinematic_viscosity) &
              * water_vapour_diffusion*leaf_width_1 * 0.5d0 * lai

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
    canopy_wind = ustar * vonkarman_1 * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength(doy,lat)

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate day length in hours and seconds

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
    aob = max(-1d0,min(1d0,sinld / cosld))

    ! estimate day length in hours and seconds and upload to module variables
    dayl_hours = 12d0 * ( 1d0 + 2d0 * asin( aob ) * pi_1 )
    dayl_seconds = dayl_hours * seconds_per_hour

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
         transmitted_fraction, & ! fraction of LW which is not incident on the canopy
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy, & ! assuming isothermal condition (W.m-2)
            trans_lw_fraction, &
        reflected_lw_fraction, &
         absorbed_lw_fraction, &
      canopy_release_fraction, & ! fraction of longwave emitted from within the canopy to ultimately be released
   canopy_absorption_from_sky, & ! canopy absorbed radiation from downward LW (W.m-2)
  canopy_absorption_from_soil, & ! canopy absorbed radiation from soil surface (W.m-2)
                  canopy_loss, & ! longwave radiation released from canopy surface (W.m-2).
                                 ! i.e. this value is released from the top and
                                 ! the bottom
       soil_incident_from_sky, &
     soil_absorption_from_sky, & ! soil absorbed radiation from sky (W.m-2)
  soil_absorption_from_canopy    ! soil absorbed radiation emitted from canopy (W.m-2)

    ! local parameters
    double precision, parameter :: clump = 1d0    & ! leaf clumping factor (1 = uniform)
                                  ,decay = -0.5d0   ! decay coefficient for incident radiation

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! First, we consider how much radiation is likely to be incident on the
    ! canopy, or put another way what fraction passes straight through the
    ! canopy?
    transmitted_fraction = exp(decay * lai * clump)

    ! second, we partition the radiation which is incident on the canopy into
    ! that which is transmitted, reflected or absorbed.

    ! Likewise we assume that the reflectance and transmittance are equal
    ! However, the non-linear interception under Beer's Law means that actual
    ! canopy transmittenace to the soil surface and reflectance back to the sky
    ! skews towards reduced transmittance at higher LAI. Both transmittance and
    ! reflectance follow linear a relationship with a common intercept
    ! NOTE: 0.02 = (1-emissivity) * 0.5
    trans_lw_fraction     = 0.02d0 * (1d0 - (lai/(lai+lai_half_lwrad_transmitted)))
    reflected_lw_fraction = 0.02d0 * (1d0 - (lai/(lai+lai_half_lwrad_reflected)))
    ! Absorption is the residual
    absorbed_lw_fraction = 1d0 - trans_lw_fraction - reflected_lw_fraction

    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky
    canopy_release_fraction = max_lw_escape * (1d0 - (max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release))

    !!!!!!!!!!
    ! Distribute longwave from sky
    !!!!!!!!!!

    ! Estimate the radiation which directly bypasses the canopy...
    soil_incident_from_sky = lwrad * transmitted_fraction
    ! ...and update the canopy intercepted radiation
    lwrad = lwrad - soil_incident_from_sky

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * absorbed_lw_fraction
    ! Long wave absorbed by soil from the sky, soil absorption assumed to be
    ! equal to emissivity
    soil_incident_from_sky = soil_incident_from_sky + (trans_lw_fraction * lwrad)
    soil_absorption_from_sky = soil_incident_from_sky * emissivity
    ! Long wave reflected directly back into sky
    sky_lwrad_Wm2 = lwrad * reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! Calculate longwave radiation coming up from the soil plus the radiation
    ! which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (soil_incident_from_sky * (1d0-emissivity))
    ! First how much directly bypasses the canopy...
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * transmitted_fraction)
    canopy_absorption_from_soil = canopy_absorption_from_soil * (1d0 - transmitted_fraction)
    ! Second, use this to estimate the longwave returning to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * trans_lw_fraction)
    ! Third, now calculate the longwave from the soil surface absorbed by the
    ! canopy
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

    !!!!!!!!!!
    ! Approximate daily impact on isothermal to net longwave
    !!!!!!!!!!

    ! Isothermal -> net radiation in SPA is largely linear, therefore we
    ! retrieve the linear correction here
    soil_lwrad_Wm2 = (soil_lwrad_Wm2 * soil_iso_to_net_coef) + soil_iso_to_net_const

  end subroutine calculate_longwave_isothermal
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_radiation_balance

    implicit none

    ! subroutine call ensures that both shortwave and longwave radiation balance
    ! are calculated at the same time but with the more readable code split
    ! between a shortwave and longwave specific subroutines.

    ! NOTE: that this code provides a daily timescale linear correction on
    ! isothermal longwave balance to net based on soil surface incident shortwave
    ! radiation

    ! declare local variables
    double precision :: delta_iso

    ! Estimate shortwave radiation balance
    call calculate_shortwave_balance
    ! Estimate isothermal long wave radiation balance
    call calculate_longwave_isothermal(meant,meant)
    ! Apply linear correction to soil surface isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = soil_iso_to_net_coef * (soil_swrad_MJday * 1d6 * dayl_seconds_1) + soil_iso_to_net_const
    soil_lwrad_Wm2 = soil_lwrad_Wm2 + delta_iso

  end subroutine calculate_radiation_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation
    ! (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

    ! local variables
    double precision :: balance                    &
                       ,transmitted_fraction       &
                       ,par,nir                    &
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
    double precision, parameter :: clump = 1d0    & ! leaf clumping factor (1 = uniform)
                                  ,decay = -0.5d0   ! decay coefficient for incident radiation

    !!!!!!!!!!
    ! Determine canopy absorption, reflectance and transmittance as function of
    ! LAI
    !!!!!!!!!!

    ! First, we consider how much radiation is likely to be incident on the
    ! canopy, or put another way what fraction passes straight through the
    ! canopy?
    transmitted_fraction = exp(decay * lai * clump)

    ! Second, of the radiation which is incident on the canopy what fractions
    ! are transmitted through, reflected from or absorbed by the canopy

    ! Canopy transmitted of PAR & NIR radiation towards the soil
    trans_par_fraction = max(0d0,max_lai_par_transmitted * lai + max_par_transmitted)
    trans_nir_fraction = max(0d0,max_lai_nir_transmitted * lai + max_nir_transmitted)
    ! Canopy reflected of near infrared and photosynthetically active radiation
    reflected_nir_fraction = 1d0 - ( (lai*max_lai_nir_reflection) &
                                    /(lai+lai_half_nir_reflection) )
    reflected_par_fraction = 1d0 - ( (lai*max_lai_par_reflection) &
                                    /(lai+lai_half_par_reflection) )
    reflected_nir_fraction = max_nir_reflected * reflected_nir_fraction
    reflected_par_fraction = max_par_reflected * reflected_par_fraction
    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = 1d0 - reflected_nir_fraction - trans_nir_fraction
    absorbed_par_fraction = 1d0 - reflected_par_fraction - trans_par_fraction

    !!!!!!!!!!
    ! Estimate canopy absorption of incoming shortwave radiation
    !!!!!!!!!!

    ! Estimate multiple use par and nir components
    par = sw_par_fraction * swrad
    nir = (1d0 - sw_par_fraction) * swrad

    ! Estimate the radiation which directly bypasses the canopy...
    trans_par_MJday = par * transmitted_fraction
    trans_nir_MJday = nir * transmitted_fraction
    ! ...and update the canopy intercepted radiation
    par = par - trans_par_MJday
    nir = nir - trans_nir_MJday

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected
    ! by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = par * absorbed_par_fraction
    canopy_nir_MJday = nir * absorbed_nir_fraction
    trans_par_MJday = trans_par_MJday + (par * trans_par_fraction)
    trans_nir_MJday = trans_nir_MJday + (nir * trans_nir_fraction)
    refl_par_MJday = par * reflected_par_fraction
    refl_nir_MJday = nir * reflected_nir_fraction

    !!!!!!!!!
    ! Estimate soil absorption of shortwave passing through the canopy
    !!!!!!!!!

    ! Then the radiation incident and ultimately absorbed by the soil surface
    ! itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday * soil_swrad_absorption
    soil_nir_MJday = trans_nir_MJday * soil_swrad_absorption
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = soil_nir_MJday + soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection / absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! calculate multiple use variables
    par = trans_par_MJday-soil_par_MJday
    nir = trans_nir_MJday-soil_nir_MJday
    ! how much of the reflected radiation directly bypasses the canopy...
    refl_par_MJday = refl_par_MJday + (par * transmitted_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir * transmitted_fraction)
    ! ...and update the canopy on this basis
    par = par * (1d0-transmitted_fraction)
    nir = nir * (1d0-transmitted_fraction)

    ! Update the canopy radiation absorption based on the reflected radiation
    ! (MJ.m-2.day-1)
    canopy_par_MJday = canopy_par_MJday + (par * absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + (nir * absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is
    ! now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + (par * trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir * trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = canopy_par_MJday + canopy_nir_MJday

    ! check energy balance
!    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
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
    water_flux = 0d0 
    ratio = 0d0 ; ratio(1) = 1d0 ; root_mass = 0d0
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    !    transpiration_resistance = (gplant * lai)**(-1d0)
    transpiration_resistance = canopy_height / (gplant * max(min_lai,lai))

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! calculate soil depth to which roots reach
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
        root_mass(1) = fine_root_biomass * 0.5d0
        ! Then quantify how much additional root is found in the top soil layer
        ! assuming that the top 25 % depth is found somewhere within the top
        ! layer
        bonus = (fine_root_biomass-root_mass(1)) &
              * (layer_thickness(1)-root_depth_50) / (root_reach - root_depth_50)
        root_mass(1) = root_mass(1) + bonus
        ! partition the remaining root biomass between the seconds and third
        ! soil layers
        if (root_reach > sum(layer_thickness(1:2))) then
            root_mass(2) = (fine_root_biomass - root_mass(1)) &
                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
            root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
        else
            root_mass(2) = fine_root_biomass - root_mass(1)
        endif

    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then

        ! Greater than 50 % of fine root biomass found in the top two soil
        ! layers. We will divide the root biomass uniformly based on volume,
        ! plus bonus for the second layer (as done above)
        root_mass(1) = fine_root_biomass * (layer_thickness(1)/root_depth_50)
        root_mass(2) = fine_root_biomass * ((root_depth_50-layer_thickness(1))/root_depth_50)
        root_mass(1:2) = root_mass(1:2) * 0.5d0

        ! determine bonus for the seconds layer
        bonus = (fine_root_biomass-sum(root_mass(1:2))) &
              * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach-root_depth_50))
        root_mass(2) = root_mass(2) + bonus
        root_mass(3) = fine_root_biomass - sum(root_mass(1:2))

    else

        ! Greater than 50 % of fine root biomass stock spans across all three
        ! layers
        root_mass(1:2) = fine_root_biomass * 0.5d0 * (layer_thickness(1:2)/root_depth_50)
        root_mass(3) = fine_root_biomass - sum(root_mass(1:2))

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
    demand = abs(minlwp)+head*canopy_height
    ! now loop through soil layers, where root is present
    do i = 1, nos_root_layers
       if (root_mass(i) > 0d0) then
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
    if (meant < 1d0) then
        water_flux(1) = 0d0
        ratio(1) = 0d0
        ratio(2:nos_root_layers) = layer_thickness(2:nos_root_layers) / sum(layer_thickness(2:nos_root_layers))
    else
        ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif

    ! calculate sum value
    sum_water_flux = sum(water_flux)
    if (sum_water_flux <= vsmall) then
        uptake_fraction = 0d0 ; uptake_fraction(1) = 1d0
    else
        ! calculate uptake fraction
        uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
    endif

    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum_water_flux

    ! finally convert transpiration flux (mmolH2O.m-2.s-1)
    ! into kgH2O.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
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
    double precision  sqrt_cd1_lai, local_lai
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
!                            ustar_Uh_max = 0.3,   & ! Maximum observed ratio of
                                                    ! (friction velocity / canopy top wind speed) (m.s-1)
        ustar_Uh_max = 1d0, ustar_Uh_min = 0.2d0, &
                                        Cw = 2d0, &  ! Characterises roughness sublayer depth (m)
                                     phi_h = 0.19314718056d0 ! Roughness sublayer influence function;

    ! describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law


    ! assign new value to min_lai to avoid max min calls
    local_lai = max(min_lai,lai)
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (1d0-((1d0-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    ! phi_h = log(Cw)-1d0+Cw**(-1d0) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl = ((1d0-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

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
                                    ,VPDfac_range_1,base_leaf_fall        &
                                    ,pot_leaf_fall,pot_leaf_growth        &
                                    ,mean_min_airt,mean_daylength,mean_vpd&
                                    ,deltaWP,Rtot,GPP_current,foliage     &
                                    ,gpp_crit_frac,GSI,leaf_fall,leaf_growth)

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
                                   ,base_leaf_fall & !
                                    ,pot_leaf_fall & !
                                  ,pot_leaf_growth

    double precision, intent(inout) :: GSI(nodays) &
                                      ,leaf_fall,leaf_growth

    ! declare local variables
    integer :: gsi_lag, m, interval
    double precision :: infi     &
                       ,tmp      &
                       ,deltaNCE &
                       ,deltaGPP &
                        ,deltaRm &
                       ,C_invest &
                         ,NCE_LL &
                       ,lai_save &
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

    ! GSI is the product of 3 limiting factors for temperature, photoperiod and
    ! vapour pressure deficit that scale linearly between 0 to 1 as a function
    ! of calibrated min and max value.
    ! Photoperiod, VPD and avgTmin are direct input

    ! Temperature limitation, then restrict to 0-1; correction for k-> oC
    itemp(current_step) = min(1d0,max(0d0,(mean_min_airt-(Tfac_min-freeze)) * Tfac_range_1))
    ! Photoperiod limitation (seconds)
    iphoto(current_step) = min(1d0,max(0d0,(mean_daylength-Photofac_min) * Photofac_range_1))
    ! VPD limitation (kPa)
    ivpd(current_step) = min(1d0,max(0d0,1d0 - ((mean_VPD-VPDfac_min) * VPDfac_range_1)))
    ! Calculate and store the GSI index
    GSI(current_step) = itemp(current_step) * ivpd(current_step) * iphoto(current_step)

    ! load lag for linear regression
    gsi_lag = gsi_lag_remembered

    !!!
    ! Estimate of change (i.e. gradient) in the GSI / NCE

    ! Determine GSI / NCE section to have linear regression applied to and
    ! determine the number of values, i.e. the interval
    if (current_step < gsi_lag) then
        if (current_step == 1) then
            gsi_history(2) = GSI(current_step)
            interval = 2
        else
            gsi_history(1:current_step) = GSI(1:current_step) 
            interval = current_step
        endif
    else 
        gsi_history(1:gsi_lag) = GSI((current_step-gsi_lag+1):current_step)
        interval = gsi_lag
    end if
    ! Now calculate the linear gradient
    gradient = linear_model_gradient(tmp_x(1:interval),gsi_history(1:interval),interval)

    ! store lag to keep fresh in memory - yes this is a hack to get around a
    ! memory problem
    gsi_lag_remembered = gsi_lag
    
    ! first assume that nothing is happening
    leaf_fall = 0d0   ! leaf turnover
    leaf_growth = 0d0 ! leaf growth

    ! Reset marginal return variables
    deltaGPP = 0d0 ; deltaRm = 0d0 ; deltaNCE = 0d0

    ! Everything else in here was needed to keep track of GSI values but
    ! ultimately if there is not labile available no growth can occur
    if (gradient > fol_turn_crit .and. deltaWP < 0d0 .and. avail_labile > 0d0) then

        ! Estimate approximate the potential leaf area increment
        ! using the Reich maintence respiration Q10 as proxy for potential
        ! metabolic activity
!        leaf_growth = pot_leaf_growth*Q10_adjustment(n)
        leaf_growth = pot_leaf_growth*GSI(current_step)
        ! calculate potential C allocation to leaves
        tmp = avail_labile * &
             (1d0-(1d0-leaf_growth)**deltat(current_step))*deltat_1(current_step)
        ! Store total investmentl, then adjust for loss to maintenance Rg
        C_invest = tmp ; tmp = tmp * one_Rg_fraction
        ! Determine the daily increment in maintenance respiration
        deltaRm = Rm_leaf_per_gC * tmp

        ! calculate new leaf area, GPP return
        lai = (foliage+tmp) * SLA
        aerodynamic_conductance = aerodynamic_conductance * (lai / lai_save)
        call calculate_stomatal_conductance(abs(deltaWP),Rtot)
        call calculate_shortwave_balance
        if (lai_save < vsmall) then
            call calculate_aerodynamic_conductance
        endif ! lai_save < vsmall
        ! calculate stomatal conductance of water
        if (stomatal_conductance > vsmall) then
             call acm_gpp_stage_1
            tmp = acm_gpp_stage_2(stomatal_conductance)
        else
            tmp = 0d0
        endif
        ! Estimate total daily GPP return from new leaf
        deltaGPP = (tmp - GPP_current)

        ! Estimate the change in net carbon export by the canopy per day,
        ! then scale the initial investment costs by leaf lifespan (days) and
        ! substract the initial investment cost
        ! i.e. gC/m2/day/(gCinvest/LL)
!        deltaNCE = (((deltaGPP - deltaRm) * leaf_life) - C_invest) / C_invest
        C_invest = C_invest / leaf_life
        deltaNCE = (((deltaGPP - deltaRm)) - C_invest) / C_invest
        ! Is the marginal return for GPP (over the mean life of leaves)
        ! less than increase in maintenance respiration and C required
        ! to growth?
        if (deltaNCE < gpp_crit_frac) then
            leaf_growth = 0d0 
        end if

    endif ! deltaWP < 0

    if (avail_labile < vsmall) then

        ! We want to lose a chunk of leaf quickly!
        leaf_fall = leaf_fall + pot_leaf_fall

    else if (gradient < fol_turn_crit .or. GSI(current_step) < vsmall) then

        ! We want to lose a chunk of leaf quickly!
        leaf_fall = (pot_leaf_fall * (1d0-GSI(current_step)))

    end if ! gradient decline / increase choice
    ! If either growing or turnover assume that some minimum turnover occurs
    if (leaf_fall == 0d0 .and. leaf_growth == 0d0) leaf_fall = base_leaf_fall

    ! restore original value back from memory
    lai = lai_save
    canopy_lwrad_Wm2 = canopy_lw_save ; soil_lwrad_Wm2 = soil_lw_save
    canopy_swrad_MJday = canopy_sw_save ; canopy_par_MJday = canopy_par_save
    soil_swrad_MJday = soil_sw_save ; stomatal_conductance = gs_save

  end subroutine calculate_leaf_dynamics
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
    ! calculations have been split.                    !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    ! local variables..
    double precision :: denominator, numerator

    arrhenious = a * exp( b * (t - 25d0) / (t + freeze) )

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
         opt_max_scaling = 0d0
    else
         dummy     = exp( log((max_val - current) / (max_val - optimum)) * kurtosis * (max_val - optimum) )
         opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
    end if

  end function opt_max_scaling
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

    ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
    soilR2 = root_resist / (root_mass*root_reach_in)
    ! Estimate the soil to plant flow of water mmolH2O/m2/s
    plant_soil_flow = demand/(transpiration_resistance + soilR2)

    ! return
    return

  end function plant_soil_flow
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
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
    !                      / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )
    ! Linear regression done as single line to reduce assignment requirements
    linear_model_gradient = ( (dble(interval)*sum(x*y)) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sum(x*x)) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !--------------------------------------------------------------------------
  !
  double precision function Rm_reich_Q10(air_temperature)

    ! Calculate Q10 temperature adjustment used in estimation of the
    ! Maintenance respiration (umolC.m-2.s-1) calculated based on modified
    ! version of the Reich et al (2008) calculation.

    ! arguments
    double precision, intent(in) :: air_temperature ! input temperature of metabolising tissue (oC)

    ! local variables
    double precision, parameter :: Q10 = 2d0,  & ! Q10 response of temperature (baseline = 20oC) ;INITIAL VALUE == 2
                                                 ! Mahecha, et al. (2010) Global Convergence in the Temperature Sensitivity
                                                 ! of Respiration at Ecosystem Level. Science 329 , 838 (2010);
                                                 ! DOI: 10.1126/science.1189587. value reported as 1.4
                         Q10_baseline = 20d0     ! Baseline temperature for Q10 ;INITIAL VALUE == 20;

    ! calculate instantaneous Q10 temperature response
    Rm_reich_Q10 = Q10**((air_temperature-Q10_baseline)*0.1d0)

    ! explicit return command
    return

  end function Rm_reich_Q10
  !
  !--------------------------------------------------------------------------
  !
  double precision function Rm_reich_N(CN_pool &
                                      ,N_exponential_response  &
                                      ,N_scaler_intercept)

    ! Calculate the nitrgen response on maintenance respiration (nmolC.g-1.s-1)
    ! calculated based on modified version of the Reich et al (2008)
    ! calculation.
    ! Note the output here is invarient for given CN ratio

    ! NOTE: Rm_tissue =  Rm_reich_Q10 * Rm_reich_N * C_POOL * 2 * 0.001
    !       umolC/m2/s = dimensionless * nmolC/g/s * gC/m * (correct g->gC) *
    !       (nmolC->umolC)

    ! arguments
    double precision, intent(in) ::        CN_pool, & ! C:N ratio for current pool (gC/gN)
                            N_exponential_response, & ! N exponential response vcoefficient (1.277/1.430)
                                N_scaler_intercept    ! N scaler (baseline) (0.915 / 1.079)

    ! local variables
    double precision, parameter :: N_g_to_mmol = 71.42857 ! i.e. (1d0/14d0)*1d3 where 14 = atomic weight of N
    double precision :: Nconc ! Nconc =mmol g-1

    ! calculate leaf maintenance respiration (nmolC.g-1.s-1)
    ! NOTE: that the coefficients in Reich et al., 2008 were calculated from
    ! log10 linearised version of the model, thus N_scaler is already in log10()
    ! scale. To remove the need of applying log10(Nconc) and 10**Rm_reich the
    ! scaler is reverted instead to the correct scale for the exponential form
    ! of the equations.

    ! calculate N concentration per g biomass.
    ! A function of C:N
    Nconc = ((CN_pool*2d0)**(-1d0)) * N_g_to_mmol

    ! leaf maintenance respiration (nmolC.g-1.s-1) at 20 oC
    Rm_reich_N = (10d0**N_scaler_intercept) * Nconc ** N_exponential_response

    ! explicit return command
    return

  end function Rm_reich_N
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
    double precision,intent(in) :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
      double precision function func( xval )
        double precision ,intent(in) :: xval
      end function func
    end interface

    ! local variables..
    integer            :: iter
    integer, parameter :: ITMAX = 10
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
end module CARBON_MODEL_MOD
