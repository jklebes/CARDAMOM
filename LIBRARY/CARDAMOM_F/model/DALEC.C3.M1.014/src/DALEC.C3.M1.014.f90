!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
! assimilate observations and ecological theory to retrieve parameters for the 
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently. 
! Copyright (C) 2024  University of Edinburgh, 
!                     Mathew Williams (mat.williams@ed.ac.uk), 
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see < https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This file contains the source code of DALEC.C3.M1
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from the following primary references:
! Sus et al., (2010), https://doi.org/10.1016/j.agee.2010.06.012.
! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
! Smallman et al., (2017), https://doi.org/10.1002/2016JG003520.
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA) and
! Oliver Sus (UoE, now at EUMETSAT, Darmstadt).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CARBON_MODEL_MOD

  implicit none

  ! make all private
  private

  ! explicit publics
  public:: CARBON_MODEL          &
           ,vsmall                &
           ,arrhenious            &
           ,acm                   &
           ,calculate_daylength   &
           ,opt_max_scaling       &
           ,freeze                &
           ,linear_model_gradient &
           ,seconds_per_day       &
           ,cica_time             &
           ,deltat_1              &
           ,DS_time               &
           ,avN_time              &
           ,doy                   &
           ,lai                   &
           ,days_per_step         &
           ,days_per_step_1       &
           ,dayl_hours            &
           ,resp_rate_temp_coeff  &
           ,nos_soil_layers       &
           ,soil_frac_clay        &
           ,soil_frac_sand        &
           ,dim_1, dim_2           &
           ,nos_trees             &
           ,nos_inputs            &
           ,leftDaughter          &
           ,rightDaughter         &
           ,nodestatus            &
           ,xbestsplit            &
           ,nodepred              &
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
  logical:: do_iWUE = .true. ! Use iWUE or WUE for stomatal optimisation
  double precision, parameter:: vsmall = tiny(0d0)*1d3 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  integer, parameter:: nos_root_layers = 2, nos_soil_layers = nos_root_layers+1
  double precision, parameter:: pi = 3.1415927d0,  &
                               pi_1 = 0.3183099d0,  & ! pi**(-1d0)
                                pi2 = 9.869604d0,   & ! pi**2d0
                             two_pi = 6.283185d0,   & ! pi*2d0
                         deg_to_rad = 0.01745329d0, & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0*deg_to_rad )
                             freeze = 273.15d0

  ! photosynthesis/respiration parameters
  double precision, parameter :: &
                        Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                      ! lost as growth respiration
                                                      ! (i.e. 0.28 .eq. xNPP)
                    one_Rg_fraction = 1d0-Rg_fraction

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                    seconds_per_day = 86400d0,        & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05      ! Inverse of seconds per day

  !!!!!!!!!
  ! Module level variables
  !!!!!!!!!

  ! hydraulic model variables
  integer:: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers):: soil_frac_clay, soil_frac_sand  ! clay and soil fractions of soil

  ! Module level variables for ACM_GPP parameters
  double precision:: ci  ! Internal CO2 concentration (ppm)

  ! Module level variables for step specific met drivers
  double precision ::  doy, & ! Day of year
                     lai_1, & ! inverse of LAI
                       lai    ! leaf area index (m2/m2)

  ! Module level variables for step specific timing and location information
  integer:: steps_per_year
  double precision ::          days_per_step, & !
                             days_per_step_1, & !
                                  dayl_hours, & ! day length in hours
                                    latitude, & ! latitude ()-90/90)
                            latitude_radians, & ! latitude in radians
                        sin_latitude_radians, & ! sin(latitude_radians)
                        cos_latitude_radians, & ! cos(latitude_radians)
                                 declination, & ! Solar declination, function of day of year
                   cosine_solar_zenith_angle    ! Cosine zenith angle of the timestep

  double precision, dimension(:), allocatable:: deltat_1, & ! inverse of decimal days
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1, &
                                                  DS_time, &
                                                 avN_time, & ! average foliar N (gN/m2)
                                                cica_time    ! Internal vs ambient CO2 concentrations

  ! variables local to this module..
  integer ::   plough_day, & ! day-of-year when field is ploughed  (default)
                  sow_day, & ! day-of-year when field is sown      (default)
              harvest_day, & ! day-of-year when field is harvested (default)
                    stmob, & ! remoblise stem C to labile (1 = on)
   turnover_labile_switch    ! begin turnover of labile C

  logical:: vernal_calcs, &  ! do vernalisation calculations?
                 ploughed, &  !
          use_seed_labile, & != .False. ! whether to use seed labile for growth
                     sown, & != .False. ! has farmer sown crop yet?
                  emerged    != .False. ! has crop emerged yet?

  double precision ::             gpp_acm, & ! gross primary productivity (gC.m-2.day-1)
                      stock_storage_organ, & ! storage organ C pool, i.e. the desired crop (gC.m-2)
                       stock_dead_foliage, & ! dead but still standing foliage (gC.m-2)
                          stock_resp_auto, & ! autotrophic respiration pool (gC.m-2)
                             stock_labile, & ! labile C pool (gC.m-2)
                            stock_foliage, & ! foliage C pool (gC.m-2)
                               stock_stem, & ! stem C pool (gC.m-2)
                              stock_roots, & ! roots C pool (gC.m-2)
                             stock_litter, & ! litter C pool (gC.m-2)
                      stock_soilOrgMatter, & ! SOM C pool (gC.m-2)
                                resp_auto, & ! autotrophic respiration (gC.m-2.d-1)
                            resp_h_litter, & ! litter heterotrophic respiration (gC.m-2.d-1)
                     resp_h_soilOrgMatter, & ! SOM heterotrophic respiration (gC.m-2.d-1)
                                      npp, & ! net primary productivity (gC.m-2.d-1)
                                nee_dalec, & ! net ecosystem exchange (gC.m-2.d-1)
                                       DS, & ! Developmental state and initial condition
                                      LCA, & ! leaf mass area (gC.m-2)
               alloc_to_storage_organ_old, & ! rolling average allocation of GPP to storage organ (gC.m-2.d-1)
                       decomposition_rate, & ! decomposition rate (frac/day)
                       frac_GPP_resp_auto, & ! fraction of GPP allocated to autotrophic carbon pool
                    turnover_rate_foliage, & ! turnover rate of foliage (frac/day)
                       turnover_rate_stem, & ! same for stem
                     turnover_rate_labile, & ! same for labile
                  turnover_rate_resp_auto, & ! same for autotrophic C pool
                   resp_cost_labile_trans, & ! labile lost to respiration per gC labile to GPP
               mineralisation_rate_litter, & ! mineralisation rate of litter
        mineralisation_rate_soilOrgMatter, & ! mineralisation rate of SOM
                                    PHUem, & ! emergance value for phenological heat units
                                      PHU, & ! phenological heat units
                                   DR_pre, & ! development rate coefficient DS 0->1
                                  DR_post, & ! development rate coefficient DS 1->2
                                     tmin, & ! min temperature for development
                                     tmax, & ! max temperature for development
                                     topt, & ! optimum temperature for development
                                   tmin_v, & ! min temperature for vernalisation
                                   tmax_v, & ! max temperature for vernalisation
                                   topt_v, & ! optimim temperature for vernalisation
                                      VDh, & ! effective vernalisation days when plants are 50 % vernalised
                                       VD, & ! count of vernalisation days
                                 RDRSHMAX, & ! maximum rate of self shading turnover (frac/day)
                                     PHCR, & ! critical value of photoperiod for development
                                     PHSC, & ! photoperiod sensitivity
                                     raso, & ! rolling average for alloc to storage organ
                                 max_raso, & ! maximum value for rolling average alloc to storage organ
                                    BM_EX, & ! biomass extracted in addition to the storage organ
                 HARVESTextracted_foliage, & ! Foliage removed by harvest activity
                    HARVESTextracted_stem, & ! Stem removed by harvest activity
            HARVESTextracted_dead_foliage, & ! Dead still standing foliage removed by harvest activity
                  HARVESTextracted_labile, & ! Labile removed by harvest activity
                    HARVESTlitter_foliage, & ! Foliage converted to litter by harvest
                       HARVESTlitter_stem, & ! Stem converted to litter by harvest
               HARVESTlitter_dead_foliage, & ! Dead standing foliage converted to litter by harvest
                  HARVESTlitter_resp_auto, & ! Autotrophic pool converted to litter by harvest
                     HARVESTlitter_labile, & ! Labile converted to litter by harvest
                       PLOUGHlitter_roots, & ! Plough induced litter generation from roots
                                       HI, & ! Harvest index, the ratio of yield to shoot C
                                    yield, & ! crop yield (gC.m-2)
                       alloc_to_resp_auto, & ! amount of carbon to allocate to autotrophic respiration pool
                      turnover_rate_roots, & ! turnover over rate of roots interpolated each time step
                                  gso_max, & !
                             max_raso_old, & !
                                 raso_old, & !
                  resp_cost_labile_to_npp, & ! respiratory cost of moving carbon..from labile to NPP pool
                  resp_cost_npp_to_labile, & ! ..from remaining NPP to labile pool
              resp_cost_foliage_to_labile, & ! ..from foliage to labile pool
                 resp_cost_stem_to_labile, & ! ..from stem to labile pool
                                resp_rate, & ! rate of respiration at given temperature
                                   Cshoot, & !
                                       DR, & !
                          fol_frac_intpol, & !
                         stem_frac_intpol, & !
                                 fP, fT, fV, & !
                        foliage_to_labile, & !
                           stem_to_labile, & !
                         root_frac_intpol, & !
                                   avtemp, & !
                   alloc_to_storage_organ, & !
                       litterfall_foliage, & !
                          litterfall_stem, & !
                         litterfall_roots, & !
                            decomposition, & !
                                npp_shoot, & !
                        alloc_from_labile, & !
                          alloc_to_labile, & !
                           alloc_to_roots, & !
                         alloc_to_foliage, & !
                            alloc_to_stem, & !
                                    RDRSH, & !
                                    RDRDV, & !
                                      RDR

  !
  ! some hardcoded crop parameters
  !

  ! defines Q10 = 2 in exponential temperature response for heterotrophic
  ! respiration
  double precision, parameter:: resp_rate_temp_coeff = 0.0693d0
  ! residue fraction of leaves left post harvest
  double precision, parameter:: lv_res = 0.1d0
  ! residue fraction of stem left post harvest
  double precision, parameter:: st_res = 0.1d0
  ! LAI above which self shading turnover occurs
  double precision, parameter:: LAICR = 4d0
  !double precision, parameter:: LAICR = 5d0
  ! allocation to storage organ relative to GPP
  double precision, parameter:: rel_gso_max = 0.35d0
  !double precision, parameter:: rel_gso_max = 0.45d0

  save

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start, finish, met, pars, deltat, nodays, lat, lai_out &
                         ,NEE_out, FLUXES, POOLS, nopars, nomet, nopools, nofluxes &
                         ,GPP_out, stock_seed_labile, DS_shoot, DS_root, fol_frac    &
                         ,stem_frac, root_frac, DS_LRLV, LRLV, DS_LRRT, LRRT)

    !
    ! The Data Assimilation Linked Ecosystem Carbon-CROP-BUCKET (DALEC_CROP) model.
    ! modified from Sus et al., (2010)
    !
    ! The Aggregated Canopy Model for Gross Primary Productivity (Williams et al., 1997)
    !
    ! This version was coded by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! Version 1.0: 15/07/2014
    ! Version 2.0: 25/04/2024-Integrated into CARDAMOM as substraction of BUCKET model

    implicit none

    ! declare input variables
    integer, intent(in):: start    &
                          ,finish   &
                          ,nopars     & ! number of paremeters in vector
                          ,nomet      & ! number of meteorological fields
                          ,nofluxes   & ! number of model fluxes
                          ,nopools    & ! number of model pools
                          ,nodays       ! number of days in simulation

    double precision, intent(in):: met(nomet, nodays)   & ! met drivers
                         ,stock_seed_labile             & ! seed carbon to get things going
                         ,deltat(nodays)                & ! time step in decimal days
                         ,pars(nopars)                  & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(:), intent(inout) ::          DS_shoot, & !
                                                               DS_root, & !
                                                              fol_frac, & !
                                                             stem_frac, & !
                                                             root_frac, & !
                                                               DS_LRLV, & !
                                                                  LRLV, & !
                                                               DS_LRRT, & !
                                                                  LRRT    !

    double precision, dimension(nodays), intent(inout):: lai_out & ! leaf area index
                                               ,GPP_out & ! Gross primary productivity
                                               ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1), nopools), intent(inout):: POOLS  ! vector of ecosystem pools

    double precision, dimension(nodays, nofluxes), intent(inout):: FLUXES  ! vector of ecosystem fluxes

    ! declare local variables
    double precision:: airt_weighting(3) &
                         ,acm_forcings(7) & ! ACM inputs (LAI+met)
                                    ,infi   ! used to calculate infinity for diagnositc

    integer:: nxp, n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY

    ! POOLS are:
    ! 1 = labile
    ! 2 = foliar
    ! 3 = root
    ! 4 = wood
    ! 5 = litter
    ! 6 = som
    ! 7 = autotrophic
    ! 8 = storage organ C
    ! 9 = dead still standing foliage

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile release
    ! 9 = alloc to storage
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som (decomposition)
    ! 16 = alloc to autotrophic pool

    ! PARAMETERS
    ! 16 values

    ! p(1) decomposition rate (frac/day)
    ! p(2) Fraction of GPP allocated to autotrophic C pool
    ! p(3) DR coef for DS (0->1)
    ! p(4) DR coef for DS (1->2)
    ! p(5) turnover rate of foliage (frac/day)
    ! p(6) Turnover rate of wood/stem (frac/day)
    ! p(7) maximum rate of foliar turnover due to self shading (frac/day)
    ! p(8) effective vernalisation days when plant is 50 % vernalised
    ! p(9) mineralisation rate of som (frac/day)
    ! p(10) mineralisation rate of litter (frac/day)
    ! p(11) = log10(avgN)
    ! p(12) = sow day
    ! p(13) = labile lost to respiration per gC labile top GPP
    ! p(14) = phenological heat units needed for emergence
    ! p(15)  ! harvest day (doy)
    ! p(16)  ! plough day (doy)
    ! p(17)  ! leaf mass area (gC.m-2)
    ! p18, p19, p20, p21, p22, p23, p24, p25 = labile, foliar, roots, stem, litter, 
    ! som, 
    ! autotrophic and storage organ pools respectively
    ! p(26)  ! min temperature for development
    ! p(27)  ! max temperature for development
    ! p(28)  ! optimum temperature for development
    ! p(29)  ! min temperature for vernalisation
    ! p(30)  ! max temperature for vernalisation
    ! p(31)  ! optimim temperature for vernalisation
    ! p(32)  ! critical value of photoperiod for development
    ! p(33)  ! photoperiod sensitivity
    ! p(34)  ! turnover rate of labile C (frac/day)
    ! p(35)  ! turnover rate of autotrophic C (frac/day)

    ! zero some values
    lai_out = 0d0; NEE_out = 0d0; GPP_out = 0d0
    ! Set some initial states for the io variables
    infi = 0d0; FLUXES = 0d0; POOLS = 0d0

    ! load some values
    acm_forcings(7) = pars(11)  ! Canopy nitrogen use efficency (gC/gNleaf/day)
    latitude_radians = lat*deg_to_rad
    sin_latitude_radians = sin(latitude_radians)
    cos_latitude_radians = cos(latitude_radians)

    ! parameters from file
    decomposition_rate                = pars(1)  ! decomposition rate (day)
    frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
    DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
    DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
    turnover_rate_foliage             = pars(5)  ! turnover_rate of foliage (day)
    turnover_rate_stem                = pars(6)  ! turnover rate of stem (day)
    RDRSHMAX                          = pars(7)  ! maximum rate of foliar turnover due to self shading (day)
    VDh                               = pars(8)  ! effective vernalisation days when plants are 50 % vernalised
    mineralisation_rate_litter        = pars(9)  ! mineralisation rate litter (day)
    mineralisation_rate_soilOrgMatter = pars(10)  ! mineralisation rate som (day)
    sow_day                           = nint(mod(pars(12), 365.25d0))  ! sow day (doy)
    resp_cost_labile_trans            = pars(13)  ! labile lost to respiration per gC labile to GPP
    PHUem                             = pars(14)  ! phenological heat units required for emergence
    harvest_day                       = nint(mod(pars(15), 365.25d0))  ! nint(mod(pars(15), 365.25))  ! harvest day (doy)
    plough_day                        = nint(mod(pars(12)-2d0, 365.25d0))  ! nint(mod(pars(16), 365.25))  ! plough day (doy)
    LCA                               = pars(17)  ! leaf mass area (gC.m-2)
    tmin                              = pars(26)-273.15d0  ! min temperature for development
    tmax                              = pars(27)-273.15d0  ! max temperature for development
    topt                              = pars(28)-273.15d0  ! optimum temperature for development
    tmin_v                            = pars(29)-273.15d0  ! min temperature for vernalisation
    tmax_v                            = pars(30)-273.15d0  ! max temperature for vernalisation
    topt_v                            = pars(31)-273.15d0  ! optimim temperature for vernalisation
    PHCR                              = pars(32)  ! critical value of photoperiod for development
    PHSC                              = pars(33)  ! photoperiod sensitivity
    turnover_rate_labile              = pars(34)  ! turnover rate labile C (day)
    turnover_rate_resp_auto           = pars(35)  ! turnover rate of autotrophic carbon for respiration (day)

    ! load stocks in first time step
    stock_labile                      = pars(18)  ! labile C
    stock_foliage                     = pars(19)  ! foliar C
    stock_roots                       = pars(20)  ! root C
    stock_stem                        = pars(21)  ! stem/wood C
    stock_litter                      = pars(22)  ! litter C
    stock_soilOrgMatter               = pars(23)  ! som C
    stock_resp_auto                   = pars(24)  ! autotrophic resp pool
    stock_storage_organ               = pars(25)  ! storage organ (i.e. desired crop)
    stock_dead_foliage                = 0d0      ! dead but still standing foliage

    ! assigning initial conditions
    POOLS(1, 1) = stock_labile
    POOLS(1, 2) = stock_foliage
    POOLS(1, 3) = stock_roots
    POOLS(1, 4) = stock_stem
    POOLS(1, 5) = stock_litter
    POOLS(1, 6) = stock_soilOrgMatter
    POOLS(1, 7) = stock_resp_auto
    POOLS(1, 8) = stock_storage_organ
    POOLS(1, 9) = stock_dead_foliage

    ! logical switches
    vernal_calcs    = .true.
    ploughed        = .false.
    sown            = .false.
    use_seed_labile = .true.
    emerged         = .false.

    ! pair incoming variables to local module levels
    ! finally set some initial conditions
    avtemp = 0d0
    yield = 0d0
    DS = -1d0
    DR = 0d0
    fV = 0d0; fT = 0d0; fP = 0d0
    alloc_to_storage_organ_old = 0d0
    PHU = 0d0; VD = 0d0
    BM_EX = 0d0; HI = 0d0
    stock_dead_foliage = 0d0
    alloc_to_labile = 0d0
    stmob = 0
    max_raso = 0d0
    raso = 0d0
    RDRDV = 0d0
    max_raso_old = 0d0
    raso_old  = 0d0
    fol_frac_intpol = 0d0
    stem_frac_intpol = 0d0
    root_frac_intpol = 0d0

    ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
    if (.not.allocated(deltat_1)) then
        allocate(deltat_1(nodays), cica_time(nodays), DS_time(nodays), &
                 avN_time(nodays))
        deltat_1 = deltat**(-1d0)
    else 
        cica_time = 0d0; DS_time = 0d0; avN_time = 0d0
    endif

    ! load some needed module level values
    lai = POOLS(1, 2)/pars(17)
    doy = ceiling(met(6, 1)-(deltat(1)*0.5d0))   ! Day of year

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! Determine day of year for time step mean
      doy = ceiling(met(6, n)-(deltat(n)*0.5d0))   ! Day of year

      ! states needed for module variables
      lai_out(n) = POOLS(n, 2)/LCA
      lai = lai_out(n)  ! leaf area index (m2/m2)

      ! Calculate solar declination for the current time step
      declination = calculate_declination(doy)
      ! calculate daylength in hours and seconds
      call calculate_daylength
      ! extract timing related values
      days_per_step = deltat(n); days_per_step_1 = deltat_1(n)

      ! Note that soil mass balance will be calculated after phenology
      ! adjustments

      ! GPP and direct allocation of GPP can only occur 
      ! if sufficient LAI to prevent numerical error
      if (lai > 1d-5) then 

          ! DS < 0.15 corresponds to the growth stage at beginning of the UK recommended period of      
          ! N fertiliser application for winter wheat (Zodocks growth stage 20) - the early tillering stage (typically mid-march to April)
          if (DS < 0.469d0) then                          
              avN = pars(36)                     
          else  ! if (DS >= 0.469d0 .and. DS <= 1.293d0) then
              ! NOTE: The slope_n parameter can be included in the MDF optimisation. 
              !       The value for this parameter has also been observed to be around-0.02.
              ! NOTE: Modified to only allow the dilution equation to dilute not enrich N content.
              !       This is to attempt to get around the N-dilultion model increasing N content 
              !       during senescence which is unrealistic.               
              avN = max(0.1d0, min(avN, (pars(37)*(POOLS(n, 2)+POOLS(n, 9))) + pars(36)))
!          else
!             ! Set LNA to 0.1 after anthesis (Zodocks growth stage 75)  
!             avN = 0.1d0                              
          end if

          ! load next met/lai values for ACM
          acm_forcings(1) = lai      ! LAI
          acm_forcings(2) = met(3, n)  ! maximum temperature (oC)
          acm_forcings(3) = met(2, n)  ! minimum temperature (oC)
          acm_forcings(4) = avN_time(n)
          acm_forcings(5) = met(5, n)  ! CO2 (ppm)
          acm_forcings(6) = met(4, n)  ! incoming short wave radiation (MJ/m2/day)

          ! GPP (gC.m-2.day-1)
          FLUXES(n, 1) = acm(acm_forcings)
          cica_time(n) = ci/met(5, n)
          
      end if
      ! Pass GPP estimate to module variable for use in crop development model
      gpp_acm = FLUXES(n, 1); GPP_out = FLUXES(n, 1)

      ! pass relevant variables into crop module memory
      avtemp = met(14, n)  ! meant

      ! calculate weighted air temperature value based on daily minimum, maximum
      ! and means. This minimises the error introduced when scaling between
      ! daily and sub-daily timesteps
      airt_weighting(1) = abs(met(3, n)-avtemp) / (met(3, n)-met(2, n))*0.5d0  ! maximum temperature weighting
      airt_weighting(2) = 0.5d0                                            ! mean temperature
      airt_weighting(3) = abs(met(2, n)-avtemp) / (met(3, n)-met(2, n))*0.5d0  ! minimum temperature weighting

      ! Heterotrophic respiration rate (Q10):  doubles with
      ! 10 degree temperature rise resprate from soil file = 0.0693
      resp_rate = 0d0
      resp_rate = resp_rate + ((0.5d0*exp( resp_rate_temp_coeff*met(3, n) )) * airt_weighting(1))
      resp_rate = resp_rate + ((0.5d0*exp( resp_rate_temp_coeff*avtemp   )) * airt_weighting(2))
      resp_rate = resp_rate + ((0.5d0*exp( resp_rate_temp_coeff*met(2, n) )) * airt_weighting(3))
      !resp_rate = 0.5*exp( resp_rate_temp_coeff*avtemp )

      ! reallocate day of year to the end of the time step for use in
      ! crop development model
      doy = met(6, n)
      ! determine development stage (DS)
      ! Note that DS must be updated here and not after management_dates. Otherwise, 
      ! there will be problems using DS_time for calculating yield:GPP calculations in MODEL_LIKELIHOOD.f90
      call development_stage(deltat(n)); DS_time(n) = DS 
      ! determine the carbon partitioning based on development stage
      call carbon_alloc_fractions(DS_shoot, DS_root, fol_frac, stem_frac, root_frac)
      ! begin carbon allocation for crops
      call calc_pools_crops(DS_LRRT, LRRT)
      ! conduct management updates at the end of the day
      call management_dates(stock_seed_labile, deltat(n))

      ! calculate the NEE (gC.m-2.d-1)
      NEE_out(n) = nee_dalec

      ! GPP (gC.m-2.d-1)
      !FLUXES(n, 1) = GPP_out(n)  ! Assigned above
      ! temprate (i.e. temperature modified rate of metabolic activity)
      FLUXES(n, 2) = resp_rate
      ! autotrophic respiration (gC.m-2.d-1)
      FLUXES(n, 3) = resp_auto+resp_cost_labile_to_npp + &
                    resp_cost_npp_to_labile+resp_cost_foliage_to_labile + &
                    resp_cost_stem_to_labile
      ! leaf production rate (gC.m-2.d-1)
      FLUXES(n, 4) = alloc_to_foliage
      ! labile production (gC.m-2.d-1)
      FLUXES(n, 5) = alloc_to_labile+foliage_to_labile+stem_to_labile
      ! root production (gC.m-2.d-1)
      FLUXES(n, 6) = alloc_to_roots
      ! stem production (gC.m-2.d-1)
      FLUXES(n, 7) = alloc_to_stem
      ! labile from NPP (gC.m-2.d-1)
      FLUXES(n, 8) = alloc_from_labile
      ! alloc to storage organ (gC.m-2.d-1)
      FLUXES(n, 9) = alloc_to_storage_organ
      ! total leaf litter production (gC.m-2.d-1)
      FLUXES(n, 10) = litterfall_foliage
      ! total stem litter production (gC.m-2.d-1)
      FLUXES(n, 11) = litterfall_stem
      ! total root litter production (gC.m-2.d-1)
      FLUXES(n, 12) = litterfall_roots
      ! respiration heterotrophic litter (gC.m-2.d-1)
      FLUXES(n, 13) = resp_h_litter
      ! respiration heterotrophic som (gC.m-2.d-1)
      FLUXES(n, 14) = resp_h_soilOrgMatter
      ! litter to som (gC.m-2.d-1)
      FLUXES(n, 15) = decomposition
      ! alloc to autotrophic pool (gC.m-2.d-1)
      FLUXES(n, 16) = alloc_to_resp_auto
      ! harvest yield (gC.m-2.d-1)
      FLUXES(n, 21) = yield
      ! C extracted in addition to yield due to harvest activity (gC.m-2.d-1)
      FLUXES(n, 22) = BM_EX
      ! Respiration from autotrophic allocation (gC.m-2.d-1)
      FLUXES(n, 23) = resp_auto
      ! Respiration from labile to foliage translocation (gC.m-2.d-1)
      FLUXES(n, 24) = resp_cost_labile_to_npp
      ! Respiration from foliage to litter translocation (gC.m-2.d-1)
      FLUXES(n, 25) = resp_cost_npp_to_labile
      ! Respiration from foliage remobilisation (gC.m-2.d-1)
      FLUXES(n, 26) = resp_cost_foliage_to_labile
      ! Respiration from stem remobilisation (gC.m-2.d-1)
      FLUXES(n, 27) = resp_cost_stem_to_labile
      ! Foliage extracted from harvest (gC.m-2.d-1)
      FLUXES(n, 28) = HARVESTextracted_foliage
      ! Stem extracted from harvest (gC.m-2.d-1)
      FLUXES(n, 29) = HARVESTextracted_stem
      ! Dead still standing foliage extracted from harvest (gC.m-2.d-1)
      FLUXES(n, 30) = HARVESTextracted_dead_foliage
      ! Labile extracted from harvest (gC.m-2.d-1)
      FLUXES(n, 31) = HARVESTextracted_labile
      ! Foliage added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n, 32) = HARVESTlitter_foliage
      ! Stem added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n, 33) = HARVESTlitter_stem
      ! Dead still standing foliage added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n, 34) = HARVESTlitter_dead_foliage
      ! Autotrophic pool added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n, 35) = HARVESTlitter_resp_auto
      ! Labile added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n, 36) = HARVESTlitter_labile
      ! Roots pool added to litter due to plough (gC.m-2.d-1)
      FLUXES(n, 37) = PLOUGHlitter_roots

      ! labile pool
      POOLS(n+1, 1) = stock_labile
      ! foliar pool
      POOLS(n+1, 2) = stock_foliage
      ! root pool
      POOLS(n+1, 3) = stock_roots
      ! wood pool
      POOLS(n+1, 4) = stock_stem
      ! litter pool
      POOLS(n+1, 5) = stock_litter
      ! som pool
      POOLS(n+1, 6) = stock_soilOrgMatter
      ! autotrophic pool
      POOLS(n+1, 7) = stock_resp_auto
      ! storage organ pool
      POOLS(n+1, 8) = stock_storage_organ
      ! dead but still standing foliage
      POOLS(n+1, 9) = stock_dead_foliage

      do nxp = 1, nopools
         if (POOLS(n+1, nxp) /= POOLS(n+1, nxp) .or. POOLS(n+1, nxp) < 0d0) then
             print*,"step",n, "POOL",nxp
             print*,"met",met(:,n)
             print*,"POOLS",POOLS(n, :)
             print*,"FLUXES",FLUXES(n, :)
             print*,"POOLS+1",POOLS(n+1, :)
             print*,"labile, foliage",stock_labile, stock_foliage
             print*,"stem, roots",stock_stem, stock_roots
             print*,"litter, som",stock_litter, stock_soilOrgMatter
             print*,"storageOrgan, auto",stock_storage_organ, stock_resp_auto
             print*,"gpp, nee",gpp_acm, nee_dalec
             print*,"Ra, Rh_som, Rh_lit",resp_auto, resp_h_soilOrgMatter, resp_h_litter
             print*,"pars",pars
             print*,"DR",DR
             print*,"DR stuff",fT, fV, fP
             print*,"leaf","stem","rootlitter",litterfall_foliage, litterfall_stem, litterfall_roots
             print*,"daylength",dayl_hours, "VD",VD, "VDh",VDh
             print*,"avtemp",avtemp
             print*,"sown",sown, "emerged",emerged
             print*,"root_frac_intpol",root_frac_intpol
             print*,"npp_shoot",npp_shoot, "npp",npp
             print*,"RDR",RDR
             !stop
         endif
      enddo

      do nxp = 1, nofluxes

            if (FLUXES(n, nxp) /= FLUXES(n, nxp) .or. FLUXES(n, nxp) < 0d0) then
                 print*,"Special: step",n, "FLUXES",nxp
                 print*,"met",met(:,n)
                 print*,"POOLS",POOLS(n, :)
                 print*,"FLUXES",FLUXES(n, :)
                 print*,"POOLS+1",POOLS(n+1, :)
                 print*,"labile, foliage",stock_labile, stock_foliage
                 print*,"stem, roots",stock_stem, stock_roots
                 print*,"litter, som",stock_litter, stock_soilOrgMatter
                 print*,"storageOrgan, auto",stock_storage_organ, stock_resp_auto
                 print*,"gpp, nee",gpp_acm, nee_dalec
                 print*,"Ra, Rh_som, Rh_lit",resp_auto, resp_h_soilOrgMatter, resp_h_litter
                 print*,"pars",pars
                 print*,"DR",DR
                 print*,"DR stuff",fT, fV, fP
                 print*,"leaf","stem","rootlitter",litterfall_foliage, litterfall_stem, litterfall_roots
                 print*,"daylength",dayl_hours, "VD",VD, "VDh",VDh
                 print*,"avtemp",avtemp
                 print*,"sown",sown, "emerged",emerged
                 print*,"root_frac_intpol",root_frac_intpol
                 print*,"npp_shoot",npp_shoot, "npp",npp
                 print*,"RDR",RDR
                 !stop
            end if
      enddo

      if (stock_labile < 0d0 .or. stock_foliage < 0d0 .or. stock_stem < 0d0 .or. &
          stock_roots < 0d0 .or. stock_litter < 0d0 .or. stock_soilOrgMatter < 0d0 .or. &
          stock_storage_organ < 0d0 .or. stock_resp_auto < 0d0 .or. &
          stock_labile /= stock_labile .or. stock_foliage /= stock_foliage .or. &
          stock_stem /= stock_stem .or. &
          stock_roots /= stock_roots .or. stock_litter /= stock_litter .or. &
          stock_soilOrgMatter /= stock_soilOrgMatter .or. &
          stock_storage_organ /= stock_storage_organ .or. &
          stock_resp_auto /= stock_resp_auto .or.  &
          gpp_acm < 0d0 .or. gpp_acm /= gpp_acm .or. resp_rate < 0d0 .or. &
          resp_rate /= resp_rate .or. decomposition < 0d0 .or. alloc_from_labile < 0d0 .or. &
          resp_cost_labile_to_npp < 0d0 .or. alloc_to_foliage < 0d0 .or. &
          alloc_to_stem < 0d0 .or. alloc_to_roots < 0d0 .or. &
          alloc_from_labile < 0d0 .or. resp_cost_labile_to_npp < 0d0) then
          print*,"stocks less than zero or NaN", n
          print*,stock_labile, stock_foliage
          print*,stock_stem, stock_roots
          print*,stock_litter, stock_soilOrgMatter
          print*,stock_storage_organ, stock_resp_auto
          print*,gpp_acm, nee_dalec
          print*,resp_auto, resp_h_soilOrgMatter, resp_h_litter
          print*,"pars",pars
          print*,"fluxes",fluxes(n, :)
          print*,"DR",DR
          print*,"DR stuff",fT, fV, fP
          print*,"leaf","stem","rootlitter",litterfall_foliage, litterfall_stem, litterfall_roots
          print*,"daylength",dayl_hours, "VD",VD, "VDh",VDh
          print*,"avtemp",avtemp
          print*,"sown",sown, "emerged",emerged
          print*,"root_frac_intpol",root_frac_intpol
          print*,"npp_shoot",npp_shoot, "npp",npp
          print*,"RDR",RDR
          !stop
      endif

    end do  ! no days loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm(acm_forcings)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in):: acm_forcings(7)  ! acm input requirements

    ! declare local variables
    double precision, parameter:: dayl_coef = 0.0156935d0, &
                              co2_comp_point = 4.22273d0,   &
                                co2_half_sat = 208.868d0,   &
                                  dayl_const = 0.0453194d0, &
                         hydraulic_temp_coef = 0.37836d0,   &
                                    lai_coef = 7.19298d0,   &
                               temp_exponent = 0.011136d0,  &
                                   lai_const = 2.1001d0,    &
                          hydraulic_exponent = 0.789798d0,  &
                                     deltaWP = -2d0,        & ! leafWP-soilWP
                                        Rtot = 1d0            ! Total hydraulic resistance
    double precision:: gc, pn, pd, pp, qq, e0, cps, nit, NUE &
                       ,trange, mint, maxt, radiation, co2, lai, doy

    ! initial values
    gc = 0d0; pp = 0d0; qq = 0d0; ci = 0d0; e0 = 0d0; cps = 0d0 

    ! load driver values to correct local vars
    lai  = acm_forcings(1)  ! leaf area index m2/m2
    maxt = acm_forcings(2)  ! mean of daily maximum temperature C
    mint = acm_forcings(3)  ! mean of daily minimum temperature C
    nit  = acm_forcings(4)  ! mean foliar nitrogen gN/m2leaf
    co2  = acm_forcings(5)  ! ppm
    radiation = acm_forcings(6)  ! MJ/m2/day
    NUE = acm_forcings(7)  ! Nitrogen use efficiency

    ! determine temperature range 
    trange = 0.5*(maxt-mint)
    ! daily canopy conductance, of CO2 or H2O? 
    gc = abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp = pn/gc 
    qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci = 0.5*(co2+qq-pp+sqrt(((co2+qq-pp)*(co2+qq-pp))-4d0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0 = lai_coef*(lai*lai)/((lai*lai)+lai_const)

    ! calculate CO2 limited rate of photosynthesis
    pd = gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps = e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm = cps*(dayl_coef*dayl_hours+dayl_const)

    return

  end function acm
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_cosine_solar_zenith_angle

    implicit none

    ! Calculate some common variables needed for the Sellers (1985)
    ! parameteric appoximation

    ! Declare local parameters

    ! Solar angle for the hour of the day, assumed in daily model to be 10 am (or 2 pm)
    ! therefore, 2 = 14-12, where 14 = timestep of day and 12 = half the number
    ! hour-angle = ( hr-12 ) * 15. * ( pi/180.d0 )
    ! of time steps per day.
    double precision, parameter:: cos_hour_angle = cos(0.5235988d0)

    ! Declare local variables

    ! Now calculate the solar-zenith-angle, as per..
    !  sin(latitude)sin(declination) + cos(latitude)cos(declination)cos(hour_angle)
    ! where hour-angle = ( hr-12 ) * 15. * ( pi/180.d0 )
    ! and from that the solar flux.
    ! eqn 2.15 in Hartmann's Global Physical Climatology...
    ! Minimum allowed value constraint from SIB3 implementation of Sellers (1985)
    cosine_solar_zenith_angle = max(0.01747d0, sin_latitude_radians*sin(declination) + &
                                               cos_latitude_radians*cos(declination) * &
                                               cos_hour_angle)
    
  end subroutine calculate_cosine_solar_zenith_angle
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_daylength

    ! Subroutine uses day of year and latitude (-90/90 degrees) as inputs, 
    ! combined with trigonomic functions to calculate day length in hours and seconds

    implicit none

    ! local variables
    double precision:: dec, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = sin_latitude_radians*sin( declination )
    cosld = cos_latitude_radians*cos( declination )
    aob = max(-1d0, min(1d0, sinld/cosld))

    ! estimate day length in hours and seconds and upload to module variables
    dayl_hours = 12d0 * ( 1d0+2d0*asin( aob ) * pi_1 )
    !dayl_seconds = dayl_hours*seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !

  !
  !-----------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  subroutine calc_pools_crops(DS_LRRT, LRRT)

    ! Allocated GPP to NPP and various carbon pools. Based  !
    ! this on physiological responses to temperature        !
    ! vernalisation, and photoperiod. Note: this code has   !
    ! been reformulated to a daily time step as default     !
    ! rather SPAcrop which operated on a hourly time step.  !

    implicit none

    ! arguments
    double precision, dimension(:), intent(inout):: DS_LRRT, & ! Development stage corresponding to loss rate of roots
                                                        LRRT    ! loss rate of roots array

    ! local variables
    double precision:: decomp_efficency

    ! turnover rate of fine roots is now equal to the
    ! loss rate of roots as a function of development stage
    ! (Penning de Vries, 1989).
    turnover_rate_roots = interpolate( DS, DS_LRRT, LRRT, 5 )

    ! if sown turn on labile/seed turnover for growth
    if ( sown ) then
        ! turnover on
        turnover_labile_switch = 1
    else
        ! turnover off
        turnover_labile_switch = 0
    endif

    ! Initialise/reset
    yield = 0d0; BM_EX = 0d0
    HARVESTextracted_foliage = 0d0
    HARVESTextracted_stem = 0d0
    HARVESTextracted_dead_foliage = 0d0
    HARVESTextracted_labile = 0d0
    HARVESTlitter_labile = 0d0
    HARVESTlitter_foliage = 0d0
    HARVESTlitter_stem = 0d0
    HARVESTlitter_dead_foliage = 0d0
    HARVESTlitter_resp_auto = 0d0
    PLOUGHlitter_roots = 0d0

    ! Determine total allocation of labile C to NPP or Ra (gC.m-2.d-1)
    alloc_from_labile = turnover_rate_labile*resp_rate*dble(turnover_labile_switch)
    alloc_from_labile = stock_labile * (1d0-(1d0-alloc_from_labile) ** days_per_step) * days_per_step_1
    ! Determine the respiratory cost of C transfer from labile pool to NPP (gC.m-2.d-1)
    resp_cost_labile_to_npp = alloc_from_labile*resp_cost_labile_trans
    ! Determine the remaining allocation of labile C to NPP (gC.m-2.d-1)
    alloc_from_labile = alloc_from_labile-resp_cost_labile_to_npp

    ! Estimate allocation (gC/m2/day) of GPP to autotrophic respiration pool
    alloc_to_resp_auto = gpp_acm*frac_GPP_resp_auto
    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation from labile
    npp = gpp_acm+alloc_from_labile-alloc_to_resp_auto

    ! Dertermine partitioning of NPP to biomass pools
    alloc_to_roots    = root_frac_intpol*npp
    ! Calculate how much NPP is left after allocation to roots
    npp_shoot         = npp-alloc_to_roots
    ! Partition the remaining NPP to foliage and stem (gC/m2/day)
    alloc_to_foliage  = fol_frac_intpol  * npp_shoot
    alloc_to_stem     = stem_frac_intpol*npp_shoot
    ! Remaining NPP is then allocated to storage organ
    alloc_to_storage_organ = max(0d0, npp_shoot-alloc_to_foliage-alloc_to_stem)

    ! Assuming allocatio to storage organ is > 0 ensure flux is limited by
    ! maximum growth rate potential, i.e. growth potential increases with size
    ! existing yield
    alloc_to_labile = 0d0; resp_cost_npp_to_labile = 0d0
    if ( alloc_to_storage_organ > 0d0 ) then
        gso_max  = ( stock_storage_organ+0.5d0 ) * rel_gso_max
        ! Restrict allocation to storage organ based on available C and potential
        ! Where does any excess go?
        alloc_to_storage_organ = min( alloc_to_storage_organ, gso_max )
        if ( sown .and. emerged) then
            ! Assuming we shown and emerged, the remainder is now allocated to labile, 
            ! less the cost of transfer
            alloc_to_labile = ( npp_shoot-alloc_to_foliage-alloc_to_stem-alloc_to_storage_organ )
            ! Having worked out the total available C, update based on cost of respiratory transfer
            resp_cost_npp_to_labile =  alloc_to_labile*resp_cost_labile_trans
            alloc_to_labile = alloc_to_labile-resp_cost_npp_to_labile
        endif
    endif

    ! set switches to (de)activate leaf, root and stem remobliization
    raso_old = raso
    ! Running Average of growth rate of Storage Organ..
    raso = ( alloc_to_storage_organ+alloc_to_storage_organ_old ) * 0.5d0
    max_raso_old = max_raso
    max_raso = max( raso, max_raso_old )
    ! Stem remobilisation triggered once running average of storage organ growth declines
    ! Second part prevents premature remobilisation
    if ( ( raso < raso_old ) .and. &
         ( alloc_to_storage_organ > ( alloc_to_storage_organ_old+0.5d0 ) ) ) then
        stmob = 1
    else
        stmob = 0
    endif

    ! Code for calculating relative death rate of leaves (RDR) as a
    !  function of shading (RDRSH) or developmental stage (RDRT).

    ! GT 0 if LAI GT 4; 0. < RDRSH < RDRSHMAX (usually ~0.03)
    RDRSH = min( RDRSHMAX, max( 0d0, RDRSHMAX * ( lai-LAICR ) / LAICR ) )
    if ( DS < 1d0 ) then
       RDRDV = 0d0
    else
       ! RDRDV dependant on DR and DS, values range typically between 0.02 <
       ! RDRDV < 0.25
!print*,"!! What RDRDV to use? !!"
!!$      RDRDV = DR /( max( 0.1, 2. - DS ) )
!!$      RDRDV = RDRDV/24. ! to get hourly senescence rate
       ! TLS: think the source of this equation needs to be found
       RDRDV = turnover_rate_foliage*min(1d0, ( 1d0 / ( ( max( 2d0-DS, 0.1d0 ) ) * 8d0 ) ) ** 2)
    ENDIF

    ! Relative leaf death rate is the maximum value of the arguments RDRSH and
    ! RDRDV
    RDR = max( RDRSH, RDRDV )

    ! remobilization of foliar C and allocation to dead leaves pool (gC.m-2.d-1)
    litterfall_foliage = stock_foliage * (1d0-(1d0-RDR) ** days_per_step) * days_per_step_1
    ! remobilization of stem C (gC.m-2.d-1)
    litterfall_stem    = stock_stem    * (1d0-(1d0-(DR*turnover_rate_stem*dble(stmob))) ** days_per_step) * days_per_step_1
    ! remobilization of fine roots C (gC.m-2.d-1)
    litterfall_roots   = stock_roots   * (1d0-(1d0-turnover_rate_roots) ** days_per_step) * days_per_step_1

    ! remobilized C to NPP (from both leaves and stems) (gC.m-2.d-1)
    ! Assume that half the foliage litter C is available for remobilisation
    foliage_to_labile   = ( litterfall_foliage*0.5d0) * ( 1d0-resp_cost_labile_trans )
    stem_to_labile      = ( litterfall_stem ) * ( 1d0-resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates) (gC.m-2.d-1)
    resp_cost_foliage_to_labile   = ( litterfall_foliage*0.5d0) * resp_cost_labile_trans
    resp_cost_stem_to_labile      = litterfall_stem*resp_cost_labile_trans

    ! litter decomposition to soil organic matter
    decomposition = stock_litter &
                  * (1d0-(1d0-(decomposition_rate*resp_rate)) ** days_per_step) &
                  * days_per_step_1

    ! heterotrophic respiration component 1: mineralisation of litter C pool (gC.m-2.d-1)
    resp_h_litter = stock_litter &
                  * (1d0-(1d0-(mineralisation_rate_litter*resp_rate)) ** days_per_step) &
                  * days_per_step_1
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool (gC.m-2.d-1)
    resp_h_soilOrgMatter = stock_soilOrgMatter &
                         * (1d0-(1d0-(mineralisation_rate_soilOrgMatter*resp_rate)) ** days_per_step) &
                         * days_per_step_1
                  
    ! Autotrophic respiration allocated pool is a special case, where the total
    ! available for the time step is accumulated first then losses are determined.
    ! Add photosynthate allocated to autotrophic respiration (gC.m-2.d-1)
    stock_resp_auto = stock_resp_auto + (alloc_to_resp_auto*days_per_step)
    ! Determine autotrophic costs as loss term. This should probably have a temperature sensitivity
    ! or be removed to work as applied in other DALEC models? (gC.m-2.d-1)
    ! NOTE: that autotrophic respiration from other sources added below
    resp_auto = stock_resp_auto * (1d0-(1d0-turnover_rate_resp_auto) ** days_per_step) * days_per_step_1
    ! Update stock
    stock_resp_auto = max(0d0, stock_resp_auto - (resp_auto*days_per_step))

    ! nee (gC.m-2.d-1)
    nee_dalec = (resp_auto+resp_cost_labile_to_npp+resp_cost_npp_to_labile + &
                 resp_cost_foliage_to_labile+resp_cost_stem_to_labile + &
                 resp_h_litter+resp_h_soilOrgMatter) - gpp_acm

    ! Recalculate Physical Carbon Pools...
    ! TLS: Note that while max(0d0, ... conditions have been applied to all pools, I hypothese that the critical
    ! restriction is applied to stock_resp_auto
    stock_foliage       = max(0d0, stock_foliage       + ((alloc_to_foliage-litterfall_foliage) * days_per_step))
    stock_stem          = max(0d0, stock_stem          + ((alloc_to_stem-litterfall_stem) * days_per_step))
    stock_storage_organ = max(0d0, stock_storage_organ + ((alloc_to_storage_organ) * days_per_step))
    stock_roots         = max(0d0, stock_roots         + ((alloc_to_roots   - litterfall_roots) * days_per_step))
    stock_litter        = max(0d0, stock_litter        + ((litterfall_roots-resp_h_litter-decomposition) * days_per_step))
    stock_soilOrgMatter = max(0d0, stock_soilOrgMatter + ((decomposition    - resp_h_soilOrgMatter) * days_per_step))
    stock_dead_foliage  = max(0d0, stock_dead_foliage  + ((litterfall_foliage*0.5d0) * days_per_step))
    stock_labile        = max(0d0, stock_labile        + ((alloc_to_labile  - alloc_from_labile - &
                                                           resp_cost_labile_to_npp+foliage_to_labile + &
                                                           stem_to_labile) * days_per_step) )

    ! When GPP is higher than seed C content, remaining seed carbon enters litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    ! Moved to the end of the mass balance, to maintain balance with fluxes which
    ! are already concurrently occuring. In this codes previous position further
    ! up the code, resulted in C being both dumped to litter and the same carbon
    ! being allocated to respiration and NPP.
    if ( ( gpp_acm > alloc_from_labile ) .and. ( use_seed_labile ) ) then
        stock_litter = stock_litter+stock_labile
        stock_labile = 0d0
        use_seed_labile = .false.
    endif

  end subroutine calc_pools_crops
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine carbon_alloc_fractions(DS_shoot, DS_root, fol_frac, stem_frac, root_frac)

    ! Determines carbon allocation fractions as a function !
    ! of developmental stage (DS).  Allocation fractions   !
    ! are from tables published in Penning de Vries (1989) !

    implicit none

    double precision, dimension(:), intent(inout) ::   DS_shoot, & !
                                                        DS_root, & !
                                                       fol_frac, & !
                                                      stem_frac, & !
                                                      root_frac    !

    ! local variables..
    double precision, dimension(:), allocatable:: frac_shoot, frac_root

    if ( sown ) then  ! after sowing

       ! loop over three crop "organs": 1) foliage 2) stems 3) root
       ! not necessary for storage organs, as all remaining C is allocated to
       ! these

       ! use different input for foliage and stem fractions, as they are
       ! relative to the total shoot (or aboveground) allocation, 
       ! root is relative to total plant (above-and belowground) allocation.

       ! leaf development stages and corresponding fractions..

       ! interpolate between PdV allocation values with reference to
       ! developmental stage (DS)..
       fol_frac_intpol = min(1d0, max(0d0, interpolate( DS, DS_shoot, fol_frac, size(DS_shoot) )))
       ! stem DS and fracs..
       stem_frac_intpol = min(1d0, max(0d0, interpolate( DS, DS_shoot, stem_frac, size(DS_shoot) )))
       ! root DS and fracs..
       root_frac_intpol = min(1d0, max(0d0, interpolate( DS, DS_root, root_frac, size(DS_root) )))

    endif  ! after crop has been sown

  end subroutine carbon_alloc_fractions
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine development_stage(days_in_step)

    ! Based on modified Wang & Engel model (Streck et al., 2003), !
    ! but with only 2 sub-phases, vegetative and reproductive     !
    ! (i.e. only two different DRmax).   O. Sus, May 2010.        !

    implicit none

    ! agruments
    double precision:: days_in_step

    ! local variables..
    double precision ::  doptmin, & ! Difference between optimum and minimum temperature
                         dmaxmin, & ! Difference between maximum and minimum temperature
                          dttmin, & ! Difference between daiy average and minimum temperatures
                       doptmin_v, & ! Difference between optimum and minimum vernalization temperatures
                       dmaxmin_v, & ! Difference between maximum and minimum vernalization temperatures
                       dttmin_v     ! Difference between daily average and minimum vernalization temperatures

    doptmin   = topt   - tmin   ! difference between optimal and minimum cardinal temperatures
    dmaxmin   = tmax   - tmin   ! difference between maximum and minimum cardinal temperatures
    dttmin    = avtemp-tmin   ! difference between daily average and minimum cardinal temperatures
    doptmin_v = topt_v-tmin_v  ! same as above, 
    dmaxmin_v = tmax_v-tmin_v !       but for vernalization
    dttmin_v  = avtemp-tmin_v  ! cardinal temperatures

    ! Calculation of developmental function values: vernalization (fV), 
    ! temperature (fT) and
    ! photoperiod (fP) these values are multiplicative factors of DRmax (maximum
    ! developmental
    ! rate), each ranging between 0 (no development) and 1 (unrestricted
    ! development).

    ! Summation of vernalization days (VD), not before sowing and only if
    ! average temperature is within min and max cardinal temperatures..
    if ( ( avtemp > tmin_v ) .and. ( avtemp < tmax_v ) .and. sown ) then
        fV = vernalization( doptmin_v, dmaxmin_v, dttmin_v, days_in_step )
    endif

    ! Only calculate temperature coefficient if avtemp lies within (tmin, tmax)
    ! range. NOTE: (doptmin+1d0) < dmaxmin added to allow for EDC search period when "not
    ! allowed" parameter sets will be tried anyway
    if ( avtemp > tmin .and. avtemp < tmax .and. (doptmin+1d0) < dmaxmin ) then
        fT = temperature_impact( doptmin, dmaxmin, dttmin )
    else
        fT = 0d0
    endif

    ! calculation of photoperiod coefficient
    fP = photoperiod_impact( PHCR, PHSC )

    if ( emerged .and. ( DS < 2d0 ) ) then   ! sum up daily DR values between emergence and maturity (DS = 2)

       if ( DS < 1d0 ) then  ! in the vegetative phase (before flowering):

          DR = DR_pre*fT*fP   ! DR is affected by temperature, photoperiod...
          if ( vernal_calcs ) DR = DR*fV ! ...and vernalization (for winter cereals)
          DS = DS + (DR*days_in_step)    ! developmental stage (DS), calculated as the sum of daily developmental rates

       else    ! in the reproductive phase (after flowering):

          DR = DR_post*fT   ! DR is affected only by temperature
          DS = DS + (DR*days_in_step)

       endif  ! vegetative or reproductive phase

    endif  ! emerged or not

  end subroutine development_stage
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine management_dates (stock_seed_labile, days_in_step)

    ! This routine should be called at the end of each day of a crops  !
    ! simulation.  It checks whether we should plough/sow/harvest, and !
    ! during the growing establishes when the crop will emerge after   !
    ! sowing, based on heat accumulation (Phenological Heat Units).    !

    implicit none

    ! arguments
    double precision, intent(in):: stock_seed_labile, days_in_step
    ! local variables
    double precision:: tmp
    logical:: plough_sanity, sow_sanity, harvest_sanity

    ! reset
    plough_sanity = .false. ; sow_sanity = .false. ; harvest_sanity = .false.
    !plough_sanity = .true. ; sow_sanity = .true. ; harvest_sanity = .true.

    ! spring crop
    if (sow_day < harvest_day .and. nint(doy) < harvest_day) sow_sanity = .true.
    if (plough_day < harvest_day .and. nint(doy) < harvest_day) plough_sanity = .true.
    if (harvest_day > sow_day) harvest_sanity = .true.
    ! winter crops
    if (sow_day > harvest_day) sow_sanity = .true.
    if (plough_day > harvest_day) plough_sanity = .true.
    if (harvest_day < plough_day .and. nint(doy) < plough_day) harvest_sanity = .true.

    if ( .not. sown ) then

        ! fresh field...

        if ( plough_sanity .and. .not.ploughed .and. nint(doy) >= plough_day ) then
            ! the field needs ploughing..
            call plough

        else if ( sow_sanity .and. nint(doy) >= sow_day ) then

            ! ensure that the field has indeed been ploughed
            if (.not.ploughed) call plough
            ! the field needs sowing..
            sown = .true.

            ! this switch controls whether the labile carbon within the seed is used
            ! for growth
            use_seed_labile = .true.
            stock_labile = stock_seed_labile

        endif  ! plough or sow?

    else

        ! crop in field..

        ! calculate when crop emerges..
        if ( .not. emerged ) then

            ! estimate emergence date based on the accumulated phenological heat
            ! units (PHU)
            ! where PHU is the (positive) heat over tmin..
            tmp = max( avtemp-tmin, 0d0 )*days_in_step
            PHU = PHU+tmp

            ! set the development stage and emergence..
            if ( PHU >= PHUem ) then
                emerged = .true.
                DS = 0d0
            else
                emerged = .false.
                DS = -1d0
            endif

        endif  ! emerged or not

        ! note that in this case harvest day has been fixed relative to the sow
        ! day
        if ( harvest_sanity .and. nint(doy) >= harvest_day) then
            ! the field needs harvesting..
            call harvest
        endif

    endif  ! sown or not

  end subroutine management_dates
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine harvest

    implicit none

    ! Declare some local variables
    double precision:: tmp, Ctotal

    ! shoot biomass..
    Cshoot = stock_foliage+stock_stem+stock_storage_organ+stock_labile
    ! Total biomass not including the storage organ
    Ctotal = stock_foliage+stock_stem+stock_roots

    ! determine harvest index..
    HI = stock_storage_organ/Cshoot

    ! the stuff we actually want from the harvest...
    yield = stock_storage_organ*days_per_step_1

    ! How much of each pool is extracted during harvest
    HARVESTextracted_foliage      = stock_foliage * ( 1d0-lv_res ) * days_per_step_1
    HARVESTextracted_stem         = stock_stem * ( 1d0-st_res ) * days_per_step_1
    HARVESTextracted_dead_foliage = stock_dead_foliage * ( 1d0-lv_res ) * days_per_step_1
    ! How much of each pool remains as litter after harvest
    HARVESTlitter_foliage      = stock_foliage*lv_res*days_per_step_1
    HARVESTlitter_stem         = stock_stem*st_res*days_per_step_1
    HARVESTlitter_dead_foliage = stock_dead_foliage*lv_res*days_per_step_1
    HARVESTlitter_resp_auto    = stock_resp_auto*days_per_step_1

    ! Labile is a special case due to being distributed within various tissues.
    ! NOTE that extracted is calculated then the litter component is estimates as residual.
    ! The time scale adjustment is applied last, rather than inline (as above)
    HARVESTextracted_labile = (stock_labile * (stock_foliage/Ctotal) *  (1d0-lv_res )) & 
                            + (stock_labile * (stock_stem/Ctotal) *  (1d0-st_res )) 
    HARVESTlitter_labile    =  stock_labile-HARVESTextracted_labile
    HARVESTextracted_labile = HARVESTextracted_labile*days_per_step_1
    HARVESTlitter_labile    = HARVESTlitter_labile*days_per_step_1

    ! the biomass that is harvested in addition to the storage-organ..
    BM_EX  = HARVESTextracted_foliage      &
           + HARVESTextracted_stem         &
           + HARVESTextracted_dead_foliage &
           + HARVESTextracted_labile

    ! what's left (will fall to the ground)..
    stock_litter  = stock_litter               &
                  + HARVESTlitter_foliage      &
                  + HARVESTlitter_stem         &
                  + HARVESTlitter_dead_foliage &
                  + HARVESTlitter_resp_auto    & 
                  + HARVESTlitter_labile

    ! empty the plant stocks..
    stock_storage_organ = 0d0
    stock_foliage       = 0d0
    stock_stem          = 0d0
    stock_dead_foliage  = 0d0
    stock_labile        = 0d0
    stock_resp_auto     = 0d0

    ! roots stay in ground and slowly decompose (until/unless the field is
    ! ploughed)

    ! reset logical variables..
    sown    = .false.
    emerged = .false.
    ploughed = .false.
    DS = -1d0; fV = 0d0; fT = 0d0; fP = 0d0; VD = 0d0

  end subroutine harvest
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function photoperiod_impact( PH_crit, PH_sens )

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

    ! arguments..
    double precision, intent(in):: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity

    photoperiod_impact = max(0d0, 1d0-exp ( - PH_Sens * ( dayl_hours-PH_crit ) ))

  end function photoperiod_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine plough

    ! this s/r will reset various carbon pools, to mimic the effect of the
    ! farmer ploughing. !

    implicit none

    ! Track roots transition to litter here, as all others have been reset already
    PLOUGHlitter_roots = stock_roots*days_per_step_1

    ! Move all plant stocks into the litter pool.
    ! ( many of these should already be empty after the harvest, )
    ! ( e.g. the stocks for labile, foliage, storage-organ stem. )
    stock_litter        = stock_litter+stock_dead_foliage &
                          + stock_foliage+stock_labile    &
                           + stock_roots+stock_stem       &
                            + stock_storage_organ

    stock_dead_foliage  = 0d0
    stock_foliage       = 0d0
    stock_labile        = 0d0
    stock_roots         = 0d0
    stock_stem          = 0d0
    stock_storage_organ = 0d0

    ! Reset the development stage & phenological heat units..
    ploughed = .true. ; DS = -1d0; PHU = 0d0
    max_raso = 0d0; raso = 0d0; max_raso_old = 0d0; raso_old = 0d0
    alloc_to_storage_organ_old = 0d0

  end subroutine plough
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
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a, b, t )

    ! The equation is simply...                        !
    !    a*exp( b * ( t-25.0 ) / ( t+273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split.                    !

    implicit none

    ! arguments..
    double precision, intent(in):: a, b, t
    double precision            :: arrhenious

    arrhenious = a*exp( b * (t-25d0) / (t+freeze) )

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_declination(doy)

    implicit none

     ! Declare arguments
     double precision, intent(in):: doy

     ! Declination calculation
     ! NOTE: 0.002739726d0 = 1/365
     !    dec = - asin( sin( 23.45d0*deg_to_rad ) * cos( 2d0*pi * ( doy+10d0 ) / 365d0 ) )
     !    dec = - asin( sin_dayl_deg_to_rad*cos( two_pi * ( doy+10d0 ) / 365d0 ) )
     calculate_declination = - asin( sin_dayl_deg_to_rad*cos( two_pi * ( doy+10d0 ) * 0.002739726d0 ) )

     ! return to user
     return

  end function calculate_declination
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, min_val, optimum, kurtosis, current )

    ! Estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis. Minimum is assumed to be at infinity
    ! (or near enough)

    implicit none

    ! arguments..
    double precision, intent(in):: max_val, min_val, optimum, kurtosis, current

    ! Code with explicit min bound
    opt_max_scaling = exp( kurtosis*log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                    * exp( kurtosis*log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    ! Sanity check, allows for overlapping parameter ranges
    if (opt_max_scaling /= opt_max_scaling) opt_max_scaling = 0d0

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x, y, interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer:: interval
    double precision, dimension(interval):: x, y

    ! declare local variables
    double precision:: sum_x, sum_y, sumsq_x, sum_product_xy

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
  !------------------------------------------------------------------
  !
  double precision function temperature_impact( doptmin, dmaxmin, dttmin )

    ! Function to determine the coefficent for  !
    ! temperature impact on developmental rate. !
    ! From Streck et al., 2003.                 !

    implicit none

    ! arguments..
    double precision, intent(in):: doptmin, dmaxmin, dttmin   ! temperature differences

    ! local variables..
    double precision:: a, nmr, dnr

    a   = log( 2.d0 ) / ( log( ( dmaxmin ) / doptmin ) )

    nmr = 2.d0 * ( ( dttmin ) ** a ) * ( doptmin**a ) - ( ( dttmin ) ** ( 2.d0*a ) )

    dnr = doptmin ** ( 2.d0*a )

    temperature_impact = nmr/dnr

  end function temperature_impact
  !
  !------------------------------------------------------------------
  !
  double precision function vernalization( doptmin_v, dmaxmin_v, dttmin_v, days_in_step )

    ! Function to determine the coefficent for vernalization !
    ! impact on developmental rate. See Streck et al., 2003. !

    implicit none

    ! arguments..
    double precision, intent(in):: dmaxmin_v, doptmin_v, dttmin_v & ! temperature differences
                                  ,days_in_step

    ! local variables..
    double precision:: a, dnr, fvn, nmr

    a   = log( 2.d0 ) / ( log( ( dmaxmin_v ) / doptmin_v ) )
    nmr = 2.d0 * ( ( dttmin_v ) ** a ) * ( doptmin_v**a ) - ( ( dttmin_v ) ** (2.d0*a ) )
    dnr = doptmin_v ** ( 2.d0*a )
    fvn = nmr/dnr

    VD = VD + (fvn*days_in_step)

    ! final output value..
    vernalization = max( 0d0, min( 1d0, ( VD**5 ) / ( ( VDh**5 ) + (VD**5 ) ) ) )

  end function vernalization
  !
  !------------------------------------------------------------------
  !
  double precision function interpolate( x, reference_x, reference_y, row )

    ! Interpolation function.                    !
    ! x is input value, interpol is output value !
    ! reference_x/y are reference input data.    !

    implicit none

    ! arguments..
    integer, intent(in)                          :: row
    double precision, intent(in)                 :: x
    double precision, dimension(row), intent(in):: reference_x, reference_y

    ! local variables..
    integer:: i

    ! provide initial value
    interpolate = -9999d0
    do i = 1, row

       if ( x .le. reference_x(1) ) then
          interpolate = reference_y(1)
          exit
       endif

       ! cycling means growth rate remains constant between DS levels
       if ( ( x .gt. reference_x(i) ) .and. ( i .lt. row ) ) cycle

       if ( x .eq. reference_x(i) ) then
          interpolate = reference_y(i)
          exit
       endif

       if ( x .lt. reference_x(i) ) then
          interpolate = reference_y(i-1) + ( x-reference_x(i-1) ) &
                       * ( reference_y(i) - reference_y(i-1) )      &
                       / ( reference_x(i) - reference_x(i-1) )
          exit
       else
          interpolate = reference_y(row)
       endif

    enddo

    ! explicit return to ser
    return

  end function interpolate
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
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
