module CARBON_MODEL_MOD

  implicit none

  ! grassland model (DALEC.M2.#, aka DALEC.16.) developed from DALEC_GSI_DFOL_FR, aka DALEC.8.
  ! ----------------------------------------------------------------------------------------------------------------
  ! POOLS:   1.labile 2.foliar 3.root                        ! PARAMETERS:
  !          4.litter 5.som                                  !
  ! ------------------------------------------               ! 1.  Decomposition rate
  ! FLUXES:  1.GPP                                           ! 2.  Fraction of GPP respired 
  ! (daily)  2.temprate                                      ! 3.  GSI sens for leaf growth
  !          3.respiration_auto                              ! 4.  NPP belowground allocation parameter
  !          4.leaf production                               ! 5.  GSI max leaf turnover 
  !          5.labile production                             ! 6.  TOR roots
  !          6.root production                               ! 7.  TOR litter
  !          7.aboveground production                        ! 8.  TOR SOM
  !          8.labile consumption -> leaves                  ! 9.  Temp factor Q10 (1.2-1.6)
  !          9.leaffall factor                               ! 10. Photosynthetic N use efficiency
  !          10.leaf litter production                       ! 11. GSI max labile turnover
  !         *11.woodlitter production                        ! 12. GSI min temperature threshold (K)
  !          12.rootlitter production                        ! 13. GSI max temperature threshold (K)
  !          13.respiration het litter                       ! 14. GSI min photoperiod threshold (sec)
  !          14.respiration het som                          ! 15. LCA - g.C.leaf_m-2
  !          15.litter2som                                   ! 16. C labile (initialization)
  !          16.labrelease factor(leaf growth)               ! 17. C foliar (initialization)
  !         *17.carbon flux due to fire                      ! 18. C roots  (initialization)
  !          18.growing season index                         ! 19. C litter (initialization)
  !          19.animal manure C soil input (per time step)   ! 20. GSI max photoperiod threshold (sec)
  !          20.animal resp co2 (per time step)              ! 21. GSI min VPD threshold (Pa) 
  !          21.animal ch4 (per time step)                   ! 22. GSI max VPD threshold (Pa)
  !          22. labile loss (per time step)                 ! 23. critical GPP for LAI increase (gC.m-2.day-1)
  !          23. foliage loss (per time step)                ! 24. GSI senstivity for leaf senescence 
  ! ------------------------------------------               ! 25. GSI - have I just left a growing state (>1)
  ! MET:     1.run day                                       ! 26. GSI - initial GSI value
  !          2.min T (C)                                     ! 27. DM min lim for grazing (kg.DM.ha-1)
  !          3.max T (C)                                     ! 28. DM min lim for cutting (kg.DM.ha-1) 
  !          4.Radiation (MJ.m-2)                            ! 29. leaf-vs-stem allocation factor
  !          5.CO2 (ppm)                                     ! 30. C SOM (initialization)
  !          6.DOY                                           ! 31. DM demand of animal weight (fraction) 
  !         *7.lagged precip                                 ! 32. Post-grazing labile loss (fraction)
  !          8.cutting/grazing :                             ! 33. Post-cut labile loss (fraction)
  !            - spatial mode = lai removed (m2.m-2)         ! 34. Minimum grazed biomass to allow grazing
  !            - field mode = LSU.ha-1                       
  !         *9.burnt area fraction                           
  !          10.21-day avg min T (K)                        
  !          11.21-day avg photoperiod (sec) 
  !          12.21-day avg VPD (Pa)         
  !         *13.Forest mgmt after clearing
  !         *14.Mean T
  ! ----------------------------------------------------------------------------------------------------------------
  ! NOTES : '*' above means not used/applicable for grasslands 
  !         1 LSU per ha = 1 cow that weighs 650kg and grazes on 1 ha of grassland 
  !         carbon = 0.475 * dry matter 
  !         1 g.C.m-2 = 1 * 0.021 t.DM.ha-1
  !         to compile this .f90 into a python shared object (.so) run: f2py -c DALEC_GRASS.f90 -m DALEC_GRASS
  ! ----------------------------------------------------------------------------------------------------------------
  !                  autotrophic      heterotrophic     loss due to     --->    manure from          
  !                  respiration      respiration       grazing/cutting         grazing livestock       
  !                       ^            ^                   ^                       |       
  !                       |            |                   |                       V      
  !                                                                                   
  !PHOTOSYNTHESIS -----> [0] -------> [0] --------------> [0] <-----------------> [0]     
  !                 GPP         NPP            NEE                   NBE                 
  ! ----------------------------------------------------------------------------------------------------------------

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
  ! This version of DALEC is derived from the following primary references:
  ! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
  ! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
  ! Subsequent modifications by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! J. F. Exbrayat (University of Edinburgh)
  ! V. Myrgiotis (UK Centre for Ecology & Hydrology)
  ! S. Zhu (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!


  ! make all private
  ! private

  ! explicit publics
  public :: CARBON_MODEL          &
           ,acm                   &
           ,linear_model_gradient &
           ,CiCa_time             &
           ,soil_frac_clay        &
           ,soil_frac_sand        &
           ,nos_soil_layers  

  ! useful technical parameters
  double precision, parameter :: vsmall = 1d-36 &!tiny(0d0)*1d6 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  ! Mathematical constants
  double precision, parameter :: pi = 3.1415927, &
                                 deg_to_rad = pi/180d0

  ! Model specific double precision variables
  double precision :: & ! GSI phenology model
                      tmp,gradient                & 
                     ,fol_turn_crit,lab_turn_crit &
                     ,gsi_history(22),just_grown  &
                     ,ci & ! Photosynthetic variables
                     ,gpppars(12)        & ! ACM inputs (LAI+met)
                     ,constants(10)      & ! parameters for ACM
                     ,foliage_frac_res   & ! Management related variables
                     ,labile_frac_res    &
                     ,roots_frac_res     &
                     ,roots_frac_death  
                                              
  ! Model specific integer variables
  integer :: gsi_lag_remembered, & 
                   two_week_lag, & ! Management related variables
                  four_week_lag     

  ! Multiple soil layer variables, these are not used in DALEC16,
  ! but declarations are needed to here ensure compilation compatability with more complex verison of DALEC
  integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
  double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand

  ! Model specific allocatable variables
  double precision, allocatable, dimension(:) :: tmp_x, tmp_m, CiCa_time

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out,NEE,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP) 

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   & 
                          ,nodays   & ! number of days in simulation
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nopools  & ! number of model pools
                          ,nofluxes   ! number of model fluxes

    double precision, intent(in) :: deltat(nodays)    & ! time step in decimal days
                                   ,lat               & ! site latitude (degrees)
                                   ,met(nomet,nodays) & ! met drivers
                                   ,pars(nopars)        ! number of parameters

    double precision, intent(out) :: lai_out(nodays) & ! leaf area index
                                    ,GPP(nodays)     & ! Gross primary productivity
                                    ,NEE(nodays)       ! net ecosystem exchange of CO2

    double precision, intent(out) :: POOLS((nodays+1),nopools) ! vector of ecosystem pools
 
    double precision, intent(out) :: FLUXES(nodays,nofluxes) ! vector of ecosystem fluxes
    
    !f2py intent(in) :: start,finish,met,pars,deltat,nodays,lat,nopars,nomet,nopools,nofluxes  

    !f2py intent(out) ::lai_out,NEE,FLUXES,POOLS,GPP    

    integer :: f,m,n,test,gsi_lag

    double precision :: gsi_lai_reduction &
                       ,tot_abg_exp       &
                       ,fol_frac,lab_frac &
                       ,f_root,NPP                

    ! Reset pools and fluxes
    POOLS = 0d0 ; FLUXES = 0d0

    ! load some values
    gpppars(4)  = 2d0  ! g N leaf_m-2
    gpppars(7)  = lat
    gpppars(9)  = -2d0 ! leafWP-soilWP
    gpppars(10) = 1d0 ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)  = pars(11) 
    constants(2)  = 0.0156935d0
    constants(3)  = 4.22273d0
    constants(4)  = 208.868d0
    constants(5)  = 0.0453194d0
    constants(6)  = 0.37836d0
    constants(7)  = 7.19298d0
    constants(8)  = 0.011136d0
    constants(9)  = 2.1001d0
    constants(10) = 0.789798d0

    ! post-removal residues and root death | 0:none 1:all
    foliage_frac_res  = 0.05d0  ! fraction of removed foliage that goes to litter
    labile_frac_res   = 0.05d0  ! fraction of removed labile that goes to litter
    roots_frac_res    = 1d0     ! fraction of roots that die which go to litter 
    roots_frac_death  = 0.01d0  ! fraction of roots that dies
    

    ! How many steps in 2 weeks
    two_week_lag = ceiling(14d0/deltat(1))
    ! How many steps in 4 weeks
    four_week_lag = ceiling(28d0/deltat(1))

    ! assigning initial conditions
    POOLS(1,1) = pars(16)
    POOLS(1,2) = pars(17)
    POOLS(1,3) = pars(18)
    POOLS(1,4) = pars(19)
    POOLS(1,5) = pars(23)

    ! Allocate dimension to time varying module level variables
    if (.not.allocated(CiCa_time)) allocate(CiCa_time(nodays))

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
               ! else we will try and work out the gradient to see what is happening
               ! to the system over all. The default assumption will be to consider
               ! the averaging period of GSI model (i.e. 21 days). If this is not
               ! possible either the time step of the system is used (if step greater
               ! than 21 days) or all available steps (if n < 21).
               m = 0 ; test = 0
               do while (test < 21)
                  m = m+1 ; test = sum(deltat((n-m):n))
                  if (m > (n-1)) then 
                      test = 21 
                  endif
               end do
               tmp_m(n) = m
           endif ! for calculating gradient
        end do ! calc daily values once
        ! allocate GSI history dimension
        gsi_lag_remembered = max(2,maxval(nint(tmp_m)))
    end if ! .not.allocated(tmp_x)
    ! assign our starting value
    gsi_history = pars(26)-1d0
    just_grown = pars(25)

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit = pars(24)-1d0
    lab_turn_crit = pars(3)-1d0

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish  

      ! Update current LAI 
      lai_out(n) = POOLS(n,2) / pars(15)

      ! load next met / lai values for ACM
      gpppars(1) = lai_out(n)   ! LAI
      gpppars(2) = met(3,n) ! max temp
      gpppars(3) = met(2,n) ! min temp
      gpppars(5) = met(5,n) ! co2
      gpppars(6) = ceiling(met(6,n)-(deltat(n)*0.5d0)) ! doy
      gpppars(8) = met(4,n) ! radiation

      ! temprate (i.e. T modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(9)*0.5d0*(met(3,n)+met(2,n)))

      ! GPP and direct allocation of GPP can only occur 
      ! if sufficient LAI to prevent numerical error
      if (lai_out(n) > vsmall) then 

          ! GPP (gC.m-2.day-1)
          FLUXES(n,1) = acm(gpppars,constants)
          CiCa_time(n) = ci / met(5,n)

          ! Allocate to autotrophic respiration (gC.m-2.day-1)
          FLUXES(n,3) = FLUXES(n,1) * pars(2)

          !! Determine direct allocation of GPP to plant tissues   
          ! Dynamic allocation to roots vs aboveground biomass after Reyes.et.al.2017 (10.1002/2017MS001022)
          ! min/max allocation to roots as fraction of NPP
          f_root = 1d0 - exp(-1d0*pars(4)*lai_out(n))
          if (f_root < 0.1d0) f_root = 0.1d0
          if (f_root > 0.7d0) f_root = 0.7d0
      
          ! allocation to roots 
          FLUXES(n,6) = (FLUXES(n,1) - FLUXES(n,3)) * f_root
          ! allocation of ABG C to leaves
          FLUXES(n,4) = (FLUXES(n,1) - FLUXES(n,3) - FLUXES(n,6)) * (1d0 - (pars(29)*(lai_out(n)/6d0)))
          ! allocation of ABG C to labile/stem using pars(29)
          ! Ostrem.et.al.2013 (10.1080/09064710.2013.819440)            
          FLUXES(n,5) = FLUXES(n,1) - FLUXES(n,3) - FLUXES(n,6) - FLUXES(n,4)

      end if ! lai_out > vsmall, i.e. can be no allocation if no GPP

      ! Determine GSI based canopy phenology
      call gsi_phenology(nodays,n,deltat(n),gsi_lag,             &
                         met(10,n),met(11,n),met(12,n),          &
                         pars(5),pars(10),pars(15),pars(30),     &
                         pars(12),pars(13),pars(14),pars(20),    & 
                         pars(21),pars(22),                      &
                         FLUXES(n,1),FLUXES(:,18),               &
                         FLUXES(n,15),FLUXES(n,16),FLUXES(n,17), &
                         FLUXES(n,8),FLUXES(n,14),               &
                         POOLS(n,1),POOLS(n,2))

      ! FLUXES WITH TIME DEPENDENCIES

      ! labile release = P_labile * (1-(1-leafgrowth)**deltat)/deltat
      FLUXES(n,7)  = POOLS(n,1)*(1d0-(1d0-FLUXES(n,14))**deltat(n))/deltat(n)
      ! leaf litter production = P_foliar * (1-(1-leaffall)**deltat)/deltat  
      FLUXES(n,9)  = POOLS(n,2)*(1d0-(1d0-FLUXES(n,8))**deltat(n))/deltat(n)
      ! root litter production = P_root * (1-(1-rootTOR)**deltat)/deltat  
      FLUXES(n,10) = POOLS(n,3)*(1d0-(1d0-pars(6))**deltat(n))/deltat(n)

      ! FLUXES WITH TEMP AND TIME DEPENDENCIES

      ! resp het litter = P_litter * (1-(1-GPP_respired*litterTOR)**deltat)/deltat  
      FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-FLUXES(n,2)*pars(7))**deltat(n))/deltat(n)
      ! resp het som = P_som * (1-(1-GPP_respired*somTOR)**deltat)/deltat
      FLUXES(n,12) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))/deltat(n)
      ! litter to som = P_litter * (1-(1-dec_rate*temprate)**deltat)/deltat
      FLUXES(n,13) = POOLS(n,4)*(1d0-(1d0-pars(1)*FLUXES(n,2))**deltat(n))/deltat(n)

      ! NEE = resp_auto + resp_het_litter + resp_het_som - GPP [i.e. '-' when CO2 sink '+' when CO2 source ]
      NEE(n) = (FLUXES(n,3) + FLUXES(n,11) + FLUXES(n,12)) - FLUXES(n,1)
      ! GPP 
      GPP(n) = FLUXES(n,1)

      ! update pools for next timestep

      ! labile pool = labile_pool[†-1] + (lab_prod - lab_cons)*deltat
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,7))*deltat(n)
      ! foliar pool = foliar_pool[†-1] + (leaf_prod + lab_prod2 - leaf_litter_prod)*deltat
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)+FLUXES(n,7)-FLUXES(n,9))*deltat(n)
      ! root pool = root_pool[†-1] + (root_prod - root_litter_prod)*deltat
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,10))*deltat(n)
      ! litter pool = litter_pool[†-1] + (leaf_litter_prod + root_litter_prod - resp_het_litter - litter2som)*deltat
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,9)+FLUXES(n,10)-FLUXES(n,11)-FLUXES(n,13))*deltat(n)
      ! som pool = som_pool[†-1] + (litter2som - resp_het_som)
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,13)-FLUXES(n,12))*deltat(n)

      ! CUTTING 
      ! ------------------------------------------------------------------------------------------------------------- ! 

      if (met(8,n) == -1d0) then
          call grass_cutting(POOLS(n+1,1),POOLS(n+1,2),POOLS(n+1,3), & 
                             POOLS(n+1,4),POOLS(n+1,5),met(6,n),met(8,n), &
                             FLUXES(:,22),FLUXES(n,25),FLUXES(n,26),FLUXES(n,27), &
                             FLUXES(n,28),FLUXES(n,29),FLUXES(n,30), &
                             n,nodays,deltat(n), & 
                             pars(28),pars(33))
      end if

      ! GRAZING 
      ! ------------------------------------------------------------------------------------------------------------- ! 

      ! Determine whether there is any remaining losses to occur given GSI driven foliar loss
      gsi_lai_reduction = FLUXES(n,9) * deltat(n) * pars(15)
      if (met(8,n)-gsi_lai_reduction > 0d0) then
          call grass_grazing(POOLS(n+1,1),POOLS(n+1,2),POOLS(n+1,3), & 
                             POOLS(n+1,4),POOLS(n+1,5),met(6,n),met(8,n)-gsi_lai_reduction, &
                             FLUXES(n,19),FLUXES(n,20),FLUXES(n,21), &
                             FLUXES(:,22),FLUXES(:,23),FLUXES(n,31),FLUXES(n,32), &
                             FLUXES(n,33),FLUXES(n,34),FLUXES(n,35),FLUXES(n,36), &
                             n,nodays,deltat(n), & 
                             pars(15),pars(27),pars(32),pars(34))
      end if                      

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm(drivers,constants)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(12) & ! acm input requirements
                                   ,constants(10) ! ACM parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, e0, dayl, cps, dec, nit &
                       ,trange, sinld, cosld,aob  &
                       ,mint,maxt,radiation,co2,lai,doy,lat &
                       ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
                       ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
                       ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! initial values
    gc = 0d0 ; pp = 0d0 ; qq = 0d0 ; ci = 0d0 ; e0 = 0d0 
    dayl = 0d0 ; cps = 0d0 ; dec = 0d0

    ! load driver values to correct local vars
    lai  = drivers(1)
    maxt = drivers(2)
    mint = drivers(3)
    nit  = drivers(4)   
    co2  = drivers(5)
    doy  = drivers(6)
    lat = drivers(7)
    radiation = drivers(8)
    deltaWP = drivers(9)
    Rtot = drivers(10)

    ! load parameters into correct local vars
    NUE = constants(1)
    dayl_coef = constants(2)
    co2_comp_point = constants(3) 
    co2_half_sat = constants(4)
    dayl_const = constants(5)
    hydraulic_temp_coef = constants(6)
    lai_coef = constants(7)
    temp_exponent = constants(8)
    lai_const = constants(9)
    hydraulic_exponent = constants(10)

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
    ! calculate day length (hours)
    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
    sinld = sin( lat*deg_to_rad ) * sin( dec )
    cosld = cos( lat*deg_to_rad ) * cos( dec )
    aob = max(-1d0,min(1d0,sinld / cosld))
    dayl = 12d0 * ( 1d0 + 2d0 * asin( aob ) / pi )

    ! calculate CO2 limited rate of photosynthesis
    pd = gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps = e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm = cps*(dayl_coef*dayl+dayl_const)

    return

  end function acm
  !
  !------------------------------------------------------------------
  !
  subroutine grass_cutting(labile,foliage,roots,litter,som,doy,lai_reduction, &
                           harvest,HARVESTextracted_labile,                   & 
                           HARVESTextracted_foliage,HARVESTextracted_roots,   &
                           HARVESTlitter_labile,HARVESTlitter_foliage,        &
                           HARVESTlitter_roots,                               &
                           timestep,nodays,step_length,                       &
                           cutting_threshold,post_cutting_labile_loss)
  
    !! Determine whether cutting has occured and the resulting impacts on the C-cycle
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: timestep, nodays
    double precision, intent(in) :: doy, &
                            step_length, &
                          lai_reduction, &
                      cutting_threshold, &
               post_cutting_labile_loss
    
    double precision, intent(inout) :: labile, &
                                      foliage, &
                                        roots, &
                                       litter, & 
                                          som, &
                      HARVESTextracted_labile, &
                     HARVESTextracted_foliage, &
                       HARVESTextracted_roots, &
                         HARVESTlitter_labile, &
                        HARVESTlitter_foliage, &
                          HARVESTlitter_roots   
                          
    double precision, dimension(nodays), intent(inout) :: harvest
                           
    ! Local variables
    double precision :: labile_loss    & 
                       ,foliar_loss    &
                       ,roots_loss     &
                       ,labile_residue &
                       ,foliar_residue &
                       ,roots_residue  
           
    ! Determine whether cutting is plausible
    ! 1) Labile+leaf C > cutting threshold 
    ! 2) LAI > 3 (note replaced here with day of year constraints?)
    ! 3) & LAI reduction = -1 & no cut in past month 
    ! TLS: Question, conditions here are very temperate centric, can these be modifed?
    ! TLS: The timestepping of the cutting assumptions need to be dynamics in code to timestep
    if ( (labile+foliage) >= cutting_threshold   & 
          .and. doy >= 91d0 .and. doy <= 304d0 & 
          ! .and. LAI(n) >= 3 & 
          .and. sum(harvest(max(1,timestep-four_week_lag):timestep)) == 0d0 ) then
                  
        ! direct C losses
        labile_loss  = labile * post_cutting_labile_loss
        foliar_loss  = foliage * 0.95d0 ! 95% of leaves lost after cutting probably 99% lost in reality 
        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death ! allocation to roots will be reduced due to reduced LAI 

        ! fraction of harvest wasted 
        labile_residue = labile_loss * labile_frac_res
        foliar_residue = foliar_loss * foliage_frac_res
        roots_residue  = roots_loss  * roots_frac_res

        ! if havest yields > 1500 kg.DM.ha-1 proceed with cut
        ! Note converted to gC/m2 equivalent assuming 47.5 % C content
        ! yields 71.25 gC/m2
        if ( ( (foliar_loss-foliar_residue)+ &
               (labile_loss-labile_residue)+ &
               (roots_loss -roots_residue ) ) >= 71.25d0 ) then
                      
            ! Assign to output the biomass extracted due to cutting
            HARVESTextracted_labile  = (labile_loss-labile_residue) / step_length
            HARVESTextracted_foliage = (foliar_loss-foliar_residue) / step_length
            HARVESTextracted_roots   = (roots_loss -roots_residue) / step_length
            ! Assign the output the biomass entering litter due to cutting
            HARVESTlitter_labile  = labile_residue / step_length
            HARVESTlitter_foliage = foliar_residue / step_length
            HARVESTlitter_roots   = roots_residue / step_length

            ! Combine to total extracted based on cutting
            harvest(timestep) = HARVESTextracted_labile  + &
                                HARVESTextracted_foliage + & 
                                HARVESTextracted_roots 

            ! update pools 
            labile  = max(0d0,labile-labile_loss)
            foliage = max(0d0,foliage-foliar_loss)
            roots   = max(0d0,roots-roots_loss)
            litter  = max(0d0,litter + (labile_residue+foliar_residue+roots_residue))
            som     = max(0d0,som)

        endif ! Determine whether the cut is actually plausible

    endif ! end cutting process   

  end subroutine grass_cutting
  !
  !------------------------------------------------------------------
  !
  subroutine grass_grazing(labile,foliage,roots,litter,som,doy,lai_reduction,       &
                           animal_manure_to_soil,animal_respiration,animal_methane, &
                           harvest,grazing,GRAZINGextracted_labile,                 & 
                           GRAZINGextracted_foliage,GRAZINGextracted_roots,         &
                           GRAZINGlitter_labile,GRAZINGlitter_foliage,              &
                           GRAZINGlitter_roots,                                     &
                           timestep,nodays,step_length,                             &
                           lca,grazing_threshold,post_grazing_labile_loss,          &
                           min_grazing_removal_threshold)
  
    !! Determine whether grazing has occured and the resulting impacts on the C-cycle
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: timestep, nodays
    double precision, intent(in) :: doy, &
                            step_length, &
                          lai_reduction, &
                                    lca, &
                      grazing_threshold, &
          min_grazing_removal_threshold, &
               post_grazing_labile_loss
    
    double precision, intent(inout) :: labile, &
                                      foliage, &
                                        roots, &
                                       litter, & 
                                          som, & 
                      GRAZINGextracted_labile, &
                     GRAZINGextracted_foliage, &
                       GRAZINGextracted_roots, &
                         GRAZINGlitter_labile, &
                        GRAZINGlitter_foliage, &
                          GRAZINGlitter_roots, &
                        animal_manure_to_soil, &
                           animal_respiration, &
                               animal_methane
                       
    double precision, dimension(nodays), intent(inout) :: harvest, &
                                                          grazing

    ! Local variables
    double precision :: labile_loss    & 
                       ,foliar_loss    &
                       ,roots_loss     &
                       ,labile_residue &
                       ,foliar_residue &
                       ,roots_residue 
        
    ! Determine whether grazing is plausible
    ! 1) An LAI reduction is specified
    ! 2) Labile+leaf C > grazing threshold 
    ! 3) No cutting in the last 2 weeks
    if (labile+foliage >= grazing_threshold .and. & 
        sum(harvest(max(1,timestep-two_week_lag):timestep)) == 0d0) then
           
        ! direct C losses
        labile_loss  = labile * post_grazing_labile_loss
        foliar_loss  = max(0d0,(lai_reduction * lca) - labile_loss)  
        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death

        ! fraction of harvest wasted 
        labile_residue = labile_loss * labile_frac_res
        foliar_residue = foliar_loss * foliage_frac_res
        roots_residue = roots_loss * roots_frac_res

        ! extracted C via grazing: if remaining AGB > pre-grazing limit DM & grazed biomass > pars(34) gC/m2/step
        if ( (labile+foliage)-foliar_loss-labile_loss >= grazing_threshold .and. & 
               (foliar_loss+labile_loss) >= min_grazing_removal_threshold*step_length ) then

            ! Assign to output the biomass extracted due to cutting
            GRAZINGextracted_labile = (labile_loss-labile_residue) / step_length
            GRAZINGextracted_foliage = (foliar_loss-foliar_residue) / step_length
            GRAZINGextracted_roots = (roots_loss -roots_residue) / step_length
            ! Assign the output the biomass entering litter due to cutting
            GRAZINGlitter_labile = labile_residue / step_length
            GRAZINGlitter_foliage = foliar_residue / step_length
            GRAZINGlitter_roots = roots_residue / step_length

            ! Combine to total extracted based on grazing
            grazing(timestep) = GRAZINGextracted_labile  + &
                                GRAZINGextracted_foliage + &
                                GRAZINGextracted_roots
               
            ! Constants used for animal C fluxes from Vertes.et.al.2019 (10.1016/B978-0-12-811050-8.00002-9)
            ! TLS: Meaning of the constants needs to be defined.

            ! animal manure-C production (gC/m2)
            animal_manure_to_soil = (grazing(timestep) * 0.32d0) / step_length
            ! animal respiration CO2-C (gC/m2)
            animal_respiration = (grazing(timestep) * 0.54d0) / step_length
            ! animal CH4-C (gC/m2)
            animal_methane = (grazing(timestep) * 0.04d0) / step_length                

            ! update pools 
            labile  = max(0d0,labile-labile_loss)
            foliage = max(0d0,foliage-foliar_loss)
            roots   = max(0d0,roots-roots_loss)
            litter  = max(0d0,litter+((animal_manure_to_soil+GRAZINGlitter_labile + &
                                       GRAZINGlitter_foliage+GRAZINGlitter_roots) * step_length))
            som     = max(0d0,som)

        endif 

        ! Determine whether to graze, assuming grazing criterior above was not met 
        ! extracted C via grazing: if not done above & postgraze AGB < pre-grazing AGB 
        if ( (grazing(timestep) == 0d0) .and. & 
             ((labile+foliage)-foliar_loss-labile_loss) <= grazing_threshold ) then
              
            ! Bulk C losses
            labile_loss  = labile * post_grazing_labile_loss
            foliar_loss  = foliage - grazing_threshold - (labile - labile_loss)
            !foliar_loss  = foliage - (grazing_threshold - labile_loss)
            roots_loss   = 0d0 ! roots * roots_frac_death

            ! fraction of harvest wasted 
            labile_residue = labile_loss * labile_frac_res
            foliar_residue = foliar_loss * foliage_frac_res
            roots_residue = roots_loss * roots_frac_res

            ! proceed if simulating this grazing will remove > ~0.5 gCm-2 from AGB
            if ((foliar_loss+labile_loss) >= min_grazing_removal_threshold) then

                ! Assign to output the biomass extracted due to cutting
                GRAZINGextracted_labile  = (labile_loss-labile_residue) / step_length
                GRAZINGextracted_foliage = (foliar_loss-foliar_residue) / step_length
                GRAZINGextracted_roots   = (roots_loss -roots_residue) / step_length
                ! Assign the output the biomass entering litter due to cutting
                GRAZINGlitter_labile  = labile_residue / step_length
                GRAZINGlitter_foliage = foliar_residue / step_length
                GRAZINGlitter_roots   = roots_residue / step_length

                ! Combine to total extracted based on grazing
                grazing(timestep) = GRAZINGextracted_labile + &
                                    GRAZINGextracted_foliage + &
                                    GRAZINGextracted_roots

                ! animal manure-C production (gC/m2)
                animal_manure_to_soil = (grazing(timestep) * 0.32d0) / step_length
                ! animal respiration CO2-C (gC/m2)
                animal_respiration = (grazing(timestep) * 0.54d0) / step_length
                ! animal CH4-C (gC/m2)
                animal_methane = (grazing(timestep) * 0.04d0) / step_length                

                ! update pools 
                labile  = max(0d0,labile-labile_loss)
                foliage = max(0d0,foliage-foliar_loss)
                roots   = max(0d0,roots-roots_loss)
                litter  = max(0d0,litter+((animal_manure_to_soil+GRAZINGlitter_labile + &
                                           GRAZINGlitter_foliage+GRAZINGlitter_roots) * step_length))
                som     = max(0d0,som)

            endif ! carry out grazing?

        endif ! Is grazing plausible

    endif ! end grazing process

  end subroutine grass_grazing
  !
  !------------------------------------------------------------------
  !
  subroutine gsi_phenology(nodays,timestep,step_length,gsi_lag, &
                           maxt_21day,dayl_21day,vpd_21day,     &
                           potential_foliar_turnover,potential_labile_turnover, &
                           lca,critical_gpp_return, &
                           Tfac_min,Tfac_max,Pfac_min,Pfac_max, & 
                           VPDfac_min,VPDfac_max,               &
                           current_gpp,gsi,     & 
                           Tfac,Pfac,VPDfac,    &
                           foliar_turnover_fraction, &
                           labile_turnover_fraction, &
                           labile,foliage)
 
    !! Description
  
    implicit none
  
    ! Declare arguments
    integer, intent(in) :: nodays, gsi_lag, timestep
    double precision, intent(in) :: step_length, &
                                     maxt_21day, & 
                                     dayl_21day, & 
                                      vpd_21day, &
                      potential_foliar_turnover, & 
                      potential_labile_turnover, & 
                              Tfac_min,Tfac_max, &
                              Pfac_min,Pfac_max, & 
                          VPDfac_min,VPDfac_max, & 
                                            lca, &
                                    current_gpp, &
                            critical_gpp_return, &
                                         labile, & 
                                        foliage
    double precision, dimension(nodays), intent(inout) :: gsi
    double precision, intent(out) :: foliar_turnover_fraction, & 
                                     labile_turnover_fraction, &
                                                         Tfac, & 
                                                         Pfac, & 
                                                       VPDfac
  
    ! Declare local arguments
    integer :: interval
    double precision :: tmp

    ! Calculate the Growing Season Index based on Jolly et al. 
    ! doi: 10.1111/j.1365-2486.2005.00930.x doi:10.1029/2010JG001545.
    ! It is the product of 3 limiting factors for temperature, photoperiod and
    ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated 
    ! min and max value. Photoperiod, VPD and avgTmin are direct input

    ! temperature limitation, then restrict to 0-1; correction for k-> oC
    ! Tfac = (met(10,n)-(pars(12)-273.15)) / (pars(13)-pars(12)) ! no need to K->C 
    Tfac = ( maxt_21day-(Tfac_min-273.15d0)) / (Tfac_max-Tfac_min )
    Tfac = min(1d0,max(0d0,Tfac))
    ! photoperiod limitation (seconds)
    Pfac = ( dayl_21day-Pfac_min) / (Pfac_max-Pfac_min )
    Pfac = min(1d0,max(0d0,Pfac))
    ! VPD limitation (Pa)
    VPDfac = 1d0 - ( (vpd_21day-VPDfac_min) / (VPDfac_max-VPDfac_min) )
    VPDfac = min(1d0,max(0d0,VPDfac))
    ! calculate and store the GSI index
    gsi(timestep) = Tfac * Pfac * VPDfac
  
    ! Determine GSI section to have linear regression applied to and
    ! determine the number of values, i.e. the interval
    if (timestep < gsi_lag) then
        if (timestep == 1) then
            gsi_history(2) = gsi(timestep)
            interval = 2
        else
            gsi_history(1:timestep) = gsi(1:timestep)
            interval = timestep
        endif
    else
        gsi_history(1:gsi_lag) = gsi((timestep-gsi_lag+1):timestep)
        interval = gsi_lag
    end if
 
    ! Now calculate the linear gradient
    gradient = linear_model_gradient(tmp_x(1:interval),gsi_history(1:interval),interval)      
    gsi_lag_remembered = gsi_lag

    ! now update foliage and labile conditions based on gradient calculations
    if (gradient < fol_turn_crit .or. gsi(timestep) == 0d0) then
        ! we are in a decending condition so foliar turnover
        foliar_turnover_fraction = potential_foliar_turnover*(1d0-gsi(timestep))
        just_grown = 0.5d0
    else if (gradient > lab_turn_crit) then
        ! we are in a assending condition so labile turnover
        labile_turnover_fraction = potential_labile_turnover*gsi(timestep)
        just_grown = 1.5d0
        ! check carbon return
        tmp = labile*(1d0-(1d0-labile_turnover_fraction)**step_length)/step_length
        tmp = (foliage+tmp)/lca
        gpppars(1) = tmp
        tmp = acm(gpppars,constants)
        ! determine if increase in LAI leads to an improvement in GPP greater
        ! than critical value, if not then no labile turnover allowed      
        if ( ((tmp - current_gpp)/current_gpp) < critical_gpp_return ) then
            labile_turnover_fraction = 0d0
        endif
    else
        ! probably we want nothing to happen, however if we are at the seasonal
        ! maximum we will consider further growth still
        if (just_grown >= 1d0) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            labile_turnover_fraction = potential_labile_turnover*gsi(timestep)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = labile*(1d0-(1d0-labile_turnover_fraction)**step_length)/step_length
            tmp = (foliage+tmp)/lca
            gpppars(1) = tmp
            tmp = acm(gpppars,constants)
            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than critical value, if not then no labile turnover allowed
            if ( ((tmp - current_gpp)/current_gpp) < critical_gpp_return ) then
                labile_turnover_fraction = 0d0
            endif
        end if ! Just grown?
    endif ! gradient choice
        
  end subroutine gsi_phenology
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
    double precision :: sum_x, sum_y!, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) )
    !&
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
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD