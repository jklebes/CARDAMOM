
module model_likelihood_module
  implicit none

  ! make all private
  private

  ! which to make open
  public :: model_likelihood, find_edc_initial_values

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: EDC
    integer :: DIAG
    integer :: PASSFAIL(100) ! allow space for 100 possible checks
    integer :: nedc ! number of edcs being assessed
  end type
  type (EDCDIAGNOSTICS), save :: EDCD

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine find_edc_initial_values
    use MCMCOPT, only: PI, MCOUT, MCO
    use cardamom_structures, only: DATAin ! will need to change due to circular dependance
    use cardamom_io, only: restart_flag
    use MHMCMC_MODULE, only: MHMCMC

    ! subroutine deals with the determination of initial parameter and initial
    ! conditions which are consistent with EDCs

    implicit none

    ! declare local variables
    integer :: n, counter_local, EDC_iter
    double precision :: PEDC, ML, ML_prior
    double precision, dimension(PI%npars+1) :: EDC_pars
    ! declare parameters
    integer, parameter :: EDC_iter_max = 1

    ! set MCMC options needed for EDC run
    MCO%APPEND = 0
    MCO%nADAPT = 1000
    MCO%fADAPT = 1d0
    MCO%nOUT = 100000
    MCO%nPRINT = 0
    MCO%nWRITE = 0
    ! the next two lines ensure that parameter inputs are either given or
    ! entered as -9999
    MCO%randparini = .true.
    MCO%returnpars = .true.
    MCO%fixedpars  = .false.

    ! Set initial priors to vector...
    PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
    ! ... and assume we need to find random parameters
    PI%parfix = 0

    ! if the prior is not missing and we have not told the edc to be random
    ! keep the value
    do n = 1, PI%npars
       if (PI%parini(n) /= -9999d0 .and. DATAin%edc_random_search < 1) PI%parfix(n) = 1
    end do ! parameter loop

    ! set the parameter step size at the beginning
    PI%stepsize = 1d0 ; PI%beta_stepsize = 0.005d0
    PI%parvar = 1d0 ; PI%Nparvar = 0d0
    PI%use_multivariate = .false.
    ! Covariance matrix cannot be set to zero therefore set initial value to a
    ! small positive value along to variance access
    PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
    do n = 1, PI%npars
       PI%covariance(n,n) = 1d0
    end do

    EDC_pars = 1d0

    ! if this is not a restart run, i.e. we do not already have a starting
    ! position we must being the EDC search procedure to find an ecologically
    ! consistent initial parameter set
    if (.not. restart_flag) then

        do EDC_iter = 1, EDC_iter_max

           ! set up edc log likelihood for MHMCMC initial run
           PEDC = -1 ; counter_local = 0
           do while (PEDC < 0)

              write(*,*)"Beginning EDC search attempt"
              ! call the MHMCMC directing to the appropriate likelihood
              call MHMCMC(EDC_MODEL_LIKELIHOOD)

              ! store the best parameters from that loop
              PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
              ! turn off random selection for initial values
              MCO%randparini = .false.

              ! call edc likelihood function to get final edc probability
              call edc_model_likelihood(PI%parini,PEDC,ML_prior)

              ! keep track of attempts
              counter_local = counter_local + 1
              ! periodically reset the initial conditions
              if (PEDC < 0d0 .and. mod(counter_local,5) == 0) then
                  PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
                  ! reset to select random starting point
                  MCO%randparini = .true.
                  ! reset the parameter step size at the beginning of each attempt
                  PI%stepsize = 1d0 ; PI%beta_stepsize = 0.005d0
                  PI%parvar = 1d0 ; PI%Nparvar = 0d0
                  ! Covariance matrix cannot be set to zero therefore set initial value to a
                  ! small positive value along to variance access
                  PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
                  PI%use_multivariate = .false.
                  do n = 1, PI%npars
                     PI%covariance(n,n) = 1d0
                  end do
              endif

           end do ! for while condition

           ! check for actual likelihood score
           call model_Likelihood(PI%parini,ML,ML_prior)

           if (EDC_pars(PI%npars+1) > 0d0 .or. (ML+ML_prior) > EDC_pars(PI%npars+1)) then

               ! either this is our first EDC starting point or a better one, so
               ! we best make a record of these
               EDC_pars(1:PI%npars) = PI%parini(1:PI%npars)
               EDC_pars(PI%npars+1) = ML+ML_prior

           end if ! EDC_pars(PI%npars+1) > 0 .or. (ML+ML_prior) > EDC_pars(PI%npars+1)

           ! unless this is the last time we want to reset to control switches
           ! to reset the parameters for new staring points
           if (EDC_iter < EDC_iter_max) then
               PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
               ! reset to select random starting point
               MCO%randparini = .true.
               ! reset the parameter step size at the beginning of each
               ! attempt
               PI%stepsize = 1d0 ; PI%beta_stepsize = 0.005d0
               PI%parvar = 1d0 ; PI%Nparvar = 0d0
               ! Covariance matrix cannot be set to zero therefore set
               ! initial value to a
               ! small positive value along to variance access
               PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
               PI%use_multivariate = .false.
               do n = 1, PI%npars
                  PI%covariance(n,n) = 1d0
               end do
           endif ! EDC_iter < EDC_iter_max

        end do ! EDC_iter

        ! now pass the best parameter set to the inital value to the main EDC
        PI%parini(1:PI%npars) = EDC_pars(1:PI%npars)

    endif ! if for restart

    ! reset so that currently saved parameters will be used
    ! starting point in main MCMC
    PI%parfix(1:PI%npars) = 0
    MCOUT%best_pars = 0d0

  end subroutine find_edc_initial_values
  !
  !------------------------------------------------------------------
  !
  subroutine edc_model_likelihood(PARS, ML_obs_out, ML_prior_out)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use CARBON_MODEL_MOD, only: carbon_model

    ! Model likelihood function specifically intended for the determination of
    ! appropriate initial parameter choices, consistent with EDCs for DALEC2 /
    ! DALEC_GSI

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS
    ! output
    double precision, intent(inout) :: ML_obs_out, ML_prior_out

    ! declare local variables
    integer ::  n
    double precision :: tot_exp, ML, EDC1, EDC2, infini

    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 1
    ML_obs_out = 0d0 ; ML_prior_out = 0d0

    ! call EDCs which can be evaluated prior to running the model
    call EDC1_CDEA_LU_FIRES(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

    ! assess post running EDCs
    call EDC2_CDEA_LU_FIRES(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                     ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                     ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                     ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

    ! calculate the likelihood
    tot_exp = sum(1d0-EDCD%PASSFAIL(1:EDCD%nedc))
!    tot_exp = 0d0
!    do n = 1, EDCD%nedc
!       tot_exp=tot_exp+(1d0-EDCD%PASSFAIL(n))
!       if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
!    end do ! checking EDCs
!    ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) then
!        print*,"Found it!" ; stop
!    endif

    ! convert to a probability
!    ML_obs_out = -0.5d0*(tot_exp*10d0)*DATAin%EDC
    ML_obs_out = -5d0*tot_exp*DATAin%EDC

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine EDC1_CDEA_LU_FIRES(PARS, npars, meantemp, meanrad, EDC1)

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (Bloom et al., 2014).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    double precision :: torfol ! yearly leaf loss fraction

    ! set initial value
    torfol = 1d0/(pars(5)*365.25d0)
    EDC1 = 1
    DIAG = EDCD%DIAG
    fauto = pars(2)
    ffol = (1d0-fauto)*pars(3)
    flab = (1d0-fauto-ffol)*pars(13)
    froot = (1d0-fauto-ffol-flab)*pars(4)
    fwood = 1d0-fauto-ffol-flab-froot
    fsom = fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(8))

    ! set all EDCs to 1 (pass)
    EDCD%nedc = 100
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > pars(8))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    endif

    ! litter2som greater than som to atm rate
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(1) < pars(9))) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    endif

    ! turnover of foliage faster than turnover of wood
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > torfol) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    end if

    ! root turnover greater than som turnover at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(7) < (pars(9)*exp(pars(10)*meantemp)))) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    endif

    ! GPP allocation to foliage and labile cannot be 5 orders of magnitude
    ! difference from GPP allocation to roots
    if ((EDC1 == 1 .or. DIAG == 1) .and. ((ffol+flab) > (5.*froot) .or. ((ffol+flab)*5.) < froot)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    endif

    ! could always add more / remove some

  end subroutine EDC1_CDEA_LU_FIRES
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_CDEA_LU_FIRES(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    ! Determines whether the dynamical contraints for the search of the initial
    ! parameters has been successful or whether or not we should abandon the
    ! current set and move on

    implicit none

    ! declare input variables
    integer, intent(in) :: npars    & ! number of model parameters
                          ,nomet    & ! number of met drivers
                          ,nofluxes & ! number of fluxes from model
                          ,nopools  & ! number of pools in model
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_LAI(nodays)               & ! LAI output from current model simulation
                                   ,M_NEE(nodays)               & ! NEE output from current model simulation
                                   ,M_GPP(nodays)               & ! GPP output from current model simulation
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    integer :: n, DIAG, no_years, y, PEDC, nn, num_EDC
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF
    double precision, dimension(:,:), allocatable :: mean_annual_pools ! FE - 29/06/2018 changed dimensions to keep
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom  & ! fraction of GPP som under eqilibrium conditions
             ,flit    ! fraction of GPP to litter under equilibrium conditions

    !JFE - 27/06/2018 newly defined variables for updated EDCs
    double precision :: FT(nofluxes), Fin(nopools), Fout(nopools)
    double precision :: fin_fout_lim, Sprox, Sprox0
    integer :: nd, fl


    ! update initial values
    DIAG = EDCD%DIAG
    EDC2 = 1
    fauto = pars(2)
    ffol = (1d0-fauto)*pars(3)
    flab = (1d0-fauto-ffol)*pars(13)
    froot = (1d0-fauto-ffol-flab)*pars(4)
    fwood = 1d0-fauto-ffol-flab-froot
    fsom = fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(8))
    flit = (froot+flab+ffol)

    ! derive mean pools
    do n = 1, nopools
       mean_pools(n) = cal_mean_pools(M_POOLS,n,nodays+1,nopools)
    end do

    !
    ! Begin EDCs here
    !

    ! EDC 6
    ! ensure ratio between Cfoilar and Croot is less than 5
    if ((EDC2 == 1 .or. DIAG == 1) .and. (mean_pools(2) > (mean_pools(3)*5.) .or. (mean_pools(2)*5.) < mean_pools(3)) ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(6) = 0
    end if

    ! EDC 7
    ! Assess growth factor over 10 years; order magnitude changes are not
    ! allowed
    no_years = int(floor(sum(deltat)/365d0))
    G = 0.1
    allocate(mean_annual_pools(no_years,nopools)) ! JFE - 29/06/2018 changed allocated dimensions

    ! generate mean annual pool values
    do n = 1, nopools
       ! rapid pool growth is not allowed, increase is restricted by G growth
       ! over N years. Rapid decay is dealth with in a later EDC
       do y = 1, no_years
          ! derive mean annual pools
          mean_annual_pools(y,n)=cal_mean_annual_pools(M_POOLS(1:nodays+1,n),y,deltat,nodays+1)
      !    print *, y,n, mean_annual_pools(y,n)
       end do ! year loop
       ! now check the growth rate
       if ((EDC2 == 1 .or. DIAG == 1) .and. &
          ((mean_annual_pools(no_years,n)/mean_annual_pools(1,n)) > (1.+G*real(no_years)))) then
          EDC2 = 0d0 ; EDCD%PASSFAIL(7+n-1) = 0
       endif
    end do ! pool loop

    ! done now so clean up
    !deallocate(mean_annual_pools)

    ! EDC 8
    ! assesses the exponential decay of each model pool

    ! JFE - EDC8 (expdecay) commented out to match EDCs in PNAS paper - 27/06/2018
    ! loop through each pool in turn
    !do n = 1, nopools
    !   if (EDC2 == 1 .or. DIAG == 1) then
    !      decay_coef=expdecay2(M_POOLS,n,deltat,nopools,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
    !      if (abs(-log(2.)/decay_coef) < (365.25*real(no_years)) .and. decay_coef < 0. ) then
    !         EDC2 = 0d0 ; EDCD%PASSFAIL(8)=0
    !      end if ! EDC conditions
    !   end if ! EDC .or. DIAG condition
    !end do ! pools loop

    ! SOM attractor - must be within a factor of 2 from Csom0
    ! eqiulibrium factor (in comparison with initial conditions)
    EQF = 2d0 ! JFE replaced 10 by 2 - 27/06/2018

    ! initialise and then calculate mean gpp values
    !meangpp=sum(M_GPP(1:nodays))/real(nodays)

    ! EDC 9 - SOM steady state within order magnitude of initial conditions - 27/06/2018
    ! JFE - EDC9 (steady-state proximity) commented outto match EDCs in PNAS paper
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(pars(9)*exp(pars(10)*meantemp))) > (pars(23)*EQF)) then
    !   EDC2 = 0d0 ; EDCD%PASSFAIL(9) = 0
    !end if
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(pars(9)*exp(pars(10)*meantemp))) < (pars(23)/EQF)) then
    !   EDC2 = 0d0 ; EDCD%PASSFAIL(9) = 0
    !endif

    ! EDC 10 - Litter steady state assumptions
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(pars(8)*exp(pars(10)*meantemp))) > (pars(22)*EQF)) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(10) = 0
    !endif
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(pars(8)*exp(pars(10)*meantemp))) < (pars(22)/EQF)) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(10) = 0
    !endif

    ! EDC 11 - Wood steady state assumptions
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fwood)/pars(6)) > (pars(21)*EQF)) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(11) = 0
    !end if
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fwood)/pars(6)) < (pars(21)/EQF)) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(11) = 0
    !endif

    ! EDC 12 - Root steady state assumptions
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*froot)/pars(7)) > (pars(20)*EQF)) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(12) = 0
    !endif
    !if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*froot)/pars(7)) < (pars(20)/EQF)) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(12) = 0
    !endif

    ! 27/06/2018 - JFE new EDC 8 to replace 9-12 to avoid dismissing simulations with early
    ! fires that could lead to exponential regrowth / decay (see Bloom et al PNAS 2016 SI)

    ! first calculate total flux for the whole simulation period
    do fl = 1, nofluxes
        FT(fl) = 0
        do nd = 1, nodays
            FT(fl) = FT(fl) + M_FLUXES(nd,fl)*deltat(nd)
        end do
    end do

    ! get total in and out for each pool
    ! labile
    Fin(1)  = FT(5)
    Fout(1) = FT(8)+FT(18)+FT(24)
    ! foliar
    Fin(2)  = FT(4)+FT(8)
    Fout(2) = FT(10)+FT(19)+FT(25)
    ! root
    Fin(3)  = FT(6)
    Fout(3) = FT(12)+FT(20)+FT(26)
    ! wood
    Fin(4)  = FT(7)
    Fout(4) = FT(11)+FT(21)+FT(27)
    ! litter
    Fin(5)  = FT(10)+FT(12)+FT(24)+FT(25)+FT(26)
    Fout(5) = FT(13)+FT(15)+FT(22)+FT(28)
    ! som
    Fin(6)  = FT(11)+FT(15)+FT(27)+FT(28)
    Fout(6) = FT(14)+FT(23)

    ! iterate to check whether Fin/Fout is within EQF limits
    do n = 1, nopools
        if (abs(log(Fin(n)/Fout(n))) > log(EQF)) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(13+n-1) = 0
        end if
    end do

    ! 27/06/2018 - JFE new EDC 9 to check steady-state proximity
    ! see Bloom et al PNAS 2016 SI eq. S3, S4 and S5
    fin_fout_lim = 0.05

    do n = 1, nopools
        Sprox  = Fin(n) / Fout(n)
        !Sprox0 = Sprox * (mean_pools(n) / M_POOLS(1,n))
        ! JFE - 29/06/2018 perform check on first year rather than first time step
        Sprox0 = Sprox * (mean_pools(n) / mean_annual_pools(1,n))
      !  print *, n, Sprox, Sprox0
        if (abs(Sprox-Sprox0) > fin_fout_lim) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(19+n) = 0
        end if
    end do

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 35 - 40 of the PASSFAIL array

    ! All pools must confirm to the prior ranges
    do n = 1, nopools
       if ((EDC2 == 1 .or. DIAG == 1) .and. (M_POOLS(1,n) > parmax(n+npars-nopools))) then
          EDC2 = 0d0 ; EDCD%PASSFAIL(35) = 0
       end if ! prior ranges conditions
    end do ! loop pools

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0. .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2 = 0d0 ; PEDC = 0 ; EDCD%PASSFAIL(35+n) = 0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

    ! finally assess how meny EDCs passed or failed
!    num_EDC = 100
!    if (DIAG == 1) then
!       do n = 1, num_EDC
!          if (EDCD%PASSFAIL(n) == 0) then
!              EDC2 = 0d0
!          endif
!       end do
!    endif

    !JFE to print EDCs
    !do n = 1, 14
   !     print*, n, EDCD%PASSFAIL(n)
   ! end do

    !JFE - 29/06/2018 now deallocate mean_annual_pools
    deallocate(mean_annual_pools)

  end subroutine EDC2_CDEA_LU_FIRES
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_pools(pools,pool_number,averaging_period,nopools)

    ! Function calculate the mean values of model pools / states across the
    ! entire simulation run

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools          & !
                          ,pool_number      & !
                          ,averaging_period   !

    double precision,dimension(averaging_period,nopools), intent (in) :: pools

    ! declare local variables
    integer :: c

    ! initial conditions
    cal_mean_pools = 0d0

    ! loop through now
    cal_mean_pools = sum(pools(1:averaging_period,pool_number))/dble(averaging_period)

    ! ensure return command issued
    return

  end function cal_mean_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_annual_pools(pools,year,interval,averaging_period)

    ! Function calculates the mean model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: year           & ! which year are we working on
                          ,averaging_period ! number of days in analysis period

    double precision, intent(in) :: pools(averaging_period) & ! input pool state variables
                                 ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    double precision :: startday, endday

    ! calculate some constants
    startday = floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday = floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual mean values
    cal_mean_annual_pools = sum(pools(startday:endday))/(endday-startday)

    ! ensure function returns
    return

  end function cal_mean_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_max_annual_pools(pools,year,interval,averaging_period)

    ! Function calculates the max model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: year            & ! which year are we working on
                          ,averaging_period  ! number of days in analysis period

    double precision, intent(in) :: pools(averaging_period) & ! input pool state variables
                                 ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday

    ! calculate some constants
    startday = floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday = floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual max values
    cal_max_annual_pools = maxval(pools(startday:endday))

    ! ensure function returns
    return

  end function cal_max_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function expdecay2(pools,interval,averaging_period)

   ! Function to calculate the exponential decay coefficients used several EDCs.
   ! We assumpe the equation Cexp= a + b*exp(c*t)

   implicit none

   ! declare input variables
   integer, intent(in) :: averaging_period ! i.e. nodays + 1

   double precision, intent(in) :: pools(averaging_period) & ! input pool state variables
                                  ,interval((averaging_period-1))      ! model time step in decimal days

   ! declare local variables
   integer :: n, aw_int
   integer, parameter :: os = 1 ! offset days
   double precision :: aw, aw_1 &
                      ,MP0   & ! mean pool (year 1 to year end-2)
                      ,MP1   & ! mean pool (year 2 to year end-1)
                      ,MP0os & ! mean pool (year 1+os to year end-2+os)
                      ,MP1os & ! mean pool (year 2+os to year end-2+os)
                      ,dcdt1 & ! gradient of exponential over time in second year
                      ,dcdt0   ! gradient of exponential over time in first year

   ! declare initial values / constants
   aw = floor(365.25d0/(sum(interval)/dble(averaging_period-1))) ! averaging window
   aw_1 = aw ** (-1d0) ; aw_int = int(aw)
   MP0 = 0d0 ; MP1 = 0d0 ; MP0os = 0d0 ; MP1os = 0d0

   ! estimate mean stock for first year
   MP0 = sum(pools(1:aw_int))
   MP0 = MP0*aw_1

   ! estimate mean stock for second year
   MP1 = sum(pools((aw_int+1):(aw_int*2)))
   MP1 = MP1*aw_1

   ! estimate mean stock for first year with offset
   MP0os = sum(pools((1+os):(aw+os)))
   MP0os = MP0os*aw_1

   ! estimate mean stock for second year with offset
   MP1os = sum(pools((aw+os+1):((aw*2)+os)))
   MP1os = MP1os*aw_1

   ! derive mean gradient ratio (dcdt1/dcdt0)
   ! where dcdt1 is the numeric gradient between n+1 and n+365+1
   ! and dcdt0 os the numeric gradient between n and n+365
   dcdt1 = MP1os-MP0os
   dcdt0 = MP1-MP0

   ! using multiple year mean to determine c
   if ((dcdt1 > 0d0 .and. dcdt0 < 0d0) .or. (dcdt1 < 0d0 .and. dcdt0 > 0d0) &
       .or. dcdt1 == 0d0 .or. dcdt0 == 0d0) then
       ! then return error values
       expdecay2 = 1d0
   else
       expdecay2 = log(dcdt1/dcdt0) / (dble(os)*(sum(interval)/dble(averaging_period-1)))
   end if

   ! ensure return
   return

  end function expdecay2
  !
  !------------------------------------------------------------------
  !
  subroutine model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call EDC1_CDEA_LU_FIRES(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)

    endif !

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call EDC2_CDEA_LU_FIRES(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                     ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                     ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                     ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + likelihood(PI%npars,PARS)

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood_p(npars,parpriors,parpriorunc,pars)
    ! function calculates the parameter based log-likelihood for the current set
    ! of parameters. This assumes that we have any actual priors / prior
    ! uncertainties to be working with. This does include initial states, as we
    ! consider them to be parameters

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars      & ! current parameter vector
                                                     ,parpriors & ! prior values for parameters
                                                     ,parpriorunc ! prior uncertainties

    ! declare local variables
    integer :: n

    ! set initial value
    likelihood_p = 0d0

    ! now loop through defined parameters for their uncertainties
    do n = 1, npars
       ! if there is actually a value
       if (parpriors(n) > -9999d0) then
           ! uncertainty provided as +/-
           likelihood_p = likelihood_p-((pars(n)-parpriors(n))/parpriorunc(n))**2
           !likelihood_p = likelihood_p-0.5d0*((pars(n)-parpriors(n))/parpriorunc(n))**2
           ! uncertainty provided as fraction of observed value
           !likelihood_p=likelihood_p-0.5d0*((pars(n)-parpriors(n))/(parpriors(n)*parpriorunc(n)))**2
           ! uncertainty provided in log scale
           !likelihood_p=likelihood_p-0.5d0*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2
       end if
    end do

    ! apply the 0.5 multiplicative which is part of the main likelihood calculation, here once.
    likelihood_p = likelihood_p * 0.5d0

    ! dont for get to return
    return

  end function likelihood_p
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood(npars,pars)
    use cardamom_structures, only: DATAin

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y
    double precision :: tot_exp, tmp_var, infini
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    likelihood = 0d0 ; infini = 0d0

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
       likelihood = likelihood-tot_exp
    endif

    ! LAI log-likelihood
    if (DATAin%nlai > 0) then
       ! loop split to allow vectorisation
       tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (DATAin%M_LAI(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp = tot_exp+(-log(infini))
         endif
       end do
       likelihood = likelihood-tot_exp
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       likelihood = likelihood-tot_exp
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cwood increment log-likelihood
    if (DATAin%nwoo > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nwoo
         dn = DATAin%woopts(n)
         ! note that division is the uncertainty
         ! tot_exp = tot_exp+(log((DATAin%M_POOLS(dn,4)-DATAin%M_POOLS(dn-365,4)) &
         !                   / DATAin%WOO(dn))/log(DATAin%WOO_unc(dn)))**2
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%M_POOLS(dn-365,4)) / DATAin%WOO_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCfol_stock
         dn = DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp = tot_exp+(log(DATAin%M_POOLS(dn,2)/DATAin%Cfol_stock(dn))/log(2.))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Cfol_stock(dn)) / DATAin%Cfol_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       no_years = int(nint(sum(DATAin%deltat)/365.25d0))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),2),y,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn = DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25d0 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_stock
         dn = DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,4)/DATAin%Cwood_stock(dn))/log(2.))**2.
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/DATAin%Cwood_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCroots_stock
         dn = DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,3)/DATAin%Croots_stock(dn))/log(2.))**2.
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Croots_stock(dn)) / DATAin%Croots_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nClit_stock
         dn = DATAin%Clit_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+((log((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
!                           *DATAin%M_POOLS(dn,5))/DATAin%Clit_stock(dn))/log(2.))**2d0
         tot_exp = tot_exp+(((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                           *(DATAin%M_POOLS(dn,5))-DATAin%Clit_stock(dn))/DATAin%Clit_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCsom_stock
         dn = DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,6)/DATAin%Csom_stock(dn))/log(2.))**2.
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,6)-DATAin%Csom_stock(dn))/DATAin%Csom_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > 0) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        likelihood = likelihood-((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2d0
    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk liklihood
    ! hear
    likelihood = likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (likelihood /= likelihood) then
       likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function likelihood
  !
  !------------------------------------------------------------------
  !
end module model_likelihood_module
