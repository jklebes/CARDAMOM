
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
    MCO%nADAPT = 500
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
    call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

    ! assess post running EDCs
    call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
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
    ML_obs_out = -0.5d0*(tot_exp*10d0)*DATAin%EDC

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine EDC1_GSI(PARS, npars, meantemp, meanrad, EDC1)

    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: sw_par_fraction, &
                                opt_max_scaling, &
                                     emissivity

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (modified from Bloom et al., 2014).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local parameters
    double precision, dimension(10), parameter :: lai = (/0.000001d0,0.03125d0,0.0625d0,0.125d0,0.25d0,0.5d0,1d0,2.5d0,5d0,10d0/)

    ! declare local variables
    integer :: n, DIAG, i
    double precision :: tmp, tmp1, tmp2 &
             ,par_trans, par_refl, nir_trans, nir_refl, lw_trans, lw_refl &
             ,fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    ! set initial value
    EDC1 = 1
    DIAG = EDCD%DIAG


    ! set all EDCs to 1 (pass)
    EDCD%nedc = 100
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    do i = 1,10

       ! Canopy transmitted of PAR & NIR radiation towards the soil
       par_trans = 1d0 - (lai(i)*pars(16)) / (lai(i)+pars(17))
       nir_trans = 1d0 - (lai(i)*pars(18)) / (lai(i)+pars(19))
       lw_trans = 1d0 - (lai(i)*pars(6)) / (lai(i)+pars(7))
       ! Canopy reflected of near infrared and photosynthetically active radiation
       nir_refl = (lai(i)*pars(8)) / (lai(i)+pars(9))
       par_refl = (lai(i)*pars(11)) / (lai(i)+pars(12))
       lw_refl = (lai(i)*pars(22)) / (lai(i)+pars(13))

       ! the transmittance and reflection of LW radiation should be less than 1
       if ((EDC1 == 1 .or. DIAG == 1) .and. lw_trans+lw_refl >= 1d0) then
           EDC1 = 0 ; EDCD%PASSFAIL(1) = 0
       endif
       ! the transmittance of NIR should always be > than PAR at all LAI values
       if ((EDC1 == 1 .or. DIAG == 1) .and. nir_trans < par_trans) then
           EDC1 = 0 ; EDCD%PASSFAIL(2) = 0
       endif
       ! reflectance and transmittance of NIR should always be < 1
       if ((EDC1 == 1 .or. DIAG == 1) .and. nir_trans + nir_refl >= 1d0 ) then
           EDC1 = 0 ; EDCD%PASSFAIL(3) = 0
       endif
       ! reflectance and transmittance of PAR should always be < 1
       if ((EDC1 == 1 .or. DIAG == 1) .and. par_trans + par_refl >= 1d0 ) then
           EDC1 = 0 ; EDCD%PASSFAIL(4) = 0
       endif
       ! reflectance should be greater for NIR than for PAR
       if ((EDC1 == 1 .or. DIAG == 1) .and. nir_refl < par_refl) then
           EDC1 = 0 ; EDCD%PASSFAIL(5) = 0
       endif
    enddo

    ! maximum temperature for photosythesis cannot be smaller than optimum
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(3) > pars(2)) then
        EDC1 = 0 ; EDCD%PASSFAIL(6) = 0
    endif

    ! assume that photosynthesis limitation at 0C should be between 10 % and 20 %
    ! of potential. Fatchi et al (2013), New Phytologist, https://doi.org/10.1111/nph.12614
    tmp = opt_max_scaling(pars(2),pars(3),pars(4),0d0)
    if ((EDC1 == 1 .or. DIAG == 1) .and. tmp > 0.20d0) then
       EDC1 = 0 ; EDCD%PASSFAIL(7) = 0
    endif

  end subroutine EDC1_GSI
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_GSI(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    use cardamom_structures, only: DATAin

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
    integer :: n, DIAG, no_years, y, PEDC, nn, num_EDC, max_location(1),i
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF &
                       ,model_living_C, target_living_C(2),hold
    double precision, dimension(:), allocatable :: mean_annual_pools,tmp
    double precision :: max_wood & !
                       ,torfol   & ! yearly average turnover
                       ,torlab   & !
                       ,fauto    & ! Fractions of GPP to autotrophic respiration
                       ,ffol     & ! Fraction of GPP to foliage
                       ,flab     & ! Fraction of GPP to labile pool
                       ,froot    & ! Fraction of GPP to root
                       ,fwood    & ! Fraction of GPP to wood
                       ,fsom     & ! fraction of GPP som under eqilibrium conditions
                       ,flit     & ! fraction of GPP to litter under equilibrium condition
                       ,delta_gsi

    ! set initial values
    DIAG=EDCD%DIAG
    EDC2=1

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! the maximum value for all fluxes must be greater than zero
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,1)) < 0.5d0) then
         EDC2 = 0 ; EDCD%PASSFAIL(8) = 0
    endif
    ! the maximum value for all fluxes must be greater than zero
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,2)) < 0.5d0) then
         EDC2 = 0 ; EDCD%PASSFAIL(9) = 0
    endif
    ! the maximum value for all fluxes must be greater than zero
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,3)) < 0.5d0) then
         EDC2 = 0 ; EDCD%PASSFAIL(10) = 0
    endif
    ! the maximum value for all fluxes must be greater than zero
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,4)) < 0.5d0) then
         EDC2 = 0 ; EDCD%PASSFAIL(11) = 0
    endif

    ! additional faults can be stored in locations 35 - 40 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0d0 .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2 = 0 ; PEDC = 0 ; EDCD%PASSFAIL(35+n) = 0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

    ! ensure FLUXES values are /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nofluxes .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_FLUXES(nn,n) /= M_FLUXES(nn,n)) then
                 EDC2 = 0 ; PEDC = 0 ; EDCD%PASSFAIL(35+nopools+n) = 0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

  end subroutine EDC2_GSI
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
        call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

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
        call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
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
           if (n >= 1 .or. n <= 19) then
               ! Gaussian uncertainties for CO2 compensation and half saturation
               ! points
               likelihood_p = likelihood_p-0.5d0*((pars(n)-parpriors(n))/parpriorunc(n))**2
           else
               likelihood_p = likelihood_p-0.5d0*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2
           endif
       end if
    end do

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
    double precision :: tot_exp, pool_dynamics, tmp_var, infini
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    likelihood = 0d0 ; infini = 0d0

    ! GPP Log-likelihood
    tot_exp = 0d0
    if (DATAin%ngpp > 0) then
        do n = 1, DATAin%ngpp
          dn = DATAin%gpppts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_GPP(dn)-DATAin%GPP(dn))/DATAin%GPP_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

    ! Evapotranspiration (kgH2O.m-2.day-1) Log-likelihood
    ! in this case transpiration only
    tot_exp = 0d0
    if (DATAin%nEvap > 0) then
        do n = 1, DATAin%nEvap
          dn = DATAin%Evappts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,2)-DATAin%Evap(dn))/DATAin%Evap_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

    ! Borrowed wood increment to provide soil evaporation for ACM recal  (kgH2O.m-2.day-1) Log-likelihood
    tot_exp = 0d0
    if (DATAin%nwoo > 0) then
        do n = 1, DATAin%nwoo
          dn = DATAin%woopts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,3)-DATAin%woo(dn))/DATAin%woo_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

    ! Borrowed Cfol_stock to provide wet canopy evaporation for ACM recal  (kgH2O.m-2.day-1) Log-likelihood
    tot_exp = 0d0
    if (DATAin%nCfol_stock > 0) then
        do n = 1, DATAin%nCfol_stock
          dn = DATAin%Cfol_stockpts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,4)-DATAin%Cfol_stock(dn))/DATAin%Cfol_stock_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

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
