
module model_likelihood_module
  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! which to make open
  public :: model_likelihood, find_edc_initial_values, &
            sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: EDC
    integer :: DIAG
    integer :: PASSFAIL(100) ! allow space for 100 possible checks
    integer :: nedc ! number of edcs being assessed
  end type
  type (EDCDIAGNOSTICS), save :: EDCD

  ! Has the model sanity check been conducted yet?
  logical :: sanity_check = .false.

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
    integer :: n, counter_local, EDC_iter, nOUT_save, nWRITE_save, nADAPT_save, append_save
    double precision :: PEDC, PEDC_prev, ML, ML_prior, P_target
    double precision, dimension(PI%npars+1) :: EDC_pars

    ! Hold for later
    nOUT_save = MCO%nOUT ; nWRITE_save = MCO%nWRITE ; nADAPT_save = MCO%nADAPT
    append_save = MCO%append

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
    MCO%fixedpars  = .true. ! TLS: changed from .false. for testing 16/12/2019

    ! Set initial priors to vector...
    PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
    ! ... and assume we need to find random parameters
    PI%parfix = 0
    ! Target likelihood allows for controlling when the MCMC will stop
    P_target = 0d0

    ! if the prior is not missing and we have not told the edc to be random
    ! keep the value
!    do n = 1, PI%npars
!       if (PI%parini(n) /= -9999d0 .and. DATAin%edc_random_search < 1) PI%parfix(n) = 1
!    end do ! parameter loop

    ! set the parameter step size at the beginning
    PI%parvar = 1d0 ; PI%Nparvar = 0d0
    PI%use_multivariate = .false.
    ! Covariance matrix cannot be set to zero therefore set initial value to a
    ! small positive value along to variance access
    PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
    do n = 1, PI%npars
       PI%covariance(n,n) = 1d0
    end do

    ! if this is not a restart run, i.e. we do not already have a starting
    ! position we must being the EDC search procedure to find an ecologically
    ! consistent initial parameter set
    if (.not. restart_flag) then

        ! set up edc log likelihood for MHMCMC initial run
        PEDC_prev = -1000d0 ; PEDC = -1d0 ; counter_local = 0
        do while (PEDC < 0d0)

           write(*,*)"Beginning EDC search attempt"
           ! call the MHMCMC directing to the appropriate likelihood
           call MHMCMC(P_target,model_likelihood,edc_model_likelihood)

           ! store the best parameters from that loop
           PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
           ! turn off random selection for initial values
           MCO%randparini = .false.
           write(*,*)"...intermediate EDC search progress check"
           ! call edc likelihood function to get final edc probability
           call edc_model_likelihood(PI%parini,PEDC,ML_prior)

           ! keep track of attempts
           counter_local = counter_local + 1
           ! periodically reset the initial conditions
           if (PEDC < 0d0 .and. PEDC <= PEDC_prev .and. counter_local > 5) then
               ! Reset the previous EDC likelihood score
               PEDC_prev = -1000d0
               ! Reset parameters back to default
               PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
               ! reset to select random starting point
               MCO%randparini = .true.
               ! reset the parameter step size at the beginning of each attempt
               PI%parvar = 1d0 ; PI%Nparvar = 0d0
               ! Covariance matrix cannot be set to zero therefore set initial value to a
               ! small positive value along to variance access
               PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
               PI%use_multivariate = .false.
               do n = 1, PI%npars
                  PI%covariance(n,n) = 1d0
               end do
           else
               PEDC_prev = PEDC
           endif

        end do ! for while condition

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

    ! Perform a more aggressive sanity check which compares the bulk difference
    ! in all fluxes and pools from multiple runs of the same parameter set
    if (.not.sanity_check) call model_sanity_check(PI%parini)

    ! call EDCs which can be evaluated prior to running the model
    call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)
!print*,"edc_model_likelihood: carbon_model done"
    ! assess post running EDCs
    call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                    ,DATAin%nodays,DATAin%deltat,DATAin%steps_per_year   &
                    ,PI%parmax,PARS,DATAin%MET &
                    ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                    ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

    ! calculate the likelihood
    tot_exp = sum(1d0-EDCD%PASSFAIL(1:EDCD%nedc))
!    tot_exp = 0d0
!    do n = 1, EDCD%nedc
!       tot_exp=tot_exp+(1d0-EDCD%PASSFAIL(n))
!       if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
!    end do ! checking EDCs
!   ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) then
!        print*,"Found it" ; stop
!    endif

    ! convert to a probability
!    ML_obs_out = -0.5d0*(tot_exp*10d0)*DATAin%EDC
    ML_obs_out = -5d0*tot_exp*DATAin%EDC

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine sub_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible for running the model,
    ! calculation of the log-likelihood on a subsample of observation
    ! for comparison assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

!    ! Debugging print statements
!    print*,"sub_model_likelihood:"

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

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
!print*,"sub_model_likelihood: carbon_model done"
    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%deltat,DATAin%steps_per_year     &
                        ,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + scale_likelihood(PI%npars,PARS)

!    ! Debugging print statements
!    print*,"sub_model_likelihood: done"

  end subroutine sub_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine sqrt_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible for running the model,
    ! calculation of the log-likelihood on a subsample of observation
    ! for comparison assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

!    ! Debugging print statements
!    print*,"sub_model_likelihood:"

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

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
!print*,"sqrt_model_likelihood: carbon_model done"
    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%deltat,DATAin%steps_per_year     &
                        ,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + sqrt_scale_likelihood(PI%npars,PARS)

!    ! Debugging print statements
!    print*,"sqrt_model_likelihood: done"

  end subroutine sqrt_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine log_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible for running the model,
    ! calculation of the log-likelihood on a subsample of observation
    ! for comparison assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

!    ! Debugging print statements
!    print*,"sub_model_likelihood:"

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

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
!print*,"log_model_likelihood: carbon_model done"
    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%deltat,DATAin%steps_per_year     &
                        ,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + log_scale_likelihood(PI%npars,PARS)

!    ! Debugging print statements
!    print*,"log_model_likelihood: done"

  end subroutine log_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine model_sanity_check(PARS)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use CARBON_MODEL_MOD, only: carbon_model

    ! Carries out multiple carbon model iterations using the same parameter set
    ! to ensure that model outputs are consistent between iterations, i.e. that
    ! the model is numerically secure. Reproducible outputs from the models is
    ! essential for successful mcmc anlaysis

    implicit none

    ! Arguments
    double precision, dimension(PI%npars), intent(in) :: PARS

    ! Local arguments
    integer :: i
    double precision, dimension((DATAin%nodays+1),DATAin%nopools) :: local_pools
    double precision, dimension(DATAin%nodays,DATAin%nofluxes) :: local_fluxes
    double precision :: pool_error, flux_error

    ! Run model

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)
!print*,"sanity_check: carbon_model done 1"
    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,local_fluxes,local_pools,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)
!print*,"sanity_check: carbon_model done 2"
    ! Compare outputs
    flux_error = sum(abs(DATAin%M_FLUXES - local_fluxes))
    pool_error = sum(abs(DATAin%M_POOLS - local_pools))
    ! If error between runs exceeds precision error then we have a problem
    if (pool_error > (tiny(0d0)*(DATAin%nopools*DATAin%nodays)) .or. &
        flux_error > (tiny(0d0)*(DATAin%nofluxes*DATAin%nodays)) .or. &
        pool_error /= pool_error .or. flux_error /= flux_error) then
        print*,"Error: multiple runs of the same parameter set indicates an error"
        print*,"Cumulative POOL error = ",pool_error
        print*,"Cumulative FLUX error = ",flux_error
        do i = 1,DATAin%nofluxes
           print*,"Sum abs error over time: flux = ",i
           print*,sum(abs(DATAin%M_FLUXES(:,i) - local_fluxes(:,i)))
        end do
        do i = 1, DATAin%nopools
           print*,"Sum abs error over time: pool = ",i
           print*,sum(abs(DATAin%M_POOLS(:,i) - local_pools(:,i)))
        end do
        stop
    end if

    ! Update the user
    print*,"Sanity check completed"

    ! Set Sanity check as completed
    sanity_check = .true.

  end subroutine model_sanity_check
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC1(PARS, npars, meantemp, meanrad, EDC1)
 
    use cardamom_structures, only: DATAin

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
    double precision :: fauto, & ! Fractions of GPP to autotrophic respiration
                         ffol, & ! Fraction of GPP to foliage
                         flab, & ! Fraction of GPP to labile pool
                        froot, & ! Fraction of GPP to root
                        fwood, & ! Fraction of GPP to wood
                         fsom, & ! fraction of GPP som under eqilibrium conditions
                temp_response, &
                       torfol    ! yearly leaf loss fraction

    ! set initial value
    EDC1 = 1
    DIAG = EDCD%DIAG

    ! yearly leaf loss fraction
    torfol = 1d0/(pars(5)*365.25d0)

    ! set all EDCs to 1 (pass)
    EDCD%nedc = 100
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! EDCs from Python code
    ! No individual daily rate of GPP should be > 25 gC/m2/day
    ! Mean annual GPP should be 500 < gpp < 2800 gC/m2/yr
    ! No individual daily rate of Reco should be > 20 gC/m2/day
    ! Mean annual Reco should be not less than 500 gC/m2/yr
    ! Mean annual Reco should be not more than 2600 gC/m2/yr
    ! Soil C shuld not change by more than 5 % over life of analysis
    ! or ( (rem[0,:]*21/float(650*0.035) > 70).any() ) # max total LSU_ha_week
    ! All specified cuts in input file are imposed

    ! calculate temperature response of decomposition processes
    temp_response = exp(pars(9)*meantemp) !Q10 
    
    ! Straight forward GSI parameter bounds
    ! Temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(13) < pars(12))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    end if
    ! Photoperiod
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(20) < pars(14))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    end if
    ! VPD - may not be needed
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(22) < pars(21))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    end if

    ! Photoperiod minimum cannot be substantially less than the observed minimum day length
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(14) < minval(DATAin%MET(11,:))-14400d0)) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    end if
    ! Photoperiod maximum cannot be greater than the observed maximum day length
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(20) > maxval(DATAin%MET(11,:)))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    end if

    ! VPD at which stress in at maximum should be no larger than max(VPDlag21) +
    ! 1500 Pa from the max VPD tolerated parameter
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(22) > maxval(DATAin%MET(12,:))+1500d0)) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
    end if

    ! Rhet from litter should be faster than Rhet from som
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(8) > pars(7) ) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
    endif

    ! Turnover of litter towards som (pars(1)) should be faster than turnover of som (pars(8))
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(8) > pars(1) ) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
    endif

    ! root turnover (pars(6)) should be greater than som turnover (pars(8)) at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(8)*temp_response) > pars(6)) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(9) = 0
    endif

    ! IMPLICIT Combustion completeness for foliage should be greater than soil
    ! IMPLICIT Combustion completeness for fol+root litter should be greater than soil

!    ! Combustion completeness for foliage should be greater than non-photosynthetic tissues
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(29) < pars(30)) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
!    endif
!    ! Combustion completeness for non-photosynthetic tissue should be greater than soil
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(30) < pars(31)) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
!    endif
!    ! Combustion completeness for foliar + fine root litter should be greater than non-photosynthetic tissue
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(32) < pars(30)) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
!    endif

    ! could always add more / remove some

  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,deltat,steps_per_year &
                        ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                        ,meantemp,EDC2)
    use cardamom_structures, only: DATAin

    ! Determines whether the dynamical contraints for the search of the initial
    ! parameters has been successful or whether or not we should abandon the
    ! current set and move on

    implicit none

    ! declare input variables
    integer, intent(in) :: npars          & ! number of model parameters
                          ,nomet          & ! number of met drivers
                          ,nofluxes       & ! number of fluxes from model
                          ,nopools        & ! number of pools in model
                          ,nodays         & ! number of days in simulation
                          ,steps_per_year

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
    integer :: n, nn, nnn, DIAG, y, PEDC, steps_per_month, nd, fl, &
               io_start, io_finish
    double precision :: infi !, EQF, etol
    double precision, dimension(nodays) :: lab_ratio
    double precision, dimension(nopools) :: jan_mean_pools, jan_first_pools, &
                                            mean_pools, Fin, Fout, &
                                            Fin_yr1, Fout_yr1
    double precision, dimension(nofluxes) :: FT, FT_yr1
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
                       ,ffol  & ! Fraction of GPP to foliage
                       ,flab  & ! Fraction of GPP to labile pool
                       ,froot & ! Fraction of GPP to root
                       ,fwood   ! Fraction of GPP to wood

    ! Steady State Attractor:
    ! Log ratio difference between inputs and outputs of the system.
    double precision, parameter :: EQF1_5 = log(1.5d0), & ! 10.0 = order magnitude; 2 = double and half
                                   EQF2 = log(2d0),   & ! 10.0 = order magnitude; 2 = double and half
                                   EQF5 = log(5d0),   &
                                   EQF10 = log(10d0), &
                                   EQF15 = log(15d0), &
                                   EQF20 = log(20d0), &
                                    etol = 0.05d0 ! 0.20d0 lots of AGB !0.10d0 global / site more data !0.05d0 global 1 or 2 AGB estimates

!    ! Debugging print statements
!    print*,"assess_EDC2: "

    ! update initial values
    DIAG = EDCD%DIAG
    EDC2 = 1
    infi = 0d0

   ! derive mean pools
    do n = 1, nopools
       mean_pools(n) = cal_mean_pools(M_POOLS,n,nodays+1,nopools)
    end do
    
    !!!!!!!!!
    ! Assess emergent ratios for ecological consistency
    !!!!!!!!!
    
    ! Ensure ratio between Cfoliar and Croot is less than 5
    if ((EDC2 == 1 .or. DIAG == 1) .and. &
        (mean_pools(2) > (mean_pools(3)*5d0) .or. (mean_pools(2)*5d0) < mean_pools(3)) ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(10) = 0
    end if

    ! We would not expect that the mean labile stock is greater than
    ! 8 % of the total ecosystem carbon stock, as we need structure to store
    ! labile.
    ! Gough et al (2009) Agricultural and Forest Meteorology. Avg 11, 12.5, 3 %
    ! (Max across species for branch, bole and coarse roots). Provides evidence that
    ! branches accumulate labile C prior to bud burst from other areas.
    ! Wurth et al (2005) Oecologia, Clab 8 % of living biomass (DM) in tropical forest
    ! Richardson et al (2013), New Phytologist, Clab 2.24 +/- 0.44 % in temperate (max = 4.2 %)
    ! Estimate the labile ratio, also used below
    lab_ratio = M_POOLS(1:nodays,1) / (M_POOLS(1:nodays,1) + M_POOLS(1:nodays,2) + M_POOLS(1:nodays,3))
    if (EDC2 == 1 .or. DIAG == 1) then
        if (maxval(lab_ratio) > 0.125d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(11) = 0
        endif
    endif ! EDC2 == 1 .or. DIAG == 1
    if (EDC2 == 1 .or. DIAG == 1) then
        if (sum(lab_ratio)/dble(nodays) > 0.08d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(12) = 0
        endif
    endif ! EDC2 == 1 .or. DIAG == 1
    
    ! GSI model should reach above 0.5 at least once during the analysis
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,18)) < 0.5d0) then 
        EDC2 = 0d0 ; EDCD%PASSFAIL(13) = 0
    end if ! EDC2 == 1 .or. DIAG == 1
    
    ! First calculate total flux for the simulation period
    io_start = (steps_per_year*2) + 1 ; io_finish = nodays
    if (DATAin%nos_years < 3) io_start = 1
    do fl = 1, nofluxes
       FT(fl) = sum(M_FLUXES(io_start:io_finish,fl)*deltat(io_start:io_finish))
       FT_yr1(fl) = sum(M_FLUXES(1:steps_per_year,fl)*deltat(1:steps_per_year))
       !FT_yr2(fl) = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),fl) &
       !                *deltat((steps_per_year+1):(steps_per_year*2)))
    end do

    ! get total in and out for each pool
    ! labile
    Fin(1)  = FT(5)
    Fout(1) = FT(7)+FT(25)+FT(28)+FT(31)+FT(34)
    Fin_yr1(1)  = FT_yr1(5)
    Fout_yr1(1) = FT_yr1(7)+FT_yr1(25)+FT_yr1(28)+FT_yr1(31)+FT_yr1(34)
    ! roots
    Fin(3)  = FT(6)
    Fout(3) = FT(10)+FT(27)+FT(30)+FT(33)+FT(36)
    Fin_yr1(3)  = FT_yr1(6)
    Fout_yr1(3) = FT_yr1(10)+FT_yr1(27)+FT_yr1(30)+FT_yr1(33)+FT_yr1(36)
    ! litter
    Fin(4)  = FT(9)+FT(10)+FT(19)+FT(28)+FT(29)+FT(30)+FT(34)+FT(35)+FT(36)
    Fout(4) = FT(11)+FT(13)
    Fin_yr1(4)  = FT_yr1(9)+FT_yr1(10)+FT_yr1(19)+FT_yr1(28)+FT_yr1(29)+FT_yr1(30)+FT_yr1(34)+FT_yr1(35)+FT_yr1(36)
    Fout_yr1(4) = FT_yr1(11)+FT_yr1(13)
    ! som
    Fin(5)  = FT(13)
    Fout(5) = FT(12)
    Fin_yr1(5)  = FT_yr1(13)
    Fout_yr1(5) = FT_yr1(12)

    if (EDC2 == 1 .or. DIAG == 1) then

        ! Labile pool
        n = 1
        ! Restrict mean rates of increase
        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(13+n-1) = 0
        end if
        ! Restrict exponential behaviour at initialisation
        if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                  abs(log(Fin(n)/Fout(n))) ) > etol ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
        end if

        ! Roots pool
        n = 3
        ! Restrict mean rates of increase
        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(13+n-1) = 0
        end if
        ! Restrict exponential behaviour at initialisation
        if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                  abs(log(Fin(n)/Fout(n))) ) > etol ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
        end if

        ! Dead pools
        do n = 4, 5
           ! Restrict rates of increase
           if (abs(log(Fin(n)/Fout(n))) > EQF2) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(13+n-1) = 0
           end if
           ! Restrict exponential behaviour at initialisation
           if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                     abs(log(Fin(n)/Fout(n))) ) > etol ) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
           end if
        end do

    end if ! EDC2 == 1 .or. DIAG == 1

    !!!!!!!!!
    ! Currently empty, probably should contain something here....
    !!!!!!!!!

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 55 - 61 of the PASSFAIL array

    ! ensure minimum pool values are >= 0, /= NaN or Inf
    if (EDC2 == 1 .or. DIAG == 1) then

       do n = 1, nopools
          if (minval(M_POOLS(1:nodays,n)) < 0d0 .or. &
              maxval(abs(M_POOLS(1:nodays,n))) == abs(log(infi)) .or. &
              minval(M_POOLS(1:nodays,n)) /= minval(M_POOLS(1:nodays,n))) then
              EDC2 = 0d0 ; EDCD%PASSFAIL(55+n) = 0
          endif
       end do

       do n = 1, nofluxes
          if (maxval(abs(M_FLUXES(:,n))) == abs(log(infi)) .or. &
              minval(M_FLUXES(:,n)) /= minval(M_FLUXES(:,n))) then
              EDC2 = 0d0 ; EDCD%PASSFAIL(55+nopools+n) = 0
          endif
       end do

    end if ! min pool assessment

!    ! Debugging print statements
!    print*,"assess_EDC2: done"

  end subroutine assess_EDC2
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
    integer :: startday, endday

    ! calculate some constants
    startday = floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday = floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual mean values
    cal_mean_annual_pools = sum(pools(startday:endday))/dble(endday-startday)

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
   MP0os = sum(pools((1+os):(aw_int+os)))
   MP0os = MP0os*aw_1

   ! estimate mean stock for second year with offset
   MP1os = sum(pools((aw_int+os+1):((aw_int*2)+os)))
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

!    ! Debugging print statements
!    print*,"model_likelihood:"

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

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
!print*,"model_likelihood: carbon_model done"
    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then
!print*,"model_likelihood: EDC2"
        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%deltat,DATAin%steps_per_year     &
                        ,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)
!print*,"model_likelihood: EDC2 done"
        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
!print*,"model_likelihood: EDC2 likelihood update"
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
!print*,"model_likelihood: update likelihood score"
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + likelihood(PI%npars,PARS)
!print*,"model_likelihood: update likelihood score done"
!    ! Debugging print statements
!    print*,"model_likelihood: done"

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood_p(npars,parpriors,parpriorunc,parpriorweight,pars)
    ! function calculates the parameter based log-likelihood for the current set
    ! of parameters. This assumes that we have any actual priors / prior
    ! uncertainties to be working with. This does include initial states, as we
    ! consider them to be parameters

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars         & ! current parameter vector
                                                     ,parpriors    & ! prior values for parameters
                                                     ,parpriorunc  & ! prior uncertainties
                                                     ,parpriorweight ! prior weighting

    ! declare local variables
    integer :: n
    double precision, dimension(npars) :: local_likelihood
!print*,"likelihood_p:"
    ! set initial value
    likelihood_p = 0d0 ; local_likelihood = 0d0

    ! now loop through defined parameters for their uncertainties
    where (parpriors > -9999) local_likelihood = parpriorweight*((pars-parpriors)/parpriorunc)**2
    likelihood_p = sum(local_likelihood) * (-0.5d0)
!print*,"likelihood_p: done"
    ! dont for get to return
    return

  end function likelihood_p
  !------------------------------------------------------------------
  !
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
    integer :: n, dn, y, s, f
    double precision :: tot_exp, tmp_var, infini, input, output, obs, model, unc
    double precision, dimension(DATAin%nodays) :: mid_state
    double precision, dimension(DATAin%steps_per_year) :: sub_time
    double precision, allocatable :: mean_annual_pools(:)

!    ! Debugging print statement
!    print*,"likelihood: "

    ! initial value
    likelihood = 0d0 ; infini = 0d0 ; mid_state = 0d0 ; sub_time = 0d0

!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       likelihood = likelihood-tot_exp
!    endif
    ! NBE Log-likelihood
    ! NBE partitioned between the mean flux and seasonal anomalies
    if (DATAin%nnbe > 0) then
!        ! Determine the mean value for model, observtion and uncertainty estimates
!        obs   = sum(DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!              / dble(DATAin%nnbe)
!        model = sum(DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+ &
!                    DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!              / dble(DATAin%nnbe)
!        unc   = sqrt(sum(DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe))**2)) &
!              / dble(DATAin%nnbe)
!        ! Update the likelihood score with the mean bias
!        likelihood = likelihood - (((model - obs) / unc) ** 2)
!        ! Determine the anomalies based on substraction of the global mean
!        ! Greater information would come from breaking this down into annual estimates
!        tot_exp = sum( (((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)-model) - &
!                        (DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe)) - obs)) / unc)**2 )
!        likelihood = likelihood-tot_exp
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%NBE(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean NBE from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%NBE(s:f)*sub_time) / sum(sub_time)
               model = sum((DATAin%M_NEE(s:f)+ &
                            DATAin%M_FLUXES(s:f,24))*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%NBE_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               likelihood = likelihood - (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               likelihood = likelihood - sum( (sub_time*(((DATAin%M_NEE(s:f)+DATAin%M_FLUXES(s:f,24)-model) - &
                                                         (DATAin%NBE(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
    endif ! nnbe > 0

!    ! GPP Log-likelihood
!    if (DATAin%ngpp > 0) then
!       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
!                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
!       likelihood = likelihood-tot_exp
!    endif
    ! GPP Log-likelihood
    ! GPP partitioned between the mean flux and seasonal anomalies
    if (DATAin%ngpp > 0) then
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%GPP(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean GPP from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%GPP(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,1)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%GPP_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               likelihood = likelihood - (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               likelihood = likelihood - sum( (sub_time*(((DATAin%M_FLUXES(s:f,1)-model) - &
                                                          (DATAin%GPP(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
    endif ! ngpp > 0

!    ! Evap Log-likelihood
!    if (DATAin%nEvap > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),29)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
!                       /DATAin%Evap_unc(DATAin%evappts(1:DATAin%nEvap)))**2)
!       likelihood = likelihood-tot_exp
!    endif

!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),24)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       likelihood = likelihood-tot_exp
!    endif
    ! Fire Log-likelihood
    ! Fire partitioned between the mean flux and seasonal anomalies
    if (DATAin%nFire > 0) then
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%Fire(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean fire from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%Fire(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,24)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%Fire_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               likelihood = likelihood - (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               likelihood = likelihood - sum( (sub_time*(((DATAin%M_FLUXES(s:f,24)-model) - &
                                                          (DATAin%Fire(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
    endif ! nFire > 0

    ! Assume physical property is best represented as the mean of value at beginning and end of times step
    if (DATAin%nlai > 0) then
       ! Create vector of (LAI_t0 + LAI_t1) * 0.5, note / pars(17) to convert foliage C to LAI
       mid_state = ( ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                   * 0.5d0 ) / pars(15)
       ! Split loop to allow vectorisation
       tot_exp = sum(((mid_state(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       ! loop split to allow vectorisation
       !tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
       !                /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (mid_state(dn) < 0d0) then
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

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       ! Create vector of (FOL_t0 + FOL_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)) &
                       -DATAin%Cfol_stock(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))&
                     / DATAin%Cfol_stock_unc(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))**2)
       ! Sum with current likelihood score
       likelihood = likelihood-tot_exp
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(DATAin%nos_years))
       ! determine the annual max for each pool
       do y = 1, DATAin%nos_years
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

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       ! Create vector of (root_t0 + root_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,3) + DATAin%M_POOLS(2:(DATAin%nodays+1),3) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)) &
                       -DATAin%Croots_stock(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))&
                     / DATAin%Croots_stock_unc(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       ! Create vector of (lit_t0 + lit_t1) * 0.5
       !mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
       !          * 0.5d0
       mid_state = sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,28)+DATAin%M_FLUXES(:,29)+ & 
                       DATAin%M_FLUXES(:,34)+DATAin%M_FLUXES(:,35))
       mid_state = (mid_state / &
                   (mid_state + sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,30)+DATAin%M_FLUXES(:,36)))) & 
                 * DATAin%M_POOLS(:,4)
       !mid_state = (sum(DATAin%M_FLUXES(:,9))/sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,10))) &
       !          * DATAin%M_POOLS(:,4)
       mid_state = (mid_state(1:DATAin%nodays) + mid_state(2:(DATAin%nodays+1))) * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Clit_stockpts(1:DATAin%nClit_stock)) &
                       -DATAin%Clit_stock(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))&
                     / DATAin%Clit_stock_unc(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       ! Create vector of (som_t0 + som_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,5) + DATAin%M_POOLS(2:(DATAin%nodays+1),5) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)) &
                       -DATAin%Csom_stock(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))&
                     / DATAin%Csom_stock_unc(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

!    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s ->
!    ! kg/m2/day)
!    if (DATAin%otherpriors(4) > -9998) then
!        tot_exp = sum(DATAin%M_FLUXES(:,29)) / sum(DATAin%MET(7,:) * 86400d0)
!        likelihood = likelihood-((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
!    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
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
  double precision function scale_likelihood(npars,pars)
    use cardamom_structures, only: DATAin

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, y, s, f
    double precision :: tot_exp, tmp_var, infini, input, output, model, obs, unc
    double precision, dimension(DATAin%nodays) :: mid_state
    double precision, dimension(DATAin%steps_per_year) :: sub_time
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    scale_likelihood = 0d0 ; infini = 0d0 ; mid_state = 0d0 ; sub_time = 0d0

!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnbe))
!    endif
    ! NBE Log-likelihood
    ! NBE partitioned between the mean flux and seasonal anomalies
    if (DATAin%nnbe > 0) then
!        ! Determine the mean value for model, observtion and uncertainty estimates
!        obs   = sum(DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!              / dble(DATAin%nnbe)
!        model = sum(DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+ &
!                    DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!              / dble(DATAin%nnbe)
!        unc   = sqrt(sum(DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe))**2)) &
!              / dble(DATAin%nnbe)
!        ! Update the likelihood score with the mean bias
!        scale_likelihood = scale_likelihood - (((model - obs) / unc) ** 2)
!        ! Determine the anomalies based on substraction of the global mean
!        ! Greater information would come from breaking this down into annual estimates
!        tot_exp = sum( (((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)-model) - &
!                        (DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe)) - obs)) / unc)**2 )
!        likelihood = likelihood-tot_exp
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%NBE(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean NBE from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%NBE(s:f)*sub_time) / sum(sub_time)
               model = sum((DATAin%M_NEE(s:f)+ &
                            DATAin%M_FLUXES(s:f,24))*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%NBE_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_NEE(s:f)+DATAin%M_FLUXES(s:f,24)-model) - &
                                                    (DATAin%NBE(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnbe))
    endif ! nnbe > 0

!    ! GPP Log-likelihood
!    if (DATAin%ngpp > 0) then
!       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
!                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%ngpp))
!    endif
    ! GPP Log-likelihood
    ! GPP partitioned between the mean flux and seasonal anomalies
    if (DATAin%ngpp > 0) then
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%GPP(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean GPP from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%GPP(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,1)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%GPP_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_FLUXES(s:f,1)-model) - &
                                                    (DATAin%GPP(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%ngpp))
    endif

!    ! Evap Log-likelihood
!    if (DATAin%nEvap > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),29)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
!                       /DATAin%Evap_unc(DATAin%Evappts(1:DATAin%nEvap)))**2)
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nEvap))
!    endif

!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),24)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nFire))
!    endif
    ! Fire Log-likelihood
    ! Fire partitioned between the mean flux and seasonal anomalies
    if (DATAin%nFire > 0) then
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%Fire(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean fire from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%Fire(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,24)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%Fire_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_FLUXES(s:f,24)-model) - &
                                                    (DATAin%Fire(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nFire))
    endif

    ! LAI log-likelihood
    ! Assume physical property is best represented as the mean of value at beginning and end of times step
    if (DATAin%nlai > 0) then
       ! Create vector of (LAI_t0 + LAI_t1) * 0.5, note / pars(17) to convert foliage C to LAI
       mid_state = ( ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0 ) / pars(15)
       ! Split loop to allow vectorisation
       tot_exp = sum(((mid_state(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       ! loop split to allow vectorisation
       !tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
       !                /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (mid_state(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp = tot_exp+(-log(infini))
         endif
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nlai))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnee))
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
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nreco))
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       ! Create vector of (FOL_t0 + FOL_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)) &
                       -DATAin%Cfol_stock(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))&
                     / DATAin%Cfol_stock_unc(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))**2)
       ! Sum with current likelihood score
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCfol_stock))
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(DATAin%nos_years))
       ! determine the annual max for each pool
       do y = 1, DATAin%nos_years
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
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCfolmax_stock))
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       ! Create vector of (root_t0 + root_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,3) + DATAin%M_POOLS(2:(DATAin%nodays+1),3) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)) &
                       -DATAin%Croots_stock(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))&
                     / DATAin%Croots_stock_unc(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))**2)
       ! Combine with existing likelihood estimate
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCroots_stock))
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       ! Create vector of (lit_t0 + lit_t1) * 0.5
       !mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
       !          * 0.5d0
       mid_state = sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,28)+DATAin%M_FLUXES(:,29)+ & 
                       DATAin%M_FLUXES(:,34)+DATAin%M_FLUXES(:,35))
       mid_state = (mid_state / &
                   (mid_state + sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,30)+DATAin%M_FLUXES(:,36)))) & 
                 * DATAin%M_POOLS(:,4)
       !mid_state = (sum(DATAin%M_FLUXES(:,9))/sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,10))) &
       !          * DATAin%M_POOLS(:,4)
       mid_state = (mid_state(1:DATAin%nodays) + mid_state(2:(DATAin%nodays+1))) * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Clit_stockpts(1:DATAin%nClit_stock)) &
                       -DATAin%Clit_stock(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))&
                     / DATAin%Clit_stock_unc(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))**2)
       ! Combine with existing likelihood estimate
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nClit_stock))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       ! Create vector of (som_t0 + som_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,5) + DATAin%M_POOLS(2:(DATAin%nodays+1),5) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)) &
                       -DATAin%Csom_stock(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))&
                     / DATAin%Csom_stock_unc(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))**2)
       ! Combine with existing likelihood estimate
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCsom_stock))
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

!    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation
!    if (DATAin%otherpriors(4) > -9998) then
!        tot_exp = sum(DATAin%M_FLUXES(:,29)) / sum(DATAin%MET(7,:) * 86400d0)
!        scale_likelihood = scale_likelihood-((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
!    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    scale_likelihood = scale_likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (scale_likelihood /= scale_likelihood) then
       scale_likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function scale_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function sqrt_scale_likelihood(npars,pars)
    use cardamom_structures, only: DATAin

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, y, s, f
    double precision :: tot_exp, tmp_var, infini, input, output, model, obs, unc
    double precision, dimension(DATAin%nodays) :: mid_state
    double precision, dimension(DATAin%steps_per_year) :: sub_time
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    sqrt_scale_likelihood = 0d0 ; infini = 0d0 ; mid_state = 0d0 ; sub_time = 0d0

!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nnbe)))
!    endif
    ! NBE Log-likelihood
    ! NBE partitioned between the mean flux and seasonal anomalies
    if (DATAin%nnbe > 0) then
!        ! Determine the mean value for model, observtion and uncertainty estimates
!        obs   = sum(DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!              / dble(DATAin%nnbe)
!        model = sum(DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+ &
!                    DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!              / dble(DATAin%nnbe)
!        unc   = sqrt(sum(DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe))**2)) &
!              / dble(DATAin%nnbe)
!        ! Update the likelihood score with the mean bias
!        sqrt_scale_likelihood = sqrt_scale_likelihood - (((model - obs) / unc) ** 2)
!        ! Determine the anomalies based on substraction of the global mean
!        ! Greater information would come from breaking this down into annual estimates
!        tot_exp = sum( (((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)-model) - &
!                        (DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe)) - obs)) / unc)**2 )
!        likelihood = likelihood-tot_exp
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%NBE(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean NBE from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%NBE(s:f)*sub_time) / sum(sub_time)
               model = sum((DATAin%M_NEE(s:f)+ &
                            DATAin%M_FLUXES(s:f,24))*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%NBE_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_NEE(s:f)+DATAin%M_FLUXES(s:f,24)-model) - &
                                                    (DATAin%NBE(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nnbe)))
    endif ! nnbe > 0

!    ! GPP Log-likelihood
!    if (DATAin%ngpp > 0) then
!       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
!                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%ngpp)))
!    endif
    ! GPP Log-likelihood
    ! GPP partitioned between the mean flux and seasonal anomalies
    if (DATAin%ngpp > 0) then
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%GPP(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean gpp from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%GPP(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,1)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%GPP_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_FLUXES(s:f,1)-model) - &
                                                    (DATAin%GPP(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%ngpp)))
    endif

!    ! Evap Log-likelihood
!    if (DATAin%nEvap > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),29)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
!                       /DATAin%Evap_unc(DATAin%Evappts(1:DATAin%nEvap)))**2)
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nEvap)))
!    endif

!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nFire)))
!    endif
    ! Fire Log-likelihood
    ! Fire partitioned between the mean flux and seasonal anomalies
    if (DATAin%nFire > 0) then
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%Fire(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean fire from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%Fire(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,24)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%Fire_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_FLUXES(s:f,24)-model) - &
                                                    (DATAin%Fire(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nFire)))
    endif

    ! LAI log-likelihood
    ! Assume physical property is best represented as the mean of value at beginning and end of times step
    if (DATAin%nlai > 0) then
       ! Create vector of (LAI_t0 + LAI_t1) * 0.5, note / pars(17) to convert foliage C to LAI
       mid_state = ( ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0 ) / pars(15)
       ! Split loop to allow vectorisation
       tot_exp = sum(((mid_state(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       ! loop split to allow vectorisation
       !tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
       !                /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (mid_state(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp = tot_exp+(-log(infini))
         endif
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nlai)))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nnee)))
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
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nreco)))
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       ! Create vector of (FOL_t0 + FOL_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)) &
                       -DATAin%Cfol_stock(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))&
                     / DATAin%Cfol_stock_unc(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))**2)
       ! Sum with current likelihood score
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCfol_stock)))
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(DATAin%nos_years))
       ! determine the annual max for each pool
       do y = 1, DATAin%nos_years
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
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCfolmax_stock)))
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       ! Create vector of (root_t0 + root_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,3) + DATAin%M_POOLS(2:(DATAin%nodays+1),3) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)) &
                       -DATAin%Croots_stock(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))&
                     / DATAin%Croots_stock_unc(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))**2)
       ! Combine with existing likelihood estimate
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCroots_stock)))
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       ! Create vector of (lit_t0 + lit_t1) * 0.5
       !mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
       !          * 0.5d0
       mid_state = sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,28)+DATAin%M_FLUXES(:,29)+ & 
                       DATAin%M_FLUXES(:,34)+DATAin%M_FLUXES(:,35))
       mid_state = (mid_state / &
                   (mid_state + sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,30)+DATAin%M_FLUXES(:,36)))) & 
                 * DATAin%M_POOLS(:,4)
       !mid_state = (sum(DATAin%M_FLUXES(:,9))/sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,10))) &
       !          * DATAin%M_POOLS(:,4)
       mid_state = (mid_state(1:DATAin%nodays) + mid_state(2:(DATAin%nodays+1))) * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Clit_stockpts(1:DATAin%nClit_stock)) &
                       -DATAin%Clit_stock(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))&
                     / DATAin%Clit_stock_unc(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))**2)
       ! Combine with existing likelihood estimate
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nClit_stock)))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       ! Create vector of (som_t0 + som_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,5) + DATAin%M_POOLS(2:(DATAin%nodays+1),5) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)) &
                       -DATAin%Csom_stock(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))&
                     / DATAin%Csom_stock_unc(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))**2)
       ! Combine with existing likelihood estimate
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCsom_stock)))
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

!    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation
!    if (DATAin%otherpriors(4) > -9998) then
!        tot_exp = sum(DATAin%M_FLUXES(:,29)) / sum(DATAin%MET(7,:) * 86400d0)
!        sqrt_scale_likelihood = sqrt_scale_likelihood-((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
!    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    sqrt_scale_likelihood = sqrt_scale_likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (sqrt_scale_likelihood /= sqrt_scale_likelihood) then
       sqrt_scale_likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function sqrt_scale_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function log_scale_likelihood(npars,pars)
    use cardamom_structures, only: DATAin

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, y, s, f
    double precision :: tot_exp, tmp_var, infini, input, output, model, obs, unc
    double precision, dimension(DATAin%nodays) :: mid_state
    double precision, dimension(DATAin%steps_per_year) :: sub_time
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    log_scale_likelihood = 0d0 ; infini = 0d0 ; mid_state = 0d0 ; sub_time = 0d0

!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nnbe))))
!    endif
    ! NBE Log-likelihood
    ! NBE partitioned between the mean flux and seasonal anomalies
    if (DATAin%nnbe > 0) then
!        ! Determine the mean value for model, observtion and uncertainty estimates
!        obs   = sum(DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!              / dble(DATAin%nnbe)
!        model = sum(DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+ &
!                    DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)) &
!              / dble(DATAin%nnbe)
!        unc   = sqrt(sum(DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe))**2)) &
!              / dble(DATAin%nnbe)
!        ! Update the likelihood score with the mean bias
!        log_scale_likelihood = log_scale_likelihood - (((model - obs) / unc) ** 2)
!        ! Determine the anomalies based on substraction of the global mean
!        ! Greater information would come from breaking this down into annual estimates
!        tot_exp = sum( (((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),24)-model) - &
!                        (DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe)) - obs)) / unc)**2 )
!        likelihood = likelihood-tot_exp
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%NBE(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean NBE from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%NBE(s:f)*sub_time) / sum(sub_time)
               model = sum((DATAin%M_NEE(s:f)+ &
                            DATAin%M_FLUXES(s:f,24))*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%NBE_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_NEE(s:f)+DATAin%M_FLUXES(s:f,24)-model) - &
                                                    (DATAin%NBE(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nnbe))))
    endif ! nnbe > 0

!    ! GPP Log-likelihood
!    if (DATAin%ngpp > 0) then
!       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
!                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%ngpp))))
!    endif
    ! GPP Log-likelihood
    ! GPP partitioned between the mean flux and seasonal anomalies
    if (DATAin%ngpp > 0) then
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%GPP(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean NBE from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%GPP(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,1)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%GPP_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_FLUXES(s:f,1)-model) - &
                                                    (DATAin%GPP(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%ngpp))))
    endif

!    ! Evap Log-likelihood
!    if (DATAin%nEvap > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),29)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
!                       /DATAin%Evap_unc(DATAin%Evappts(1:DATAin%nEvap)))**2)
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nEvap))))
!    endif

!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nFire))))
!    endif
    ! Fire Log-likelihood
    ! Fire partitioned between the mean flux and seasonal anomalies
    if (DATAin%nFire > 0) then
        ! Reset variable
        tot_exp = 0d0
        ! Loop through each year
        do y = 1, DATAin%nos_years
           ! Reset selection variable
           sub_time = 0d0
           ! Determine the start and finish of the current year of interest
           s = ((DATAin%steps_per_year*(y-1))+1) ; f = (DATAin%steps_per_year*y)
           where (DATAin%Fire(s:f) > -9998d0) sub_time = 1d0
           if (sum(sub_time) > 0d0) then
               ! Determine the current years mean Fire from observations and model
               ! assuming we select only time steps in the model time series which have
               ! an estimate in the observations
               obs   = sum(DATAin%Fire(s:f)*sub_time) / sum(sub_time)
               model = sum(DATAin%M_FLUXES(s:f,24)*sub_time) / sum(sub_time)
               unc   = sqrt(sum((sub_time*DATAin%Fire_unc(s:f))**2)) / sum(sub_time)
               ! Update the likelihood score with the mean bias
               tot_exp = tot_exp + (((model - obs) / unc) ** 2)
               ! Determine the anomalies based on substraction of the annual mean
               ! Greater information would come from breaking this down into annual estimates
               tot_exp = tot_exp + sum( (sub_time*(((DATAin%M_FLUXES(s:f,24)-model) - &
                                                    (DATAin%Fire(s:f) - obs)) / unc))**2 )
           end if ! sum(sub_time) > 0
        end do ! loop years
        ! Update the likelihood score with the anomaly estimates
        log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nFire))))
    endif

    ! LAI log-likelihood
    ! Assume physical property is best represented as the mean of value at beginning and end of times step
    if (DATAin%nlai > 0) then
       ! Create vector of (LAI_t0 + LAI_t1) * 0.5, note / pars(17) to convert foliage C to LAI
       mid_state = ( ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0 ) / pars(15)
       ! Split loop to allow vectorisation
       tot_exp = sum(((mid_state(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       ! loop split to allow vectorisation
       !tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
       !                /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (mid_state(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp = tot_exp+(-log(infini))
         endif
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nlai))))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nnee))))
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
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nreco))))
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       ! Create vector of (FOL_t0 + FOL_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)) &
                       -DATAin%Cfol_stock(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))&
                     / DATAin%Cfol_stock_unc(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))**2)
       ! Sum with current likelihood score
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nCfol_stock))))
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(DATAin%nos_years))
       ! determine the annual max for each pool
       do y = 1, DATAin%nos_years
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
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nCfolmax_stock))))
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       ! Create vector of (root_t0 + root_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,3) + DATAin%M_POOLS(2:(DATAin%nodays+1),3) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)) &
                       -DATAin%Croots_stock(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))&
                     / DATAin%Croots_stock_unc(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))**2)
       ! Combine with existing likelihood estimate
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nCroots_stock))))
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       ! Create vector of (lit_t0 + lit_t1) * 0.5
       !mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
       !          * 0.5d0
       mid_state = sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,28)+DATAin%M_FLUXES(:,29)+ & 
                       DATAin%M_FLUXES(:,34)+DATAin%M_FLUXES(:,35))
       mid_state = (mid_state / &
                   (mid_state + sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,30)+DATAin%M_FLUXES(:,36)))) & 
                 * DATAin%M_POOLS(:,4)
       !mid_state = (sum(DATAin%M_FLUXES(:,9))/sum(DATAin%M_FLUXES(:,9)+DATAin%M_FLUXES(:,10))) &
       !          * DATAin%M_POOLS(:,4)
       mid_state = (mid_state(1:DATAin%nodays) + mid_state(2:(DATAin%nodays+1))) * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Clit_stockpts(1:DATAin%nClit_stock)) &
                       -DATAin%Clit_stock(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))&
                     / DATAin%Clit_stock_unc(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))**2)
       ! Combine with existing likelihood estimate
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nClit_stock))))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       ! Create vector of (som_t0 + som_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,5) + DATAin%M_POOLS(2:(DATAin%nodays+1),5) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)) &
                       -DATAin%Csom_stock(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))&
                     / DATAin%Csom_stock_unc(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))**2)
       ! Combine with existing likelihood estimate
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0+log(dble(DATAin%nCsom_stock))))
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

    ! Initial soil water prior. The actual prior is linked to a fraction of
    ! field capacity so here is were that soil water at t=1
    ! is actually assessed against an observation
!    if (DATAin%otherpriors(1) > -9998) then
!        tot_exp = (DATAin%M_POOLS(1,7) * 1d-3) / layer_thickness(1) ! convert mm -> m3/m3
!        log_scale_likelihood = log_scale_likelihood-((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2
!    end if

!    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation
!    if (DATAin%otherpriors(4) > -9998) then
!        tot_exp = sum(DATAin%M_FLUXES(:,29)) / sum(DATAin%MET(7,:) * 86400d0)
!        log_scale_likelihood = log_scale_likelihood-((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
!    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    log_scale_likelihood = log_scale_likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (log_scale_likelihood /= log_scale_likelihood) then
       log_scale_likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function log_scale_likelihood
  !
  !------------------------------------------------------------------
  !
end module model_likelihood_module
