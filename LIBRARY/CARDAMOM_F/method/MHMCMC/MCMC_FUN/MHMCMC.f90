module MHMCMC_MODULE

! module contains all subroutine and functions relevant specifically to the
! MHMCMC method. The choice of EDC, likelihood and model are made else where and
! are thus contains within a seperate module

implicit none

! make all private
private

! specify what can be seen
public :: MHMCMC, par_minstepsize,par_initstepsize

! declare any module level variables needed

! related to random number generator
integer :: uniform, unif_length
double precision, allocatable, dimension(:) :: uniform_random_vector
! MHMCMC step size
double precision, parameter :: par_minstepsize = 0.001d0 & ! 0.0005
                              ,par_maxstepsize = 0.01d0  &
                              ,par_initstepsize = 0.005d0
! Delayed Rejection flag
double precision :: DR_scaler = 1d0 &
                   ,opt_scaling ! scd = 2.381204 the optimal scaling parameter
                                ! for MCMC search, when applied to  multivariate proposal.
                                ! NOTE 1: 2.38 / sqrt(npars) sometimes used when applied to the Cholesky
                                ! factor. NOTE 2: 2.381204 ** 2 = 5.670132
double precision, parameter :: beta = 0.05d0 !, delta = 0.5d0
! Is current proposal multivariate or not?
logical :: multivariate_proposal = .false.!, do_DR = .true.
integer, parameter :: N_before_mv = 10d0 !3

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine MHMCMC (P_target,model_likelihood_option)
    use MCMCOPT, only: PI, MCO, MCOUT, COUNTERS
    use math_functions, only: randn, random_uniform, log_par2nor, log_nor2par, par2nor, nor2par &
                             ,dsymv
    use cardamom_io, only: write_parameters,write_variances,write_covariance_matrix &
                          ,write_covariance_info,restart_flag
    use cardamom_structures, only: DATAin
    use model_likelihood_module, only: model_likelihood

    implicit none

    !/* ***********INPUTS************
    ! *
    ! * MODEL_LIKELYHOOD: A function wholly responsible for
    ! * (a) running the model given the DATA and parameters,
    ! * (b) comparing it to observations,and
    ! * (c) returning  the (log) likelihood.
    ! * The function will be run as MODEL_LIKELIHOOD(DATA,PI,PARS);
    ! * To facilitate this, ALL data can be
    ! * passed to the MHMCMC function as a structure (in order to avoid
    ! * repeated read/write computational time).
    ! *
    ! * DATA: All data needed for the MODEL_LIKELYHOOD. It can include
    ! * drivers, observations, etc.
    ! *
    ! * PARINFO: This structure contains information on
    ! * (a) pmin, pmax:      parameter ranges (compulsory)
    ! * (b) initpars:        parameter starting values (optional/recommended).
    ! * (c) npars:           number of pars
    ! *
    ! * MCO: This structure contains option values for the MCMC run.
    ! * These will be set to default values if empty. Options include:
    ! * (a) number of runs
    ! * (b) filename for writing file with results
    ! * (c) step adaptation frequency
    ! * (d) initial step size
    ! * */
    !
    !/* **************OUTPUTS*************
    ! *
    ! * RESULTS FILE: File includes (a) results (b) likelihood and (c) final step
    ! size
    ! *
    ! * */

    ! declare interface for the model likelihood function.
    ! NOTE that the inputted MODEL_LIKEIHOOD_OPTION could be multiple subroutines,
    ! the interface allows for making the requirements of this explicit
    interface
      subroutine model_likelihood_option(param_vector, ML_obs_out, ML_prior_out)
        use cardamom_structures, only: DATAin, emulator_pars
        use MCMCOPT, only: PI
        use CARBON_MODEL_MOD, only: carbon_model
           implicit none
           ! declare input variables
           double precision, dimension(PI%npars), intent(inout) :: param_vector
           ! output
           double precision,intent(inout) :: ML_obs_out, ML_prior_out
      end subroutine model_likelihood_option
    end interface

    ! Remaining Arguments
    double precision, intent(in) :: P_target

    ! declare any local variables
    type ( counters ) :: N
    double precision, dimension(PI%npars) :: deltaPARS2_iC &
                                            ,deltaPARS0_iC &
                                            ,norPARS0      & ! normalised parameter values for current state
                                            ,norPARS       & ! normalised parameter values for current proposal
                                            ,norPARS2      & ! normalised parameters for the DR proposal
                                            ,PARS0         & ! parameter values for current state
                                            ,PARS          & ! parameter values for current proposal
                                            ,PARS2         & ! parameters for the DR proposal
                                            ,BESTPARS        ! best set of parameters so far

    double precision, dimension(PI%npars,MCO%nADAPT) :: PARSALL ! All accepted normalised parameters since previous step adaption
    double precision :: infini &
                       ,target_P &
                       ,burn_in_period &
                       ,crit1  & ! random numbers log(0->1) used to accept / reject
                       ,AM_likelihood &
                       ,DR_likelihood &
                       ,AM_vs_DR_P &
                       ,DR_vs_current &
                       ,DR_vs_current_pars &
                       ,outputP0      &
                       ,outputP0prior  &
                       ,DR_P, DR_Pprior & ! likelihood scores from Delayed Rejection step
                       ,Pmax, P0prior, Pprior & ! as below but for priors only
                       ,P0 & ! previously accepted observation based log-likelihood
                       ,P    ! current observation based log-likelihood
    integer :: i

    ! initial values
    uniform = 1 
    P = -1d0 ; Pprior = -1d0
    N%ACC = 0d0 ; N%ACC_first = 0d0 ; N%ITER = 0d0 ; N%ACC_beta = 0d0
    N%ACCLOC = 0d0 ; N%ACCRATE = 0d0 ; N%ACCRATE_GLOBAL = 0d0
    N%ACCLOC_beta = 0d0 ; N%ITER_beta = 0d0 ; N%ACCRATE_beta = 0d0

    ! Determine how long we will continue to adapt our proposal covariance
    ! matrix and use of Delayed Rejection
    burn_in_period = MCO%fADAPT*dble(MCO%nOUT)

    ! See step() for relevant references.
    ! scd = 2.381204 the optimal scaling parameter for MCMC search, when applied
    ! to multivariate proposal.
    ! NOTE 1: 2.38 / sqrt(npars) sometimes used when applied to the Cholesky factor
    ! NOTE 2: 2.381204 ** 2 = 5.670132
    opt_scaling = 5.670132d0 / dble(PI%npars)

    ! calculate initial vector of uniform random values
    unif_length = MCO%nADAPT * 5
    allocate(uniform_random_vector(unif_length))
    call random_uniform(uniform_random_vector,unif_length)

    ! assume this is a restart which we want to load previous values
!    if (restart_flag) then
!!        N%ITER = dble(MCO%nWRITE * accepted_so_far) ! number of iterations is directly linked to the number output
!!        N%ACC = ceiling(accept_rate * N%ITER)
!        ! Load the total number of iterations rather than starting from one
!        N%ITER = MCOUT%nos_iterations
!        N%ACC = MCOUT%acceptance_rate
!        N%ACC_first = N%ACC ! conservative approch to assume all proposals are first accepted
!!        N%ACC_beta = ceiling(N%ACC * beta)
!    endif

    ! add something here to delete previous files if wanted later
    if (MCO%APPEND == 0 .and. MCO%nWRITE > 0) then
        write(*,*) "Oooops have requested that existing files be deleted but you have not finished the code to do so...."
    end if

    ! start random sampling if MCO%randparini set
    PI%parfix = 0d0
    do i = 1, PI%npars
       ! parfix = 1 stay at prior value, parfix = 0 randomly search
       if (MCO%fixedpars .and. PI%parini(i) /= -9999d0) PI%parfix(i) = 1d0 
       ! only assign random parameters if (a) randparini == .true. or (b) PI$parini(n) == -9999)
       if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
           call nor2par(1,uniform_random_vector(uniform),PI%parmin(i),PI%parmax(i),PI%parini(i))
!           call log_nor2par(1,uniform_random_vector(uniform),PI%parmin(i),PI%parmax(i),PI%parini(i))
           uniform = uniform + 1
       end if
       ! write(*,*) parameter values to screen
       write(*,*) "p",i,"=",PI%parini(i)
    end do ! for PI%npar loop

    ! Inform the user
    write(*,*) "Have loaded / randomly assigned PI%parini - now begin the MHMCMC"

    ! Initialise the prior and best pars vectors
    PARS0(1:PI%npars) = PI%parini(1:PI%npars)
    BESTPARS(1:PI%npars) = PI%parini(1:PI%npars)

    ! Track normalised initial proposal
    do i = 1, PI%npars
       call par2nor(1,PARS0(i),PI%parmin(i),PI%parmax(i),norPARS0(i))
!       call log_par2nor(1,PARS0(i),PI%parmin(i),PI%parmax(i),norPARS0(i))
    end do

    ! calculate the initial probability / log likelihood.
    ! NOTE: passing P0 -> P is needed during the EDC searching phase where we
    ! could read an EDC consistent parameter set in the first instance
    call model_likelihood_option(PI%parini,P0,P0prior) ; P = P0 ; Pprior = P0prior
    write(*,*) "Starting likelihood = ",P0,"+",P0prior
    Pmax = P0 + P0prior

    ! checks whether the EDCs (combined with P0 not P0prior) have been met in the initial parameter set
    infini = 0d0
    if (P0 == log(infini)) then
        write(*,*) "WARNING! P0 = ",P0," - MHMCMC may get stuck, if so please check initial conditins"
    endif

    ! Begin the main MHMCMC loop
!    do while (N%ITER < MCO%nOUT .and. (Pmax < P_target .or. MCO%nWRITE > 0))
    do while (N%ITER < MCO%nOUT .and. Pmax < P_target)

       ! take a step in parameter space
       call step(N,PARS0,PARS,norPARS0,norPARS,.false.)

       ! if parameter proposal in bounds check the model
       if (minval(norPARS) > 0d0 .and. maxval(norPARS) < 1d0) then

           ! calculate the model likelihood
           call model_likelihood_option(PARS, P, Pprior)

           ! accept or reject, draw uniform distribution (0,1)
           crit1 = log(uniform_random_vector(uniform))
           uniform = uniform + 1
           ! if we are near to the end re-generate some more values
           if (uniform >= unif_length) then
               ! calculate new vector of uniform random values
               call random_uniform(uniform_random_vector,unif_length)
               ! and reset uniform counter
               uniform = 1
           endif

           ! determine accept or reject the current proposal
           AM_likelihood = ( (P-P0) + (Pprior-P0prior) )

       else

           ! proposal out of parameter bounds, set likelihoods to ensure
           ! rejection
           AM_likelihood = log(infini) ; P = AM_likelihood ; Pprior = P
           crit1 = 0d0

       end if ! in bound

       if ( AM_likelihood > crit1) then

           ! Store accepted parameter proposals
           ! keep record of all parameters accepted since step adaption
           PARSALL(1:PI%npars,(nint(N%ACCLOC)+1)) = norPARS(1:PI%npars)
           PARS0(1:PI%npars) = PARS(1:PI%npars)
           norPARS0(1:PI%npars) = norPARS(1:PI%npars)
           ! specifically store the best parameter set
           if ((P+Pprior) >= Pmax) then
               BESTPARS = PARS ; Pmax = P+Pprior
           endif
           N%ACCLOC = N%ACCLOC + 1d0
           if (.not.multivariate_proposal) then
               ! accepted from beta step proposal
!               N%ACCLOC_beta = N%ACCLOC_beta + 1d0
           else
               ! accepted first proposal from multivarite
               N%ACC_first = N%ACC_first + 1d0
           end if
           P0 = P ; P0prior = Pprior

!       else ! we have rejected - apply Delayed Rejection approach to enhance local searching
!
!           if (do_DR .and. multivariate_proposal) then
!
!                DR_scaler = 0.01d0
!
!                ! take a step in parameter space
!                call step(N,PARS0,PARS2,norPARS0,norPARS2,.true.)
!
!                if (minval(norPARS2) > 0d0 .and. maxval(norPARS2) < 1d0) then
!
!                    ! calculate the model likelihood
!                    call model_likelihood_option(PARS2, DR_P, DR_Pprior)
!
!                    ! Determine whether we accept or reject the DR proposal
!
!                    ! First what is the likelihood ratio between the 1st propsal (AM) and the
!                    ! second propsal (DR), note we are setting maximum weighting towards AM proposal
!                    ! if DR is worse than AM
!                    if (P+Pprior == log(infini)) then
!                        ! set to zero probability
!                        AM_vs_DR_P = 0d0
!                    else
!                        ! extimate in log-likelihood and convert to probability, may be zero frequently
!                        AM_vs_DR_P = min(1d0,exp((P-DR_P) + (Pprior-DR_Pprior)))
!                    endif
!                    ! Second estimate the liklihood of the DR propsal compared to the P0 and P0prior
!                    DR_vs_current = (DR_P-P0) + (DR_Pprior-P0prior)
!                    ! Third estimate likelihood based on both the current (just rejected) proposal
!                    ! and current position in chain assuming mulitvariate gaussian errors between
!                    call dsymv("U",PI%npars,1d0,PI%iC,PI%npars,norPARS2-norPARS,1,0d0,deltaPARS2_iC,1)
!                    call dsymv("U",PI%npars,1d0,PI%iC,PI%npars,norPARS0-norPARS,1,0d0,deltaPARS0_iC,1)
!                    DR_vs_current_pars = -0.5d0 * &
!                                        (sum(deltaPARS2_iC * (norPARS2 - norPARS)) - &
!                                         sum(deltaPARS0_iC * (norPARS0 - norPARS)) )
!                    ! Finally calculate the final likelihood of the DR proposal vs AM and the existing parameter set
!                    DR_likelihood = DR_vs_current + DR_vs_current_pars + log((1d0-AM_vs_DR_P)/(1d0-exp(AM_likelihood)))
!
!                    ! accept or reject, draw uniform distribution (0,1)
!                    crit1 = log(uniform_random_vector(uniform))
!                    uniform = uniform + 1
!                    ! if we are near to the end re-generate some more values
!                    if (uniform >= unif_length) then
!                       ! calculate new vector of uniform random values
!                       call random_uniform(uniform_random_vector,unif_length)
!                       ! and reset uniform counter
!                       uniform = 1
!                    endif
!
!                else
!
!                    ! parameters are out of bounds
!                    DR_likelihood = log(infini) ; DR_P = DR_likelihood ; DR_Pprior = DR_P
!                    crit1 = 0d0
!
!                endif ! inbounds?
!
!                ! now do we accept or reject the DR step as a function of both PARS0 and PARS
!                if ( DR_likelihood > crit1 ) then
!
!                   ! accepted DR proposal
!
!                   ! Store accepted parameter solutions
!                   ! keep record of all parameters accepted since step adaption
!                   PARSALL(1:PI%npars,(nint(N%ACCLOC)+1)) = norPARS2(1:PI%npars)
!                   PARS0(1:PI%npars) = PARS2(1:PI%npars)
!                   norPARS0(1:PI%npars) = norPARS2(1:PI%npars)
!                   ! specifically store the best parameter set
!                   if ((DR_P+DR_Pprior) >= Pmax) then
!                       BESTPARS = PARS2 ; Pmax = DR_P+DR_Pprior
!                   endif
!                   N%ACCLOC = N%ACCLOC + 1d0
!                   P0 = DR_P ; P0prior = DR_Pprior
!
!                end if
!
!                ! return DR scaler back to 1d0
!                DR_scaler = 1d0
!
!            endif ! multivariate_proposal

       endif ! accept or reject condition

       ! count iteration whether the current proposal is accepted or rejected
       N%ITER = N%ITER + 1

       if (MCO%nWRITE > 0 .and. mod(nint(N%ITER),MCO%nWRITE) == 0) then
           ! calculate the likelhooh for the actual uncertainties - this avoid
           ! issues with different phases of the MCMC which may use sub-samples
           ! of observations or inflated uncertainties to aid parameter
           ! searching
           call model_likelihood(PARS0, outputP0, outputP0prior)
           call write_variances(PI%parvar,PI%npars,N%ACCRATE) ! should this be the global rate?
           call write_parameters(PARS0,(outputP0+outputP0prior),PI%npars)
           call write_covariance_matrix(PI%covariance,PI%npars,.false.)
           call write_covariance_info(PI%mean_par,PI%Nparvar,PI%npars)
       end if ! write or not to write

       ! time to adapt?
       if (mod(nint(N%ITER),MCO%nADAPT) == 0) then

           ! Total accepted values
           N%ACC = N%ACC + N%ACCLOC
!           N%ACC_beta = N%ACC_beta + N%ACCLOC_beta
           ! Calculate global acceptance rate
           N%ACCRATE_GLOBAL = N%ACC / N%ITER

           ! Estimate acceptance rate for the multivariate set at first proposal only
!           if (burn_in_period > N%ITER .or. (N%ACC_first / N%ITER) < 0.05d0) then
!               do_DR = .true.
!           else
!               do_DR = .false.
!           endif

           ! First, we can only do a covariance adaption if we have accepted some
           ! parameters in the last period
!           if (N%ACCLOC > 0d0) then

               ! Calculate local acceptance rate (i.e. since last adapt)
               N%ACCRATE = N%ACCLOC / dble(MCO%nADAPT)
!               N%ACCRATE_beta = N%ACCLOC_beta / N%ITER_beta

               ! Second, are we still in the adaption phase?
               if (burn_in_period > N%ITER .or. (N%ACC_first / N%ITER) < 0.05d0 .or. .not.PI%use_multivariate) then

                   ! Once covariance matrix has been created just update based on a
                   ! single parameter set from each period.
                   if (PI%cov) then
                       N%ACCLOC = 1d0 ; PARSALL(1:PI%npars,nint(N%ACCLOC)) = norPARS0(1:PI%npars)
                   else if (N%ACCLOC > 3d0) then
                       PARSALL(1:PI%npars,2) = PARSALL(1:PI%npars,ceiling(N%ACCLOC*0.5d0))
                       PARSALL(1:PI%npars,3) = PARSALL(1:PI%npars,N%ACCLOC)
                       N%ACCLOC = 3d0
                   endif

                   ! adapt the covariance matrix for multivariate proposal
                   call adapt_step_size(PARSALL,N)

               end if !  have enough parameter been accepted

!           else ! any newly accepted parameters?

               ! no new parameters accepted...

               !...so acceptance rates are zero but no need to update N%ACC or
               ! N%ACC_beta
               N%ACCRATE = 0d0 !; N%ACCRATE_beta = 0d0
              
!           end if ! have any parameters been accepted in the last period?

           ! resets to local counters
           N%ACCLOC = 0d0 !; N%ACCLOC_beta = 0d0

       end if ! time to adapt?

       ! Should I be write(*,*)ing to screen or not?
       if (MCO%nPRINT > 0 .and. (mod(nint(N%ITER),MCO%nPRINT) == 0)) then
           write(*,*)"Using multivariate sampling = ",PI%use_multivariate
           write(*,*)"Total proposal = ",N%ITER," out of ",MCO%nOUT
           write(*,*)"Total accepted = ",N%ACC
           write(*,*)"Overall acceptance rate    = ",N%ACC / N%ITER
!           write(*,*)"Overall beta accept rate   = ",N%ACC_beta / N%ITER_beta
           write(*,*)"Local   acceptance rate    = ",N%ACCRATE
!           write(*,*)"Local beta acceptance rate = ",N%ACCRATE_beta
           write(*,*)"Current obs   = ",P0,"proposed = ",P," log-likelihood"
           write(*,*)"Current prior = ",P0prior,"proposed = ",Pprior," log-likelihood"
           write(*,*)"Maximum likelihood = ",Pmax
           ! NOTE: that -infinity in current obs only indicates failure of EDCs
           ! but -infinity in both obs and parameter likelihood scores indicates
           ! that proposed parameters are out of bounds
       end if ! write(*,*) to screen or not

    end do ! while conditions

    ! write out final covariance matrix for the analysis
    if (MCO%nWRITE > 0) call write_covariance_matrix(PI%covariance,PI%npars,.false.)

    ! record the best single set of parameters
    MCOUT%best_pars(1:PI%npars) = BESTPARS(1:PI%npars)
    ! set the initial parameter set the final one accepted
    PI%parini(1:PI%npars) = PARS0(1:PI%npars)
    ! record how many iterations were taken to complete
    MCOUT%nos_iterations = MCOUT%nos_iterations + N%ITER
    ! set flag MCMC completed
    MCOUT%complete = 1
    ! tidy up
    deallocate(uniform_random_vector)

    ! completed MHMCMC loop
    write(*,*)"MHMCMC loop completed"
    write(*,*)"Overall final acceptance rate = ",N%ACC / N%ITER
    write(*,*)"Best log-likelihood = ",Pmax

  end subroutine MHMCMC
  !
  !------------------------------------------------------------------
  !
  subroutine adapt_step_size(PARSALL,N)
    use cardamom_io, only: write_covariance_matrix, write_covariance_info
    use MCMCOPT, only: MCO, PI, COUNTERS
    use math_functions, only: nor2par, par2nor, log_nor2par, log_par2nor, &
                              cholesky_factor, std, covariance_matrix, &
                              increment_covariance_matrix, inverse_matrix

    ! Update the multivariate propsal distribution.
    ! Ensure that this subroutine is only called if at least 1 parameter propsal
    ! has been accepted in the last adaption period.

    implicit none

    ! declare input types
    type ( counters ), intent(inout) :: N


    ! declare inputs variables
    double precision, intent(in) :: PARSALL(PI%npars,MCO%nADAPT) ! collection of recently accepted normalised parameter combinations

    ! declare local variables
    integer p, i, info ! counters
    double precision, dimension(PI%npars,PI%npars) :: cholesky, cov_backup
    double precision, dimension(PI%npars) :: mean_par_backup
    double precision :: Nparvar_backup    

    ! if we have a covariance matrix then we want to update it, if not then we need to create one
    if (PI%cov) then

        ! Increment the variance-covariance matrix with new accepted parameter sets
        ! NOTE: that this also increments the total accepted counter (PI%Nparvar)
        cov_backup = PI%covariance ; mean_par_backup = PI%mean_par ; Nparvar_backup = PI%Nparvar
        call increment_covariance_matrix(PARSALL(1:PI%npars,1:nint(N%ACCLOC)),PI%mean_par,PI%npars &
                                        ,PI%Nparvar,nint(N%ACCLOC),PI%covariance)

        ! Calculate the cholesky factor as this includes a determination of
        ! whether the covariance matrix is positive definite.
        cholesky = PI%covariance
        call cholesky_factor( PI%npars, cholesky, info )
        ! If the updated covariance matrix is not positive definite we should
        ! reject the update in favour of the existing matrix
        if (info == 0) then
            ! Set multivariate sampling to true
            PI%use_multivariate = .true.
            ! Re-calculate inverse matrix as it may now be used,
            ! direct inversion of covariance matrix based on Doolittle LU
            ! factorization
 !           if (do_DR) call inverse_matrix( PI%npars, PI%covariance, PI%iC )
        else 
            ! The current addition of a parameter leads to a matrix which is not
            ! positive definite. If we previously had a matrix which is positive
            ! definite then we should reject totally the new matrix, if not we
            ! should keep it and accumulate the information
            if (PI%use_multivariate) then
                ! return original matrix to place
                PI%covariance = cov_backup
                PI%mean_par = mean_par_backup
                PI%Nparvar = Nparvar_backup
            else
                ! Keep accumulating the information
                PI%use_multivariate = .false.
            end if
        endif

    else ! PI%cov == .true.

        ! we have not yet created a covariance matrix based on accepted
        ! parameters. Assuming we have some then create one...
        if (N%ACCLOC > 2d0) then

            ! estimate covariance matrix
            call covariance_matrix(PARSALL(1:PI%npars,1:nint(N%ACCLOC)),PI%mean_par,PI%npars,nint(N%ACCLOC),PI%covariance)
            PI%cov = .true. ; PI%Nparvar = N%ACCLOC

            ! Calculate the cholesky factor as this includes a determination of
            ! whether the covariance matrix is positive definite.
            cholesky = PI%covariance
            call cholesky_factor ( PI%npars, cholesky, info )
            ! If not positive definite then we should not use the multivariat
            ! step at this time.
            if (info /= 0) then
                ! Keep accumulating information until positive definite matrix
                ! calculated
                PI%use_multivariate = .false.
            else
                ! Positive definite found straight away - might as well use it!
                PI%use_multivariate = .true.
                ! Re-calculate inverse matrix as it may now be used,
                ! direct inversion of covariance matrix based on Doolittle LU
                ! factorization
!                if (do_DR) call inverse_matrix( PI%npars, PI%covariance, PI%iC )
            endif

            ! write out first covariance matrix, this will be compared with the final covariance matrix
            if (MCO%nWRITE > 0) then
                call write_covariance_matrix(PI%covariance,PI%npars,.true.)
                call write_covariance_info(PI%mean_par,PI%Nparvar,PI%npars)
            endif

        end if ! N%ACCLOC > 2

    end if ! PI%cov == .true.

    ! adjust step size by local variance
    do p = 1, PI%npars
       ! calculate standards deviation (variability) for the local
       ! window extracted from the variance
       PI%parvar(p) = PI%covariance(p,p)
    end do ! p

    return

  end subroutine adapt_step_size
  !
  !------------------------------------------------------------------
  !
  subroutine step(N,pars0,pars,norpars0,norpars,DR_step)
    use math_functions, only: par2nor, nor2par, log_par2nor, log_nor2par, &
                              random_normal, random_multivariate, random_uniform
    use MCMCOPT, only: PI, COUNTERS

    ! carries out the next step to parameters in the MCMC search

    implicit none

    ! declare input variables
    type ( counters ), intent(inout) :: N
    logical, intent(in) :: DR_step
    double precision, dimension(PI%npars), intent(inout) :: norpars0 & ! normalised current parameters
                                                           ,norpars  & ! normalised proposal
                                                           ,pars0    & ! current parameters
                                                           ,pars       ! proposal

    ! declare local variables
    integer :: p
    double precision :: rn(PI%npars), mu(PI%npars), rn2(PI%npars) &
                       ,tmp, delta_scaler(PI%npars)

    ! reset values
    mu = 0d0

    ! Begin sampling parameter space, first estimate multivariate random number
    ! Multivariate sample around a mean of zero

    ! if we are near to the end re-generate some more values
    if (uniform > (unif_length-3)) then
        ! calculate new vector of uniform random values
        call random_uniform(uniform_random_vector,unif_length)
        ! and reset uniform counter
        uniform = 1
    endif
    ! Draw from uniform distribution to determine whether multivariate step will be used
    tmp = uniform_random_vector(uniform) ; uniform = uniform + 1

    ! Increment step size via different proposal based on whether we have a sufficiently development covariance matrix
    ! Splitting step calculation based on number of parameter vectors accepted
    ! is linked to the need build a covarianc matrix prior to multivariate
    ! sampling.
    ! See Roberts and Rosenthal, Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367, 2009.
!    if (PI%use_multivariate .and. nint(PI%Nparvar) > N_before_mv*PI%npars .and. (tmp > beta .or. DR_step)) then
    if ((PI%use_multivariate .and. nint(PI%Nparvar) > N_before_mv*PI%npars)) then

        ! Is this step a multivariate proposal or not
        multivariate_proposal = .true.

        ! Draw from multivariate random distribution
        ! NOTE: if covariance matrix provided is not positive definite
        !       a sample form normal distribution is returned
!        call random_multivariate(PI%npars, 1, uniform, unif_length, &
!                                 uniform_random_vector, PI%covariance*DR_scaler, mu, rn)
        call random_multivariate(PI%npars, 1, uniform, unif_length, &
                                 uniform_random_vector, PI%covariance, mu, rn)
        ! Sample random normal distribution (mean = 0, sd = 1)
        do p = 1, PI%npars
           call random_normal(uniform,unif_length,uniform_random_vector,rn2(p))
        end do !

        ! Estimate the step to be applied to the current parameter vector to
        ! create the new proposal. scd = a scaling parameter linking searching
        ! stepping to the number of parameters being retrieved by the analysis.
        ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
        ! and references therein.
!        norpars = norpars0 + (rn * opt_scaling)
!        if (tmp > delta) then
            ! Do Beta Gaussian step (par_minstepsize)
            norpars = norpars0 + (rn * opt_scaling * (1d0-beta)) + (par_minstepsize * rn2 * beta)
!        else
!            ! Do Delta Gaussian step (par_maxstepsize)
!            norpars = norpars0 + (rn * opt_scaling * (1d0-beta)) + (par_maxstepsize * rn2 * beta)
!        endif

    else ! nint(PI%Nparvar) > N_before_mv*PI%npars .and. tmp > beta

       ! is this step a multivariate proposal or not
       multivariate_proposal = .false.
       ! track information on whether this is a beta step proposal
!       N%ITER_beta = N%ITER_beta + 1d0

       ! Sample random normal distribution (mean = 0, sd = 1)
       do p = 1, PI%npars
          call random_normal(uniform,unif_length,uniform_random_vector,rn2(p))
       end do !

       ! Estimate the step to be applied to the current parameter vector to
       ! create the new proposal. The Beta and Delta steps are intended to
       ! provide a Gaussian proposal but of a small (beta) or large (delta) size 
       ! Draw from uniform distribution to determine which step size is used
!       tmp = uniform_random_vector(uniform) ; uniform = uniform + 1
!       if (tmp < delta) then
!           ! Delta Step
!           norpars = norpars0 + (par_maxstepsize*rn2)
!       else
!           ! Beta Step
           norpars = norpars0 + (par_minstepsize*rn2)
!       endif

    end if ! nint(PI%Nparvar) > 10*PI%npars .and. tmp > beta

    ! reverse normalisation on the new parameter step
    do p = 1, PI%npars
       call nor2par(1,norpars(p),PI%parmin(p),PI%parmax(p),pars(p))
!       call log_nor2par(1,norpars(p),PI%parmin(p),PI%parmax(p),pars(p))
    end do

  end subroutine step
  !
  !------------------------------------------------------------------
  !
end module MHMCMC_module
