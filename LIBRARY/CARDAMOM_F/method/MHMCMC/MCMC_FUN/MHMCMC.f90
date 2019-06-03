module MHMCMC_MODULE

! module contains all subroutine and functions relevant specifically to the
! MHMCMC method. The choice of EDC, likelihood and model are made else where and
! are thus contains within a seperate module

implicit none

! make all private
private

! specify what can be seen
public :: MHMCMC

! declare any module level variables needed

! related to random number generator
integer :: uniform
double precision, allocatable, dimension(:) :: uniform_random_vector
! MHMCMC step size
double precision :: minstepsize = 0.0001d0

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine MHMCMC (model_likelihood_option,PI,MCO,MCOUT)
!TLS    use, intrinsic :: ieee_arithmetic
    use MCMCOPT, only: MCMC_OUTPUT, MCMC_OPTIONS, PARAMETER_INFO, COUNTERS
    use math_functions, only: randn, random_uniform, par2nor, nor2par
    use cardamom_io, only: write_results,restart_flag,accepted_so_far
    implicit none

    !/* ***********INPUTS************
    ! *
    ! * MODEL_LIKELYHOOD: A function wholly responsible for
    ! * (a) running the model given the DATA and parameters,
    ! * (b) comparing it to observations,and
    ! * (c) returning  the (log) likelyhood.
    ! * The function will be run as MODEL_LIKELYHOOD(DATA,PI,PARS);
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
    ! * RESULTS FILE: File includes (a) results (b) likelyhood and (c) final step
    ! size
    ! *
    ! * */

    ! declare interface for the model likelihood function.
    ! NOTE that the inputted MODEL_LIKEIHOOD_OPTION could be multiple subroutines,
    ! the interface allows for making the requirements of this explicit
    interface
      subroutine model_likelihood_option(param_info, param_vector, ML_obs_out, ML_prior_out)
        use cardamom_structures, only: DATAin, emulator_pars
        use MCMCOPT, only: PARAMETER_INFO
        use CARBON_MODEL_MOD, only: carbon_model
           implicit none
           ! declare input variables
           type ( parameter_info ), intent(inout) :: param_info
           double precision, dimension(param_info%npars), intent(inout) :: param_vector
           ! output
           double precision,intent(inout) :: ML_obs_out, ML_prior_out
      end subroutine model_likelihood_option
    end interface

    ! declare input variables
    type ( parameter_info ), intent(inout) :: PI
    type ( mcmc_options ), intent(inout) :: MCO
    type ( mcmc_output ), intent(inout) :: MCOUT

    ! declare any local variables
    type ( counters ) :: N
    double precision, dimension(PI%npars) :: PARS0    & ! parameter values in initial / reference iteration
                                            ,PARS     & ! parameter values in current interation
                                            ,BESTPARS   ! best set of parameters so far

    double precision, dimension(PI%npars,MCO%nADAPT) :: PARSALL ! All accepted parameters since previous step adaption
    double precision :: infini &
                       ,crit1,crit2 & ! random numbers log(0->1) used to accept / reject
                       ,Pmax, P0prior, Pprior & ! as below but for priors only
                       ,P0 & ! previously accepted observation based log-likelihood
                       ,P    ! current observation based log-likelihood
    integer :: i

    ! initial values
    P = -1d0 ; Pprior = -1d0
    uniform = 1
    N%ACC = 0 ; N%ITER = 0 ; N%ACCEDC = 0
    N%ACCLOC = 0 ; N%ACCRATE = 0d0
    N%beta_step = .false. ; N%Nbeta = 0 ; N%ACCLOC_beta = 0 ; N%ACCRATE_beta = 0d0

    ! calculate initial vector of uniform random values
    allocate(uniform_random_vector(MCO%nOUT))
    call random_uniform(uniform_random_vector,size(uniform_random_vector))
    ! assume this is a restart which we want to load previous values
    if (restart_flag) then
        N%ACC = accepted_so_far
        N%ACCEDC = N%ACC
        N%ITER = nint(dble(N%ACC) / 0.23d0) ! approximation of current iteration when restarting
    endif

    ! add something here to delete previous files if wanted later
    if (MCO%APPEND == 0 .and. MCO%nWRITE > 0) then
       write(*,*) "Oooops have requested that existing files be deleted but you have not finished the code to do so...."
    end if

    ! start random sampling if MCO%randparini set
    do i = 1, PI%npars
       ! parfix = 1 stay at prior value, parfix = 0 randomly search
       if (.not.MCO%fixedpars) PI%parfix(i) = 0d0 ! also impact STEP function
       ! only assign random parameters if (a) randparini == 1 or (b) PI$parini(n) == -9999)
       if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
           call nor2par(1,uniform_random_vector(uniform),PI%parmin(i),PI%parmax(i),PI%parini(i))
!           PI%parini(i) = nor2par(randn(0),PI%parmin(i),PI%parmax(i))
           uniform = uniform + 1
       end if
       ! write(*,*) parameter values to screen
       write(*,*) "p",i,"=",PI%parini(i)
    end do ! for PI%npar loop

    ! inform the user
    write(*,*) "Have loaded / randomly assigned PI%parini - now begin the MHMCMC"

    ! initialise the prior and best pars vectors
    PARS0(1:PI%npars) = PI%parini(1:PI%npars)
    BESTPARS(1:PI%npars) = PI%parini(1:PI%npars)

    ! calculate the initial probability / log likelihood.
    ! NOTE: passing P0 -> P is needed during the EDC searching phase where we
    ! could read an EDC consistent parameter set in the first instance
    call model_likelihood_option(PI, PI%parini,P0,P0prior) ; P = P0 ; Pprior = P0prior
    write(*,*) "Starting likelihood = ",P0,"+",P0prior
    Pmax = P0+P0prior

    ! checks whether the EDCs (combined with P0 not P0prior) have been met in the initial parameter set
    infini = 0d0
    if (P0 == log(infini)) then
        write(*,*) "WARNING! P0 = ",P0," - MHMCMC may get stuck, if so please check initial conditins"
    endif

    ! begin the main MHMCMC loop
!    do while (N%ACC < MCO%nOUT .and. (P < 0d0 .or. MCO%nWRITE > 0))
!    do while (N%ITER < MCO%nOUT .and. (P < 0d0 .or. MCO%nWRITE > 0))
    do while (N%ACCEDC < MCO%nOUT .and. (P < 0d0 .or. MCO%nWRITE > 0))

       ! take a step in parameter space
       call step(PARS0,PARS,PI,N)
       ! calculate the model likelihood
       call model_likelihood_option(PI, PARS, P, Pprior)

       ! accept or reject, draw uniform distribution (0,1)
       crit1 = log(uniform_random_vector(uniform)) !crit=log(randn(0))
       uniform = uniform + 1
       crit2 = log(uniform_random_vector(uniform))
       uniform = uniform + 1
       ! if we are near to the end re-generate some more values
       if (uniform >= size(uniform_random_vector)-4) then
           ! calculate new vector of uniform random values
           call random_uniform(uniform_random_vector,size(uniform_random_vector))
           ! and reset uniform counter
           uniform = 1
       endif

       ! determine accept or reject the current proposal
       if ((P-P0) > crit1 .and. (Pprior-P0prior) > crit2) then
!       if ( ( (P+Pprior)-(P0+P0prior) ) >= crit1) then

          ! Store accepted parameter solutions
          ! keep record of all parameters accepted since step adaption
          PARSALL(1:PI%npars,(N%ACCLOC+1)) = PARS(1:PI%npars)
          PARS0(1:PI%npars) = PARS(1:PI%npars)
          ! specifically store the best parameter set
          if ((P+Pprior) >= Pmax) then
              BESTPARS = PARS ; Pmax = P+Pprior
          endif
          ! keep track of how many accepted solutions (global and local)
          if (N%beta_step) N%ACCLOC_beta = N%ACCLOC_beta + 1
          N%ACC = N%ACC + 1 ; N%ACCLOC = N%ACCLOC + 1
          P0 = P ; P0prior = Pprior

          ! write out parameter, log-likelihood and step if appropriate
!          if (MCO%nWRITE > 0 .and. mod(N%ACC,MCO%nWRITE) == 0) then
!             call write_results(PARS0,(P0+P0prior),PI)
!          end if ! write or not to write

       endif ! accept or reject condition

       ! write out parameter, log-likelihood and step if appropriate
!       if (MCO%nWRITE > 0 .and. mod(N%ITER,MCO%nWRITE) == 0) then
!          call write_results(PARS0,(P0+P0prior),PI)
!       end if ! write or not to write

       ! count iteration whether accepted or rejected
       N%ITER = N%ITER + 1
       ! count how many EDC compatible iterations there have been
       if ((P+Pprior) > log(infini)) then
           N%ACCEDC = N%ACCEDC + 1
           if (MCO%nWRITE > 0 .and. mod(N%ACCEDC,MCO%nWRITE) == 0) then
              call write_results(PARS0,(P0+P0prior),PI)
           end if ! write or not to write
       endif

       ! time to adapt?
       if (mod(N%ITER,MCO%nADAPT) == 0) then
           ! work out local acceptance rate (i.e. since last adapt)
           N%ACCRATE = dble(N%ACCLOC) / dble(MCO%nADAPT)
           N%ACCRATE_beta = dble(N%ACCLOC_beta) / dble(N%Nbeta)
           ! have few enough parameters been accepted to consider adapting
           if ((MCO%fADAPT*dble(MCO%nOUT)) > dble(N%ACC)) then
               call adapt_step_size(PARSALL,PI,N,MCO)
           end if !  have enough parameter been accepted
           ! resets to local counter
           N%ACCLOC = 0 ; N%ACCLOC_beta = 0 ; N%Nbeta = 0
       end if ! time to adapt?

       ! should I be write(*,*)ing to screen or not?
       if (MCO%nPRINT > 0 .and. (mod(N%ITER,MCO%nPRINT) == 0)) then
           write(*,*)"Total accepted = ",N%ACC," out of ",MCO%nOUT
           write(*,*)"Overall acceptance rate  = ",dble(N%ACC) / dble(N%ITER)
           write(*,*)"Local   acceptance rate  =",N%ACCRATE
           write(*,*)"Current maximum stepsize =",maxval(PI%stepsize)
           write(*,*)"Current obs   log-likelihood = ",P0
           write(*,*)"Current prior log-likelihood = ",P0prior
       end if ! write(*,*) to screen or not

    end do ! while conditions

    ! fill MCOUT details
    MCOUT%best_pars(1:PI%npars) = BESTPARS(1:PI%npars)
    ! set flag MCMC completed
    MCOUT%complete = 1
    ! tidy up
    deallocate(uniform_random_vector)

    ! completed MHMCMC loop
    write(*,*)"MHMCMC loop completed"
    write(*,*)"Final acceptance rate = ",dble(N%ACC) / dble(N%ITER)

  end subroutine MHMCMC
  !
  !------------------------------------------------------------------
  !
  subroutine adapt_step_size(PARSALL,PI,N,MCO)
    use MCMCOPT, only: MCMC_OPTIONS, COUNTERS, PARAMETER_INFO
    use math_functions, only: nor2par, par2nor, &
                              std, covariance_matrix, increment_covariance_matrix

    implicit none

    ! declare input types
    type ( parameter_info ), intent(inout) :: PI
    type ( counters ), intent(inout) :: N
    type ( mcmc_options ), intent(inout) :: MCO

    ! declare inputs variables
    double precision, intent(in) :: PARSALL(PI%npars,MCO%nADAPT) ! collection of recently accepted parameter combinations

    ! declare local variables
    integer p,i ! counters
    double precision, dimension(PI%npars,MCO%nADAPT) :: norparvec ! normaised parameter values
    double precision, dimension(PI%npars,PI%npars) :: local_covariance
    double precision :: dble_accloc &
                       ,norparstd     ! normalised parameter value standard deviation
    double precision, parameter :: fac = 2d0,   & ! factor used to determine whether to reduce step size
                                 fac_1 = 0.5d0, & ! ...and its inverse
                              adaptfac = 1.25d0, & ! fraction applied to reduce / increase step size...
                            adaptfac_1 = 0.80d0, & ! ...and its inverse
                         sqrt_adaptfac = 1.118033988749d0

    ! Multiple use conversion
    dble_accloc = dble(N%ACCLOC)
    ! calculate adaption factor, the scaling potential reduces the further along the analysis we are
!    adaptfac=(1d0-((dble(N%ACC)/dble(MCO%nOUT))*0.5d0))*0.001d0+1d0
    ! determine local acceptance rate
    N%ACCRATE = dble_accloc/dble(MCO%nADAPT)
    ! calculate new minimum stepsize
    minstepsize = min(0.01,10000d0/dble(N%iter))
!    print*,"...minstep = ",minstepsize," ACCRATE = ",N%ACCRATE
!    print*," minval(step) = ",minval(PI%stepsize)!," maxval(step) = ",maxval(PI%stepsize)
!    print*," minval(pstd) = ",minval(PI%parstd)," maxval(pstd) = ",maxval(PI%parstd)

    ! default stepsize increment
!    if (N%ACCLOC > 0 .and. nint(PI%Nparstd) < 10*PI%npars) then
!        if (N%ACCRATE < 0.23d0) then
!            ! make step size smaller
!            PI%stepsize = PI%stepsize * adaptfac_1
!        else !if (N%ACCRATE > 0.44d0) then
!            ! make step size bigger
!            PI%stepsize = PI%stepsize * adaptfac
!        end if ! conditional if acceptance rate low or high
!    end if ! N%ACCLOC > 0

    ! default stepsize increment
    if (N%ACCLOC_beta > 0) then
        if (N%ACCRATE_beta < 0.23d0) then
            ! make step size smaller
            PI%stepsize = PI%stepsize * adaptfac_1
        else !if (N%ACCRATE > 0.44d0) then
            ! make step size bigger
            PI%stepsize = PI%stepsize * adaptfac
        end if ! conditional if acceptance rate low or high
    end if ! N%ACCLOC > 0

    ! Next do dimension / parameter specific adjustments
    ! this is the adaptive part (Bloom & Williams 2015)

    if (N%ACCLOC > 0) then

        ! normalise all parameters from current vector
        do p = 1, PI%npars
           call par2nor(N%ACCLOC,PARSALL(p,1:N%ACCLOC),PI%parmin(p),PI%parmax(p),norparvec(p,1:N%ACCLOC))
        end do

        ! if we have a covariance matrix then we want to update it, if not then we need to create one
        if (PI%cov) then

            ! Increment the variance-covariance matrix with new accepted parameter sets
            ! NOTE: that this also increments the total accepted counter (PI%Nparstd)
            call increment_covariance_matrix(norparvec(1:PI%npars,1:N%ACCLOC),PI%mean_par,PI%npars &
                                            ,PI%Nparstd,N%ACCLOC,PI%covariance)

        else ! PI%cov == .true.

            ! we have not yet created a covariance matrix based on accepted
            ! parameters. Assuming we have some then create one...
            if (N%ACCLOC > 3) then
                ! estimate covariance matrix
                call covariance_matrix(norparvec(1:PI%npars,1:N%ACCLOC),PI%mean_par,PI%npars,N%ACCLOC,local_covariance)
                PI%cov = .true. ; PI%Nparstd = dble_accloc
                ! replace the initial matrix with current estimate
                PI%covariance = local_covariance
            end if ! N%ACCLOC > 3

        end if ! PI%cov == .true.

        ! adjust step size by local standard deviation
        do p = 1, PI%npars
           ! calculate standards deviation (variability) for the local window
           ! extracted from the variance
           PI%parstd(p) = sqrt(PI%covariance(p,p))
        end do ! p

    endif ! if N%ACCLOC > 0

!    ! Next do dimension / parameter specific adjustments
!    ! this is the adaptive part (Bloom & Williams 2015)
!    if (N%ACCLOC > 3 .and. N%ACCRATE < 0.23d0) then
!        do p = 1, PI%npars
!!           do i = 1, N%ACCLOC
!!              norparvec(i) = par2nor(PARSALL((PI%npars*(i-1))+p),PI%parmin(p),PI%parmax(p))
!!           end do
!           call par2nor(N%ACCLOC,PARSALL(p,1:N%ACCLOC),PI%parmin(p),PI%parmax(p),norparvec(1:N%ACCLOC))
!           ! calculate standards deviation (variability) for the local window
!           norparstd = std(norparvec,N%ACCLOC)
!           ! weight the current SD with the existing history
!           PI%parstd(p) = (PI%parstd(p) * PI%Nparstd) + (norparstd * dble(N%ACCLOC))
!           PI%parstd(p) = PI%parstd(p) / (PI%Nparstd + dble(N%ACCLOC))
!           ! if stepsize is smaller than the variability in accepted parameters,
!           ! and acceptance rate is low, this indicates that we are stuck in a
!           ! poor local minima and should moderate the reduction in step size
!!           if (PI%stepsize(p) < norparstd*fac_1) then
!           if (PI%stepsize(p) < PI%parstd(p)*fac_1) then
!               PI%stepsize(p) = PI%stepsize(p)*(sqrt_adaptfac)
!           endif
!        end do ! p
!        ! update for next iterations
!        PI%Nparstd = PI%Nparstd + dble(N%ACCLOC)
!    endif ! if N%ACCLOC > 3

    !!!!!!!
    ! carry out final checks
    !!!!!!!

    ! step size can't be greater than 1
!    where (PI%stepsize > 1d0) PI%stepsize = PI%stepsize / adaptfac
    ! if stepsize below minimum allowed value increase
    where (PI%stepsize < minstepsize) PI%stepsize = PI%stepsize * adaptfac
!    if (minval(PI%stepsize) < minstepsize) PI%stepsize = PI%stepsize * adaptfac
    ! if stepsize still below minimum allowed value then set to minimum
    !where (PI%stepsize < minstepsize) PI%stepsize = minstepsize

  end subroutine adapt_step_size
  !
  !------------------------------------------------------------------
  !
  subroutine step(pars0,pars,PI,N)
    use math_functions, only: par2nor, nor2par, &
                              random_normal, random_multivariate, random_uniform
    use MCMCOPT, only: PARAMETER_INFO, COUNTERS

    ! carries out the next step to parameters in the MCMC search

    implicit none

    ! declare input variables
    type ( parameter_info ), intent(in) :: PI
    type ( counters ), intent(inout) :: N
    double precision, dimension(PI%npars), intent(inout) :: pars0 & ! current set of parameters
                                                           ,pars    ! next set of parameters

    ! declare local variables
    integer :: p, fp
    !double precision  :: npar(1), npar_new(1), rn(PI%npars), mu(PI%npars)
    double precision  :: npar(PI%npars), npar_new(PI%npars), stepping(PI%npars) &
                        ,rn(PI%npars), mu(PI%npars), rn2(PI%npars), scd, Id, tmp
    double precision, parameter :: beta = 0.05d0

    ! reset values
    npar = 0d0 ; npar_new = 0d0 ; rn = 0d0 ; mu = 0d0
    Id = sqrt(dble(PI%npars)) ; scd = 2.381204d0 / Id

!    ! Now iterate through the parameters updating them in turn
!    do n = 1, PI%npars
!       fp = 0
!       ! normalise parameters first
!       call par2nor(1,pars0(n),PI%parmin(n),PI%parmax(n),npar(1))
!       ! then apply step...
!       do while (fp == 0)
!          ! get a normally distributed random number
!          call random_normal(uniform, size(uniform_random_vector), uniform_random_vector, rn)
!          ! apply to our current parameter value
!!          npar_new(1) = npar(1) + (rn*PI%stepsize(n)*PI%parstd(n)*(1d0-PI%parfix(n)))
!          npar_new(1) = npar(1) + (rn(n)*PI%stepsize(n)*(1d0-PI%parfix(n)))
!!          npar_new(1) = min(1d0,max(0d0,npar_new(1)))
!          ! ensure the new parameter value is contrained between 0 and 1
!          if (npar_new(1) >= 0d0 .and. npar_new(1) <= 1d0) then
!              fp = 1
!              call nor2par(1,npar_new(1),PI%parmin(n),PI%parmax(n),pars(n))
!          end if
!       end do ! while conditions
!    end do ! parameter loop

    ! Begin sampling parameter space, first estimate multivariate random number
    ! Multivariate sample around a mean of zero

    ! Normalise parameters first
    do p = 1, PI%npars
       call par2nor(1,pars0(p),PI%parmin(p),PI%parmax(p),npar(p))
    end do !

    tmp = uniform_random_vector(uniform)
    uniform = uniform + 1
    ! if we are near to the end re-generate some more values
    if (uniform >= size(uniform_random_vector)-4) then
        ! calculate new vector of uniform random values
        call random_uniform(uniform_random_vector,size(uniform_random_vector))
        ! and reset uniform counter
        uniform = 1
    endif

    ! Increment step size via different proposal based on whether we have a sufficiently development covariance matrix
    ! Splitting step calculation based on number of parameter vectors accepted
    ! is linked to the need build a covarianc matrix prior to multivariate
    ! sampling.
    ! See Roberts and Rosenthal, Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367, 2009.
    if (nint(PI%Nparstd) > 10*PI%npars .and. tmp > beta) then

        fp = 0 ; N%beta_step = .false.
        do while (fp == 0)

           ! Get random normal mean = 0, sd = 1
           do p = 1, PI%npars
              call random_normal(uniform, size(uniform_random_vector),uniform_random_vector, rn2(p))
           end do !
           ! Draw from multivariate random distribution
           ! NOTE: if covariance matrix provided is not positive definite
           !       a sample form normal distribution is returned
           call random_multivariate(PI%npars, 1, uniform, size(uniform_random_vector), &
                                    uniform_random_vector, PI%covariance, mu, rn)

           ! Estimate the step to be applied to the current parameter vector to
           ! create the new proposal. scd = a scaling parameter linking searching
           ! stepping to the number of parameters being retrieved by the analysis.
           ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
           ! and references therein.
           stepping = rn*(1d0-PI%parfix)*scd
           npar_new = npar + stepping
           ! ensure the new parameter value is contrained between 0 and 1
           if (minval(npar_new) > 0d0 .and. maxval(npar_new) < 1d0) fp = 1

        end do ! while conditions

    else ! nint(PI%Nparstd) > 2*PI%npars

        fp = 0 ; N%beta_step = .true. ; N%Nbeta = N%Nbeta + 1
        do while (fp == 0)

          ! Get random normal mean = 0, sd = 1
          do p = 1, PI%npars
             call random_normal(uniform, size(uniform_random_vector),uniform_random_vector, rn2(p))
          end do !

          ! Estimate the step to be applied to the current parameter vector to
          ! create the new proposal. scd = a scaling parameter linking searching
          ! stepping to the number of parameters being retrieved by the analysis.
          ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
          ! and references therein.
          stepping = (PI%stepsize*rn2/Id)*(1d0-PI%parfix)
          npar_new = npar + stepping
          ! ensure the new parameter value is contrained between 0 and 1
          if (minval(npar_new) > 0d0 .and. maxval(npar_new) < 1d0) fp = 1

        end do ! while conditions

    end if ! nint(PI%Nparstd) > 2*PI%npars

    ! reverse normalisation on the new parameter step
    do p = 1, PI%npars
       call nor2par(1,npar_new(p),PI%parmin(p),PI%parmax(p),pars(p))
    end do

! SHOULD THE EDCS BE USED TO APPLY A COST ON THE OUTPUTS RATHER THAN BLOCKAGE?

  end subroutine step
  !
  !------------------------------------------------------------------
  !
end module MHMCMC_module
