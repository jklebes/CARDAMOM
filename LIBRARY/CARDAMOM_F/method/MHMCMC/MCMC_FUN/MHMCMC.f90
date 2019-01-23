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

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine MHMCMC (model_likelihood_option,PI,MCO,MCOUT)
!TLS    use, intrinsic :: ieee_arithmetic
    use MCMCOPT, only: MCMC_OUTPUT, MCMC_OPTIONS, PARAMETER_INFO, COUNTERS
    use math_functions, only: randn, random_uniform
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
!TLS        use, intrinsic :: ieee_arithmetic
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

    double precision, dimension(PI%npars*MCO%nADAPT) :: PARSALL ! All accepted parameters since previous step adaption
    double precision :: infini &
                       ,crit1,crit2 & ! random numbers log(0->1) used to accept / reject
                       ,P0prior, Pprior & ! as below but for priors only
                       ,P0 & ! previously accepted observation based log-likelihood
                       ,P    ! current observation based log-likelihood
    integer :: i

    ! initial values
    P = -1d0 ; Pprior = -1d0
    uniform = 1
    N%ACC = 0 ; N%ITER = 0
    N%ACCLOC = 0 ; N%ACCRATE = 0d0

    ! calculate initial vector of uniform random values
    allocate(uniform_random_vector(MCO%nOUT))
    call random_uniform(uniform_random_vector,size(uniform_random_vector))
    ! assume this is a restart which we want to load previous values
    if (restart_flag) then
        N%ACC = accepted_so_far
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
           PI%parini(i) = nor2par(uniform_random_vector(uniform),PI%parmin(i),PI%parmax(i))
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
    write(*,*) "Starting likelihood = ",(P0+P0prior)

    ! checks whether the EDCs (combined with P0 not P0prior) have been met in the initial parameter set
    infini = 0d0
    if (P0 == log(infini)) then
        write(*,*) "WARNING! P0 = ",P0," - MHMCMC may get stuck, if so please check initial conditins"
    endif

    ! begin the main MHMCMC loop
    do while (N%ACC < MCO%nOUT .and. (P < 0d0 .or. MCO%nWRITE > 0))

       ! take a step in parameter space
       call step(PARS0,PARS,PI)
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

       ! determine accept or reject, should this criterion be > crit + ? (or
       ! other limit, such as p = 0.05 -> -3)
!print*,(P-P0)
       if ((P-P0) >= crit1 .and. (Pprior-P0prior) >= crit2) then

          ! store accepted parameter solutions
          ! keep record of all parameters accepted since step adaption
          PARSALL(((N%ACCLOC*PI%npars)+1):((N%ACCLOC*PI%npars)+PI%npars)) = PARS(i)
          PARS0(1:PI%npars) = PARS(1:PI%npars)
          ! specifically store the best parameter set
!          if (P > P0) BESTPARS = PARS
          if ((P+Pprior) > (P0+P0prior)) BESTPARS = PARS
          ! keep track of how many accepted solutions (global and local)
          N%ACC = N%ACC + 1 ; N%ACCLOC = N%ACCLOC + 1 
          P0 = P ; P0prior = Pprior

          ! write out parameter, log-likelihood and step if appropriate
          if (MCO%nWRITE > 0 .and. mod(N%ACC,MCO%nWRITE) == 0) then
             call write_results(PARS,(P+Pprior),PI)
          end if ! write or not to write

       endif ! accept or reject condition

       ! count iteration whether accepted or rejected
       N%ITER = N%ITER + 1

       ! time to adapt?
       if (mod(N%ITER,MCO%nADAPT) == 0) then
           ! work out local acceptance rate (i.e. since last adapt)
           N%ACCRATE = dble(N%ACCLOC) / dble(MCO%nADAPT)
!print*,"acceptance = ",N%ACCRATE
           ! have few enough parameters been accepted to consider adapting
           if ((MCO%fADAPT*dble(MCO%nOUT)) > dble(N%ACC)) then
               call adapt_step_size(PARSALL,PI,N,MCO)
           end if !  have enough parameter been accepted
           ! resets to local counter
           N%ACCLOC = 0
       end if ! time to adapt?

       ! should I be write(*,*)ing to screen or not?
       if (MCO%nPRINT > 0 .and. (mod(N%ITER,MCO%nPRINT) == 0)) then
           write(*,*)"Total accepted = ",N%ACC," out of ",MCO%nOUT
           write(*,*)"Local acceptance rate ",N%ACCRATE*100
           write(*,*)"Current observation log-likelihood = ",P0
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
    write(*,*)"MHMCMC loop completed, next please..."

  end subroutine MHMCMC
  !
  !------------------------------------------------------------------
  !
  subroutine adapt_step_size(PARSALL,PI,N,MCO)
    use MCMCOPT, only: MCMC_OPTIONS, COUNTERS, PARAMETER_INFO
    use math_functions, only: std

    implicit none

    ! declare input types
    type ( parameter_info ), intent(inout) :: PI
    type ( counters ), intent(inout) :: N
    type ( mcmc_options ), intent(inout) :: MCO

    ! declare inputs variables
    double precision, intent(in) :: PARSALL(MCO%nADAPT*PI%npars) ! collection of recently accepted parameter combinations

    ! declare local variables
    integer p,i ! counters
    double precision minstepsize         & ! minimum step size
                    ,norparstd           & ! normalised parameter value standard deviation
                    ,norparvec(N%ACCLOC)   ! normaised parameter values
    double precision, parameter :: fac = 2d0,   & ! factor used to determine whether to reduce step size
                                 fac_1 = 0.5d0, & ! ...and its inverse
                              adaptfac = 1.5d0, & ! fraction applied to reduce / increase step size...
                            adaptfac_1 = 0.6666667d0, & ! ...and its inverse
                         sqrt_adaptfac = 1.224745d0


    ! calculate constants
    minstepsize = 10000d0/dble(N%ITER)
    if (minstepsize > 0.01d0) minstepsize = 0.01d0
!    minstepsize = 100d0/dble(N%ACC) ! should this be linked to number of accepted parameter?
!    if (minstepsize > 0.01d0) minstepsize = 0.01d0

    ! determine local acceptance rate
    N%ACCRATE = dble(N%ACCLOC)/dble(MCO%nADAPT)
!print*,"...................................",minstepsize,N%ACCRATE
    ! default stepsize increment
    if (N%ACCLOC > 0 .and. N%ACCRATE < 0.23d0) then
        ! make step size smaller
        PI%stepsize = PI%stepsize * adaptfac_1
    else if (N%ACCRATE > 0.44d0) then
        ! make step size bigger
        PI%stepsize = PI%stepsize * adaptfac
    end if ! conditional if acceptance rate low or high

    ! Next do dimension / parameter specific adjustments
    ! this is the adaptive part (Bloom & Williams 2015)
    ! NOTE: original value was > 3, however this result in a biased estimate of
    ! the standard deviation to a lower value.
    if (N%ACCLOC > 5 .and. N%ACCRATE < 0.23d0) then
        do p = 1, PI%npars
           do i = 1, N%ACCLOC
              norparvec(i) = par2nor(PARSALL((PI%npars*(i-1))+p),PI%parmin(p),PI%parmax(p))
           end do ! i
           ! calculate standard deviation (variability)
           norparstd = std(norparvec,N%ACCLOC)
           ! if stepsize is smaller than the variability in accepted parameters,
           ! and acceptance rate is low, this indicates that we are stuck in a
           ! poor local minima and should moderate the reduction in step size
           if (PI%stepsize(p) < norparstd*fac_1) then
               PI%stepsize(p) = PI%stepsize(p)*(sqrt_adaptfac)
           endif
        end do ! p
    endif ! if N%ACCLOC > 10

    ! keep track of how long we have been stuck somewhere
    ! if (N%ACCLOC == 0) N%ACCLOC_ZEROS = N%ACCLOC_ZEROS + 1

    !!!!!!!
    ! carry out final checks
    !!!!!!!

    ! step size can't be greater than 1
    where (PI%stepsize > 1d0) PI%stepsize = PI%stepsize * adaptfac_1
    ! if stepsize below minimum allowed value increase
    where (PI%stepsize < minstepsize) PI%stepsize = PI%stepsize * adaptfac
    ! if stepsize still below minimum allowed value then set to minimum
    where (PI%stepsize < minstepsize) PI%stepsize = minstepsize

  end subroutine adapt_step_size
  !
  !------------------------------------------------------------------
  !
  double precision function par2nor(initial_par,min_par,max_par)

    ! functions to normalised log parameter values and return them back to
    ! unlogged / un-normalised value

    ! converting parameters on log scale between 0-1 for min/max values
    implicit none
    double precision initial_par, min_par, max_par

    if (max_par > 0d0 .and. min_par < 0d0) then
        ! then normalise without logs and we cross zero
        par2nor = (initial_par-min_par)/(max_par-min_par)
    else
        par2nor = log(initial_par/min_par)/log(max_par/min_par)
    end if

    ! explicit return
    return

  end function par2nor
  !
  !---------------------and vise versa ------------------------------
  !
  double precision function nor2par(initial_par,min_par,max_par)

    ! Converting values back from normalised (0-1) to 'real' numbers

    implicit none
    double precision initial_par, min_par, max_par

    ! determine whether we used log normalisation or not
    if ( max_par > 0d0 .and. min_par < 0d0) then
        ! ...then un-normalise without logs as we cross zero and logs wont work
        nor2par = min_par+(max_par-min_par)*initial_par
    else
        ! ...then un-normalised assuming logs-normalisation was used
        nor2par = min_par*((max_par/min_par)**initial_par)
    endif

    ! explicit return
    return

  end function nor2par
  !
  !------------------------------------------------------------------
  !
  subroutine step (pars0,pars,PI)
    use math_functions, only: randn, random_normal
    use MCMCOPT, only: PARAMETER_INFO

    ! carries out the next step to parameters in the MCMC search

    implicit none

    ! declare input variables
    type ( parameter_info ), intent(in) :: PI
    double precision, dimension(PI%npars), intent(inout) :: pars0 & ! current set of parameters
                                                 ,pars    ! next set of parameters

    ! declare local variables
    integer :: n, fp
    double precision  :: npar, npar_new, rn

    ! default values
    npar = 0d0 ; rn = 0d0

    ! begin sampling parameter space
    do n = 1, PI%npars
      fp = 0
      ! normalise parameters first
      npar = par2nor(pars0(n),PI%parmin(n),PI%parmax(n))
      ! then apply step...
      do while (fp == 0)
         ! get a normally distributed random number
!         rn = randn(1)
         call random_normal(uniform, size(uniform_random_vector), uniform_random_vector, rn)
         ! apply to our current parameter value
         npar_new = npar + (rn*PI%stepsize(n)*(1d0-PI%parfix(n)))
!         npar_new = min(1d0,max(0d0,npar_new))
         ! ensure the new parameter value is contrained between 0 and 1
         if (npar_new > 0d0 .and. npar_new < 1d0) then
            fp = 1
            pars(n) = nor2par(npar_new,PI%parmin(n),PI%parmax(n))
         end if
      end do ! while conditions
    end do ! parameter loop

  end subroutine step
  !
  !------------------------------------------------------------------
  !
end module MHMCMC_module
