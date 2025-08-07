module samplers_shared
implicit none
!> A collection of info about the model's parameters
!> npars and min, max bounds as two arrays
!> Closely related to model fct; model fct must take this number
!> and type of paramters
!> Unlike in previous versions, intended to be intent(in) only
!> TODO we could bundle this in a type with the function:
!> the fortran quasi-object
   type PARINFO
      integer:: npars
      double precision, allocatable, dimension(:):: parmin, parmax, paradj
      logical, allocatable, dimension(:):: fix_pars
      logical, allocatable, dimension(:):: parfix
   end type PARINFO

!> A collection of input options to the MCMC sampler run
!> contains default values
   type MCMC_options
      integer:: MAXITER = 10000  ! overall steps, if convergence not reached
      integer:: nadapt = 1000  ! steps per "local" sampling period, between adaptation steps
      integer:: Nchains = 1  ! consider setting OMP env to something compatible
      integer:: nwrite = 1000
      integer:: nprint = 1000
      integer:: nout
      real:: P_target  
      !! termination criterion-a loglikelihood to stop at (optional)
!> file names
      character(350):: outfile = "parout.txt"
      character(350):: stepfile = "stepout.txt"
      character(350) ::  covfile = "covout.txt"
      character(350):: covifile = "covinfoout.txt"
      logical:: append
      real:: fadapt  ! TODO fraction adapt-move to outside
      logical:: randparini
      logical:: returnpars  
      !! a variable that is never used and has no effect, needs deleting in all model likelihood files
      logical:: restart = .false.
      !! is it a restart ?
      logical:: fixedpars  
      !! never used-yet-for array of flags to hold some parameters constant
! Adaptive
!> setting for adaptive AP-MCMC step size
      double precision:: par_minstepsize = 0.001d0 & ! 0.0005 -> 0.001 -> 0.01 -> 0.1 -> 0.005
                                           , par_maxstepsize = 0.01d0 &
                                                               , par_initstepsize = 0.005d0
      double precision:: beta = 0.05d0  ! weighting for gaussian step in multivariate proposals
!> Optimal scaling variable for parameter searching
      double precision:: opt_scaling_const = 2.381204**2  ! scd = 2.381204 the optimal scaling parameter
      ! for MCMC search, when applied to  multivariate proposal.
      ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky
      ! factor. NOTE 2: 2.381204**2 = 5.670132
      double precision:: N_before_mv = 10d0
!! step
!> Is current proposal multivariate or not?
      logical:: multivariate_proposal = .false.
      logical:: use_multivariate
   end type MCMC_OPTIONS

!> Collection of info for output of the sampling run
!> , can also be passed to next run to continue from the last state
!> Note output is mainly via file writing
   type MCMC_OUTPUT
      double precision:: bestll
      !! best (maximum) loglikelihood value found so far
      double precision:: ll 
      !! latest loglikelihood value
      double precision, allocatable, dimension(:):: bestpars
      !! best (loglikelihood-maximizing) parameter values found so far
      double precision, allocatable, dimension(:):: pars
      !! latest parameter values
      double precision:: acceptance_rate
      !! acceptance rate
      logical:: complete
      !! Did the main loop finish?
      integer:: nos_iterations
      !! number main loop interations run so far
!stats collection:
      integer:: Nparvar, Nparvar_local
      !! Number of states in history that have gone into running
      !! mean, variance, and covariance calculation.
      double precision, allocatable, dimension(:):: parvar
      !! variance
      double precision, allocatable, dimension(:):: meanpar
      !! mean of parameters
      double precision, allocatable, dimension(:, :):: covariance
      !! covariance matrix measured during sampling
      logical:: cov = .false.
      !! Does the covariance matrix exist yet?
      logical:: use_multivariate
      !! Are we in the later simulation phase where step size depend
      !! on covariance matrix ?  i.e. after enough data has been observed
      !! that a useful covariance matrix exists
      logical:: multivariate_proposal
      !! TODO general setting to ever use multivariate or not?
   end type MCMC_OUTPUT


contains

!! utilities

   logical function is_infinity(ll)
!! Check whether a loglikelihood is infinity/likelihood is zero
!! usually signalling hard reject of the state due to boundary conditions and physical constraints
      double precision, intent(in):: ll
! undo the log
      is_infinity = (ll <= log(epsilon(1d0)))  ! check approx zero
   end function

!! core sampler math

!> Accept or reject
!> First argument is new/proposed log(!)likelihood, 
!> second is old log likelihood
!> return logical
   logical function metropolis_choice(new_loglikelihood, old_loglikelihood)
      double precision, intent(in):: new_loglikelihood, old_loglikelihood
      double precision:: r  ! draw random number 0 to 1
      call random_number(r)
! TODO add optional pregen random
! l1/l2 > r  <=> logl1-logl2 > log(r)
      metropolis_choice = ((new_loglikelihood-old_loglikelihood) > log(r))
   end function

!!!! routines for random initialization

   subroutine init_pars_random(PI, pars0, fix_pars_flag, uniform_random_vector)
      use random_uniform, only: UNIF_VECTOR, next_random_uniform
      use samplers_math, only: log_nor2par
      implicit none
      type(PARINFO), intent(in):: PI  ! give number, bounds of params
      double precision, dimension(PI%npars), intent(inout):: pars0  ! return random initial values-nonnormalized
      logical, dimension(PI%npars), optional, intent(in):: fix_pars_flag  ! flags .true. to keep inidividual pars
      logical, dimension(PI%npars):: fix_pars_flag_  ! internal version
      type(UNIF_VECTOR), optional:: uniform_random_vector  ! object supplying pre-generated randoms 0 to 1  ! TODO make optional
      integer:: i

      if (.not. (present(fix_pars_flag))) then
         fix_pars_flag_ = .false.
      else
         fix_pars_flag_ = fix_pars_flag
      end if

      do i = 1, PI%npars
         ! parfix = 1 stay at prior value, parfix = 0 randomly search
         ! TODO condense this horrible logic to restart and fix_pars_flag somewhere outside of this fct
         !if (MCO%fixedpars .and. PI%parini(i) /= -9999d0) PI%parfix(i) = 1d0
         ! only assign random parameters if (a) randparini == .true. or (b) PI$parini(n) == -9999)
         !if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
!
         !TODO vectorize
         if (.not. (fix_pars_flag_(i))) then
            ! make sure to give it a random number vector unique to the chain
            ! scale each to the parameter's range
            pars0(i) = log_nor2par(uniform_random_vector%next_random_uniform(), PI%parmin(i), PI%parmax(i), PI%paradj(i))
         end if

      end do  ! for PI%npar loop
      !TODO test
      !TODO merge, optional uniform_random_vector arg
   end subroutine

   subroutine init_latin_square(PI, pars0, n_chains)  ! TODO
      use samplers_math, only: nor2par
      implicit none
      type(PARINFO), intent(in):: PI  ! give number, bounds, and potentially current value of params
      integer, intent(in):: n_chains
      double precision, dimension(PI%npars, n_chains), intent(out):: pars0  ! return initial values-nonnormalized
      double precision, dimension(PI%npars, n_chains):: points
      integer:: i, j
      ! generate N initial points on the (0..1)^N space in latin hypercube distribution...
      ! Must be run for all chains at once outside of parallel regions

      ! convert to real parameter values space
      do j = 1, n_chains
         do i = 1, PI%npars
            pars0(i, j) = nor2par(points(i, j), PI%parmin(i), PI%parmax(i))
         end do
      end do

   end subroutine

   pure logical function bounds_check(PI, PARS)
      type(PARINFO), intent(in):: PI
      double precision, dimension(:), intent(in):: PARS

      ! given real-space params ... convert to lognorm space, check if all in (0, 1)
      ! or check directly against real boundary values
      bounds_check = all((PARS > PI%parmin) .and. (PARS < PI%parmax))

   end function

   subroutine number_filenames(outfile, stepfile, covfile, covifile, chainid)
      character(len=*), intent(inout):: outfile, stepfile, covfile, covifile
      integer, intent(in):: chainid
      character(4):: chainid_str
      !! internal char version of chainid number, for filenames
         ! internal write to convert int -> str
         write (chainid_str, '(i0)') chainid
         ! append number to file names
         outfile = trim(outfile)//"_"//trim(chainid_str)
         stepfile = trim(stepfile)//"_"//trim(chainid_str)
         covfile = trim(covfile)//"_"//trim(chainid_str)
         covifile = trim(covifile)//"_"//trim(chainid_str)
   end subroutine

end module
