module DEMCz

   !-
   ! On normalized parameters:
   ! All calculation, statistics, and writing is done on "raw" values, 
   ! NOT parameters normalized to (0, 1).
   !
   ! How to use:
   !  Create an object of type PARINFO : number and bounds of parameters
   !  and DEMCzOPT-options, containing as many or few of the fields as needed, the rest default to the
   !                 default values in type definiton here
   !  Create an object of type MCMC_OUTPUT to write reults to.
   !  Create a double precision function loglikelihood taking a vector of npars (same as in PARINFO) parameters and
   !                 returning rel:: loglikelihood.
   !  (not implemented yet) Optionally set OMP_NUM_THREADS
   !  Call subroutine run_DEMCz(fct, parinfo, demczopt, mcmcout)
   !!!!
   use samplers_shared, only: PARINFO, bounds_check
   use samplers_math, only: log_nor2par
   use random_uniform, only: UNIF_VECTOR
   use cardamom_MHMCMC, only: MCMC_OUTPUT, MCMC_options  
   use samplers_io, only: io_buffer_space, initialize_buffers, open_output_files
   use OMP_LIB

   implicit none

   public

   !> A collection of input options to the DEMCz sampler run
   !> contains default values
   type, extends(MCMC_options):: DEMCZOPT
! DEMcz algorithm parameters
      double precision:: differential_weight = 0.8  
!! differential weight gamma, [0, 2]
      double precision:: crossover_probability = 0.9  
!! crossover probability CR, [0, 1]
   end type DEMCzOPT

contains

   !> Main DEMCz sampler subroutine
   !> Write all OMP parallelization at this top level only
   !> IN: model_loglikelihood a function parameters (array) -> loglikelihood (real)
   !> IN: model_loglikelihood_write: an alternative function of same type for writing to file
   !>                                (optional, defaults to using model_loglikelihood )
   !> IN: PI type(ParInfo) collection of parameter bounds
   !> IN: OPT type(DEMCz) collection of sampling options
   !> OUT: MCOUT type(DEMCzOUT) collection of results
   !> Also writes history to file/output stream and progress to console.
   subroutine run_DEMCz(model_likelihood, PI, MCO, MCOUT_list, model_likelihood_write_in, restart, nchains)
      implicit none

      !! input and output structs
      type(PARINFO), intent(in):: PI
      !!PARINFO struct from model giving number, bounds of parameters
      type(DEMCZOPT), intent(inout):: MCO
      !! struct of options for the run, shared between all threads
      type(MCMC_OUTPUT), dimension(:), allocatable, intent(inout):: MCOUT_list  ! array of MCOUT objects
      !! Array of MCMC_OUTPUT structs for each thread's results

      logical, optional, intent(in):: restart
      !! is it a restart ? (i.e. start from data in MCOUT instead of initializing new), optional, default .false.
      logical:: restart_
      !! internal restart flag, equal to optional input flag 'restart' or .false.
      integer, optional:: nchains
      !! number chains optional, default 1

      !> Matrix X, (npars x nchains), holding current state of the n chains
      double precision, allocatable, dimension(:, :):: PARS_current
      double precision, dimension(PI%npars ):: norPARS
      ! and their current loglikelihood values
      double precision, allocatable, dimension(:):: l0
      ! and their best likelihood values and best pars so far
      double precision, allocatable, dimension(:):: l_best
      double precision, allocatable, dimension(:, :):: PARS_best

      !> history Matrix Z, (npars x (nchains*maxiter))
      double precision, allocatable, dimension(:, :):: PARS_history

      double precision:: differential_weight
      type(UNIF_VECTOR), allocatable, dimension(:):: random_uniform_vectors

      integer:: npars, MAXITER
      integer:: P_target
      integer:: seed
      integer:: i, j, k, len_history  ! counters

      !> the function to maximize.
      !> Completely agnostic, samples any functions vector -> double
      !> Typically a loglikelihood evaluation of a model against observation data
      !> given the inputted parameter values.
      interface
         subroutine model_likelihood(param_vector, n, ML, id)
            implicit none
            double precision, dimension(n), intent(inout):: param_vector  ! intent(in), inout for compatibility with R via C
            integer, intent(in):: n
            double precision, intent(out):: ML
            integer, intent(in), optional:: id
         end subroutine model_likelihood
      end interface

      !> optionally  give a second function with same shape as model_likelihood, 
      !> for writing to file.  model_likelihood_write_in is the input arg, which may not be present.
      procedure(model_likelihood), optional:: model_likelihood_write_in
      !> A second function with same shape as model_likelihood, 
      !> for writing to file.  Internal variable equal to model_likelihood_write_in if present
      !> or (default) same as main model_likelihood function
      procedure(model_likelihood), pointer:: model_likelihood_write

      ! Argument processing  !!!!!!!!!!!!!!!

      if (.not. present(restart)) then
         restart_ = .false.
      else
         restart_ = restart
      end if

      if (.not. present(nchains)) then
         MCO%nchains = 3  ! at least 3 are mandatory for this method to work
      else
      MCO%nchains = nchains
      end if

      ! process function arguments
      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default, if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if


      ! Extract from types
      differential_weight = MCO%differential_weight

      npars = PI%npars
      MAXITER = MCO%nout
      P_target = MCO%P_target

      ! prepare outputs
      if (.not.allocated(MCOUT_list)) allocate(MCOUT_list(mco%nchains))

      ! Allocate arrays

      allocate (PARS_current(npars, mco%Nchains))
      ! but we need these to persist between parallel regions
      allocate (l0(mco%nchains))
      allocate (l_best(mco%nchains))
      allocate (PARS_best(npars, mco%nchains))
      allocate (PARS_history(npars, MAXITER*mco%nchains))

      allocate (random_uniform_vectors(mco%nchains))
      ! where we are in filling in the history matrix so far
      len_history = 0

      !!! Initial state

!$OMP PARALLEL DO private(norpars)
      do j = 1, mco%nchains
         ! choose initial values
         ! TODO better function for initial state : latin square
         if (.not. restart_) then 
         call init_random(npars, norpars)
         pars_current(:,j) = log_nor2par(npars, norpars, PI%parmin, PI%parmax, pi%paradj)
         write(*,*) "randomized", pars_current(:,j)

      else
         pars_current(:,j) =mcout_list(j)%pars
         write(*,*) "did not randomize", pars_current(:,j)
      end if
         ! also set the loglikelihoof of the state generated
         call model_likelihood(PARS_current(:, j), npars, l0(j), j)

         ! potential burnin steps
         ! ... TODO

         ! write first values to history matrix
         PARS_history(:, j) = PARS_current(:, j)
      ! Initialize pregenerated random numbers, if using-local to this chain
      seed = irand()  ! TODO record later  ! TODO always the same ?
      write(*,*) "seed", seed
        call random_uniform_vectors(j)%initialize_random(seed)
      end do
!$OMP END PARALLEL DO

      len_history = len_history+mco%nchains

      !!! Main "time" loop

      do i = 2, MAXITER/mco%nadapt+2

         ! evolve each chain independently for nsteps (nsteps = K in ter Braak & Vrugt)
!$OMP PARALLEL DO 
         do j = 1, mco%nchains
            do k = 1, mco%nadapt
               call step_chain(PARS_current(:, j), l0(j), model_likelihood, &
                               PI, PARS_history, len_history, differential_weight, j, random_uniform_vectors(j))
         
              if (mod(i*mco%nadapt+k, MCO%nprint) == 0) then
                write(*,*) "thread", j
                write(*,*) "loglikelihood", l0(j)
                write(*,*) "pars", pars_current(:,j)
              endif 
            end do
            ! write the chain's state after nadapt steps to Z
            PARS_history(:, len_history+j) = PARS_current(:, j)
         end do
!$OMP END PARALLEL DO !!Barrier implicit

         ! increment length M (filled so far) of Z
         len_history = len_history+mco%nchains

         ! check convergence

         ! Reorder for best chains ?  Then write to Z later.
         
      end do

      !!! Write out
      write (*, *) l0

   end subroutine

   !> Initialize the chain's state with random values from
   !> parameter ranges.
   !> Works with normalized values : returns a number between 0 and 1
   !> For each parameter
   subroutine init_random(npars, norpars)
      integer, intent(in):: npars
      double precision, dimension(:), intent(out):: norpars
      integer:: i
      do i = 1, npars
         call random_number(norpars(i))
      end do
      ! also output the loglikelihood of the state generated

   end subroutine

   !> Evolve the state of one chain for 1 step
   subroutine step_chain(X_i, l0, model_likelihood, &
                         PI, PARS_history, len_history, differential_weight, &
                         thread_id, random_uniform_vector)
      use samplers_shared, only: metropolis_choice
      type(PARINFO), intent(in):: PI
      integer, intent(in):: thread_id  ! thread_id-to pass to model evaluation in cases where it matters
      double precision, dimension(:), intent(inout):: X_i  ! current state of the chain; normalized values of all pars
      double precision, dimension(PI%npars):: previous_vector, proposed_vector  ! internal: save previous state, proposed new state
      double precision, dimension(:, :), intent(in):: PARS_history  ! the matrix Z so far, to read 2 rows from
      integer, intent(in):: len_history  ! length to which Z is filled
      double precision, intent(inout):: l0  ! likelihood of previous accepted params
      type(UNIF_VECTOR), intent(inout):: random_uniform_vector
      double precision             :: l  ! likelihood of proposed values

      integer:: R1, R2  ! indices of 2 random rows
      double precision:: rand
      double precision:: differential_weight
      integer:: i  ! counters

      ! the function, vector of normalized par values -> loglikelihood

      interface
         subroutine model_likelihood(param_vector, n, ML, id)
            implicit none
            double precision, dimension(n), intent(inout):: param_vector
            integer, intent(in):: n
            double precision, intent(out):: ML
            integer, intent(in), optional:: id
         end subroutine model_likelihood
      end interface

      R1 = random_int(len_history)
      R2 = random_int(len_history)
      do while (R1 == R2)  ! should not equal R1 ("without replacement")
         R2 = random_int(len_history)
      end do
      previous_vector = X_i
      call step(proposed_vector, previous_vector, PARS_history(:, R1), PARS_history(:, R2), differential_weight, &
         & random_uniform_vector, PI%npars)
      if (bounds_check(PI, proposed_vector)) then
      call model_likelihood(proposed_vector, PI%npars, l, thread_id)
      !write(*,*) "loglikelihood proposed", l
      !write(*,*) "loglikelihood previous", l0
      !write(*,*) "pars proposed", proposed_vector
      !write(*,*) "pars previous", previous_vector
      !write(*,*) "accept", metropolis_choice(l, l0)
      if (metropolis_choice(l, l0)) then
         X_i = proposed_vector
         l0 = l
      end if
      endif
   end subroutine

   ! Generate new proposed state from currect state and history
   ! ter Braak & Vrugt eq 2
   subroutine step(vout, v1, v2, v3, differential_weight, random_uniform_vector, npars)
      use samplers_math, only: random_normal
      integer, intent(in):: npars
      double precision, dimension(:), intent(out):: vout
      double precision, dimension(:), intent(in):: v1, v2, v3
      type(UNIF_VECTOR), intent(inout):: random_uniform_vector
      double precision:: differential_weight
      double precision:: rn(npars)
      integer:: p
      ! get differential_weight, corssover_probability from module data
      do p = 1, npars
         call random_normal(random_uniform_vector, rn(p))
      end do
      vout = v1+differential_weight*(v2-v3)  + .0001*rn*
   end subroutine

   integer function random_int(N)
      integer, intent(in):: N
      double precision:: r
      call random_number(r)
      random_int = floor(N*r) + 1
   end function

end module DEMCz
