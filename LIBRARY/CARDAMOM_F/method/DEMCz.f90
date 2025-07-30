module DEMCz
   !-
   ! Differential evolution sampler (with history z), see ter Braak & Vrugt, Stat Comput 2008, 
   ! "Differential Evolution Markov Chain with snooker updater and fewer chains".
   !
   !  jklebes 2024-2025
   !
   !  See also R BayesianSamplers DEMCz ... we implement a similar algorithm so that behavior can be 
   ! compared.
   !
   !  This sampler should function similarly to other cardamom-samplers in terms of 
   !  invocation, printing and file writing behavior.
   !
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
   !-
   use samplers_shared, only: PARINFO, bounds_check, init_pars_random
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
   !> IN: model_loglikelihoodi: a subroutine parameters (array) -> loglikelihood (real), to maximize
   !> IN: model_loglikelihood_write: an alternative function of same type for writing to file
   !>                                (optional, defaults to using model_loglikelihood )
   !> IN: PI type(ParInfo) collection of parameter bounds
   !> IN: OPT type(DEMCz) collection of sampling options
   !> OUT: MCOUT_list type(DEMCzOUT) collection of results, one object for each chain
   !> IN: restart (optional, default .false.) : .true. -> restart from state in MCOU_list
   !>                          .false. -> start from random initial positions
   !> IN : nchains (optional, default 3) : number chains in swarm.
   !> Also writes history to file/output stream and progress to console.
   subroutine run_DEMCz(model_likelihood, PI, MCO, MCOUT_list, model_likelihood_write_in, restart, nchains)
      use samplers_shared, only: metropolis_choice
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
      double precision             :: l  
        !! likelihood of proposed values (private on each thread)
      double precision, dimension(PI%npars):: proposed_vector
        !! proposed values (private on each thread)
      double precision, dimension(PI%npars):: proposed_vector_real

      !> history Matrix Z, (npars x (nchains*maxiter))
      double precision, allocatable, dimension(:, :):: PARS_history

      double precision:: differential_weight
      type(UNIF_VECTOR), allocatable, dimension(:):: random_uniform_vectors

      integer:: npars, MAXITER
      integer:: P_target
      integer:: seed
      integer, dimension(:), allocatable:: ACC
      !! acceptance counter for each thread
      integer, dimension(:), allocatable:: ACCLOC
      !! acceptance counter for each thread, local to each nadapt phase 
      integer:: i, j, k, ITER, len_history  ! counters
      integer:: R1, R2 
      !!random indices in history

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
      allocate (ACC(mco%nchains))
      allocate (ACCLOC(mco%nchains))
      allocate (PARS_best(npars, mco%nchains))
      allocate (PARS_history(npars, MAXITER*mco%nchains))

      allocate (random_uniform_vectors(mco%nchains))
      ! where we are in filling in the history matrix so far
      len_history = 0

      !!! Initial state

!$OMP PARALLEL DO private(norpars)
      do j = 1, mco%nchains
      ! Initialize pregenerated random numbers, if using-local to this chain
      seed = irand()  ! TODO record later  ! TODO always the same ?
        call random_uniform_vectors(j)%initialize_random(seed)
         ! choose initial values
         ! TODO better function for initial state : latin square
         if (.not. restart_) then 
         call init_pars_random(PI, pars_current(:,j), PI%fix_pars, random_uniform_vectors(j))
      else
         pars_current(:,j) =  mcout_list(j)%pars
      end if
         ! also set the loglikelihoof of the state generated
         call model_likelihood(PARS_current(:, j), npars, l0(j), j)
         l_best(j) = l0(j)
         PARS_best(:,j) = PARS_current(:,j)

         write(*,*) "start", j, l0(j), PARS_current(:,j)
         ! potential burnin steps
         ! ... TODO

         ! write first values to history matrix
         PARS_history(:, j) = PARS_current(:, j)
      end do
!$OMP END PARALLEL DO

      len_history = len_history+mco%nchains

      !!! Main "time" loop

      do i = 2, MAXITER/mco%nadapt+2

         ! evolve each chain independently for nsteps (nsteps = K in ter Braak & Vrugt)
!$OMP PARALLEL DO private(ITER, R1, R2, l, proposed_vector, proposed_vector_real)
         do j = 1, mco%nchains
            do k = 1, mco%nadapt
                     ITER = i*mco%nadapt+k
!               call step_chain(PARS_current(:, j), l0(j), model_likelihood, &
!                               PI, PARS_history, len_history, differential_weight, j, random_uniform_vectors(j))

                     R1 = random_int(len_history)
                     R2 = random_int(len_history)
                     do while (R1 == R2)  ! should not equal R1 ("without replacement")
                        R2 = random_int(len_history)
                     end do
                     call step_real(proposed_vector, PARS_current(:,j), PARS_history(:, R1), PARS_history(:, R2), differential_weight, &
                        & random_uniform_vectors(j), PI)
                     if (bounds_check(PI, proposed_vector)) then
                     call model_likelihood(proposed_vector, PI%npars, l, j)
                     if (metropolis_choice(l, l0(j))) then
                        PARS_current(:,j) = proposed_vector
                        l0(j) = l
                        ACCLOC(j) = ACCLOC(j) + 1
                        ! check for new best
                        if (l > l_best(j)) then
                           l_best(j) = l
                           pars_best(:,j) = proposed_vector
                        endif
                     !else 
                        ! else-no accept 
                     end if
                     endif
                        

                    if (mod(ITER, MCO%nprint) == 0) then
                     write (*, *) "Chain ", j, "of", mco%nchains 
                     write (*, *) "Total proposal = ", ITER, " out of ", MAXITER
                     write (*, *) "Total accepted = ", ACC(j)
                     !write (*, *) "Overall acceptance rate    = ", dble(ACC)/dble(ITER)
                     !write (*, *) "Local   acceptance rate    = ", ACCRATE
                     write (*, *) "Current obs   = ", l0(j), "proposed = ", l, " log-likelihood"
                     write (*, *) "Maximum likelihood = ", l_best(j)
                    endif 
                  end do
                  ACC(j) = ACC(j) + ACCLOC(j)
                  ! write the chain's state after nadapt steps to Z
                  PARS_history(:, len_history+j) = PARS_current(:, j)
               end do
      !$OMP END PARALLEL DO !!Barrier implicit

         ! increment length M (filled so far) of Z
         len_history = len_history+mco%nchains

         ! check convergence

         ! Reorder for best chains ?  Then write to Z later.
         
      end do
      
      ! Final summary output from each chain individually
      !$OMP PARALLEL DO
      do j = 1, mco%nchains
      ! completed 
      write (*, *) "chain ",j, ": DEMCZ completed"
      write (*, *) "Overall acceptance rate (approx)  = ", dble(ACC(j))/dble(MAXITER)
      !write (*, *) "Final local acceptance rate = ", ACCRATE
      write (*, *) "Best log-likelihood = ", l_best(j)
      write (*, *) "Best parameters = ", pars_best(:,j)
      end do
      !$OMP END PARALLEL DO

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


   !> Thin wrapper on "step" which takes three vectors in the space of raw parameter values; 
   !> they are converted to lognormalized parameters (0, 1) before being lineraly combined in 
   !> the core "step" function.  After cardamom-MHMCMC and 
   !> sampling is observed to be better when this is done on lognormalized parameters.
   subroutine step_real(vout, v1, v2, v3, differential_weight, random_uniform_vector, PI)
      use samplers_math, only: log_nor2par_scalar, log_par2nor_scalar
      type(PARINFO), intent(in):: PI
      double precision, dimension(PI%npars), intent(out):: vout
      !! proposed step on raw parameter space
      double precision, dimension(PI%npars):: vout_lognorm
      !! proposed step on lognormed parameter space
      double precision, dimension(PI%npars), intent(in):: v1, v2, v3
      !! input: current and two random states from history on parameter space
      double precision, dimension(PI%npars):: v1_lognorm, v2_lognorm, v3_lognorm
      !! input: current and two random states from history on lognorm space
      type(UNIF_VECTOR), intent(inout):: random_uniform_vector
      double precision:: differential_weight
      v1_lognorm = log_par2nor_scalar(v1, PI%parmin, PI%parmax, PI%paradj) 
      v2_lognorm = log_par2nor_scalar(v2, PI%parmin, PI%parmax, PI%paradj) 
      v3_lognorm = log_par2nor_scalar(v3, PI%parmin, PI%parmax, PI%paradj) 
      call step(vout_lognorm,  log_par2nor_scalar(v1, PI%parmin, PI%parmax, PI%paradj), &
                 log_par2nor_scalar(v2, PI%parmin, PI%parmax, PI%paradj), &
                log_par2nor_scalar(v3, PI%parmin, PI%parmax, PI%paradj), &
               differential_weight, random_uniform_vector, PI%npars)
      vout = log_nor2par_scalar(vout_lognorm, PI%parmin, PI%parmax, PI%paradj)
   end subroutine

   !> Generate new proposed state from currect state and history
   !> ter Braak & Vrugt eq 2
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
      vout = v1+differential_weight*(v2-v3)  + .000001*rn
   end subroutine

   integer function random_int(N)
      integer, intent(in):: N
      double precision:: r
      call random_number(r)
      random_int = floor(N*r) + 1
   end function

end module DEMCz
