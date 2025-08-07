module cardamom_main_utils
   implicit none
   public

contains

   subroutine initialize_stats(MCOUT, npars)
      use cardamom_MHMCMC, only: MCMC_OUTPUT
      integer, intent(in):: npars
      type(MCMC_OUTPUT), intent(inout):: MCOUT
      allocate (MCOUT%covariance(npars, npars), MCOUT%parvar(npars), MCOUT%meanpar(npars))
      call reset_stats(MCOUT, npars)
   end subroutine
   subroutine reset_stats(MCOUT, npars)
      use cardamom_MHMCMC, only: MCMC_OUTPUT
      type(MCMC_OUTPUT), intent(inout):: MCOUT
      integer, intent(in):: npars
      integer:: n
      MCOUT%parvar = 1d0; MCOUT%Nparvar = 0d0
      ! Covariance matrix cannot be set to zero therefore set initial
      ! value to a small positive value along to variance access
      MCOUT%covariance = 0d0; MCOUT%meanpar = 0d0; MCOUT%cov = .false.
      MCOUT%use_multivariate = .false.
      do n = 1, npars
         MCOUT%covariance(n, n) = 1d0
      end do
   end subroutine

   !
   !------------------------------------------------------------------
   !
   subroutine find_edc_initial_values(MCO, MCOUT_list, nchains)
    !! subroutine deals with the determination of initial parameter and initial
    !! conditions which are consistent with EDCs
    !! pre-loop, Run MCMC sampler with modified likelihood fct
      use model_shared, only: PI
      use cardamom_io, only: restart_flag
      use cardamom_MHMCMC, only: MCMC_OUTPUT, run_mcmc, run_parallel_mcmc, MCMC_options
      !use model_likelihood_module, only: model_likelihood, &
      !sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood  ! to replace soon with wrappers
      use model_likelihood_wrapper  ! TODO next refactoring step
      use cardamom_structures, only: DATAin

      implicit none

      ! declare local variables
      integer, intent(in):: nchains
      type(MCMC_OUTPUT), dimension(:), allocatable, intent(inout):: MCOUT_list
      type(MCMC_OUTPUT), dimension(:), allocatable:: MCOUT_list_tmp
      type(mcmc_OPTIONS), intent(out):: MCO
      integer:: i, counter_local(nchains), nOUT_save, nWRITE_save, nADAPT_save
      integer:: success_count
      logical:: append_save
      logical:: restart(nchains)
      double precision:: ll
      double precision:: PEDC(nchains), PEDC_prev(nchains), P_target
      double precision, dimension(PI%npars):: parini  ! local variable, or array

      allocate (MCOUT_list_tmp(nchains))

      ! Hold for later
      nOUT_save = MCO%nOUT; nWRITE_save = MCO%nWRITE; nADAPT_save = MCO%nADAPT
      append_save = MCO%append

      ! set MCMC options needed for EDC run
      ! same for all chains
      MCO%APPEND = .false.
      MCO%nADAPT = 500
      MCO%fADAPT = 1d0
      MCO%nOUT = 100000
      MCO%nPRINT = 0
      MCO%nWRITE = 0
      ! the next two lines ensure that parameter inputs are either given or
      ! entered as-9999
      MCO%randparini = .true.
      MCO%returnpars = .true.
      MCO%fixedpars = .true. ! TLS: changed from .false. for testing 16/12/2019

      ! Set initial priors to vector...
      ! TODO array for multichain
      parini = DATAin%parpriors(1:PI%npars)
      ! ... and assume we need to find random parameters
      ! Target likelihood allows for controlling when the MCMC will stop
      P_target = 0d0

      ! if the prior is not missing and we have not told the edc to be random
      ! keep the value
!    do n = 1, PI%npars
!       if (PI%parini(n) /= -9999d0 .and. DATAin%edc_random_search < 1) PI%parfix(n) = 1
!    end do  ! parameter loop

      do i = 1, nchains
         MCOUT_list_tmp(i)%PARS = parini
         call initialize_stats(MCOUT_list_tmp(i), PI%npars)
      end do

      ! TODO sort out omp loop over this
      ! if this is not a restart run, i.e. we do not already have a starting
      ! position we must being the EDC search procedure to find an ecologically
      ! consistent initial parameter set
      restart = .false.
      if (.not. restart_flag) then  ! TODO outside-only run this fct if not restart

         ! set up edc log likelihood for MHMCMC initial run
         PEDC_prev = -1000d0; PEDC = -1d0; counter_local = 0
         success_count = 0

            write (*, *) "Beginning EDC search attempt "
            write (*, *) nchains, "chains working ... "
            MCO%randparini = .false.
            ! TODO limit number to number of available hardware threads
            !$omp parallel do private(ll) 
            do i = 1, nchains
             do while (success_count < nchains)!(PEDC < 0d0)
               ! call the MHMCMC directing to the appropriate likelihood function
               MCO%restart = restart(i)
            call run_mcmc(edc_model_likelihood_fct, PI, MCO, MCOUT_list_tmp(i), model_likelihood_fct, chainid = i)
               restart(i) = .true.

               ! turn off random selection for initial values
               write (*, *) "...intermediate EDC search progress check"

               ! store the best parameters from that loop
               PEDC(i) = MCOUT_list_tmp(i)%bestll
               write (*, *) "Found best loglikelihood", PEDC(i)
               ! if any chains's MCOUT object is success (reached loglikelihood = 0), 
               ! copy it to MCOUT_list
               !$omp critical
               if (MCOUT_list_tmp(i)%ll >= (0d0-10*epsilon(1.0d0))) then
                  success_count = success_count+1
                  write (*, *) "Found ", success_count, "EDC-compatible starting points (", i, ")"
                  if (success_count <= nchains) then
                     MCOUT_list(success_count) = MCOUT_list_tmp(i)
                  end if
                  call reset_stats(MCOUT_list_tmp(i), PI%npars)
                  MCO%randparini = .true. ! TODO problem
                  PEDC_prev(i) = -1000d0
                  MCOUT_list_tmp(i)%PARS = DATAin%parpriors(1:PI%npars)
                  restart(i) = .false.
                  counter_local(i) = 0

               end if
               !$omp end critical

               ! keep track of attempts
               counter_local(i) = counter_local(i) + 1
               ! periodically reset the initial conditions
               if (PEDC(i) < 0d0 .and. PEDC(i) <= PEDC_prev(i) .and. counter_local(i) > 5) then
                  ! Reset the previous EDC likelihood score
                  PEDC_prev(i) = -1000d0
                  ! Reset parameters back to default
                  MCOUT_list_tmp(i)%PARS = DATAin%parpriors(1:PI%npars)
                  ! reset to select random starting point
                  MCO%randparini = .true. ! TODO problem
                  ! reset the parameter step size at the beginning of each attempt
                  call reset_stats(MCOUT_list_tmp(i), PI%npars)
                  counter_local(i) = 0
                  restart(i) = .false.
                  write (*, *) "resetting to initial"
               else
                  PEDC_prev(i) = PEDC(i)
               end if
               end do  ! for while condition
            end do
            !$omp end parallel do


      end if  ! if for restart

      ! reset so that currently saved parameters will be used
      ! starting point in main MCMC
      ! PI%parfix(1:PI%npars) = 0  ! TODO
      !MCOUT%bestpars = 0d0

      end subroutine
end module
