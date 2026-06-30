module cardamom_main_stresstest

contains
   subroutine run_stresstest()
      use cardamom_io, only: initialize, &
                             read_options, open_output_files, &
                             check_for_existing_output_files, restart_flag, &
                             update_for_restart_simulation, write_covariance_matrix, &
                             close_output_files, write_covariance_info
      use MHMCMC_module, only: MHMCMC, par_minstepsize, par_initstepsize, N_before_mv
      use MHMCMC_StressTests, only: StressTest_likelihood, StressTest_sublikelihood, prepare_for_stress_test
      use model_likelihood_wrapper  ! TODO next step

 !!!!!!!!!!!
      ! Authorship contributions
      !
      ! This code is based on the original C verion of the University of Edinburgh
      ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
      ! All code translation into Fortran, integration into the University of
      ! Edinburgh CARDAMOM code and subsequent modifications by:
      ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
      ! J. F. Exbrayat (University of Edinburgh)
      ! D. T. Milodowski (d.t.milodowski@ed.ac.uk, University of Edinburgh)
      ! See function/subroutine specific comments for exceptions and contributors
 !!!!!!!!!!!

      ! Created: Anthony A. Bloom
      ! Major modification history:
      ! Version 1.0: C language MHMCMC and io created by Anthony A. Bloom
      ! Version 1.1: Translated into Fortran, including some generic function replacement by T. L. Smallman
      ! Version 1.2: Analysis restart capacity added by T.L. Smallman
      ! Version 1.3: MHMCMC updated to APMCMC by T.L. Smallman with advice from Anthony A. Bloom
      !            : Following Haario et al., (2001, 2006) and Roberts & Rosenthal (2009).
      ! Version 1.4: Pre-APMCMC phase using normalised likelihoods added by T. L. Smallman
      !            : Pre-APMCMC allows for rapidly moving towards observations from very bad starting points.
      ! Version 1.5: Added options for varied scaling approaches for the cost function (e.g. division by sample size (i.e. n), sqrt(n), 1+log(n))
      ! Specific citations for developments included in the code.

      ! This is the main subroutine for the CARDAMOM framework. The specific model
      ! method combinations are achieved through case specific compilation of the
      ! case while maintaining strict consistent io formats to allow for these
      ! combinations

      ! Command line inputs are:
      ! 1) file in
      ! 2) file out
      ! 3) integer number of solutions requested
      ! 4) print-to-screen frequency
      ! 5) write-to-file frequency
      ! 6) 0/1 flag to use normalised log-likelihood pre-mcmc
      ! 7) Flag to select the cost function normalisation approach

      implicit none

      ! declare local variables
      character(350):: infile, outfile, solution_wanted_char, freq_print_char, &
                       freq_write_char, do_inflate_char, cost_func_scaling_char
      integer:: solution_wanted, freq_print, freq_write, time1, time2, time3, n, &
                nOUT_save, do_inflate_dble, cost_func_scaling_dble
      logical:: do_inflate = .false.

      ! user update
      write (*, *) "Beginning read of the command line"

      ! read user options from the command line
      call get_command_argument(1, infile)
      call get_command_argument(2, outfile)
      call get_command_argument(3, solution_wanted_char)
      call get_command_argument(4, freq_print_char)
      call get_command_argument(5, freq_write_char)
      call get_command_argument(6, do_inflate_char)
      call get_command_argument(7, cost_func_scaling_char)

      ! now convert relevant ones to integeter
      ! Note: that I10 is the maximum allowed with the default integer kind (kind = 4).
      !       allow for larger number of iterations integer (kind  = 8) is needed throughout the code.
      read (solution_wanted_char, '(I10)') solution_wanted
      read (freq_print_char, '(I10)') freq_print
      read (freq_write_char, '(I10)') freq_write
      read (do_inflate_char, '(I10)') do_inflate_dble
      read (cost_func_scaling_char, '(I10)') cost_func_scaling_dble

      ! Assign inflate logical condition
      if (do_inflate_dble == 1) do_inflate = .true.
      ! Sanity check of the cost_function_scaling_char
      if (cost_func_scaling_dble > -1 .and. cost_func_scaling_dble < 4 .and. freq_write > 0) then
         ! All is well
      else
         ! All is not well-complain
         print *, "ERROR: Command line argument to specify the cost function or write to file frequency is incorrect."
         print *, "Command line should have 7 arguments (in addition to the cardamom.exe)."
         print *, "These are: "
         print *, "1) input file path."
         print *, "2) output file path-note that PARS, STEP, COV, COVINFO will be appended to this name path outfile."
         print *, "3) No. of parameter proposals to make."
         print *, "4) Iteration freq. for printing to screen (main MCMC phase only)."
         print *, "5) Iteration freq. for writing results to files."
         print *, "6) do sample size normalisation phase - "
         print *, "  0 = FALSE"
         print *, "  1 = TRUE"
         print *, "7) Select likelihood cost function - "
         print *, "  0 = no scaling"
         print *, "  1 = scaling by sample size (n)"
         print *, "  2 = scaling by sqrt(n)"
         print *, "  3 = scaling by log(n)"
         stop
      end if

      ! user update
      write (*, *) "Command line options read, moving on now"

      ! seed the random number generator
      ! determine unique (sort of) seed value; based on system time
      call system_clock(time1, time2, time3)
      ! set seed value outside of the function, idum must be a negative number
      idum = dble(time1 + time2 + time3)
      call rnstrt(nint(idum))

      ! Determine whether ot not we are doing a real analysis or running a stress trest
      if (trim(infile) == "StressTest") then
         ! call special functions to prepare for stress test
         call prepare_for_stress_test(infile, outfile)
      else
         write (*, *) "ERROR expected an infile called StressTest"
         STOP 1
      end if

      ! load module variables needed for restart check
      ! NOTE: THIS MUST HAPPEN BEFORE CHECKING FOR RESTART
      call read_options(solution_wanted, freq_print, freq_write, outfile)
      ! check whether this is a restart?
      call check_for_existing_output_files(PI%npars, MCO%nOUT, MCO%nWRITE, MCO%sub_fraction, &
                                           MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile)
      ! Initialise MCMC output, possibly a bit of a redundent subroutine...
      call initialise_mcmc_output  ! TODO now happens in run_mcmc if not restart
      ! Open the relevant output files TODO now happens in run_mcmc
      call open_output_files(MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile)

      ! Initialise counters used to track the output of parameter sets
      !TODO now each chain has its own
      io_space%io_buffer_count = 0
      io_space%io_buffer = min(1000, max(10, (MCO%nOUT/MCO%nWRITE)/10))

      ! Allocate variables used in io buffering,
      ! these could probably be moved to a more sensible place within cardamom_io.f90 DONE
      allocate (io_space%variance_buffer(PI%npars, io_space%io_buffer), &
                io_space%meanpars_buffer(PI%npars, io_space%io_buffer), &
                io_space%pars_buffer(PI%npars, io_space%io_buffer), &
                io_space%prob_buffer(io_space%io_buffer), &
                io_space%nsample_buffer(io_space%io_buffer), &
                io_space%accept_rate_buffer(io_space%io_buffer))

      ! Report which model ID we are using
      write (*, *) "Running model version ", DATAin%ID

      ! Check whether we are doing a stress test again
      if (DATAin%ID < 0) then
         ! We are doing a stress test
         write (*, *) "Carrying out a stress test analysis"
         write (*, *) "Any existing files will be ignored"

         ! Reset interations counter
         MCOUT%nos_iterations = 0
         ! Ensure that we use a random starting point
         MCO%randparini = .true.
         MCO%returnpars = .true.
         MCO%fixedpars = .false.
         restart_flag = .false.

         ! Do we do the initial MCMC period where we normalise the likelihood by
         ! number of observations
         ! This process allows for very bad starting points to more easily move
         ! towards the general area of the observatons.
         if (MCOUT%nos_iterations < (MCO%nOUT*MCO%sub_fraction) .and. do_inflate) then

            ! Having found an EDC compliant parameter vector, we want to do a MCMC
            ! search on inflated uncertainties. This inflation search allows us to
            ! more easily move towards higher likelihoods, where the inflation
            ! allows easier movement through parameter space.

            ! The inflated search phase will make use of three stages over which
            ! the
            ! inflation will be reduced. Phase 1 used half of the allocated
            ! iterations for the inflation while the second and third is

            ! Set flag to indicate this phase has occurred and make a record of the
            ! total iterations to be attempted
            MCO%sub_sample_complete = .true.; nOUT_save = MCO%nOUT

            ! Report to the user
            write (*, *) "Beginning parameter search on sample size normalised likelihoods"

            MCO%nOUT = nint(dble(nOUT_save)*MCO%sub_fraction) - MCOUT%nos_iterations
            write (*, *) "Nos iterations to be proposed = ", MCO%nOUT
            MCO%fADAPT = 1d0 !; MCO%nADAPT = 1000
            call MHMCMC(1d0, StressTest_likelihood, StressTest_sublikelihood)
            ! Use the best parameter set as the starting point for the next stage
! REALLY NOT SURE I SHOULD BE DOING THIS-SHOULD BE PROGRESSING FROM THE LAST ACCEPTED PARAMETER SET?
            PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
            ! Leave parameter and covariance structures as they come out form the
            ! sub-sample-but reset the number of samples used in the update
            ! weighting
            if (PI%cov .and. PI%use_multivariate) then
               PI%Nparvar = (N_before_mv*dble(PI%npars)) + 1d0
            else
               ! reset the parameter step size at the beginning of each attempt
               PI%parvar = 1d0; PI%Nparvar = 0d0
               ! Covariance matrix cannot be set to zero therefore set initial
               ! value to a small positive value along to variance access
               PI%covariance = 0d0; PI%mean_par = 0d0; PI%cov = .false.
               PI%use_multivariate = .false.
               do n = 1, PI%npars
                  PI%covariance(n, n) = 1d0
               end do
            end if  ! do we need a new covariance matrix or can we use the existing one?

            ! Assume that sub-sampling process, if completed, will use 10 % of the
            ! simulation time therefore we want to adjust the output frequency to
            ! correct for this
            MCO%nOUT = max(1, nOUT_save - MCOUT%nos_iterations)
            ! by pass read_options file for StressTest special case
            MCO%append = 1
            MCO%nADAPT = 1000
            MCO%fADAPT = 0.5d0
            MCO%randparini = .false.
            MCO%returnpars = .false.
            MCO%fixedpars = .true.

         end if  ! restart flag

         ! Let the user know how many more we will propose
         write (*, *) "Nos iterations to be proposed = ", MCO%nOUT
         ! Call the AP-MCMC
         call MHMCMC(1d0, StressTest_likelihood, StressTest_likelihood)
         ! Tell the user the best parameter set
         print *, "Best parameters = ", MCOUT%best_pars

      else  ! We are not doing a stress test
         write (*, *) "ERROR missing data for stress test"
         STOP 1
      end if  ! stress test or not

      ! tidy up by closing all files
      call close_output_files

      ! Final message to the user
      write (*, *) "==========================================================="
      write (*, *) "==== CARDAMOM analysis for the current chain completed ===="
      write (*, *) "==========================================================="
      write (*, *) "=========================Honestly=========================="

   end subroutine run_stresstest
end module cardamom_main_stresstest
