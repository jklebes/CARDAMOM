!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to
! assimilate observations and ecological theory to retrieve parameters for the
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently.
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk),
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see < https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This is the main subroutine for the CARDAMOM framework. The specific model
! method combinations are achieved through case specific compilation of the
! case while maintaining strict consistent io formats to allow for these
! combinations
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! D. T. Milodowski (d.t.milodowski@ed.ac.uk, University of Edinburgh)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cardamom_framework
   use cardamom_MHMCMC, only: MCMC_OUTPUT, MCMC_OPTIONS, run_mcmc, run_parallel_mcmc
   use model_shared, only: PI, initialize_carbon_model
   use samplers_shared, only: init_infinity ! should be in Utils ?
   use cardamom_structures, only: DATAin
   use cardamom_io, only: initialize, &
                          read_options, &
                          update_obs_scaling_normal, update_obs_scaling_nsamples, &
                          update_obs_scaling_sqrt_nsamples, update_obs_scaling_log_nsamples
   use samplers_io, only: open_output_files, &
                          check_for_existing_output_files, &
                          update_for_restart_simulation, &
                          write_covariance_matrix, &
                          close_output_files, write_covariance_info
   !use MHMCMC_StressTests, only: StressTest_likelihood_fct, StressTest_sublikelihood_fct, prepare_for_stress_test
   use model_likelihood_wrapper, only: model_likelihood_fct, edc_model_likelihood_fct, scaled_model_likelihood_fct
   use cardamom_main_utils

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

   implicit none(type, external)

   ! declare local variables
   character(350):: infile, outfile, solution_wanted_char, freq_print_char, &
                    freq_write_char, do_inflate_char, cost_func_scaling_char
   integer:: solution_wanted, freq_print, freq_write, time1, time2, time3, &
             nOUT_save, do_inflate_dble, cost_func_scaling_dble
   logical:: do_inflate = .false.
   logical:: sub_sample_complete = .false.
   double precision:: sub_fraction = 0.2d0
    !! run this percentage of simulation with variant function
   type(MCMC_OUTPUT), dimension(:), allocatable:: MCOUT_list
    !! array of output objects from each thread
   type(MCMC_OPTIONS):: MCO
     !! options for sampler
   logical:: restart

   ! TODO not to hardcode, from command line argument
   integer:: nchains = 3
   integer:: i

   call init_infinity()

   allocate (MCOUT_list(nchains))

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

   ! TODO get n seeds, here or smoewhere else
   ! TODO make sure seeds are saved
   ! seed the random number generator
   ! determine unique (sort of) seed value; based on system time
   ! call system_clock(time1, time2, time3)
   ! set seed value outside of the function, idum must be a negative number
   !idum = dble(time1+time2+time3)
   !call rnstrt(nint(idum))

   ! Determine whether ot not we are doing a real analysis or running a stress trest
   if (trim(infile) == "StressTest") then
      !call run_stresstest()
      ! call prepare_for_stress_test(infile, outfile)  ! sets cardamom_structures:: DATAin
      stop
   end if

   ! read input data file
   call initialize(infile) ! = initialize_parinfo, read_check_binary_data, initialize_model
   ! sets cardamom_structures:: DATAin

   call initialize_carbon_model(nchains)
   do i = 1, nchains
      call initialize_stats(MCOUT_list(i), PI%npars)
   end do

   ! having filled PI%npars from model file, we can allocate stats array in MCOUT

   ! load module variables needed for restart check
   ! NOTE: THIS MUST HAPPEN BEFORE CHECKING FOR RESTART
   call read_options(solution_wanted, freq_print, freq_write, outfile, MCO)

   ! check whether this is a restart from aborted simulation
   MCO%restart = .true. !to gather results onto-false as soon as some file is not found
   do i = 1, nchains
      call check_for_existing_output_files(PI%npars, MCO, sub_fraction, i, restart)
      MCO%restart = MCO%restart .and. restart
      !if all nchains files were found, read them to get a starting point
      if (MCO%restart) then
         call update_for_restart_simulation(MCO, MCOUT_list(i), PI%npars)
      end if
   end do

   ! Report which model ID we are using
   write (*, *) "Running model version ", DATAin%ID

   if (.not. MCO%restart) then
      ! Begin search for initial conditions
      write (*, *) "Beginning search for initial parameter conditions"
      ! Determine initial values, this requires using the AP-MCMC
      call find_edc_initial_values(MCO, MCOUT_list, nchains)
      ! Having done EDC search phase, flag to start the next phase from this state
      MCO%fixedpars = .true.
      do i = 1, nchains
         MCOUT_list(i)%nos_iterations = 0
      end do
   end if

   do i = 1, nchains

      ! Reset the MCMC parameters for the next stage
      call read_options(solution_wanted, freq_print, freq_write, outfile, MCO)

      ! Reset stepsize and covariance for main DRAM-MCMC
      call reset_stats(MCOUT_list(i), PI%npars)

   end do

   ! sub-sampling phase, first sub_fraction% of the simulation with variant loglikelihood
   if (DATAin%total_obs > 0 .and. MCOUT_list(1)%nos_iterations < (MCO%nOUT*sub_fraction) .and. do_inflate) then

      ! Having found an EDC compliant parameter vector, we want to do a MCMC
      ! search on inflated uncertainties. This inflation search allows us to
      ! more easily move towards higher likelihoods, where the inflation
      ! allows easier movement through parameter space.

      ! The inflated search phase will make use of three stages over which the
      ! inflation will be reduced. Phase 1 used half of the allocated
      ! iterations for the inflation while the second and third is

      ! Set flag to indicate this phase has occurred and make a record of the
      ! total iterations to be attempted
      sub_sample_complete = .true.; nOUT_save = MCO%nOUT

      ! Report to the user
      write (*, *) "Beginning parameter search on sample size normalised likelihoods"

      MCO%nOUT = nint(dble(nout_save)*sub_fraction)
      MCO%fADAPT = 1d0
      MCO%fixedpars = .true. ! start from end points of EDC phase
      ! Second phase, run Mcmc with sub scaling
      call update_obs_scaling_nsamples
      call run_parallel_mcmc(scaled_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains=nchains)
      MCO%fixedpars = .true.
      do i = 1, nchains
         ! Use the best parameter set as the starting point for the next stage
         MCOUT_list(i)%pars(1:PI%npars) = MCOUT_list(i)%bestpars(1:PI%npars)

         ! Leave parameter and covariance structures as they come out form the
         ! sub-sample-but reset the number of samples used in the update
         ! weighting
         if (MCOUT_list(i)%cov .and. MCOUT_list(i)%use_multivariate) then
            ! TODO check this branch is happening
            write (*, *) "in this branch"
            MCOUT_list(i)%Nparvar = MCO%N_before_mv*PI%npars + 1
            write (*, *) "Set Nparvar to", MCOUT_list(i)%Nparvar
         else
            ! reset the parameter step size at the beginning of each attempt
            call reset_stats(MCOUT_list(i), PI%npars)
         end if  ! do we need a new covariance matrix or can we use the existing one?
      end do

      MCOUT_list(i)%nos_iterations = 0

   end if

   ! Restore module variables needed for the run-these components could be split
   ! into two subroutines to avoid double calling of file name creation
   ! components.
   call read_options(solution_wanted, freq_print, freq_write, outfile, MCO)

   MCO%nOUT = nout_save*(1 - sub_fraction)  !number of steps in final phase
   MCO%fixedpars = .true.

   ! Update the user
   write (*, *) "Beginning parameter search in real likelihoods"
   write (*, *) "Nos iterations to be proposed = ", MCO%nOUT - MCOUT_list(1)%nos_iterations

   ! Call the main MCMC
   ! The specific normalisation of the cost function is determined here.
   ! But to avoid getting through the EDC do_inflate sections before finding
   ! out that the cost_function_scaling has not been set correctly,
   ! ensure code after the command line read (above) has been correctly maintained
   if (cost_func_scaling_dble == 0) then
      call update_obs_scaling_normal
   else if (cost_func_scaling_dble == 1) then
      call update_obs_scaling_nsamples
   else if (cost_func_scaling_dble == 2) then
      call update_obs_scaling_sqrt_nsamples
   else if (cost_func_scaling_dble == 3) then
      call update_obs_scaling_log_nsamples
   end if  ! cost_func_scaling_dble ==
   do i = 1, nchains
   end do
   call run_parallel_mcmc(scaled_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains=nchains)

   ! Let the user know we are done
   write (*, *) "AP-MCMC done now, moving on ..."

   ! tidy up by closing all files
   do i = 1, nchains
      call close_output_files(i)
   end do

   ! Final message to the user
   write (*, *) "==========================================================="
   write (*, *) "==== CARDAMOM analysis for the current chain completed ===="
   write (*, *) "==========================================================="
   write (*, *) "=========================Honestly=========================="
contains

end program cardamom_framework
