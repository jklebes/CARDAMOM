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
 use math_functions, only:  rnstrt, idum  ! TODO redo random seeds
 use MHMCMC, only: MCMC_OUTPUT, MCMC_OPTIONS, run_mcmc, run_parallel_mcmc 
 use model_shared, only: PI
 use cardamom_structures, only: DATAin 
 use cardamom_io, only: initialize, &
                        read_options, & 
                        restart_flag,   &
                        update_for_restart_simulation 
 use samplers_io, only:  open_output_files, &
                        check_for_existing_output_files,  &
                        write_covariance_matrix, &
                        close_output_files, write_covariance_info
 !use MHMCMC_module, only: MHMCMC, par_minstepsize, par_initstepsize, N_before_mv
 use MHMCMC_StressTests, only: StressTest_likelihood_fct, StressTest_sublikelihood_fct, prepare_for_stress_test
 !use model_likelihood_module, only: model_likelihood, &
 !   sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood  ! to replace soon with wrappers
 use model_likelihood_wrapper, only: model_likelihood_fct, log_model_likelihood_fct, sqrt_model_likelihood_fct, &
     sub_model_likelihood_fct, edc_model_likelihood_fct
 use CARBON_MODEL_MOD, only: initialize_carbon_model
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

 implicit none

 ! declare local variables
 character(350):: infile, outfile, solution_wanted_char, freq_print_char, &
                   freq_write_char, do_inflate_char, cost_func_scaling_char
 integer:: solution_wanted, freq_print, freq_write, time1, time2, time3, n, &
            nOUT_save, do_inflate_dble, cost_func_scaling_dble
 logical:: do_inflate = .false.
 logical:: sub_sample_complete = .false.
 double precision:: sub_fraction = 0.2d0
 double precision:: ll
 !double precision:: idum  ! TODO redo seeds
 type(MCMC_OUTPUT):: MCOUT
 type(MCMC_OUTPUT), dimension(:), allocatable:: MCOUT_list  ! for parallel-could keep single here and make interface
 type(MCMC_OPTIONS):: MCO

 ! TODO not to hardcode, from command line argument
 integer:: nchains = 4
 integer:: i

 allocate(MCOUT_list(nchains))

 ! user update
 write(*,*)"Beginning read of the command line"

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
 read(solution_wanted_char, '(I10)') solution_wanted
 read(freq_print_char, '(I10)') freq_print
 read(freq_write_char, '(I10)') freq_write
 read(do_inflate_char, '(I10)') do_inflate_dble
 read(cost_func_scaling_char, '(I10)') cost_func_scaling_dble

 ! Assign inflate logical condition
 if (do_inflate_dble == 1) do_inflate = .true.
 ! Sanity check of the cost_function_scaling_char
 if (cost_func_scaling_dble > -1 .and. cost_func_scaling_dble < 4 .and. freq_write > 0) then
     ! All is well
 else
     ! All is not well-complain
     print*, "ERROR: Command line argument to specify the cost function or write to file frequency is incorrect."
     print*, "Command line should have 7 arguments (in addition to the cardamom.exe)."
     print*, "These are: "
     print*, "1) input file path."
     print*, "2) output file path-note that PARS, STEP, COV, COVINFO will be appended to this name path outfile."
     print*, "3) No. of parameter proposals to make."
     print*, "4) Iteration freq. for printing to screen (main MCMC phase only)."
     print*, "5) Iteration freq. for writing results to files."
     print*, "6) do sample size normalisation phase - "
     print*, "  0 = FALSE"
     print*, "  1 = TRUE"
     print*, "7) Select likelihood cost function - "
     print*, "  0 = no scaling"
     print*, "  1 = scaling by sample size (n)"
     print*, "  2 = scaling by sqrt(n)"
     print*, "  3 = scaling by log(n)"
     stop
 end if

 ! user update
 write(*,*)"Command line options read, moving on now"

 ! TODO get n seeds, here or smoewhere else
 ! TODO make sure seeds are saved 
 ! seed the random number generator
 ! determine unique (sort of) seed value; based on system time
 call system_clock(time1, time2, time3)
 ! set seed value outside of the function, idum must be a negative number
 idum = dble(time1+time2+time3)
 call rnstrt(nint(idum))

 ! Determine whether ot not we are doing a real analysis or running a stress trest
 if (trim(infile) == "StressTest") then
     ! call special functions to prepare for stress test
     call prepare_for_stress_test(infile, outfile)  ! sets cardamom_structures:: DATAin
 else
    call initialize(infile) ! = initialize_parinfo, read_check_binary_data, initialize_model  ! sets cardamom_structures:: DATAin
    call initialize_carbon_model(nchains)
    call initialize_stats(MCOUT, PI%npars)
 end if
 ! having filled PI%npars from model file, we can allocate stats array in MCOUT

 ! load module variables needed for restart check
 ! NOTE: THIS MUST HAPPEN BEFORE CHECKING FOR RESTART
 call read_options(solution_wanted, freq_print, freq_write, outfile, MCO, MCOUT)
 ! check whether this is a restart?
 ! PI lives in model_shared and its info can be read after call to initiialize_model
 ! TODO not sure about MCO at this point
 do i = 1, nchains
 call check_for_existing_output_files(PI%npars, MCO%nOUT, MCO%nWRITE, sub_fraction, &
                                      MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile, i)
 end do
 ! Initialise MCMC output, possibly a bit of a redundent subroutine...
 !call initialise_mcmc_output  ! TODO now happens in run_mcmc if not restart
 ! Open the relevant output files TODO now happens in run_mcmc
 !call open_output_files(MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile)

 ! Initialise counters used to track the output of parameter sets
 !TODO now each chain has its own 
 !io_space%io_buffer_count = 0
 ! io_space%io_buffer = min(1000, max(10, (MCO%nOUT/MCO%nWRITE) / 10))

 ! Allocate variables used in io buffering, 
 ! these could probably be moved to a more sensible place within cardamom_io.f90 DONE
 !allocate(io_space%variance_buffer(PI%npars, io_space%io_buffer), &
 !         io_space%meanpars_buffer(PI%npars, io_space%io_buffer), &
 !         io_space%pars_buffer(PI%npars, io_space%io_buffer), &
 !         io_space%prob_buffer(io_space%io_buffer), &
 !         io_space%nsample_buffer(io_space%io_buffer), &
 !         io_space%accept_rate_buffer(io_space%io_buffer))

 ! Report which model ID we are using
 write(*,*) "Running model version ", DATAin%ID  ! TODO where does DATAin live and where did it get filled


 ! Check whether we are doing a stress test again
 if (DATAin%ID < 0) then
     !TODO move to testing

     ! We are doing a stress test
     write(*,*)"Carrying out a stress test analysis"
     write(*,*)"Any existing files will be ignored"

     ! Reset interations counter
     MCOUT%nos_iterations = 0
     ! Ensure that we use a random starting point
     ! MCO%randparini = .true.
     ! MCO%fixedpars  = .false.
     ! restart_flag = .false. ! default all random if no further arguments passed to run_mcmc

     ! Do we do the initial MCMC period where we normalise the likelihood by
     ! number of observations
     ! This process allows for very bad starting points to more easily move
     ! towards the general area of the observatons.
     if (MCOUT%nos_iterations < (MCO%nOUT*sub_fraction) .and. do_inflate) then

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
         sub_sample_complete = .true. ; nOUT_save = MCO%nOUT

         ! Report to the user
         write(*,*)"Beginning parameter search on sample size normalised likelihoods"

         MCO%nOUT = nint(dble(nOUT_save) * sub_fraction) - MCOUT%nos_iterations
         write(*,*)"Nos iterations to be proposed = ",MCO%nOUT
         MCO%fADAPT = 1d0 !; MCO%nADAPT = 1000
         !call run_mcmc(1d0, StressTest_likelihood, StressTest_sublikelihood)

         allocate(MCOUT_list(nchains))
         MCOUT_list(1) = MCOUT
         call run_parallel_mcmc(stresstest_sublikelihood_fct, PI, MCO, MCOUT_list, stresstest_likelihood_fct, nchains = nchains)
         MCOUT = MCOUT_list(1)

         ! Use the best parameter set as the starting point for the next stage
! REALLY NOT SURE I SHOULD BE DOING THIS-SHOULD BE PROGRESSING FROM THE LAST ACCEPTED PARAMETER SET?
         MCOUT%pars = MCOUT%bestpars
         ! Leave parameter and covariance structures as they come out form the
         ! sub-sample-but reset the number of samples used in the update
         ! weighting
         if (MCOUT%cov .and. MCOUT%use_multivariate) then
             MCOUT%Nparvar = (MCO%N_before_mv*dble(PI%npars)) + 1d0 
         else
             call reset_stats(MCOUT, PI%npars)
             ! reset the parameter step size at the beginning of each attempt  ! TODO where does this comment come from, to do?
         endif  ! do we need a new covariance matrix or can we use the existing one?

         ! Assume that sub-sampling process, if completed, will use 10 % of the
         ! simulation time therefore we want to adjust the output frequency to
         ! correct for this
         MCO%nOUT = max(1, nOUT_save-MCOUT%nos_iterations)
         ! by pass read_options file for StressTest special case
         MCO%append = .true.
         MCO%nADAPT = 1000
         MCO%fADAPT = 0.5d0
         MCO%randparini = .false.
         MCO%fixedpars  = .true.

     end if  ! restart flag

     ! Let the user know how many more we will propose
     write(*,*)"Nos iterations to be proposed = ",MCO%nOUT
     ! Call the AP-MCMC
     MCOUT_list(1) = MCOUT
     call run_parallel_mcmc(stresstest_likelihood_fct, PI, MCO, MCOUT_list, stresstest_likelihood_fct, nchains = nchains)
     MCOUT = MCOUT_list(1)
     ! Tell the user the best parameter set
     print*,"Best parameters = ",MCOUT%bestpars

 else  ! We are not doing a stress test

     ! Begin search for initial conditions
     write(*,*) "Beginning search for initial parameter conditions"
     ! Determine initial values, this requires using the AP-MCMC
     call find_edc_initial_values(MCO, MCOUT_list, nchains) 

     ! Reset the iterations counter-if not then the wrong number of iterations will be attempted
     do i = 1, nchains
     !MCOUT_list(i)%pars = (/8.2541992449544470E-004,  0.29269005698013606,       0.10918820581034874,       0.65410628727175912, &
     !1.4116626415828890,        6.1425101753895630E-004,   3.2985120921817839E-003,   1.4904751343118651E-004, &
     !3.1947208234022437E-005, &
     !6.6555453603566697E-002,   62.881210637480095,        551.42881422649646,        3.8930261063292386E-002, 18.225050185374187, &
     !1161.4281598899104,        32.225636460595346,        98.522275278559022,       3.6863718165541406, 38.956182582098513, &
     !92.906961634604770,        240.14375534689168,        172.37962589728303  ,      5107.9915424005903,     0.17472418107258553, &
     !0.22608381258586158,        2283.2151160970793,        9.7709099175705525,        4.1058161174047723E-002, &
     !5.4968083425296944E-002,   4.8002121930932840E-002,   1.1374800014070333E-002,  0.17644136331124649/)
     MCOUT_list(i)%bestpars = MCOUT_list(i)%pars
     MCOUT_list(i)%nos_iterations = 0
                  write(*,*) "with sub_model_likelihood 2"
                  call sub_model_likelihood_fct(MCOUT_LIST(i)%pars, PI%npars, ll, i)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 4"
                      write(*,*) i, MCOUT_LIST(i)%pars 
                      stop 1 
                end if

     ! Reset the MCMC parameters for the next stage
     call read_options(solution_wanted, freq_print, freq_write, outfile, MCO, MCOUT_list(i))
                  write(*,*) "with sub_model_likelihood 3"
                  call sub_model_likelihood_fct(MCOUT_LIST(i)%pars, PI%npars, ll, i)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 5"
                      write(*,*) i, MCOUT_LIST(i)%pars 
                      stop 1 
                end if

     ! Reset stepsize and covariance for main DRAM-MCMC
     call reset_stats(MCOUT_list(i), PI%npars)
                  write(*,*) "with sub_model_likelihood 4"
                  call sub_model_likelihood_fct(MCOUT_LIST(i)%pars, PI%npars, ll, i)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 6"
                      write(*,*) i, MCOUT_LIST(i)%pars 
                      stop 1 
                end if

     if (restart_flag) then
         ! Restarting an old one
         print*, "beginning restart simulation"
         ! now begin update of model timing variables and parameter values if this is a
         ! restart. NOTE that this include information determining the number of
         ! iterations already completed...
         call update_for_restart_simulation(MCO, MCOUT_list(i))
     else
         ! Brand new analysis
         !print*,"writing initial covariance matrix"
         ! write out first covariance matrix, this will be compared with the final covariance matrix
         !if (MCO%nWRITE > 0) then
         !    call write_covariance_matrix(mcout_list(i)%covariance, PI%npars, .true., i)
         !    call write_covariance_info(mcout_list(i)%meanpar, mcout_list(i)%Nparvar, PI%npars, i)
         !endif
         !...so the reset for nos_iterations must only occur when not a restart run
         MCOUT_list(i)%nos_iterations = 0
     endif  ! restart run or not
     end do

     ! Do we do the initial MCMC period where we normalise the likelihood by number of observations
     ! This process allows for very bad starting points to more easily move towards the general area of the observatons.
     if (DATAin%total_obs > 0 .and. MCOUT%nos_iterations < (MCO%nOUT*sub_fraction) .and. do_inflate) then

         ! Having found an EDC compliant parameter vector, we want to do a MCMC
         ! search on inflated uncertainties. This inflation search allows us to
         ! more easily move towards higher likelihoods, where the inflation
         ! allows easier movement through parameter space.

         ! The inflated search phase will make use of three stages over which the
         ! inflation will be reduced. Phase 1 used half of the allocated
         ! iterations for the inflation while the second and third is

         ! Set flag to indicate this phase has occurred and make a record of the
         ! total iterations to be attempted
         sub_sample_complete = .true. ; nOUT_save = MCO%nOUT

         ! Report to the user
         write(*,*)"Beginning parameter search on sample size normalised likelihoods"

         MCO%nOUT = nint(dble(nOUT_save) * sub_fraction) - MCOUT%nos_iterations
         write(*,*)"Nos iterations to be proposed = ",MCO%nOUT
         MCO%fADAPT = 1d0 !; MCO%nADAPT = 1000
         MCO%nwrite = 0
         MCO%nprint = 1000
         do i = 1, nchains
         write(*,*) i, "MCOUT_list(i)%pars", MCOUT_list(i)%pars
                  write(*,*) "with sub_model_likelihood 5"
                  call sub_model_likelihood_fct(MCOUT_LIST(i)%pars, PI%npars, ll, i)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 7"
                      write(*,*) i, MCOUT_LIST(i)%pars 
                      stop 1 
                end if
         end do
         call run_parallel_mcmc(sub_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains = nchains, restart=.true.)
         !call run_mcmc(1d0, model_likelihood, sub_model_likelihood)
         ! call MHMCMC(PI, MCO, model_likelihood, sub_model_likelihood)
         ! Use the best parameter set as the starting point for the next stage
         MCOUT%pars(1:PI%npars) = MCOUT%bestpars(1:PI%npars)
         MCO%fixedpars  = .true.
         ! Leave parameter and covariance structures as they come out form the
         ! sub-sample-but reset the number of samples used in the update
         ! weighting
         if (MCOUT%cov .and. MCOUT%use_multivariate) then
             MCOUT%Nparvar = (MCO%N_before_mv*dble(PI%npars)) + 1d0
         else
             ! reset the parameter step size at the beginning of each attempt
             ! TODO fct for this
             MCOUT%parvar = 1d0; MCOUT%Nparvar = 0d0
             ! Covariance matrix cannot be set to zero therefore set initial
             ! value to a small positive value along to variance access
             MCOUT%covariance = 0d0; MCOUT%meanpar = 0d0; MCOUT%cov = .false.
             MCOUT%use_multivariate = .false.
             do n = 1, PI%npars
                MCOUT%covariance(n, n) = 1d0
             end do
         endif  ! do we need a new covariance matrix or can we use the existing one?

     end if  ! restart flag

     ! Restore module variables needed for the run-these components could be split
     ! into two subroutines to avoid double calling of file name creation
     ! components.
     call read_options(solution_wanted, freq_print, freq_write, outfile, MCO, MCOUT)

     ! Update the user
     write(*,*)"Beginning parameter search in real likelihoods"
     write(*,*)"Nos iterations to be proposed = ",MCO%nOUT

     ! Call the main MCMC
     ! The specific normalisation of the cost function is determined here.
     ! But to avoid getting through the EDC do_inflate sections before finding
     ! out that the cost_function_scaling has not been set correctly, 
     ! ensure code after the command line read (above) has been correctly maintained
     if (cost_func_scaling_dble == 0) then
        ! Caution: order of loglikelihood function arguments is being switched so that the first 
        ! function in the (maybe scaled) one to do the sampling calculation with, second optional 
        ! function argument is the one for writing only
         call run_parallel_mcmc(model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains = nchains, restart=.true.)
     else if (cost_func_scaling_dble == 1) then
         call run_parallel_mcmc(sub_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains = nchains, restart=.true.)
     else if (cost_func_scaling_dble == 2) then
         call run_parallel_mcmc(sqrt_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains = nchains, restart=.true.)
     else if (cost_func_scaling_dble == 3) then
         call run_parallel_mcmc(log_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains = nchains, restart=.true.)
     !else if (cost_func_scaling_dble == 4) then
     !    call MHMCMC(1d0, model_likelihood, log_model_likelihood_dtm)
     end if  ! cost_func_scaling_dble == 
     
     ! Let the user know we are done
     write(*,*)"AP-MCMC done now, moving on ..."

 end if  ! stress test or not

 ! tidy up by closing all files
 do i = 1, nchains
 call close_output_files(i)
 end do

 ! Final message to the user
 write(*,*)"==========================================================="
 write(*,*)"==== CARDAMOM analysis for the current chain completed ===="
 write(*,*)"==========================================================="
 write(*,*)"=========================Honestly=========================="
  contains 


  !
  !------------------------------------------------------------------
  !
  subroutine find_edc_initial_values(MCO, MCOUT_list, nchains)
    !! subroutine deals with the determination of initial parameter and initial
    !! conditions which are consistent with EDCs
    !! pre-loop, Run MCMC sampler with modified likelihood fct  
    use model_shared, only: PI
    use MHMCMC, only: MCMC_OUTPUT, MCMC_OPTIONS, MCSTATS, run_mcmc, run_parallel_mcmc
    !use model_likelihood_module, only: model_likelihood, &
    !sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood  ! to replace soon with wrappers
    use model_likelihood_wrapper  ! TODO next refactoring step
    use cardamom_structures, only: DATAin 


    implicit none

    ! declare local variables
    integer, intent(in):: nchains
    type(MCMC_OUTPUT), dimension(:), allocatable, intent(inout):: MCOUT_list
    type(MCMC_OUTPUT), dimension(:), allocatable:: MCOUT_list_tmp
    type(MCMC_OPTIONS), intent(out):: MCO
    integer:: n, i, counter_local(nchains), EDC_iter, nOUT_save, nWRITE_save, nADAPT_save, j
    integer:: success_count
    logical:: append_save
    logical:: restart(nchains)
    double precision:: ll
    double precision:: PEDC(nchains), PEDC_prev(nchains), ML, ML_prior, P_target
    double precision, dimension(PI%npars+1):: EDC_pars
    double precision, dimension(PI%npars):: parini  ! local variable, or array

    allocate(MCOUT_list_tmp(nchains))

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
    MCO%fixedpars  = .true. ! TLS: changed from .false. for testing 16/12/2019

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
        do while (success_count < nchains)!(PEDC < 0d0)

           write(*,*)"Beginning EDC search attempt "
           write(*,*)  nchains, "chains working ... "
           MCO%randparini = .false.
        !$omp parallel do private(ll)
           do i = 1, nchains
           ! call the MHMCMC directing to the appropriate likelihood function
           call run_mcmc(edc_model_likelihood_fct, PI, MCO, MCOUT_list_tmp(i), model_likelihood_fct, restart = restart(i), chainid = i)
           restart(i) = .true.

           ! turn off random selection for initial values
           write(*,*)"...intermediate EDC search progress check"

           ! store the best parameters from that loop
           PEDC(i) = MCOUT_list_tmp(i)%bestll 
           write(*,*) "Found best loglikelihood",PEDC(i) 
           ! if any chains's MCOUT object is success (reached loglikelihood = 0), 
           ! copy it to MCOUT_list
                  !$omp critical
              if (MCOUT_list_tmp(i)%ll >= (0d0-10*epsilon(1.0d0)) ) then
                  success_count = success_count+1
                  write(*,*) "Found ", success_count, "EDC-compatible starting points (", i, ")"
                  write(*,*) "with sub_model_likelihood"
                  call sub_model_likelihood_fct(MCOUT_LIST_tmp(i)%pars, PI%npars, ll, i)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 8"
                      write(*,*) i, MCOUT_LIST_tmp(i)%pars 
                      stop 1 
                end if
                  if (success_count <= nchains) then  
                  MCOUT_list(success_count) = MCOUT_list_tmp(i)
                  endif
                  call reset_stats(MCOUT_list_tmp(i), PI%npars)
                  MCO%randparini = .true. ! TODO problem
                  PEDC_prev(i) = -1000d0
                  MCOUT_list_tmp(i)%PARS = DATAin%parpriors(1:PI%npars)
                  restart(i) = .false.
                  counter_local(i) = 0

                  write(*,*) "with edc_model_likelihood"
                  write(*,*) MCOUT_LIST(success_count)%pars, PI%npars, ll, i
                  call edc_model_likelihood_fct(MCOUT_LIST(success_count)%pars, PI%npars, ll, i)
                  write(*,*) "with sub_model_likelihood"
                  call sub_model_likelihood_fct(MCOUT_LIST(success_count)%pars, PI%npars, ll, i)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 10"
                      write(*,*) i, MCOUT_LIST(success_count)%pars 
                      stop 1 
              endif 

          end if
                  !$omp end critical

           ! keep track of attempts
           counter_local(i) = counter_local(i) +1 
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
               write(*,*) "resetting to initial"
           else
               PEDC_prev(i) = PEDC(i)
           endif
        end do
        !$omp end parallel do


        end do  ! for while condition

        do i = 1, nchains
        do j = 1, nchains

                 write(*,*) "check", i, j
                  write(*,*) "with edc_model_likelihood"
                  write(*,*) MCOUT_LIST(i)%pars, PI%npars, ll, i
                  call edc_model_likelihood_fct(MCOUT_LIST(i)%pars, PI%npars, ll, j)
                  write(*,*) ll
                  write(*,*) "with sub_model_likelihood"
                  call sub_model_likelihood_fct(MCOUT_LIST(i)%pars, PI%npars, ll, j)
                  write(*,*) ll
                  if (ll < -9999999) then
                      write (*,*) "Infinity check 9"
                      write(*,*) j, MCOUT_LIST(i)%pars 
                end if
                end do
        end do

    endif  ! if for restart

    ! reset so that currently saved parameters will be used
    ! starting point in main MCMC
    ! PI%parfix(1:PI%npars) = 0  ! TODO
    !MCOUT%bestpars = 0d0

  end subroutine find_edc_initial_values


end program cardamom_framework
