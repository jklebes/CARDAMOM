program cardamom_framework
   use math_functions, only: rnstrt, idum  ! TODO redo random seeds
   use DEMCz, only: demczOPT, run_demcz
   use cardamom_MHMCMC, only: MCMC_OUTPUT, mcmc_options
   use model_shared, only: PI, initialize_carbon_model
   use cardamom_structures, only: DATAin
   use cardamom_io, only: initialize, &
                          read_options, &
                          restart_flag, &
                          update_for_restart_simulation, &
                        update_obs_scaling_normal, update_obs_scaling_nsamples, &
                        update_obs_scaling_sqrt_nsamples, update_obs_scaling_log_nsamples
   use samplers_io, only: open_output_files, &
                          check_for_existing_output_files, &
                          write_covariance_matrix, &
                          close_output_files, write_covariance_info
   use MHMCMC_StressTests, only: StressTest_likelihood_fct, StressTest_sublikelihood_fct, prepare_for_stress_test
   use model_likelihood_wrapper, only: model_likelihood_fct, edc_model_likelihood_fct, scaled_model_likelihood_fct
   use cardamom_main_utils

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
   double precision:: ll
   type(MCMC_OUTPUT):: MCOUT
   type(MCMC_OUTPUT), dimension(:), allocatable:: MCOUT_list  ! for parallel-could keep single here and make interface
   type(mcmc_options):: edc_MCO
   type(DEMCZOPT):: MCO

   ! TODO not to hardcode, from command line argument
   integer:: nchains = 46
   integer:: i

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
   do i = 1, nchains
      call initialize_stats(MCOUT_list(i), PI%npars)
   end do
   end if
   ! having filled PI%npars from model file, we can allocate stats array in MCOUT

   ! load module variables needed for restart check
   ! NOTE: THIS MUST HAPPEN BEFORE CHECKING FOR RESTART
   call read_options(solution_wanted, freq_print, freq_write, outfile, MCO, MCOUT)
   ! check whether this is a restart?
   ! PI lives in model_shared and its info can be read after call to initiialize_model
   ! TODO not sure about MCO at this point
   !do i = 1, nchains
   !   call check_for_existing_output_files(PI%npars, MCO%nOUT, MCO%nWRITE, sub_fraction, &
   !                                        MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile, i)
   !end do
   ! Initialise MCMC output, possibly a bit of a redundent subroutine...
   !call initialise_mcmc_output  ! TODO now happens in run_mcmc if not restart
   ! Open the relevant output files TODO now happens in run_mcmc
   !call open_output_files(MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile)

   ! Report which model ID we are using
   write (*, *) "Running model version ", DATAin%ID  ! TODO where does DATAin live and where did it get filled

   ! Check whether we are doing a stress test again
   if (DATAin%ID < 0) then
      !TODO move to testing

      ! We are doing a stress test
      write (*, *) "Carrying out a stress test analysis"
      write (*, *) "Any existing files will be ignored"

      ! Reset interations counter
      MCOUT%nos_iterations = 0
      ! Ensure that we use a random starting point
      ! MCO%randparini = .true.
      ! MCO%fixedpars  = .false.
      ! restart_flag = .false. ! default all random if no further arguments passed to run_mcmc

      ! Let the user know how many more we will propose
      write (*, *) "Nos iterations to be proposed = ", MCO%nOUT
      ! Call stresstest sampling
      write(*,*) nchains
      call run_demcz(stresstest_likelihood_fct, PI, MCO, MCOUT_list, stresstest_likelihood_fct, nchains = nchains)
      ! Tell the user the best parameter set
      print *, "Best parameters = ", MCOUT%bestpars

   else  ! We are not doing a stress test

      ! Begin search for initial conditions
      write (*, *) "Beginning search for initial parameter conditions"
      ! Determine initial values, this requires using the AP-MCMC
      call find_edc_initial_values(edc_MCO, MCOUT_list, nchains)

      ! Reset the iterations counter-if not then the wrong number of iterations will be attempted
      mco%mcmc_options = edc_mco
      do i = 1, nchains

         ! Reset the MCMC parameters for the next stage
         call read_options(solution_wanted, freq_print, freq_write, outfile, MCO, MCOUT_list(i))

         ! Reset stepsize and covariance for main DRAM-MCMC
         call reset_stats(MCOUT_list(i), PI%npars)

         if (restart_flag) then
            ! Restarting an old one
            print *, "beginning restart simulation"
            ! now begin update of model timing variables and parameter values if this is a
            ! restart. NOTE that this include information determining the number of
            ! iterations already completed...
            call update_for_restart_simulation(MCO, MCOUT_list(i))
         else
            ! Brand new analysis
            !print*,"writing initial covariance matrix"
            ! write out first covariance matrix, this will be compared with the final covariance matrix
            !if (MCO%nWRITE > 0) then %TODO problem because files not opened yet?
                !call write_covariance_matrix(mcout_list(i)%covariance, PI%npars, .true., i)
                !call write_covariance_info(mcout_list(i)%meanpar, mcout_list(i)%Nparvar, PI%npars, i)
            !endif
            !...so the reset for nos_iterations must only occur when not a restart run
            MCOUT_list(i)%nos_iterations = 0
         end if  ! restart run or not
      end do


      ! Restore module variables needed for the run-these components could be split
      ! into two subroutines to avoid double calling of file name creation
      ! components.
      call read_options(solution_wanted, freq_print, freq_write, outfile, MCO, MCOUT)
      ! Since they all get the same nout setting, assume all sub_model phase simulation 
      ! were same length  ! TODO

      ! Update the user
      write (*, *) "Beginning parameter search in real likelihoods"
      write (*, *) "Nos iterations to be proposed = ", MCO%nOUT, MCOUT_list(i)%nos_iterations


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
      ! TODO scaled model likelihood fct
      call run_demcz(scaled_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains = nchains, restart=.true.)

      ! Let the user know we are done
      write (*, *) "AP-MCMC done now, moving on ..."

   end if  ! stress test or not

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

   !
   !------------------------------------------------------------------

end program cardamom_framework
