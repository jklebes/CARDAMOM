
program cardamom_framework

 use math_functions, only: idum, randn, rnstrt, inverse_matrix
 use MCMCOPT, only: MCO, MCOUT, PI, initialise_mcmc_output
 use cardamom_structures, only: DATAin
 use cardamom_io, only: read_pari_data, read_options, open_output_files, &
                        check_for_existing_output_files,restart_flag,   &
                        update_for_restart_simulation, write_covariance_matrix, &
                        close_output_files, write_covariance_info
 use MHMCMC_module, only: MHMCMC, par_minstepsize, par_initstepsize
 use model_likelihood_module, only: model_likelihood, & 
                                    find_edc_initial_values, &
                                    sub_model_likelihood

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

 implicit none

 ! declare local variables
 character(350) :: infile, outfile, solution_wanted_char, freq_print_char, freq_write_char
 integer :: solution_wanted, freq_print, freq_write, time1, time2, time3, i, n, nOUT_save
 logical :: do_inflate 

 ! user update
 write(*,*)"Beginning read of the command line"

 ! read user options from the command line
 call get_command_argument(1 ,infile)
 call get_command_argument(2 ,outfile)
 call get_command_argument(3 ,solution_wanted_char)
 call get_command_argument(4 ,freq_print_char)
 call get_command_argument(5 ,freq_write_char)

 ! now convert relevant ones to integeter
 read(solution_wanted_char,'(I10)') solution_wanted
 read(freq_print_char,'(I10)') freq_print
 read(freq_write_char,'(I10)') freq_write

 ! user update
 write(*,*)"Command line options read, moving on now"

 ! seed the random number generator
 ! determine unique (sort of) seed value; based on system time
 call system_clock(time1,time2,time3)
 ! set seed value outside of the function, idum must be a negative number
 idum = dble(time1+time2+time3)
 call rnstrt(nint(idum))

 ! read input data (DATAin located in module)
 call read_pari_data(infile)

 ! load module variables needed for restart check 
 ! NOTE: THIS MUST HAPPEN BEFORE CHECKING FOR RESTART
 call read_options(solution_wanted,freq_print,freq_write,outfile)
 ! check whether this is a restart?
 call check_for_existing_output_files(PI%npars,MCO%nOUT,MCO%nWRITE,MCO%sub_fraction, &
                                      MCO%outfile,MCO%stepfile,MCO%covfile,MCO%covifile)
 ! Initialise MCMC output, possibly a bit of a redundent subroutine...
 call initialise_mcmc_output
 ! Open the relevant output files
 call open_output_files(MCO%outfile,MCO%stepfile,MCO%covfile,MCO%covifile)

 ! Report which model ID we are using
 write(*,*) "Running model version ", DATAin%ID

 ! Begin search for initial conditions
 write(*,*) "Beginning search for initial parameter conditions"
 ! Determine initial values, this requires using the MHMCMC
 call find_edc_initial_values
 ! Reset the MCMC parameters for the next stage
 call read_options(solution_wanted,freq_print,freq_write,outfile)

 ! Reset stepsize and covariance for main DRAM-MCMC
 PI%Nparvar = 0d0 ; PI%parvar = 0d0
 PI%covariance = 0d0 ; PI%mean_par = 0d0 
 PI%cov = .false. ; PI%use_multivariate = .false.
 do n = 1, PI%npars
    PI%covariance(n,n) = 1d0
 end do

 if (restart_flag) then
     print*, "beginning restart simulation"
     ! now begin update of model timing variables and parameter values if this is a
     ! restart. NOTE that this include information determining the number of
     ! iterations already completed...
     call update_for_restart_simulation
 else
     print*,"writing initial covariance matrix"
     ! write out first covariance matrix, this will be compared with the final covariance matrix
     if (MCO%nWRITE > 0) then
         call write_covariance_matrix(PI%covariance,PI%npars,.true.)
         call write_covariance_info(PI%mean_par,PI%Nparvar,PI%npars)
     endif
     !...so the reset for nos_iterations must only occur when not a restart run
     MCOUT%nos_iterations = 0
 endif

 do_inflate = .true.
! if (DATAin%total_obs > 0 .and. .not.restart_flag .and. do_inflate) then
 if (DATAin%total_obs > 0 .and. MCOUT%nos_iterations < (MCO%nOUT*MCO%sub_fraction) .and. do_inflate) then

     ! Having found an EDC compliant parameter vector, we want to do a MCMC
     ! search on inflated uncertainties. This inflation search allows us to
     ! more easily move towards higher likelihoods, where the inflation
     ! allows easier movement through parameter space.

     ! The inflated search phase will make use of three stages over which the
     ! inflation will be reduced. Phase 1 used half of the allocated
     ! iterations for the inflation while the second and third is

     ! Set flag to indicate this phase has occurred and make a record of the
     ! total iterations to be attempted
     MCO%sub_sample_complete = .true. ; nOUT_save = MCO%nOUT
     ! Estimate the maximum number of iterations to be used in the analysis
!     MCO%nOUT = nint(dble(nOUT_save) * MCO%sub_fraction * (1d0/3d0))

     ! Report to the user
     write(*,*)"Beginning parameter search on inflated likelihoods"

     MCO%nOUT = nint(dble(nOUT_save) * MCO%sub_fraction) - MCOUT%nos_iterations
     call MHMCMC(1d0,sub_model_likelihood)
     ! Leave parameter and covariance structures as they come out form the
     ! sub-sample - but reset the number of samples used in the update
     ! weighting
     if (PI%cov .and. PI%use_multivariate) then
!         PI%Nparvar = 1d0
     else
         ! reset the parameter step size at the beginning of each attempt
         PI%parvar = 1d0 ; PI%Nparvar = 0d0
         ! Covariance matrix cannot be set to zero therefore set initial
         ! value to a small positive value along to variance access
         PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
         PI%use_multivariate = .false.
         do n = 1, PI%npars
            PI%covariance(n,n) = 1d0
         end do
     endif

     ! Estimate the maximum number of iterations to be used in the inflation analysis
!     MCO%nOUT = nint(dble(nOUT_save) * MCO%sub_fraction * (1d0/3d0))
     ! set the factor by which the uncertainty will be inflated by
!     MCO%inflation_factor = 5d0
!     do i = 1, 3
!         ! Report to the user
!         write(*,*)"Current inflation factor = ",MCO%inflation_factor
!
!         ! Run a short MCMC
!         call MHMCMC(-dble(DATAin%total_obs)/MCO%inflation_factor,sub_model_likelihood)
!         call MHMCMC(1d0,sub_model_likelihood)
!         ! use the final accepted parameter set as the starting point (or best?)
!!         PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)        
!         ! Leave parameter and covariance structures as they come out form the
!         ! sub-sample - but reset the number of samples used in the update
!         ! weighting
!         if (PI%cov .and. PI%use_multivariate) then
!             PI%Nparvar = 1d0
!         else
!             ! reset the parameter step size at the beginning of each attempt
!             PI%parvar = 1d0 ; PI%Nparvar = 0d0
!             ! Covariance matrix cannot be set to zero therefore set initial
!             ! value to a
!             ! small positive value along to variance access
!             PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
!             PI%use_multivariate = .false.
!             do n = 1, PI%npars
!                PI%covariance(n,n) = 1d0
!             end do
!         endif
!
!         ! update the factor by which the uncertainty will be inflated by
!         if (i == 1) MCO%inflation_factor = 2d0 !MCO%inflation_factor - 1d0
!         if (i == 2) MCO%inflation_factor = 1.5d0
!
!     end do ! Inflation loop

     ! Check whether we should be running the large inflation stage
!     if (MCOUT%nos_iterations < nOUT_save*MCO%sub_fraction*(1d0/3d0)) then
!
!         ! Report to the user
!         write(*,*)"Beginning parameter search on large inflated likelihoods"
!
!         ! set the factor by which the uncertainty will be inflated by
!         MCO%inflation_factor = 4.5d0
!         call MHMCMC(-dble(DATAin%total_obs*0.5d0)/MCO%inflation_factor,sub_model_likelihood)
!         ! store the best parameters from that loop
!         PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
!         ! Leave parameter and covariance structures as they come out form the
!         ! sub-sample - but reset the number of samples used in the update
!         ! weighting
!         if (PI%cov .and. PI%use_multivariate) then
!             PI%Nparvar = 1d0
!        ! else
!             ! reset the parameter step size at the beginning of each attempt
!             PI%parvar = 1d0 ; PI%Nparvar = 0d0
!             ! Covariance matrix cannot be set to zero therefore set initial
!             ! value to a
!             ! small positive value along to variance access
!             PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
!             PI%use_multivariate = .false.
!             do n = 1, PI%npars
!                PI%covariance(n,n) = 1d0
!             end do
!         endif
!
!     end if ! MCOUT%nos_iterations > nOUT_save*MCO%sub_fraction*(1d0/3d0)
!
!     ! Check whether we should be running the intermediate inflation stage
!     if (MCOUT%nos_iterations < nOUT_save*MCO%sub_fraction*(2d0/3d0)) then
! 
!         ! Report to the user
!         write(*,*)"Beginning parameter search on intermediate inflated likelihoods"
!         MCO%inflation_factor = 3.0d0
!         call MHMCMC(-dble(DATAin%total_obs*0.5d0)/MCO%inflation_factor,sub_model_likelihood)
!         ! store the best parameters from that loop
!         PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
!         ! Leave parameter and covariance structures as they come out form the
!         ! sub-sample - but reset the number of samples used in the update
!         ! weighting
!         if (PI%cov .and. PI%use_multivariate) then
!             PI%Nparvar = 1d0
!         else
!             ! reset the parameter step size at the beginning of each attempt
!             PI%parvar = 1d0 ; PI%Nparvar = 0d0
!             ! Covariance matrix cannot be set to zero therefore set initial
!             ! value to a small positive value along to variance access
!             PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
!             PI%use_multivariate = .false.
!             do n = 1, PI%npars
!                PI%covariance(n,n) = 1d0
!             end do
!         endif ! PI%cov .and. PI%use_multivariate
!
!     end if ! MCOUT%nos_iterations > nOUT_save*MCO%sub_fraction*(2d0/3d0)
!
!     ! Well we can assume that if we are in this section that we must at least
!     ! want to carry out the small inflation...
!
!     ! Report to the user
!     write(*,*)"Beginning parameter search on small inflated likelihoods"
!     MCO%inflation_factor = 1.5d0
!     call MHMCMC(-dble(DATAin%total_obs*0.5d0)/MCO%inflation_factor,sub_model_likelihood)
!     ! store the best parameters from that loop
!     PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
!     ! Leave parameter and covariance structures as they come out form the
!     ! sub-sample - but reset the number of samples used in the update
!     ! weighting
!     if (PI%cov .and. PI%use_multivariate) then
!         PI%Nparvar = 1d0
!     else
!         ! reset the parameter step size at the beginning of each attempt
!         PI%parvar = 1d0 ; PI%Nparvar = 0d0
!         ! Covariance matrix cannot be set to zero therefore set initial
!         ! value to a small positive value along to variance access
!         PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
!         PI%use_multivariate = .false.
!         do n = 1, PI%npars
!            PI%covariance(n,n) = 1d0
!         end do
!     endif ! PI%cov .and. PI%use_multivariate
!
 end if ! restart flag

 ! Restore module variables needed for the run - these components could be split
 ! into two subroutines to avoid double calling of file name creation
 ! components.
 call read_options(solution_wanted,freq_print,freq_write,outfile)

 ! Update the user
 write(*,*)"Beginning parameter search in real likelihoods"
 write(*,*)"Nos iterations = ",MCO%nOUT
 ! Call the main MCMC
 call MHMCMC(1d0,model_likelihood)
 ! Let the user know we are done
 write(*,*)"MHMCMC done now, moving on ..."

 ! tidy up by closing all files
 call close_output_files

end program cardamom_framework
