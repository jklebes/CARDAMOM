
program cardamom_framework

 use math_functions, only: idum, randn, rnstrt
 use MCMCOPT, only: MCO, MCOUT, PI, initialise_mcmc_output
 use cardamom_structures, only: DATAin
 use cardamom_io, only: read_pari_data, read_options, open_output_files, &
                        check_for_existing_output_files,restart_flag,   &
                        update_for_restart_simulation, write_covariance_matrix, &
                        close_output_files
 use MHMCMC_module, only: MHMCMC, par_minstepsize, par_initstepsize
 use model_likelihood_module, only: model_likelihood, find_edc_initial_values

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
 integer :: solution_wanted, freq_print, freq_write, time1, time2, time3, i

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
 call read_options(solution_wanted,freq_print,freq_write,outfile)
 ! check whether this is a restart?
 call check_for_existing_output_files(PI%npars,MCO%outfile,MCO%stepfile,MCO%covfile,MCO%covifile)
 ! Initialise MCMC output, possibly a bit of a redundent subroutine...
 call initialise_mcmc_output

 write(*,*) "Running model version ", DATAin%ID

 ! Begin search for initial conditions
 write(*,*) "Beginning search for initial parameter conditions"
 ! Determine initial values, this requires using the MHMCMC
 call find_edc_initial_values
 ! Reset stepsize and covariance for main DRAM-MCMC
 PI%stepsize = 1d0 ; PI%beta_stepsize = par_initstepsize !par_minstepsize
 PI%parstd = 1d0 ; PI%Nparstd = 0d0 
 PI%covariance = 0d0 ; PI%mean_par = 0d0 
 PI%cov = .false. ; PI%use_multivariate = .false.
 do i = 1, PI%npars
    PI%covariance(i,i) = 1d0
 end do

 ! Restore module variables needed for the run - these components could be split
 ! into two subroutines to avoid double calling of file name creation
 ! components...
 call read_options(solution_wanted,freq_print,freq_write,outfile)

 ! Open the relevant output files
 call open_output_files(MCO%outfile,MCO%stepfile,MCO%covfile,MCO%covifile)

 if (restart_flag) then
     print*, "beginning restart simulation"
     ! now begin update of model timing variables and parameter values if this is a
     ! restart
     call update_for_restart_simulation
! else
!     ! write out first covariance matrix, this will be compared with the final covariance matrix
!     if (MCO%nWRITE > 0) call write_covariance_matrix
 endif

 ! update the user
 write(*,*)"Beginning MHMCMC for real..."
 ! call the main MCMC
 call MHMCMC(model_likelihood)
 write(*,*)"MHMCMC done now, moving on ..."

 ! tidy up by closing all files
 call close_output_files

end program cardamom_framework
