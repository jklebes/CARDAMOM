
module MHMCMC_StressTests

  ! Module contains a number of diagnostic tests used to ensure that the MCMC
  ! is able to retrieve a known distribution of parameters for simple models.
  ! Created: 12/05/2021, T. L. Smallman (t.l.smallman@ed.ac.uk)
  ! Version history:
  ! Version 1: Estimating pi and radius of a circle

  implicit none

  ! Assume all contents private unless explicitly states
  private

  ! Explicit statement of public variables or functions
  public :: prepare_for_stress_test, StressTest_likelihood

  ! Declare any module level variables

  ! Stress Test 1 - estimate parameters for a circle
  ! Parameter 1 = pi, parameter 2 = radius
  double precision, parameter :: circle_par_1 = 3.141d0, &
                                 circle_par_2 = 2d0, &
                                 circle_obs_mean = circle_par_1 * circle_par_2 ** 2d0, &
                                 circle_obs_unc  = 0.2d0

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine circle(pars,nopars,area)

    implicit none

    ! Subroutine estimates the area of a circle using two parameters
    ! which define pi and the radius.

    ! arguments
    integer, intent(in) :: nopars
    double precision, dimension(nopars), intent(in) :: pars
    double precision, intent(out) :: area

    ! Determine the area of the circle for the current parameters
!    area = pars(1) * pars(2) ** 2d0
    area = pars(1) * circle_par_2 ** 2d0

  end subroutine circle
  !
  !--------------------------------------------------------------------
  !
  subroutine circle_parameter_prior_ranges
    use MCMCOPT, only: PI

    ! define the parameter prior ranges

    implicit none

    ! Pi
    PI%parmin(1) = 1d0
    PI%parmax(1) = 5d0

    ! Radius
!    PI%parmin(2) =  0.1d0
!    PI%parmax(2) = 10.0d0

  end subroutine circle_parameter_prior_ranges
  !
  !--------------------------------------------------------------------
  !
  subroutine prepare_for_stress_test(infile,outfile)
    use MCMCOPT, only: PI, MCO
    use cardamom_structures, only: DATAin

    ! Function by-passes the main CARDAMOM i/o code to allow
    ! for a non-standard operation of the model stress test

    implicit none

    ! Arguments
    character(350), intent(inout) :: infile, outfile

    ! local variables
    integer :: i

    ! Set internal parameters in the absence of an input file
    ! allocate the default run information

    if (outfile == "Circle") then
        ! ID = -1 StressTest - Circle
        DATAin%ID = -1
        DATAin%nodays = 1
        DATAin%nomet = 1
        DATAin%noobs = 1
        DATAin%nopools = 1
        DATAin%nopars = 1!2
        DATAin%nofluxes = 1
    else
        print*,"Valid Stress Test has not been specified"
        stop
    end if

    ! Now we have used the infile to determine that this is going to be stress test,
    ! and the specific one has been determined from the outfile,
    ! we will now overwrite the outfile to give a default output location
    outfile = "stress_test_output_"

    ! need to allocate memory to the model output variables
    allocate(DATAin%M_FLUXES(DATAin%nodays,DATAin%nofluxes)&
            ,DATAin%M_POOLS((DATAin%nodays+1),DATAin%nopools))

    ! alert the user
    write(*,*)"Created fields for model output"

    ! Begin allocating parameter info
    PI%npars = DATAin%nopars
    allocate(PI%parmin(PI%npars),PI%parmax(PI%npars),PI%parini(PI%npars) &
            ,PI%parfix(PI%npars),PI%parvar(PI%npars),PI%paradj(PI%npars) &
            ,PI%covariance(PI%npars,PI%npars),PI%mean_par(PI%npars) &
            ,PI%iC(PI%npars,PI%npars))

    ! force zero
    PI%parmin = 0d0 ; PI%parmax = 0d0 ; PI%parini = 0d0
    PI%parfix = 0d0 ; PI%parvar = 0d0 ; PI%paradj = 0d0
    PI%covariance = 0d0 ; PI%iC = 0d0

    ! load parameter max/min information
    if (DATAin%ID == -1) then
        call circle_parameter_prior_ranges
    end if

    ! For log-normalisation procedure, no parameter can be <=0.
    ! To facilitate easy of setting parameter ranges to real values
    ! we here instead calculate the adjustment need to ensure positive only values
    where (PI%parmin <= 0d0) PI%paradj = abs(PI%parmin) + 1d0

    ! defining initial MHMCMC stepsize and standard deviation
    PI%parvar = 1d0 ; PI%Nparvar = 0d0
    ! Covariance matrix cannot be set to zero therefore set initial value to a
    ! small positive value along to variance access
    PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false. ; PI%use_multivariate = .false.
    do i = 1, PI%npars
       PI%covariance(i,i) = 1d0
    end do

  end subroutine prepare_for_stress_test
  !
  !------------------------------------------------------------------
  !
  subroutine StressTest_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood

    ! local variables
    double precision :: output
    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0

    ! run the circle model
    call circle(PARS,DATAin%nopars,output)
    ! Estimate the likelihood score
    ML_obs_out = -0.5d0 * (((output - circle_obs_mean) / circle_obs_unc) ** 2)

  end subroutine StressTest_likelihood
  !
  !--------------------------------------------------------------------
  !
end module ! MHMCMC_StressTests
