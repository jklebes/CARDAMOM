module MODEL_PARAMETERS

  implicit none

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine pars_info
    use MCMCOPT, only: PI

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or
    ! possibly should go into an alternate file which can be read in.
    ! This may
    ! improve the usability when it comes to reading these information
    ! in for
    ! different PFTs

    implicit none

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    ! Decomposition litter -> som (day-1)
    PI%parmin(1) = 0.00001d0
    PI%parmax(1) = 0.01d0

    ! Fraction of GPP respired as autotrophic
    PI%parmin(2) = 0.2d0
    PI%parmax(2) = 0.8d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3) = 0.01d0
    PI%parmax(3) = 0.5d0

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.01d0
    PI%parmax(4) = 1d0

    ! Leaf Lifespan (yr)
    ! Wright et al. (2004)
    ! 55 - 2922 days
    PI%parmin(5) = 0.15d0 
    PI%parmax(5) = 8d0

    ! TOR wood (2.7 - 109 years)
    PI%parmin(6) = 0.000025d0
    PI%parmax(6) = 0.001d0

    ! TOR roots (0.27 - 27 years)
    PI%parmin(7) = 0.0001d0
    PI%parmax(7) = 0.01d0

    ! TOR litter (0.27 - 27 years)
    PI%parmin(8) = 0.0001d0
    PI%parmax(8) = 0.01d0

    ! TOR som (2.7-2737 years)
    PI%parmin(9) = 0.000001d0
    PI%parmax(9) = 0.001d0

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10) = 0.018d0
    PI%parmax(10) = 0.08d0

    ! Canopy Efficiency (gC/m2leaf/day)
    PI%parmin(11) = 5d0
    PI%parmax(11) = 50d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(12) = 5d0
    PI%parmax(12) = 200d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C foliar
    PI%parmin(13) = 1d0 
    PI%parmax(13) = 2000d0

    ! C roots
    PI%parmin(14) = 20d0
    PI%parmax(14) = 2000d0

    ! C_wood
    PI%parmin(15) = 1d0
    PI%parmax(15) = 100000d0

    ! C litter
    PI%parmin(16) = 1d0
    PI%parmax(16) = 2000d0

    ! C_som
    PI%parmin(17) = 1d0
    PI%parmax(17) = 200000d0 

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
