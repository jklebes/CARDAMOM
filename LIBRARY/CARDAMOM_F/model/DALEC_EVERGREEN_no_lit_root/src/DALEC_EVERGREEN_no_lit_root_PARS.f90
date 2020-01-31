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

    ! Fraction of GPP respired as autotrophic 
    PI%parmin(1) = 0.20d0
    PI%parmax(1) = 0.80d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(2) = 0.01d0
    PI%parmax(2) = 0.5d0

    ! Leaf Lifespan (yr) 
    ! Wright et al. (2004) 55 - 2922 days
    PI%parmin(3) = 0.15d0 
    PI%parmax(3) = 8d0

    ! TOR wood + root 
    PI%parmin(4) = 0.000009d0 ! 304  years
    PI%parmax(4) = 0.01d0     ! 0.27 years

    ! TOR litter + som
    PI%parmin(5) = 2.737851d-06 ! 1000 years 
    PI%parmax(5) = 0.01d0       ! 0.13 years

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(6) = 0.018d0
    PI%parmax(6) = 0.06d0

    ! Canopy Efficiency (gC/m2leaf/day)
    PI%parmin(7) = 5d0
    PI%parmax(7) = 50d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(8) = 5d0
    PI%parmax(8) = 200d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C foliar
    PI%parmin(9) = 1d0 
    PI%parmax(9) = 2000d0

    ! C root + wood
    PI%parmin(10) = 1d0
    PI%parmax(10) = 100000d0

    ! C_som
    PI%parmin(11) = 1d0
    PI%parmax(11) = 200000d0 

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
