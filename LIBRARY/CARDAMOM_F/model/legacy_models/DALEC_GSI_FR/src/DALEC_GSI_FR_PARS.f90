module MODEL_PARAMETERS

  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

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
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or
    ! possibly should go into an alternate file which can be read in.
    ! This may
    ! improve the usability when it comes to reading these information
    ! in for
    ! different PFTs

    implicit none

!    PI%npars=29 ! dont forget to change in cardamom_io.f90

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    ! Decomposition rate
    PI%parmin(1)=0.00001
    PI%parmax(1)=0.01

    ! Fraction of GPP respired
    PI%parmin(2)=0.3
    PI%parmax(2)=0.7

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3)=0.01
    PI%parmax(3)=0.5

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4)=0.01
    PI%parmax(4)=1.

    ! GSI max leaf turnover
    PI%parmin(5)=0.000001
    PI%parmax(5)=0.2

    ! TOR wood* - 1% loss per year value
    PI%parmin(6)=0.00001
    PI%parmax(6)=0.001

    ! TOR roots
    PI%parmin(7)=0.0001
    PI%parmax(7)=0.01

    ! TOR litter
    PI%parmin(8)=0.0001
    PI%parmax(8)=0.01

    ! TOR SOM
    PI%parmin(9)=0.0000001
    PI%parmax(9)=0.001

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10)=0.018
    PI%parmax(10)=0.08

    ! log10 avg foliar N (gN.m-2) ! Canopy Efficiency
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=10d0 ! -0.50 ! 10d0
    PI%parmax(11)=100d0! 1d0  ! 100d0

    ! GSI max labile turnover
    PI%parmin(12)=0.000001
    PI%parmax(12)=0.2

    ! Fraction to Clab*/
    PI%parmin(13)=0.01
    PI%parmax(13)=0.5

    ! GSI min temperature threshold (oC)
    PI%parmin(14)=225d0
    PI%parmax(14)=330d0

    ! GSI max temperature threshold (oC)
    PI%parmin(15)=225d0
    PI%parmax(15)=330d0

    ! GSI min photoperiod threshold (sec)
    PI%parmin(16)=21600d0 ! 6 hours
    PI%parmax(16)=64800d0 ! 18 hours

    ! LMA
    ! Kattge et al. 2011,
    PI%parmin(17)=10d0 ! 10.
    PI%parmax(17)=200d0 ! 400.

    ! GSI max photoperiod threshold (sec)
    PI%parmin(24)=21600d0 ! 6 hours
    PI%parmax(24)=64800d0 ! 18 hours

    ! GSI min VPD threshold (Pa)
    PI%parmin(25)=1d0
    PI%parmax(25)=5500d0

    ! GSI max VPD threshold (Pa)
    PI%parmin(26)=1d0
    PI%parmax(26)=5500d0

    ! critical GPP for LAI increase (gC.m-2.day-1)
    PI%parmin(27)=1e-8
    PI%parmax(27)=0.5

    ! fraction of Cwood which is branch
    PI%parmin(28)=0.05
    PI%parmax(28)=0.65

    ! fraction of Cwood which is coarse root
    PI%parmin(29)=0.15
    PI%parmax(29)=0.45


    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile
    PI%parmin(18)=1d0
    PI%parmax(18)=1000d0

    ! C foliar
    PI%parmin(19)=1d0
    PI%parmax(19)=1000d0

    ! C roots
    PI%parmin(20)=1d0
    PI%parmax(20)=1000d0

    ! C_wood
    PI%parmin(21)=1d0
    PI%parmax(21)=20000d0

    ! C litter
    PI%parmin(22)=1d0
    PI%parmax(22)=10000d0

    ! C_som
    PI%parmin(23)=100d0
    PI%parmax(23)=200000d0

    !
    ! Replanting pools values
    !

    ! C labile
    PI%parmin(30)=1.0
    PI%parmax(30)=1000.0

    ! C foliar
    PI%parmin(31)=1.0
    PI%parmax(31)=1000.0

    ! C roots
    PI%parmin(32)=1.0
    PI%parmax(32)=1000.0

    ! C_wood
    PI%parmin(33)=1.0
    PI%parmax(33)=1000.0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
