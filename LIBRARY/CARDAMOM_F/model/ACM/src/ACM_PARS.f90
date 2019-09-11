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
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

!    PI%npars=12 ! dont forget to change in cardamom_io.f90

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    !
    ! Metabolic photosynthesis
    !

    ! Nitrogen use efficiency (gC/gN per m2 at optimum temperature)
    ! Derived from Vcmax reported in Wullschleger (1993), Journal of
    ! Experimental Botany, Vol 44, No. 262, pp. 907-920.
    ! More constrain prior range of 3-30 gC/gN/m2/day with a mean of 11.197440 +/-1.32313 Kattge et al., (2011)
    PI%parmin(1)=03d0 
    PI%parmax(1)=40d0 

    ! max temperature for photosynthesis (oC)
    PI%parmin(2)=45d0
    PI%parmax(2)=60d0
    ! optimum temperature for photosynthesis (oC)
    PI%parmin(3)=20d0
    PI%parmax(3)=40d0
    ! kurtosis of photosynthesis temperature response
    PI%parmin(4)=0.10d0
    PI%parmax(4)=0.30d0

    ! light limited photosynthesis

    ! maximum canopy quantum yield (gC/MJ_PAR/m2/day)
    PI%parmin(5)=1d0  !7.19298-(0.9*7.19298)
    PI%parmax(5)=7d0  !7.19298+(0.9*7.19298)

    !
    ! Canopy longwave radiation transmittance
    !

    ! max reduction in longwave transmitted by canopy (fraction)
    PI%parmin(6)=0.5d0
    PI%parmax(6)=1d0
    ! LAI at which 50 % of reduction has occured (m2/m2)
    PI%parmin(7)=0.01d0
    PI%parmax(7)=2.5d0

    !
    ! Canopy NIR shortwave radiation reflectance
    !

    ! maximum reflectance radiation (fraction)
    PI%parmin(8)=0.10d0
    PI%parmax(8)=0.99d0
    ! LAI at which radiation absorption is at half saturation (m2/m2)
    ! SPA's radiative transfer scheme absorbs 50 % of shortwave at lai == ~1.65
    PI%parmin(9)=0.1d0
    PI%parmax(9)=2.0d0

    !
    ! canopy conductance (gc) drivers
    !

    ! leafWP-soilWP (MPa); actual SPA parameter = 2+gplant
    PI%parmin(10)=-2.01d0 !-4.0
    PI%parmax(10)=-1.99d0 !-1.0

    !
    ! Canopy PAR shortwave radiation reflectance
    !

    ! maximum reflectance
    PI%parmin(11)=0.10d0
    PI%parmax(11)=0.99d0
    ! LAI at which radiation absorption is at half saturation
    PI%parmin(12)=0.1d0
    PI%parmax(12)=2.0d0

    !
    ! Longwave reflectance back to sky
    !

    ! lai at which longwave returned to sky by canopy is half saturation
    PI%parmin(13)=0.75d0
    PI%parmax(13)=1.5d0

    !
    ! GPP / transpiration optimisation
    !

    ! iWUE (gC/m2leaf/s/mmolH2Ogs)
    ! Actual value used in SPA is 8.4e-8 (Williams et al., 1996)
    ! Other reported values are 9.0e-08 -> 1.8e-07 (Bonan et al., 2014)
    PI%parmin(14)=1d-10  
    PI%parmax(14)=0.0001d0

    !
    ! Soil shortwave radiation absorption
    !

    ! soil sw radiation absorption (fraction)
    PI%parmin(15)=0.85d0
    PI%parmax(15)=0.99d0

    !
    ! Canopy shortwave transmittance
    !

    ! max reduction in par transmitted by canopy (fraction)
    PI%parmin(16)=0.5d0
    PI%parmax(16)=1d0
    ! LAI at which 50 % occurs
    PI%parmin(17)=0.1d0
    PI%parmax(17)=2.0d0

    ! max reduction in nir transmitted by canopy (fraction)
    PI%parmin(18)=0.5d0
    PI%parmax(18)=1d0
    ! LAI at which 50 % occurs
    PI%parmin(19)=0.1d0
    PI%parmax(19)=2.0d0

    !
    ! Canopy longwave escape from canopy
    !

    ! max fraction of longwave release from canopy (0.20:0.30)
    PI%parmin(20)=0.01d0
    PI%parmax(20)=0.99d0

    ! lai adjustment for long wave release from canopy (1:3)
    PI%parmin(21)=0.01d0
    PI%parmax(21)=3.0d0

    ! max fraction longwave reflected back to sky by canopy
    PI%parmin(22)=0.01d0
    PI%parmax(22)=1d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
