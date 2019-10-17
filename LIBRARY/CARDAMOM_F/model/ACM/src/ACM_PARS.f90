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
    PI%parmin(1) = 03d0 
    PI%parmax(1) = 40d0 

    ! max temperature for photosynthesis (oC)
    PI%parmin(2) = 45d0
    PI%parmax(2) = 60d0
    ! optimum temperature for photosynthesis (oC)
    PI%parmin(3) = 20d0
    PI%parmax(3) = 40d0
    ! kurtosis of photosynthesis temperature response
    PI%parmin(4) = 0.10d0
    PI%parmax(4) = 0.25d0

    ! light limited photosynthesis

    ! maximum canopy quantum yield (gC/MJ_PAR/m2/day)
    PI%parmin(5) = 1d0  !7.19298-(0.9*7.19298)
    PI%parmax(5) = 7d0  !7.19298+(0.9*7.19298)

    !
    ! Canopy longwave radiation transmittance / reflectance coefficients
    !

    ! LAI at which LW transmittance to soil at 50 %
    PI%parmin(6) = 0.001d0 !0.01d0
    PI%parmax(6) = 0.2d0   !1.00d0
    ! LAI at which LW reflectance to sky at 50 %
    PI%parmin(7) = 0.001d0 !0.01d0
    PI%parmax(7) = 0.2d0   !2.5d0 

    !
    ! Canopy NIR shortwave radiation reflectance of intercepted radiation
    !

    ! maximum reflectance radiation (fraction)
    PI%parmin(8) = 0.10d0
    PI%parmax(8) = 0.50d0 
    ! LAI at which radiation absorption is at half saturation (m2/m2)
    ! SPA's radiative transfer scheme absorbs 50 % of shortwave at (lai == ~0.8)
    PI%parmin(9) = 0.01d0
    PI%parmax(9) = 2.00d0

    !
    ! canopy conductance (gc) drivers
    !

    ! leafWP-soilWP (MPa); actual SPA parameter = 2+gplant
    PI%parmin(10) = -4d0 !-2.01d0 !-4.0
    PI%parmax(10) = -1d0 !-1.99d0 !-1.0

    !
    ! Canopy PAR shortwave radiation reflectance of intercepted radiation
    !

    ! maximum reflectance
    PI%parmin(11) = 0.10d0
    PI%parmax(11) = 0.50d0 
    ! LAI at which radiation absorption is at half saturation (lai == ~ 0.8)
    PI%parmin(12) = 0.01d0
    PI%parmax(12) = 2.00d0

    !
    ! Longwave emitted by canopy
    !

    ! Maximum amount of LW emitted by canopy which escapes in one direction
    PI%parmin(13) = 0.25d0
    PI%parmax(13) = 0.75d0

    !
    ! GPP / transpiration optimisation
    !

    ! iWUE (gC/m2leaf/s/mmolH2Ogs)
    ! Actual value used in SPA is 8.4e-8 (Williams et al., 1996)
    ! Other reported values are 9.0e-08 -> 1.8e-07 (Bonan et al., 2014)
    PI%parmin(14) = 1d-9  
    PI%parmax(14) = 1e-5

    !
    ! Soil shortwave radiation absorption
    !

    ! soil sw radiation absorption (fraction)
    PI%parmin(15) = 0.75d0
    PI%parmax(15) = 0.99d0

    !
    ! Canopy shortwave transmittance of intercepted radiation
    !

    ! max reduction in par transmitted by canopy (fraction)
    PI%parmin(16) = -0.04d0 
    PI%parmax(16) = -0.01d0

    ! Max reduction in nir transmitted by canopy (fraction)
    PI%parmin(17) = -0.04d0 
    PI%parmax(17) = -0.01d0

    !
    ! Canopy longwave escape from canopy
    !

    ! max fractional reduction of longwave release from canopy 
    PI%parmin(18) = 0.750d0
    PI%parmax(18) = 0.999d0

    ! lai adjustment for long wave release from canopy 
    PI%parmin(19) = 2.0d0
    PI%parmax(19) = 4.0d0

    ! 
    ! Linear correction for soil isothermal to net radiation
    !

    ! The magnitude of adjustment between soil isothermal net long wave radiation is found to have
    ! a linear relationship with absorbed short-wave radiation. 

    ! Coefficient linking isothermal->net adjustment and absorbed SW
    ! SPA based prior is -1.015 (SE +/- 0.0008)
    PI%parmin(20) = -1.6d0
    PI%parmax(20) = -0.3d0

    ! Constant relating isothermal->net adjustment soil radiation (W/m2)
    ! SPA based prior is -1.84 (SE +/- 0.078)
    PI%parmin(21) =-2.5d0 
    PI%parmax(21) =-1.0d0

    ! 
    ! Max canopy intercepted PAR and NIR transmittance / reflectance
    !

    ! Max PAR transmittance
    PI%parmin(22) = 0.10d0
    PI%parmax(22) = 0.40d0

    ! Max NIR transmittance
    PI%parmin(23) = 0.10d0
    PI%parmax(23) = 0.40d0

    ! Max PAR reflectance
    PI%parmin(24) = 0.20d0
    PI%parmax(24) = 0.60d0

    ! Max NIR reflectance
    PI%parmin(25) = 0.20d0
    PI%parmax(25) = 0.60d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
