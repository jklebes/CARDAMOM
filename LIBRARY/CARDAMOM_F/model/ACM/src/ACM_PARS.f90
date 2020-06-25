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
    PI%parmin(1) = 01d0
    PI%parmax(1) = 40d0

    ! max temperature for photosynthesis (oC)
    PI%parmin(2) = 45d0
    PI%parmax(2) = 65d0
    ! optimum temperature for photosynthesis (oC)
    PI%parmin(3) = 20d0
    PI%parmax(3) = 40d0
    ! kurtosis of photosynthesis temperature response
    PI%parmin(4) = 0.10d0
    PI%parmax(4) = 0.25d0

    ! light limited photosynthesis

    ! maximum canopy quantum yield (gC/MJ_PAR/m2/day)
    ! SPA apparent canopy quantum yield 3.2
    ! Observational constraints = ~ 1

    PI%parmin(5) = 2d0   !7.19298-(0.9*7.19298)
    PI%parmax(5) = 5d0   !7.19298+(0.9*7.19298)

    !
    ! canopy conductance (gc) drivers
    !

    ! leafWP-soilWP (MPa); actual SPA parameter = 2+gplant
    PI%parmin(6) = -4d0 !-2.01d0 !-4.0
    PI%parmax(6) = -1d0 !-1.99d0 !-1.0

    !
    ! Longwave emitted by canopy
    !

    ! Maximum amount of LW emitted by canopy which escapes in one direction
    PI%parmin(7) = 0.25d0
    PI%parmax(7) = 0.75d0

    !
    ! GPP / transpiration optimisation
    !

    ! iWUE (gC/m2leaf/s/mmolH2Ogs)
    ! Actual value used in SPA is 8.4e-8 (Williams et al., 1996)
    ! Other reported values are 9.0e-08 -> 1.8e-07 (Bonan et al., 2014)
    PI%parmin(8) = 1d-9
    PI%parmax(8) = 1e-6

    !
    ! Soil shortwave radiation absorption
    !

    ! soil sw radiation absorption (fraction)
    PI%parmin(9) = 0.75d0
    PI%parmax(9) = 0.99d0

    !
    ! Canopy longwave escape from canopy
    !

    ! max fractional reduction of longwave release from canopy
    PI%parmin(10) = 0.750d0
    PI%parmax(10) = 0.999d0

    ! lai adjustment for long wave release from canopy
    PI%parmin(11) = 2.0d0
    PI%parmax(11) = 4.0d0

    !
    ! Linear correction for soil isothermal to net radiation
    !

    ! The magnitude of adjustment between soil isothermal net long wave radiation is found to have
    ! a linear relationship with absorbed short-wave radiation.
    ! NOTE: in constrast to the canopy relationship (positive) the soil one has two distinct
    ! negative relationships divided based on whether the soil is wet or not.
    ! This is not accounted for here

    ! Coefficient linking isothermal->net adjustment and absorbed SW
    ! SPA based prior is -0.015 (SE +/- 0.00077)
    PI%parmin(12) = -0.1d0
    PI%parmax(12) =  0.0d0

    ! Constant relating isothermal->net adjustment soil radiation (W/m2)
    ! SPA based prior is -1.84 W/m2 (SE +/- 0.078)
    PI%parmin(13) =-1.5d0
    PI%parmax(13) = 1.5d0

    !
    ! Max canopy intercepted PAR and NIR transmittance / reflectance
    !

    ! Max PAR transmittance
    ! SPA value = 0.16 (Sitka Spruce 0.16, grass/crop ~ 0.38)
    PI%parmin(14) = 0.10d0
    PI%parmax(14) = 0.40d0

    ! Max NIR transmittance
    ! SPA value = 0.26
    PI%parmin(15) = 0.10d0
    PI%parmax(15) = 0.40d0

    ! Max PAR reflectance
    ! SPA value = 0.16 (Sitka Spruce 0.07, grass/crop ~ 0.11)
    PI%parmin(16) = 0.10d0
    PI%parmax(16) = 0.60d0

    ! Max NIR reflectance
    ! SPA value = 0.43 (Sitka Spruce 0.16, grass/crop ~ 0.38)
    PI%parmin(17) = 0.10d0
    PI%parmax(17) = 0.60d0

    !
    ! Linear correction for canopy isothermal to net radiation
    !

    ! The magnitude of adjustment between canopy isothermal net long wave radiation is found to have
    ! a linear relationship (R2 = 0.74) with absorbed short-wave radiation.

    ! Coefficient linking isothermal->net adjustment and absorbed SW
    ! SPA based prior is 0.0688 (SE +/- 0.00045)
    PI%parmin(18) = 0.0d0
    PI%parmax(18) = 0.2d0

    ! Constant relating isothermal->net adjustment canopy radiation (W/m2)
    ! SPA based prior is 0.888 W/m2 (SE +/- 0.0316)
    PI%parmin(19) = 0d0
    PI%parmax(19) = 1.5d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
