!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to
! assimilate observations and ecological theory to retrieve parameters for the
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently.
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk),
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains all subroutine and functions relevant to determining the log-likelihood
! of DALEC.A1.C1.D2.F2.H1.P1 as a function of observations and ecological dynamical constraints.
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module model_likelihood_module
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

  ! which to make open
  public :: model_likelihood, scaled_model_likelihood, find_edc_initial_values

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: EDC
    integer :: DIAG
    integer :: PASSFAIL(150) ! allow space for 150 possible checks
    integer :: nedc ! number of edcs being assessed
  end type
  type (EDCDIAGNOSTICS), save :: EDCD

  ! Has the model sanity check been conducted yet?
  logical :: sanity_check = .false.

  contains
  !
  !------------------------------------------------------------------
  !


  subroutine assess_EDC1(PARS, npars, meantemp, meanrad, EDC1)

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (Bloom & Williams 2015).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    double precision :: torfol ! yearly leaf loss fraction

    ! set initial value
    EDC1 = 1
    DIAG = EDCD%DIAG

    ! estimate GPP allocation fractions
    fauto = pars(2)
    ffol = (1d0-fauto)*pars(3)
    flab = (1d0-fauto-ffol)*pars(13)
    froot = (1d0-fauto-ffol-flab)*pars(4)
    fwood = 1d0-fauto-ffol-flab-froot
    fsom = fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(8))

    ! yearly leaf loss fraction
    torfol = 1d0/(pars(5)*365.25d0)

    ! set all EDCs to 1 (pass)
    EDCD%nedc = 150
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(8) < pars(9))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    endif

    ! litter2som greater than som to atm rate
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(1) < pars(9))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    endif

    ! turnover of foliage faster than turnover of wood
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > torfol) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    end if

    ! root turnover greater than som turnover at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(7) < (pars(9)*exp(pars(10)*meantemp)))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    endif

    ! Initial leaf area index should not be larger than ~10 m2/m2
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(19)/pars(17)) > 10d0) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    endif

    !! GPP allocation to foliage and labile cannot be 5 orders of magnitude
    !! difference from GPP allocation to roots
    !if ((EDC1 == 1 .or. DIAG == 1) .and. ((ffol+flab) > (5d0*froot) .or. ((ffol+flab)*5d0) < froot)) then
    !    EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    !endif

    ! IMPLICIT Combustion completeness for foliage should be greater than soil
    ! IMPLICIT Combustion completeness for fol+root litter should be greater than soil

    ! Combustion completeness for foliage should be greater than non-photosynthetic tissues
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(28) < pars(29)) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
    endif
    ! Combustion completeness for non-photosynthetic tissue should be greater than soil
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(29) < pars(30)) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
    endif
    ! Combustion completeness for foliar + fine root litter should be greater foliage
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(31) < pars(28)) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
    endif

    ! could always add more / remove some

  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,nodiags,deltat,steps_per_year &
                        ,parmax,pars,met,M_POOLS,M_FLUXES,M_DIAGS &
                        ,meantemp,EDC2)

    use cardamom_structures, only: DATAin

    ! Determines whether the dynamical contraints for the search of the initial
    ! parameters has been successful or whether or not we should abandon the
    ! current set and move on

    implicit none

    ! declare input variables
    integer, intent(in) :: npars          & ! number of model parameters
                          ,nomet          & ! number of met drivers
                          ,nofluxes       & ! number of fluxes from model
                          ,nopools        & ! number of pools in model
                          ,nodays         & ! number of days in simulation
                          ,nodiags        & ! number of diagnostic variables
                          ,steps_per_year

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,M_DIAGS(nodays,nodiags)     & ! time varying diagnostics from current model simulation model
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    integer :: n, nn, nnn, DIAG, y, steps_per_month, nd, fl, fs, &
               io_start, io_finish
    double precision :: infi, tmp, tmp1, tmp2, &!, EQF, etol
                        jan_sd_lai, jan_mean_lai, jan_first_lai
    !double precision, dimension(nodays) :: tmp1, tmp2
    double precision, dimension(nopools) :: jan_mean_pools, jan_first_pools, &
                                            mean_pools, Fin, Fout, Rm, Rs, &
                                            Fin_yr1, Fout_yr1, Fin_yr2, Fout_yr2
    double precision, dimension(nofluxes) :: FT, FT_yr1, FT_yr2

    ! Steady State Attractor:
    ! Log ratio difference between inputs and outputs of the system.
    double precision, parameter :: EQF1_5 = log(1.5d0), & ! 10.0 = order magnitude; 2 = double and half
                                   EQF2 = log(2d0),   & ! 10.0 = order magnitude; 2 = double and half
                                   EQF5 = log(5d0),   &
                                   EQF10 = log(10d0), &
                                   EQF15 = log(15d0), &
                                   EQF20 = log(20d0), &
                                  C_etol = 0.10d0       ! 0.20d0 lots of AGB !0.10d0 global / site more data !0.05d0 global 1 or 2 AGB estimates

!    ! Debugging print statements
!    print*,"assess_EDC2: "

    ! update initial values
    DIAG = EDCD%DIAG
    EDC2 = 1
    infi = 0d0

    ! derive mean pools for first year
    do n = 1, nopools
       mean_pools(n) = cal_mean_pools(M_POOLS(1:steps_per_year,n),steps_per_year)
    end do

    ! number of time steps per month
    steps_per_month = ceiling(dble(steps_per_year) * 0.08333333d0)

    ! First calculate total flux for the simulation period
    io_start = (steps_per_year*2) + 1 ; io_finish = nodays
    if (DATAin%nos_years < 3) io_start = 1
    do fl = 1, nofluxes
!       FT(fl) = sum(M_FLUXES(1:nodays,fl)*deltat(1:nodays))
       FT(fl) = sum(M_FLUXES(io_start:io_finish,fl)*deltat(io_start:io_finish))
       FT_yr1(fl) = sum(M_FLUXES(1:steps_per_year,fl)*deltat(1:steps_per_year))
       !FT_yr2(fl) = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),fl) &
                       !*deltat((steps_per_year+1):(steps_per_year*2)))
    end do

    ! get total in and out for each pool
    ! labile
    Fin(1)  = FT(5)
    Fout(1) = FT(8)+FT(18)+FT(24)+FT(31)+FT(37)
    Fin_yr1(1)  = FT_yr1(5)
    Fout_yr1(1) = FT_yr1(8)+FT_yr1(18)+FT_yr1(24)+FT_yr1(31)+FT_yr1(37)
!    Fin_yr2(1)  = FT_yr2(5)
!    Fout_yr2(1) = FT_yr2(8)+FT_yr2(18)+FT_yr2(24)+FT_yr2(31)+FT_yr2(37)
    ! foliar
    Fin(2)  = FT(4)+FT(8)
    Fout(2) = FT(10)+FT(19)+FT(25)+FT(32)+FT(38)
    Fin_yr1(2)  = FT_yr1(4)+FT_yr1(8)
    Fout_yr1(2) = FT_yr1(10)+FT_yr1(19)+FT_yr1(25)+FT_yr1(32)+FT_yr1(38)
!    Fin_yr2(2)  = FT_yr2(4)+FT_yr2(8)
!    Fout_yr2(2) = FT_yr2(10)+FT_yr2(19)+FT_yr2(25)+FT_yr2(32)+FT_yr2(38)
    ! root
    Fin(3)  = FT(6)
    Fout(3) = FT(12)+FT(20)+FT(26)+FT(33)+FT(39)
    Fin_yr1(3)  = FT_yr1(6)
    Fout_yr1(3) = FT_yr1(12)+FT_yr1(20)+FT_yr1(26)+FT_yr1(33)+FT_yr1(39)
!    Fin_yr2(3)  = FT_yr2(6)
!    Fout_yr2(3) = FT_yr2(12)+FT_yr2(20)+FT_yr2(26)+FT_yr2(33)+FT_yr2(39)
    ! wood
    Fin(4)  = FT(7)
    Fout(4) = FT(11)+FT(21)+FT(27)+FT(34)+FT(40)
    Fin_yr1(4)  = FT_yr1(7)
    Fout_yr1(4) = FT_yr1(11)+FT_yr1(21)+FT_yr1(27)+FT_yr1(34)+FT_yr1(40)
!    Fin_yr2(4)  = FT_yr2(7)
!    Fout_yr2(4) = FT_yr2(11)+FT_yr2(21)+FT_yr2(27)+FT_yr2(34)+FT_yr2(40)
    ! litter
    Fin(5)  = FT(10)+FT(12)+FT(24)+FT(25)+FT(26)
    Fout(5) = FT(13)+FT(15)+FT(22)+FT(28)+FT(35)
    Fin_yr1(5)  = FT_yr1(10)+FT_yr1(12)+FT_yr1(24)+FT_yr1(25)+FT_yr1(26)
    Fout_yr1(5) = FT_yr1(13)+FT_yr1(15)+FT_yr1(22)+FT_yr1(28)+FT_yr1(35)
!    Fin_yr2(5)  = FT_yr2(10)+FT_yr2(12)+FT_yr2(24)+FT_yr2(25)+FT_yr2(26)
!    Fout_yr2(5) = FT_yr2(13)+FT_yr2(15)+FT_yr2(22)+FT_yr2(28)+FT_yr2(35)
    ! som
    Fin(6)  = FT(11)+FT(15)+FT(27)+FT(28)
    Fout(6) = FT(14)+FT(23)+FT(36)
    Fin_yr1(6)  = FT_yr1(11)+FT_yr1(15)+FT_yr1(27)+FT_yr1(28)
    Fout_yr1(6) = FT_yr1(14)+FT_yr1(23)+FT_yr1(36)
!    Fin_yr2(6)  = FT_yr2(11)+FT_yr2(15)+FT_yr2(27)+FT_yr2(28)
!    Fout_yr2(6) = FT_yr2(14)+FT_yr2(23)+FT_yr2(36)

    !
    ! Begin EDCs here
    !

!    ! ensure ratio between Cfoliar and Croot is less than 5
!    if ((EDC2 == 1 .or. DIAG == 1) .and. &
!        (mean_pools(2) > (mean_pools(3)*5d0) .or. (mean_pools(2)*5d0) < mean_pools(3)) ) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(9) = 0
!    end if

    if (DATAin%nos_years > 1) then
        ! Determine the mean and standard deviation of January LAIs
        jan_sd_lai = 0d0 ; jan_mean_lai = 0d0 ; jan_first_lai = 0d0 ! reset
        jan_first_lai = M_DIAGS(1,1) ! First January LAI
        ! Initially sum each January from each year
        do y = 1, DATAin%nos_years
           nn = 1 + (steps_per_year * (y - 1))
           jan_mean_lai = jan_mean_lai + M_DIAGS(nn,1)
        end do
        ! Calculate the mean
        jan_mean_lai = jan_mean_lai / dble(DATAin%nos_years)
        ! Calculate the standard deviation now
        do y = 1, DATAin%nos_years
           nn = 1 + (steps_per_year * (y - 1))
           jan_sd_lai = jan_sd_lai + (jan_mean_lai - M_DIAGS(nn,1))**2d0
        end do
        jan_sd_lai = sqrt(jan_sd_lai / (dble(DATAin%nos_years - 1)))
        if ((EDC2 == 1 .or. DIAG == 1) .and. &
            abs(jan_first_lai-jan_mean_lai) > (jan_sd_lai*2d0) .and. abs(jan_first_lai-jan_mean_lai) > 0.01d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(9) = 0
        end if
    end if ! nos_years > 1

    ! EDC just for DALEC_CDEA_ACM2_BUCKET due to complications linked to
    ! the empirical phenology but mechanistic hydrology / photosynthesis
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_DIAGS(1:nodays,1)) > 10d0 ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(11) = 0
    end if

    ! Equilibrium factor (in comparison with initial conditions)
!    EQF = 10d0 ! TLS 06/11/2019 !10d0 ! JFE replaced 10 by 2 - 27/06/2018
    ! Pool exponential decay tolerance
!    etol = 0.3d0 !0.1d0

    ! first calculate total flux for the whole simulation period
!    do fl = 1, nofluxes
!        FT(fl) = 0
!        do nd = 1, nodays
!            FT(fl) = FT(fl) + M_FLUXES(nd,fl)*deltat(nd)
!        end do
!    end do

    ! Iterate through C pools to determine whether they have their ratio of
    ! input and outputs are outside of steady state approximation.
    ! See Bloom et al., 2016 PNAS for details

!    ! iterate to check whether Fin/Fout is within EQF limits
!    Rm = Fin/Fout
!    Rs = Rm * (jan_mean_pools / jan_first_pools)
!    do n = 1, nopools-1
!       ! Restrict rates of increase
!       if ((EDC2 == 1 .or. DIAG == 1) .and. abs(log(Rm(n))) > log(EQF10)) then
!           EDC2 = 0d0 ; EDCD%PASSFAIL(13+n-1) = 0
!       end if
!       ! Restrict exponential decay
!       if ((EDC2 == 1 .or. DIAG == 1) .and. abs(Rs(n)-Rm(n)) > 0.1d0) then
!           EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!       end if
!    end do

     ! What are in effect the potential growth rates are modulated by the current
     ! fixed temperature sub-model used in the model. This means that the parameterised
     ! potential rates might never be achievable even if plausible. Thus the maximum
     ! parameter bound for the potential growth rates need to be increased. These EDCs
     ! prevent an emergent growth rate that is unrealistic. Here we assume that tissue
     ! growth for foliage, wood and roots cannot be greater than 10 gC/m2/day
     if ((EDC2 == 1 .or. DIAG == 1)) then
         ! Foliage
         if (maxval(M_FLUXES(:,4) + M_FLUXES(:,8)) > 10d0) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(14) = 0
         end if
         ! Fine roots
         if (maxval(M_FLUXES(:,6)) > 10d0) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(15) = 0
         end if
         ! Wood
         if (maxval(M_FLUXES(:,7)) > 10d0) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(16) = 0
         end if
     end if

    ! Average growth rates for foliage and fine roots cannot be 5 orders of magnitude different
    if ((EDC2 == 1 .or. DIAG == 1) .and. (FT(4)+FT(8)) > (5d0*FT(6))) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(17) = 0
    endif
    ! Average growth rates for foliage and fine roots cannot be 5 orders of magnitude different
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((FT(4)+FT(8))*5d0) < FT(6)) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(18) = 0
    endif

    if (EDC2 == 1 .or. DIAG == 1) then

!        ! Living pools
!        do n = 1, 3
!           ! Restrict mean rates of increase
!           if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!               EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!           end if
!           ! Restrict rates from deviating unrealistically from the mean
!           if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
!                     abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
!               EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
!           end if
!        end do
!        ! Foliage pool, note that in CDEA EDCs Fin has already been multiplied by time step
!        n = 2
!        ! Restrict mean rates of increase
!        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
!        ! Restrict exponential behaviour at initialisation
!        if (abs(abs(log(Fin_yr1(n)/Fout_yr1(n))) - abs(log(Fin_yr2(n)/Fout_yr2(n)))) > C_etol) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
        ! Fine root pool, note that in CDEA EDCs Fin has already been multiplied by time step
        n = 3
!        ! Restrict mean rates of increase
!        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
        ! Restrict rates from deviating unrealistically from the mean
        if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                  abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
        end if
!        ! Wood pool hack, note that in CDEA EDCs Fin has already been multiplied by time step
!        n = 4
!        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
!        if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
!                  abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
!        end if
        ! Dead pools
        do n = 5, 6
           ! Restrict rates of increase
           if (abs(log(Fin(n)/Fout(n))) > EQF1_5) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
           end if
           ! Restrict rates from deviating unrealistically from the mean
           if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                     abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
           end if
        end do

    end if ! EDC2 == 1 .or. DIAG == 1

    ! Ensure that the mean transit time of foliage and the LCA are consistent with the
    ! leaf economic spectrum (LES).
    ! LL (months) ~ LMA (gm2) R2 = 0.42 from
    ! Wright et al., (2004), doi: https://doi.org/10.1038/nature02403
    ! Onoda et al., (2017), doi: https://doi.org/10.1111/nph.14496
    if (EDC2 == 1 .or. DIAG == 1) then
        ! Assume that the MTT(nat,fire) foliage should be within the uncertainty bounds of the LES
        ! Mean equation LL(months) = 0.0031 * LMA**1.71, coefficient 95CI = 1.62,1.82
        ! Estimating the MTT, converting from days to years using 1/365.25 = 0.002737851
        ! 0.08333333 converts months to years for the LES equation.
        tmp = sum(M_POOLS(:,2)) / dble(nodays)
        tmp1 = sum(M_FLUXES(:,10)+M_FLUXES(:,19)+M_FLUXES(:,25)) / dble(nodays)
        tmp = (tmp / tmp1) * 0.002737851d0
        ! determine the lower and upper bound of the LES .
        ! not for the upper bound, do not allow a value less than 1 years
        tmp1 = 0.08333333d0*(0.0031d0*(pars(17)*2.083333d0)**1.62d0)
        tmp2 = max(1d0,0.08333333d0*(0.0031d0*(pars(17)*2.083333d0)**1.82d0))
        if (tmp < tmp1) then
            ! The current leaf lifespan is shorter than expected
            EDC2 = 0d0 ; EDCD%PASSFAIL(45) = 0
        endif
        if (tmp > tmp2) then
            ! The current leaf life span is longer than expected
            EDC2 = 0d0 ; EDCD%PASSFAIL(46) = 0
        endif
    endif ! EDC2 == 1 .or. DIAG == 1

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 55 - 61 of the PASSFAIL array

    ! ensure minimum pool values are >= 0, /= NaN or Inf
    if (EDC2 == 1 .or. DIAG == 1) then

       do n = 1, nopools
          if (minval(M_POOLS(1:nodays,n)) < 0d0 .or. &
              maxval(abs(M_POOLS(1:nodays,n))) == abs(log(infi)) .or. &
              minval(M_POOLS(1:nodays,n)) /= minval(M_POOLS(1:nodays,n))) then
              EDC2 = 0d0 ; EDCD%PASSFAIL(55+n) = 0
          endif
       end do

       do n = 1, nofluxes
          if (maxval(abs(M_FLUXES(:,n))) == abs(log(infi)) .or. &
              minval(M_FLUXES(:,n)) /= minval(M_FLUXES(:,n))) then
              EDC2 = 0d0 ; EDCD%PASSFAIL(55+nopools+n) = 0
          endif
       end do

    end if ! min pool assessment

!    ! Debugging print statements
!    print*,"assess_EDC2: done"

  end subroutine assess_EDC2



  !
  !------------------------------------------------------------------
  !
  subroutine calc_obs_likelihoods(ML_obs_out)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: sw_par_fraction, top_soil_depth

    ! Subroutine to control the calculation of the observation related
    ! log-likelihoods and accounting for the various scalings

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    ! local variable
    double precision, dimension(DATAin%nodays) :: mod

    !
    ! Do diagnostics (DIAGS)
    !

    ! Calculate log-likelihood for fraction of absorbed radiation
    if (DATAin%nfAPAR > 0) then
        mod = DATAin%M_DIAGS(1:DATAin%nodays,3) & ! APAR
            / (DATAin%met(4,1:DATAin%nodays)*sw_par_fraction)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nfAPAR,DATAin%fAPARpts,DATAin%fAPAR,DATAin%fAPAR_unc,DATAin%fAPAR_lag, &
                                             1d0,mod)
    endif ! nfAPAR > 0
    ! Calculate log-likelihood for leaf area index
    if (DATAin%nlai > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nlai,DATAin%laipts,DATAin%LAI,DATAin%LAI_unc,DATAin%LAI_lag, &
                                             1d0,DATAin%M_DIAGS(1:DATAin%nodays,1))
    end if

    !
    ! Do pools (POOLS)
    !

    ! Calculate log-likelihood for foliage stocks
    if (DATAin%nCfol_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCfol_stock,DATAin%Cfol_stockpts, &
                                             DATAin%Cfol_stock,DATAin%Cfol_stock_unc,DATAin%Cfol_stock_lag, &
                                             1d0,DATAin%M_POOLS(1:DATAin%nodays,2))
    endif ! nCfol_stock > 0
    ! Calculate log-likelihood for fine root stocks
    if (DATAin%nCroots_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCroots_stock,DATAin%Croots_stockpts, &
                                             DATAin%Croots_stock,DATAin%Croots_stock_unc,DATAin%Croots_stock_lag, &
                                             1d0,DATAin%M_POOLS(1:DATAin%nodays,3))
    endif ! nCroots_stock > 0
    ! Calculate log-likelihood for total wood stocks
    if (DATAin%nCwood_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_stock,DATAin%Cwood_stockpts, &
                                             DATAin%Cwood_stock,DATAin%Cwood_stock_unc,DATAin%Cwood_stock_lag, &
                                             1d0,DATAin%M_POOLS(1:DATAin%nodays,4))
    endif ! nCwood_stock > 0
    ! Calculate log-likelihood for foliage litter stocks
    if (DATAin%nClit_stock > 0) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = (sum(DATAin%M_FLUXES(1:DATAin%nodays,10))/sum(DATAin%M_FLUXES(1:DATAin%nodays,10)+DATAin%M_FLUXES(1:DATAin%nodays,12))) &
                 * DATAin%M_POOLS(1:DATAin%nodays,5)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nClit_stock,DATAin%Clit_stockpts, &
                                             DATAin%Clit_stock,DATAin%Clit_stock_unc,DATAin%Clit_stock_lag, &
                                             1d0,mod)
    endif ! nClit_stock > 0
    ! Calculate log-likelihood for soil organic matter stocks
    if (DATAin%nCsom_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCsom_stock,DATAin%Csom_stockpts, &
                                             DATAin%Csom_stock,DATAin%Csom_stock_unc,DATAin%Csom_stock_lag, &
                                             1d0,DATAin%M_POOLS(1:DATAin%nodays,6))
    endif ! nCsom_stock > 0

    !
    ! Do fluxes (FLUXES)
    !

    ! Calculate log-likelihood for evapotranspiration
    if (DATAin%nEvap > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nEvap,DATAin%Evappts,DATAin%Evap,DATAin%Evap_unc,DATAin%Evap_lag, &
                                             1d0,DATAin%M_FLUXES(1:DATAin%nodays,29))
    endif ! nEvap > 0
    ! Calculate log-likelihood for fire
    if (DATAin%nFire > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nFire,DATAin%Firepts,DATAin%Fire,DATAin%Fire_unc,DATAin%Fire_lag, &
                                             1d0,DATAin%M_FLUXES(1:DATAin%nodays,17))
    endif ! nFire > 0
    ! Calculate log-likelihood for gross primary production
    if (DATAin%ngpp > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%ngpp,DATAin%gpppts,DATAin%GPP,DATAin%GPP_unc,DATAin%GPP_lag, &
                                             1d0,DATAin%M_FLUXES(1:DATAin%nodays,1))
    endif ! ngpp > 0
    ! Calculate log-likelihood for harvest
    if (DATAin%nharvest > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nharvest,DATAin%harvestpts,DATAin%harvest,DATAin%harvest_unc,DATAin%harvest_lag, &
                                             1d0,DATAin%M_FLUXES(1:DATAin%nodays,30))
    endif ! nharvest > 0
    ! Calculate log-likelihood for net biome productivity
    if (DATAin%nnbe > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + DATAin%M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + DATAin%M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            + DATAin%M_FLUXES(1:DATAin%nodays,17) & ! Fire
            - DATAin%M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnbe,DATAin%nbepts,DATAin%NBE,DATAin%NBE_unc,DATAin%NBE_lag, &
                                             1d0,mod)
    endif ! nnbe > 0
    ! Calculate log-likelihood for net ecosystem exchange of CO2
    if (DATAin%nnee > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + DATAin%M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + DATAin%M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            - DATAin%M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnee,DATAin%neepts,DATAin%NEE,DATAin%NEE_unc,DATAin%NEE_lag, &
                                             1d0,mod)
    endif ! nnee > 0
    ! Calculate log-likelihood for ecosystem respiration
    if (DATAin%nreco > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + DATAin%M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + DATAin%M_FLUXES(1:DATAin%nodays,14)   ! Rhet som
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nreco,DATAin%recopts,DATAin%Reco,DATAin%Reco_unc,DATAin%Reco_lag, &
                                             1d0,mod)
    endif ! nreco > 0
    ! Calculate log-likelihood for total wood net increment
    if (DATAin%nCwood_inc > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,7) - DATAin%M_FLUXES(1:DATAin%nodays,11)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_inc,DATAin%Cwood_incpts, &
                                             DATAin%Cwood_inc,DATAin%Cwood_inc_unc,DATAin%Cwood_inc_lag, &
                                             1d0,mod)
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood gross increment
    if (DATAin%nCwood_growth > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_growth,DATAin%Cwood_growthpts, &
                                             DATAin%Cwood_growth,DATAin%Cwood_growth_unc,DATAin%Cwood_growth_lag, &
                                             1d0,DATAin%M_FLUXES(1:DATAin%nodays,7))
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood mortality
    if (DATAin%nCwood_mortality > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_mortality,DATAin%Cwood_mortalitypts, &
                                             DATAin%Cwood_mortality,DATAin%Cwood_mortality_unc,DATAin%Cwood_mortality_lag, &
                                             1d0,DATAin%M_FLUXES(1:DATAin%nodays,11))
    endif ! nCwood_mortality > 0

    return

  end subroutine calc_obs_likelihoods
  !
  !------------------------------------------------------------------
  !
  subroutine calc_scaled_obs_likelihoods(ML_obs_out)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: sw_par_fraction, top_soil_depth

    ! Subroutine to control the calculation of the observation related
    ! log-likelihoods and accounting for the various scalings

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    ! local variable
    double precision, dimension(DATAin%nodays) :: mod

    !
    ! Do diagnostics (DIAGS)
    !

    ! Calculate log-likelihood for fraction of absorbed radiation
    if (DATAin%nfAPAR > 0) then
        mod = DATAin%M_DIAGS(1:DATAin%nodays,3) & ! APAR
            / (DATAin%met(4,1:DATAin%nodays)*sw_par_fraction)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nfAPAR,DATAin%fAPARpts,DATAin%fAPAR,DATAin%fAPAR_unc,DATAin%fAPAR_lag, &
                                             DATAin%fAPAR_scaling,mod)
    endif ! nfAPAR > 0
    ! Calculate log-likelihood for leaf area index
    if (DATAin%nlai > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nlai,DATAin%laipts,DATAin%LAI,DATAin%LAI_unc,DATAin%LAI_lag, &
                                             DATAin%LAI_scaling,DATAin%M_DIAGS(1:DATAin%nodays,1))
    end if ! nlai > 0

    !
    ! Do pools (POOLS)
    !

    ! Calculate log-likelihood for foliage stocks
    if (DATAin%nCfol_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCfol_stock,DATAin%Cfol_stockpts, &
                                             DATAin%Cfol_stock,DATAin%Cfol_stock_unc,DATAin%Cfol_stock_lag, &
                                             DATAin%Cfol_stock_scaling,DATAin%M_POOLS(1:DATAin%nodays,2))
    endif ! nCfol_stock > 0
    ! Calculate log-likelihood for fine root stocks
    if (DATAin%nCroots_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCroots_stock,DATAin%Croots_stockpts, &
                                             DATAin%Croots_stock,DATAin%Croots_stock_unc,DATAin%Croots_stock_lag, &
                                             DATAin%Croots_stock_scaling,DATAin%M_POOLS(1:DATAin%nodays,3))
    endif ! nCroots_stock > 0
    ! Calculate log-likelihood for total wood stocks
    if (DATAin%nCwood_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_stock,DATAin%Cwood_stockpts, &
                                             DATAin%Cwood_stock,DATAin%Cwood_stock_unc,DATAin%Cwood_stock_lag, &
                                             DATAin%Cwood_stock_scaling,DATAin%M_POOLS(1:DATAin%nodays,4))
    endif ! nCwood_stock > 0
    ! Calculate log-likelihood for foliage litter stocks
    if (DATAin%nClit_stock > 0) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = (sum(DATAin%M_FLUXES(1:DATAin%nodays,10))/sum(DATAin%M_FLUXES(1:DATAin%nodays,10)+DATAin%M_FLUXES(1:DATAin%nodays,12))) &
                 * DATAin%M_POOLS(1:DATAin%nodays,5)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nClit_stock,DATAin%Clit_stockpts, &
                                             DATAin%Clit_stock,DATAin%Clit_stock_unc,DATAin%Clit_stock_lag, &
                                             DATAin%Clit_stock_scaling,mod)
    endif ! nClit_stock > 0
    ! Calculate log-likelihood for soil organic matter stocks
    if (DATAin%nCsom_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCsom_stock,DATAin%Csom_stockpts, &
                                             DATAin%Csom_stock,DATAin%Csom_stock_unc,DATAin%Csom_stock_lag, &
                                             DATAin%Csom_stock_scaling,DATAin%M_POOLS(1:DATAin%nodays,6))
    endif ! nCsom_stock > 0

    !
    ! Do fluxes (FLUXES)
    !

    ! Calculate log-likelihood for evapotranspiration
    if (DATAin%nEvap > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nEvap,DATAin%Evappts,DATAin%Evap,DATAin%Evap_unc,DATAin%Evap_lag, &
                                             DATAin%Evap_scaling,DATAin%M_FLUXES(1:DATAin%nodays,29))
    endif ! nEvap > 0
    ! Calculate log-likelihood for fire
    if (DATAin%nFire > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nFire,DATAin%Firepts,DATAin%Fire,DATAin%Fire_unc,DATAin%Fire_lag, &
                                             DATAin%Fire_scaling,DATAin%M_FLUXES(1:DATAin%nodays,17))
    endif ! nFire > 0
    ! Calculate log-likelihood for gross primary production
    if (DATAin%ngpp > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%ngpp,DATAin%gpppts,DATAin%GPP,DATAin%GPP_unc,DATAin%GPP_lag, &
                                             DATAin%GPP_scaling,DATAin%M_FLUXES(1:DATAin%nodays,1))
    endif ! ngpp > 0
    ! Calculate log-likelihood for harvest
    if (DATAin%nharvest > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nharvest,DATAin%harvestpts,DATAin%harvest,DATAin%harvest_unc,DATAin%harvest_lag, &
                                             DATAin%harvest_scaling,DATAin%M_FLUXES(1:DATAin%nodays,30))
    endif ! nharvest > 0
    ! Calculate log-likelihood for net biome productivity
    if (DATAin%nnbe > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + DATAin%M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + DATAin%M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            + DATAin%M_FLUXES(1:DATAin%nodays,17) & ! Fire
            - DATAin%M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnbe,DATAin%nbepts,DATAin%NBE,DATAin%NBE_unc,DATAin%NBE_lag, &
                                             DATAin%NBE_scaling,mod)
    endif ! nnbe > 0
    ! Calculate log-likelihood for net ecosystem exchange of CO2
    if (DATAin%nnee > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + DATAin%M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + DATAin%M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            - DATAin%M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnee,DATAin%neepts,DATAin%NEE,DATAin%NEE_unc,DATAin%NEE_lag, &
                                             DATAin%NEE_scaling,mod)
    endif ! nnee > 0
    ! Calculate log-likelihood for ecosystem respiration
    if (DATAin%nreco > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,3)  & ! Rauto
            + DATAin%M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + DATAin%M_FLUXES(1:DATAin%nodays,14)   ! Rhet som
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nreco,DATAin%recopts,DATAin%Reco,DATAin%Reco_unc,DATAin%Reco_lag, &
                                             DATAin%Reco_scaling,mod)
    endif ! nreco > 0
    ! Calculate log-likelihood for total wood net increment
    if (DATAin%nCwood_inc > 0) then
        mod = DATAin%M_FLUXES(1:DATAin%nodays,7) - DATAin%M_FLUXES(1:DATAin%nodays,11)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_inc,DATAin%Cwood_incpts, &
                                             DATAin%Cwood_inc,DATAin%Cwood_inc_unc,DATAin%Cwood_inc_lag, &
                                             DATAin%Cwood_inc_scaling,mod)
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood gross increment
    if (DATAin%nCwood_growth > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_growth,DATAin%Cwood_growthpts, &
                                             DATAin%Cwood_growth,DATAin%Cwood_growth_unc,DATAin%Cwood_growth_lag, &
                                             DATAin%Cwood_growth_scaling,DATAin%M_FLUXES(1:DATAin%nodays,7))
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood mortality
    if (DATAin%nCwood_mortality > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_mortality,DATAin%Cwood_mortalitypts, &
                                             DATAin%Cwood_mortality,DATAin%Cwood_mortality_unc,DATAin%Cwood_mortality_lag, &
                                             DATAin%Cwood_mortality_scaling,DATAin%M_FLUXES(1:DATAin%nodays,11))
    endif ! nCwood_mortality > 0

    return

  end subroutine calc_scaled_obs_likelihoods
  !
  !------------------------------------------------------------------
  !
  subroutine calc_other_likelihoods(ML_obs_out)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: top_soil_depth

    ! Subroutine to control the calculation of the 'other priors'
    ! log-likelihoods. These are typically derived variables /
    ! emergent functions averaged over the life of the anaylsis

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    ! local variable
    integer :: dummy_nodays = 1, dummy_noobs = 1
    integer, dimension(1) :: dummy_pts = 1, dummy_lag = 0
    double precision, dimension(1) :: mod
    double precision :: dummy_scaling = 1d0

    ! ...OTHERPRIOR(1)...

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(2) > -9998) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = sum(DATAin%M_FLUXES(1:DATAin%nodays,3)) / sum(DATAin%M_FLUXES(1:DATAin%nodays,1)) ! sum(Rauto) / sum(GPP)
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(2)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(2),DATAin%otherpriorunc(2),dummy_lag,dummy_scaling,mod))
    end if

    ! ...OTHERPRIOR(3)...

    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s -> kg/m2/day)
    if (DATAin%otherpriors(4) > -9998) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = sum(DATAin%M_FLUXES(1:DATAin%nodays,29)) / sum(DATAin%MET(7,1:DATAin%nodays) * 86400d0)
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(4)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(4),DATAin%otherpriorunc(4),dummy_lag,dummy_scaling,mod))
    end if

    ! Estimate the biological steady state attractor on the wood pool.
    ! NOTE: this arrangement explicitly neglects the impact of disturbance on
    ! residence time (i.e. no fire and biomass removal)
    if (DATAin%otherpriors(5) > -9998) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = (sum(DATAin%M_FLUXES(1:DATAin%nodays,7))/dble(DATAin%nodays)) &
            * ( (sum(DATAin%M_POOLS(1:DATAin%nodays,4) / (DATAin%M_FLUXES(1:DATAin%nodays,11)+DATAin%M_FLUXES(1:DATAin%nodays,25)))) / dble(DATAin%nodays))
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(5)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(5),DATAin%otherpriorunc(5),dummy_lag,dummy_scaling,mod))
    end if

    return

  end subroutine calc_other_likelihoods

!
!------------------------------------------------------------------
!
end module model_likelihood_module
