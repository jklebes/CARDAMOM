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

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains uniform prior parameter information for the DALEC.D1.F2 model.
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module MODEL_PARAMETERS
use samplers_shared, only: PARINFO

  implicit none

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine pars_info(PI)
    

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

    type(PARINFO), intent(inout):: PI

    PI%npars = 22
    if (.not. allocated(PI%parmin)) allocate(PI%parmin(PI%npars))
    if (.not. allocated(PI%parmax)) allocate(PI%parmax(PI%npars))

    ! Decomposition of litter to som (fraction / day-1)
    ! Note is modified by exponential temperature function (p10)
    PI%parmin(1) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(1) = 0.02d0      ! 0.13 years at 0oC

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
    PI%parmax(5) = 6d0 ! 8d0

    ! TOR wood* - 1% loss per year value
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! Turnover fraction of roots
    ! Gill and Jackson (2000), New Phytol., 147, 13–31
    ! Fig. 6 turnover by diameter class
    PI%parmin(7) = 0.001368925d0 ! 2    years !0.0006844627d0 ! 4 years
    PI%parmax(7) = 0.01d0        ! 0.27 years   

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(9) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC
!    PI%parmin(9) = 0.0000001d0 ! 27378.0 years at 0oC
!    PI%parmax(9) = 0.001d0     !     2.7 years at 0oC

    ! Temp factor* = Q10 = 1.2-2.2
    PI%parmin(10) = 0.019d0
    PI%parmax(10) = 0.08d0

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a prior of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    PI%parmin(11) = 1.64d0 !5d0
    PI%parmax(11) = 42d0 !50d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(12) = 5d0
    PI%parmax(12) = 200d0

    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(18) = 0.1d0
    PI%parmax(18) = 1d0
    ! Combustion completeness factor for foliage
    PI%parmin(19) = 0.01d0
    PI%parmax(19) = 0.99d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(20) = 0.01d0
    PI%parmax(20) = 0.99d0
    ! Combustion completeness factor for soil
    PI%parmin(21) = 0.001d0
    PI%parmax(21) = 0.1d0
    ! Combustion completeness factor for foliage + fine root litter
    PI%parmin(22) = 0.01d0
    PI%parmax(22) = 0.99d0

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
