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
! Subroutine to allow direct interface between ACM-GPP-ET and the R code
!
! Author: T. Luke Smallman (02/05/2024)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine racm(output_dim,met,pars,out_var,lat,nopars,nomet &
               ,nofluxes,nopools,pft,pft_specific,nodays,deltat &
               ,nos_iter,soil_frac_clay_in,soil_frac_sand_in)

  use CARBON_MODEL_MOD, only: CARBON_MODEL &
                             ,soil_frac_clay, soil_frac_sand &
                             ,nos_soil_layers, wSWP_time &
                             ,gs_demand_supply_ratio, canopy_par_MJday_time &
                             ,gs_total_canopy, gb_total_canopy, cica_time

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of parameters in vector
                        ,output_dim     & !
                        ,pft            & ! plant functional type
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(in) :: met(nomet,nodays)     & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)           & ! number of parameters
                       ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                       ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,lat                               ! site latitude (degrees)

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days


  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var

  ! local variables
  ! vector of ecosystem pools
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  integer i
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! update soil parameters
  soil_frac_clay = soil_frac_clay_in
  soil_frac_sand = soil_frac_sand_in

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do

  do i = 1, nos_iter

     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP)
!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/out.csv", &
!         status='replace',action='readwrite' )
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     ! now allocate the output the our 'output' variable
     out_var(i,1:nodays,1)  = lai
     out_var(i,1:nodays,2)  = GPP
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,2) ! Evap (kg.m-2.day-1)
     if (output_dim > 3) then
         out_var(i,1:nodays,4) = FLUXES(1:nodays,3) ! soil evap (kg.m-2.day-1)
         out_var(i,1:nodays,5) = wSWP_time(1:nodays)
         out_var(i,1:nodays,6) = FLUXES(1:nodays,4) ! wet canopy evap (kg.m-2.day-1)
         out_var(i,1:nodays,7) = gs_demand_supply_ratio(1:nodays)
         out_var(i,1:nodays,8) = gs_total_canopy(1:nodays)
         out_var(i,1:nodays,9) = canopy_par_MJday_time(1:nodays)
         out_var(i,1:nodays,10) = gb_total_canopy(1:nodays)
         out_var(i,1:nodays,11) = cica_time(1:nodays)
     endif

  end do ! nos_iter loop

  ! return back to the subroutine then
  return

end subroutine racm
