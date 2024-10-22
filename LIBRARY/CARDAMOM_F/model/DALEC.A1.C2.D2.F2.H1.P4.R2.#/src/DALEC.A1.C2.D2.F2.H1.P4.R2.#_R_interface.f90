

subroutine rdalec10(output_dim,aNPP_dim,MTT_dim,SS_dim,fire_dim &
                   ,met,pars &
                   ,out_var1,out_var2,out_var3,out_var4,out_var5 &
                   ,lat,nopars,nomet,nofluxes,nopools,pft &
                   ,nodays,nos_years,deltat,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, itemp, ivpd, iphoto,                               &
                              nos_soil_layers,                 &
                              harvest_residue_to_litter, harvest_residue_to_woodlitter,        &
                              harvest_residue_to_som,                                          &
                              harvest_residue_labile, harvest_residue_foliar,                  &
                              harvest_residue_roots, harvest_residue_wood,                     &
                              harvest_extracted_woodlitter, harvest_extracted_som,             &
                              harvest_extracted_labile, harvest_extracted_foliar,              &
                              harvest_extracted_roots, harvest_extracted_wood,                 &
                              harvest_extracted_litter,                                        &
                              fire_emiss_labile, fire_emiss_foliar, fire_emiss_roots,          &
                              fire_emiss_wood, fire_emiss_litter, fire_emiss_woodlitter,       &
                              fire_emiss_som, fire_litter_labile, fire_litter_foliar,          &
                              fire_litter_roots, fire_litter_wood, fire_litter_litter,         &
                              fire_litter_woodlitter, fire_litter_som, fire_residue_to_litter, &
                              fire_residue_to_woodlitter,fire_residue_to_som,                  &
                              gs_demand_supply_ratio, cica_time, Rg_from_labile,               &
                              gs_total_canopy, gb_total_canopy, canopy_par_MJday_time,         &
                              root_depth_time, Rm_from_labile

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
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,MTT_dim        &
                        ,SS_dim         &
                        ,fire_dim       &
                        ,pft            & ! plant functional type
                        ,nos_iter       & !
                        ,nos_years        &
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var1
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var2  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var3  ! Steady State (gC/m2)
  double precision, intent(out), dimension(nos_iter,output_dim) :: out_var4 ! Long term mean of out_var1
  double precision, intent(out), dimension(nos_iter,nos_years,output_dim) :: out_var5 ! Mean annual of out_var1

  ! local variables
  integer :: i, y, y_s, y_e, steps_per_year
  integer, dimension(nodays) :: pool_hak
  double precision, dimension(nos_iter) :: woodlitter_to_som_frac
  ! vector of ecosystem pools
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: resid_fol
  integer, dimension(nodays) :: hak ! variable to determine number of NaN
  double precision :: sumNPP
  double precision, dimension(nodays) :: tmp, tmp1 &
                                        ,lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0,total_time = 0
!real :: Rtot_track_time = 0, aero_time = 0 , soilwater_time = 0 , acm_et_time = 0 , Rm_time = 0
!call cpu_time(done)
!print*,"time taken per iter",(done-begin) / real(nos_iter)

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0
  out_var1 = 0d0 ; out_var2 = 0d0 ; out_var3 = 0d0

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do
  ! number of time steps per year
  steps_per_year = nodays/nos_years

  ! begin iterations
  do i = 1, nos_iter
     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP)
!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/out.csv", &
!         status='replace',action='readwrite' )
!write(666,*)"deltat",deltat
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     !
     ! Allocate the output the our 'output' variable
     !

     ! C ecosystem fluxes (gC/m2/day)
     out_var1(i,1:nodays,1)  = FLUXES(1:nodays,1)       ! GPP (gC/m2/day)
     out_var1(i,1:nodays,2)  = FLUXES(1:nodays,3)       ! Rauto (gC/m2/day)
     out_var1(i,1:nodays,3)  = Rg_from_labile(1:nodays) ! Rgrowth_leaf (gC/m2/day)
     out_var1(i,1:nodays,4)  = FLUXES(1:nodays,13)      ! Rhet_litter (gC/m2/day)
     out_var1(i,1:nodays,5)  = FLUXES(1:nodays,14)      ! Rhet_som (gC/m2/day)
     out_var1(i,1:nodays,6)  = FLUXES(1:nodays,4)       ! Rhet_woodlitter (gC/m2/day)
     out_var1(i,1:nodays,7)  = FLUXES(1:nodays,17)      ! Total fire (gC/m2/day)
     out_var1(i,1:nodays,8)  = FLUXES(1:nodays,21)      ! harvested material (gC/m2/day)
     ! C internal fluxes (gC/m2/day)
     out_var1(i,1:nodays,9)  = FLUXES(1:nodays,5)       ! allocation to labile (gC/m2/day)
     out_var1(i,1:nodays,10) = FLUXES(1:nodays,6)       ! allocation to fine roots (gC/m2/day)
     out_var1(i,1:nodays,11) = FLUXES(1:nodays,7)       ! allocation to wood (gC/m2/day)
     out_var1(i,1:nodays,12) = FLUXES(1:nodays,8)       ! labile to foliage (gC/m2/day)
     out_var1(i,1:nodays,13) = FLUXES(1:nodays,10)      ! foliage scenesence (gC/m2/day)
     out_var1(i,1:nodays,14) = FLUXES(1:nodays,12)      ! fine root turnover (gC/m2/day)
     out_var1(i,1:nodays,15) = FLUXES(1:nodays,11)      ! wood turnover (gC /m2/day)
     out_var1(i,1:nodays,16) = FLUXES(1:nodays,15)      ! Decomp_litter (gC/m2/day)
     out_var1(i,1:nodays,17) = FLUXES(1:nodays,20)      ! Decomp_woodlitter (gC/m2/day)
     ! C disturbance fluxes (gC/m2/day)
     out_var1(i,1:nodays,18) = fire_emiss_labile            ! fire emission from labile (gC/m2/day)
     out_var1(i,1:nodays,19) = fire_litter_labile           ! fire induced litter from labile (gC/m2/day)
     out_var1(i,1:nodays,20) = fire_emiss_foliar            ! fire emission from foliage (gC/m2/day)
     out_var1(i,1:nodays,21) = fire_litter_foliar           ! fire induced litter from foliage (gC/m2/day)
     out_var1(i,1:nodays,22) = fire_emiss_roots             ! fire emission from fine roots (gC/m2/day)
     out_var1(i,1:nodays,23) = fire_litter_roots            ! fire induced litter from fine roots (gC/m2/day)
     out_var1(i,1:nodays,24) = fire_emiss_wood              ! fire emission from wood (gC/m2/day)
     out_var1(i,1:nodays,25) = fire_litter_wood             ! fire induced litter from wood (gC/m2/day)
     out_var1(i,1:nodays,26) = fire_emiss_litter            ! fire emission from litter (gC/m2/day)
     out_var1(i,1:nodays,27) = fire_litter_litter           ! fire induced litter from litter (gC/m2/day)
     out_var1(i,1:nodays,28) = fire_emiss_woodlitter        ! fire emission from wood litter (gC/m2/day)
     out_var1(i,1:nodays,29) = fire_litter_woodlitter       ! fire induced litter from wood litter (gC/m2/day)
     out_var1(i,1:nodays,30) = fire_emiss_som               ! fire emission from som (gC/m2/day)
     out_var1(i,1:nodays,31) = harvest_extracted_labile     ! harvest extracted from labile (gC/m2/day)
     out_var1(i,1:nodays,32) = harvest_extracted_foliar     ! harvest extracted from foliage (gC/m2/day)
     out_var1(i,1:nodays,33) = harvest_extracted_roots      ! harvest extracted from fine roots (gC/m2/day)
     out_var1(i,1:nodays,34) = harvest_extracted_wood       ! harvest extracted from wood (gC/m2/day)
     out_var1(i,1:nodays,35) = harvest_extracted_litter     ! harvest extracted from litter (gC/m2/day)
     out_var1(i,1:nodays,36) = harvest_extracted_woodlitter ! harvest extracted from wood litter (gC/m2/day)
     out_var1(i,1:nodays,37) = harvest_extracted_som        ! harvest extracted from som (gC/m2/day)
     out_var1(i,1:nodays,38) = harvest_residue_labile       ! harvest litter / residue from labile (gC/m2/day)
     out_var1(i,1:nodays,39) = harvest_residue_foliar       ! harvest litter / residue from foliage (gC/m2/day)
     out_var1(i,1:nodays,40) = harvest_residue_roots        ! harvest litter / residue from fine roots (gC/m2/day)
     out_var1(i,1:nodays,41) = harvest_residue_wood         ! harvest litter / residue from wood (gC/m2/day)
     ! C pools (gC/m2)
     out_var1(i,1:nodays,42) = POOLS(1:nodays,1)        ! labile (gC/m2)
     out_var1(i,1:nodays,43) = POOLS(1:nodays,2)        ! foliage (gC/m2)
     out_var1(i,1:nodays,44) = POOLS(1:nodays,3)        ! fine root (gC/m2)
     out_var1(i,1:nodays,45) = POOLS(1:nodays,4)        ! wood (gC/m2)
     out_var1(i,1:nodays,46) = POOLS(1:nodays,5)        ! litter (gC/m2)
     out_var1(i,1:nodays,47) = POOLS(1:nodays,7)        ! wood litter (gC/m2)
     out_var1(i,1:nodays,48) = POOLS(1:nodays,6)        ! som (gC/m2)
     ! Canopy (phenology) properties
     out_var1(i,1:nodays,49) = lai                      ! LAI (m2/m2)
     out_var1(i,1:nodays,50) = FLUXES(1:nodays,18)      ! GSI value (0-1)
     out_var1(i,1:nodays,51) = itemp(1:nodays)          ! GSI temp component (0-1)
     out_var1(i,1:nodays,52) = iphoto(1:nodays)         ! GSI photoperiod component (0-1)
     out_var1(i,1:nodays,53) = ivpd(1:nodays)           ! GSI vpd component (0-1)
     ! Photosynthesis / C~water coupling related
     out_var1(i,1:nodays,54) = gs_demand_supply_ratio   ! ratio of evaporative demand over supply
     out_var1(i,1:nodays,55) = gs_total_canopy          ! stomatal conductance (mmolH2O/m2ground/day)
     out_var1(i,1:nodays,56) = canopy_par_MJday_time    ! Canopy absorbed PAR (MJ/m2ground/day)
     out_var1(i,1:nodays,57) = gb_total_canopy          ! boundary conductance (mmolH2O/m2ground/day)
     out_var1(i,1:nodays,58) = cica_time                ! ratio of leaf internal to external CO2
     ! misc
     out_var1(i,1:nodays,59) = root_depth_time          ! rooting depth (m)

     !
     ! Calculate long-term mean of out_var1
     !
     
     ! Loop across each variable
     do v = 1, output_dim
        ! Calculate mean value
        out_var4(i,v) = sum(out_var1(i,1:nodays,v)) / dble(nodays)
     end do

     !
     ! Calculate the mean annual of out_var1
     !

     ! Calculate mean annual
     s = 1 ; e = steps_per_year
     do a = 1, nos_years
        do v = 1, output_dim
           out_var5(i,a,v) = sum(out_var1(i,s:e,v)) / dble(steps_per_year)
        end do
        ! Iterate counters
        s = s + steps_per_year ; e = s + steps_per_year - 1
     end do
     
     !!!
     ! Estimate residence time information
     !!!

     ! Labile
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,1) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,8) + Rg_from_labile + &
                    fire_emiss_labile + fire_litter_labile + &
                    harvest_extracted_labile + harvest_residue_labile ) / POOLS(1:nodays,1))
     end where
     out_var2(i,1) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,2) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,10) + fire_emiss_foliar + &
                    fire_litter_foliar + harvest_extracted_foliar + &
                    harvest_residue_foliar ) / POOLS(1:nodays,2))
     end where
     out_var2(i,2) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Fine roots
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,3) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,12) + fire_emiss_roots + &
                    fire_litter_roots + harvest_extracted_roots + &
                    harvest_residue_roots ) / POOLS(1:nodays,3))
     end where
     out_var2(i,3) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Wood
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,4) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,11) + fire_emiss_wood + &
                    fire_litter_wood + harvest_extracted_wood + &
                    harvest_residue_wood ) / POOLS(1:nodays,4))
     end where
     out_var2(i,4) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage + fine root litter
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,5) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,13) + FLUXES(1:nodays,15) + &
                    fire_emiss_litter + fire_litter_litter + &
                    harvest_extracted_litter) / POOLS(1:nodays,5))
     end where
     out_var2(i,5) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Wood litter
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,7) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,20) + FLUXES(1:nodays,4) + &
                    fire_emiss_woodlitter + fire_litter_woodlitter + &
                    harvest_extracted_woodlitter) / POOLS(1:nodays,7))
            tmp1 = FLUXES(1:nodays,20) / POOLS(1:nodays,7)
     end where
     out_var2(i,6) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Keep track of the fraction of wood litter transfer to som, this value is needed for the steady state estimation
     woodlitter_to_som_frac(i) = sum(tmp1) / dble(nodays-sum(pool_hak))
     ! Soil
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,6) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,14) + fire_emiss_som + harvest_extracted_som) &
                  / POOLS(1:nodays,6))
     end where
     out_var2(i,7) = sum(tmp) / dble(nodays-sum(pool_hak))

     !
     ! Estimate pool inputs needed for steady state calculation
     !

     ! Once the canopy has closes the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs
     out_var3(i,1) = sum(FLUXES(:,5)) ! Labile
     out_var3(i,2) = sum(FLUXES(:,8)) ! Foliage
     out_var3(i,3) = sum(FLUXES(:,6)) ! Fine root
     out_var3(i,4) = sum(FLUXES(:,7)) ! Wood
     out_var3(i,5) = sum(FLUXES(:,10)+FLUXES(:,12)+ &
                         fire_residue_to_litter + harvest_residue_to_litter) ! litter (foliage + roots)
     ! While foliar and fine root litter can be reasonably estimated directly (above),
     ! wood litter and soil C inputs are still changing as the wood pool is not in steady state.
     ! Therefore, at this point we can account for disturbance inputs but NOT wood.
     ! The wood input is estimated later based on the steady state its steady state estimate
     out_var3(i,6) = sum(fire_residue_to_woodlitter + harvest_residue_to_woodlitter) ! woodlitter
     out_var3(i,7) = sum(FLUXES(:,15)+fire_residue_to_som) ! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var2 = (out_var2*365.25d0)**(-1d0) ! iter,(lab,fol,root,wood,lit,woodlitter,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var3 = (out_var3 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the labile, foliar, fine root, wood steady and fine litter states.
  out_var3(:,1:5) =  out_var3(:,1:5) * out_var2(:,1:5) ! multiply by residence time in years
  ! Using the wood SS estimate (gC/m2) the steady state input to the wood litter pool...
  out_var3(:,6) = (out_var3(:,6) + (out_var3(:,4) / out_var2(:,4))) * out_var2(:,6)
  ! ...which is then in turn used to update the soil pool
  ! NOTE: that because not all wood litter
  out_var3(:,7) = (out_var3(:,7) + ((out_var3(:,6) / out_var2(:,6))*woodlitter_to_som_frac) ) * out_var2(:,7)

  ! deallocate harvested variable
  deallocate(itemp,ivpd,iphoto)

  ! return back to the subroutine then
  return

end subroutine rdalec10
!
!------------------------------------------------------------------
!
