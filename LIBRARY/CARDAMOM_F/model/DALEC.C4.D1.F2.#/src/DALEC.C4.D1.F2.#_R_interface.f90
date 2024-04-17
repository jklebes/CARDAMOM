

subroutine rdalec12(output_dim,MTT_dim,SS_dim &
                   ,met,pars &
                   ,out_var1,out_var2,out_var3,out_var4,out_var5 &
                   ,lat,nopars,nomet &
                   ,nofluxes,nopools,nodays,nos_years,deltat &
                   ,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL

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
                        ,output_dim     & ! number of outputted variables
                        ,MTT_dim        & ! number of pools mean transit time estimates
                        ,SS_dim         & ! number of pools the steady state will be output for
                        ,nos_iter       & ! number of iterations
                        ,nos_years        & ! number of years simulated
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                             ,pars(nopars,nos_iter)   & ! number of parameters
                             ,lat                       ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var1
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var2  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var3  ! Steady State (gC/m2)
  double precision, intent(out), dimension(nos_iter,output_dim) :: out_var4 ! Long term mean of out_var1
  double precision, intent(out), dimension(nos_iter,nos_years,output_dim) :: out_var5 ! Mean annual of out_var1

  ! local variables
  ! vector of ecosystem pools
  integer :: i, y, y_s, y_e, nos_years, steps_per_year
  integer, dimension(nodays) :: pool_hak
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: tmp &
                                        ,lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2s

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
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     !
     ! Allocate the output the our 'output' variable
     !

     ! C ecosystem fluxes (gC/m2/day)
     out_var1(i,1:nodays,1)  = FLUXES(1:nodays,1)       ! GPP (gC/m2/day)
     out_var1(i,1:nodays,2)  = FLUXES(1:nodays,3)       ! Rauto (gC/m2/day)
     out_var1(i,1:nodays,3)  = FLUXES(1:nodays,13)      ! Rhet_litter + som (gC/m2/day)
     out_var1(i,1:nodays,4)  = FLUXES(1:nodays,17)      ! Total fire (gC/m2/day)
     out_var1(i,1:nodays,5)  = FLUXES(1:nodays,23)      ! total harvested material (gC/m2/day)
     ! C internal fluxes (gC/m2/day)
     out_var1(i,1:nodays,6)  = FLUXES(1:nodays,4)       ! allocation to foliage (gC/m2/day)
     out_var1(i,1:nodays,7)  = FLUXES(1:nodays,6)       ! allocation to fine roots + wood (gC/m2/day)
     out_var1(i,1:nodays,8)  = FLUXES(1:nodays,10)      ! foliage to litter (gC/m2/day)
     out_var1(i,1:nodays,9)  = FLUXES(1:nodays,11)      ! wood + root to som (gC /m2/day)
     ! C disturbance fluxes (gC/m2/day)
     out_var1(i,1:nodays,10) = FLUXES(1:nodays,18)      ! fire emission from foliage (gC/m2/day)
     out_var1(i,1:nodays,11) = FLUXES(1:nodays,21)      ! fire induced litter from foliage (gC/m2/day)
     out_var1(i,1:nodays,12) = FLUXES(1:nodays,19)      ! fire emission from root+wood (gC/m2/day)
     out_var1(i,1:nodays,13) = FLUXES(1:nodays,22)      ! fire induced litter from root+wood (gC/m2/day)
     out_var1(i,1:nodays,14) = FLUXES(1:nodays,20)      ! fire emission from dom (gC/m2/day)
     out_var1(i,1:nodays,15) = FLUXES(1:nodays,24)      ! harvest extracted from foliage (gC/m2/day)
     out_var1(i,1:nodays,16) = FLUXES(1:nodays,25)      ! harvest extracted from root+wood (gC/m2/day)
     out_var1(i,1:nodays,17) = FLUXES(1:nodays,26)      ! harvest extracted from som (gC/m2/day)
     out_var1(i,1:nodays,18) = FLUXES(1:nodays,27)      ! harvest litter / residue from foliage (gC/m2/day)
     out_var1(i,1:nodays,19) = FLUXES(1:nodays,28)      ! harvest litter / residue from root+wood (gC/m2/day)
     ! C pools (gC/m2)
     out_var1(i,1:nodays,20) = POOLS(1:nodays,1)        ! foliage (gC/m2)
     out_var1(i,1:nodays,21) = POOLS(1:nodays,2)        ! fine root + wood (gC/m2)
     out_var1(i,1:nodays,22) = POOLS(1:nodays,3)        ! dom (gC/m2)
     ! Canopy (phenology) properties
     out_var1(i,1:nodays,23) = lai                      ! LAI (m2/m2)

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

     ! Foliage
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,1) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,10) + &
                    FLUXES(1:nodays,18) + FLUXES(1:nodays,21) + &
                    FLUXES(1:nodays,24) + FLUXES(1:nodays,27)) / POOLS(1:nodays,1))
     end where
     out_var2(i,1) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Fine roots + wood
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,2) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,11) + &
                    FLUXES(1:nodays,19) + FLUXES(1:nodays,22) + &
                    FLUXES(1:nodays,25) + FLUXES(1:nodays,28)) / POOLS(1:nodays,2))
     end where
     out_var2(i,2) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage + fine root litter + som
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,3) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,13) + FLUXES(1:nodays,20) + &
                    FLUXES(1:nodays,26)) / POOLS(1:nodays,3))
     end where
     out_var2(i,3) = sum(tmp) / dble(nodays-sum(pool_hak))

     !
     ! Estimate pool inputs needed for steady state calculation
     !

     ! Once the canopy has closed the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs
     out_var3(i,1) = sum(FLUXES(:,4)) ! Foliage
     out_var3(i,2) = sum(FLUXES(:,6)) ! Fine root + wood
     ! While foliar and fine root litter can be reasonably estimated directly (above),
     ! soil C inputs are still changing as the wood pool is not in steady state.
     ! Therefore, at this point we can account for disturbance inputs but NOT wood.
     ! The wood input is estimated later based on the steady state its steady state estimate
     out_var3(i,3) = sum(FLUXES(:,10)+FLUXES(:,11)+ &
                         FLUXES(:,21)+FLUXES(:,22)+ &
                         FLUXES(:,27)+FLUXES(:,28)) ! dom
  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var2 = (out_var2*365.25d0)**(-1d0) ! iter,(lab,fol,root,wood,lit,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var3 = (out_var3 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the labile, foliar, fine root, wood and litter steady states.
  out_var3(:,1:2) = out_var3(:,1:2) * out_var2(:,1:2) ! multiply by residence time in years
  ! Using the wood SS estimate (gC/m2) the steady state input to the som litter pool...
  out_var3(:,3) = (out_var3(:,3) + (out_var3(:,2) / out_var2(:,2))) * out_var2(:,3)

  ! return back to the subroutine then
  return

end subroutine rdalec12
