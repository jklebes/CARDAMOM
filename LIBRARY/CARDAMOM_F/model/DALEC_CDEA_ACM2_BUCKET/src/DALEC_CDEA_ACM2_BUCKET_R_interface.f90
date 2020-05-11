

subroutine rdaleccdeaacm2bucket(output_dim,aNPP_dim,MTT_dim,SS_dim,met,pars,out_var,out_var2,out_var3,out_var4 &
                               ,lat,nopars,nomet &
                               ,nofluxes,nopools,pft,pft_specific,nodays,deltat &
                               ,nos_iter,soil_frac_clay_in,soil_frac_sand_in)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, extracted_C, wSWP_time &
                             ,soil_frac_clay, soil_frac_sand, nos_soil_layers &
                             ,gs_demand_supply_ratio &
                             ,gs_total_canopy, gb_total_canopy, canopy_par_MJday_time

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,pft            & ! plant functional type
                        ,output_dim     & !
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,MTT_dim        &
                        ,SS_dim         &
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                  ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                  ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var3
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var4

  ! local variables
  ! vector of ecosystem pools
  integer i
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision :: sumNPP,airt_adj
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! update settings
  if (allocated(extracted_C)) deallocate(extracted_C)
  allocate(extracted_C(nodays))

  ! update soil parameters
  soil_frac_clay = soil_frac_clay_in
  soil_frac_sand = soil_frac_clay_in

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do

  ! begin iterations
  do i = 1, nos_iter
     ! reset harvest variable
     extracted_C = 0d0
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
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,3) ! auto resp
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) ! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,4) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,6) ! som
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) & ! common pools
                              + POOLS(1:nodays,4) + POOLS(1:nodays,5) + POOLS(1:nodays,6)
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,5) ! litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     out_var(i,1:nodays,13) = extracted_C(1:nodays) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,17) ! Fire value

     out_var(i,1:nodays,18) = FLUXES(1:nodays,29) ! Evapotranspiration (kgH2O.m-2.day-1)
     out_var(i,1:nodays,19) = POOLS(1:nodays,7)   ! rootwater (kgH2O.m-2.10cmdepth)
     out_var(i,1:nodays,20) = wSWP_time(1:nodays) ! Weighted Soil Water Potential (MPa)
     out_var(i,1:nodays,21) = gs_demand_supply_ratio ! ratio of evaporative demand over supply
     out_var(i,1:nodays,22) = gs_total_canopy(1:nodays)
     out_var(i,1:nodays,23) = canopy_par_MJday_time(1:nodays)
     out_var(i,1:nodays,24) = gb_total_canopy(1:nodays)

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     sumNPP = sum(FLUXES(1:nodays,1)*(1-pars(2,i))) ! GPP * (1-Ra) fraction
     airt_adj = sum(met(3,1:nodays)) / dble(nodays)
     airt_adj = exp(pars(10,i)*airt_adj)
     out_var2(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) / sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) / sumNPP ! fine root
     out_var2(i,3) = sum(FLUXES(1:nodays,7)) / sumNPP ! wood
     ! Mean transit time
     out_var3(i,1) = pars(5,i) ! leaf life span (years)
     out_var3(i,2) = (pars(7,i)*365.25d0) ** (-1d0) ! root residence time (years)
     out_var3(i,3) = (pars(6,i)*365.25d0) ** (-1d0) ! wood residence time (years)
     out_var3(i,4) = ((pars(1,i) + pars(8,i)) * 365.25d0 * airt_adj) ** (-1d0) ! litter
     out_var3(i,5) = (pars(9,i) * 365.25d0 * airt_adj) ** (-1d0) ! som

     ! Calculate the mean inputs to each pool, needed for steady state calculation
     out_var4(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) ! fol
     out_var4(i,2) = sum(FLUXES(1:nodays,6)) ! root
     out_var4(i,3) = sum(FLUXES(1:nodays,7)) ! wood
     out_var4(i,4) = sum(FLUXES(1:nodays,10)+FLUXES(1:nodays,12)) ! lit
     out_var4(i,5) = sum(FLUXES(1:nodays,15))! som

  end do ! nos_iter loop

  ! Steady state gC/m2
  out_var4 = (out_var4 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  out_var4 = out_var4 * out_var3     ! multiply by residence time in years

  ! deallocate harvested variable
  deallocate(extracted_C)

  ! return back to the subroutine then
  return

end subroutine rdaleccdeaacm2bucket
