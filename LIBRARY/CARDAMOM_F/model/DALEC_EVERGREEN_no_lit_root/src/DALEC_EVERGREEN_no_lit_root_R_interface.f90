

subroutine rdalecevergreennolitroot(output_dim,aNPP_dim,met,pars,out_var,out_var2,lat,nopars,nomet &
                                   ,nofluxes,nopools,pft,pft_specific,nodays,deltat &
                                   ,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, extracted_C

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,pft            & ! plant functional type
                        ,output_dim     & !
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2

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
  allocate(extracted_C(nodays+1))

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do

  ! begin iterations
  do i = 1, nos_iter
     ! reset harvest variable
     extracted_C=0.
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
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,3)  ! auto resp
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) ! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,3) ! litter + som
     out_var(i,1:nodays,7)  = 0d0  
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) ! Biomass
     out_var(i,1:nodays,9)  = POOLS(1:nodays,2) ! root+wood
     out_var(i,1:nodays,10) = 0d0
     out_var(i,1:nodays,11) = 0d0  
     out_var(i,1:nodays,12) = POOLS(1:nodays,1) ! foliage
     out_var(i,1:nodays,13) = extracted_C(1:nodays) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,17) ! Fire value

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     sumNPP = sum(FLUXES(1:nodays,1)*(1d0-pars(1,i))) ! GPP * (1-Ra) fraction
     airt_adj = sum(met(3,1:nodays)) / dble(nodays)
     airt_adj = exp(pars(6,i)*airt_adj)
     out_var2(i,1) = sum(FLUXES(1:nodays,4)) / sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) / sumNPP ! wood + root
     out_var2(i,3) = 0d0
     out_var2(i,4) = pars(3,i) ! leaf life span (years)
     out_var2(i,5) = (pars(4,i)*365.25d0) ** (-1d0) ! wood+root residence time (years)
     out_var2(i,6) = 0d0
     out_var2(i,7) = (pars(5,i) * 365.25d0 * airt_adj) ** (-1d0) ! som
     out_var2(i,8) = 0d0

  end do ! nos_iter loop

  ! deallocate harvested variable
  deallocate(extracted_C)

  ! return back to the subroutine then
  return

end subroutine rdalecevergreennolitroot
