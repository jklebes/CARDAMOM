
subroutine rdaleccdealufires(output_dim,aNPP_dim,MTT_dim,SS_dim,met,pars &
                            ,out_var,out_var2,out_var3,out_var4,out_var5 &
                            ,lat,nopars,nomet &
                            ,nofluxes,nopools,pft,pft_specific,nodays,noyears,deltat &
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
                        ,MTT_dim        & ! Mean Transit time dimension
                        ,SS_dim         & ! Steady State dimension
                        ,pft_specific   & ! Crop model switch - kept for code consistency with other versions
                        ,nos_iter       & ! Number of parameter vectors
                        ,noyears        & ! Number of years to be simulated
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  ! Declare output variables
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var3
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var4
  double precision, intent(out), dimension(nos_iter,MTT_dim,noyears) :: out_var5

  ! Local variables
  ! vector of ecosystem pools
  integer :: i, y, y_s, y_e, nos_years, steps_per_year
  integer, dimension(nodays) :: fol_hak, root_hak, wood_hak, lit_hak, som_hak
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision :: sumNPP, airt_adj
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE & ! net ecosystem exchange of CO2
                                 ,fol_filter &
                                ,root_filter &
                                ,wood_filter &
                                 ,lit_filter &
                                 ,som_filter

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! Create biomass harvest variable
  if (allocated(extracted_C)) deallocate(extracted_C)
  allocate(extracted_C(nodays))

  ! Generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do
  ! Number of years in analysis
  nos_years = nint(sum(deltat)/365.25d0)
  ! number of time steps per year
  Steps_per_year = nodays/nos_years

  ! begin iterations
  do i = 1, nos_iter

     ! reset harvest variable
     extracted_C = 0d0
     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP)

     ! now allocate the output the our 'output' variable
     out_var(i,1:nodays,1)  = lai
     out_var(i,1:nodays,2)  = GPP
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,3) ! auto resp
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) ! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,4) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,6) ! som
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) & ! common pools
                            + POOLS(1:nodays,4) !+ POOLS(1:nodays,5) + POOLS(1:nodays,6)
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,5) ! litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     out_var(i,1:nodays,13) = extracted_C(1:nodays) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,17) ! Fire value

     ! Calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     ! Not strictly needed for this version of DALEC but allows code consistency with more complex versions
     sumNPP = sum(FLUXES(1:nodays,1)*(1d0-pars(2,i))) ! GPP * (1-Ra) fraction
     airt_adj = sum(met(3,1:nodays)) / dble(nodays)
     airt_adj = exp(pars(10,i)*airt_adj)
     out_var2(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) / sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) / sumNPP ! fine root
     out_var2(i,3) = sum(FLUXES(1:nodays,7)) / sumNPP ! wood

     ! Estimate C pool Mean transit (or residence) time
     out_var3(i,1) = pars(5,i) ! leaf life span (years)
     out_var3(i,2) = pars(7,i) ! root residence time (/day)
     out_var3(i,3) = pars(6,i) ! wood residence time (/day)
     out_var3(i,4) = (pars(1,i) + pars(8,i)) * airt_adj ! litter (/day)
     out_var3(i,5) = pars(9,i) * airt_adj ! som (/day)

     ! Determine locations of zeros in pools to correct turnover calculation
     ! Foliage
     fol_hak = 0 ; fol_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,2) == 0) ! protection against NaN from division by zero
           fol_hak = 1 ; fol_filter(1:nodays) = 0d0
     end where
     ! Fine roots
     root_hak = 0 ; root_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,3) == 0) ! protection against NaN from division by zero
           root_hak = 1 ; root_filter(1:nodays) = 0d0
     end where
     ! Wood
     wood_hak = 0 ; wood_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,4) == 0) ! protection against NaN from division by zero
            wood_hak = 1 ; wood_filter(1:nodays) = 0d0
     end where
     ! Fol+root litter
     lit_hak = 0 ; lit_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,5) == 0) ! protection against NaN from division by zero
            lit_hak = 1 ; lit_filter(1:nodays) = 0d0
     end where
     ! Soil
     som_hak = 0 ; som_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,6) == 0) ! protection against NaN from division by zero
           som_hak = 1 ; som_filter(1:nodays) = 0d0
     end where

     ! Estimate the MTT for each year individually, accounting for disturbance factors
     do y = 1, nos_years
        ! Estimate time steps covered by this year
        y_s = 1 + (steps_per_year * (y-1)) ; y_e = steps_per_year * y
        ! Foliage
        out_var5(i,1,y) = sum( ((FLUXES(y_s:y_e,10)+FLUXES(y_s:y_e,19)+FLUXES(y_s:y_e,25)) &
                                / POOLS(y_s:y_e,2)) * fol_filter(y_s:y_e)) / dble(steps_per_year-sum(fol_hak(y_s:y_e)))
        ! Fine roots
        out_var5(i,2,y) = sum( ((FLUXES(y_s:y_e,12)+FLUXES(y_s:y_e,20)+FLUXES(y_s:y_e,26)) &
                                / POOLS(y_s:y_e,3)) * root_filter(y_s:y_e)) / dble(steps_per_year-sum(root_hak(y_s:y_e)))
        ! Wood
        out_var5(i,3,y) = sum( ((FLUXES(y_s:y_e,11)+FLUXES(y_s:y_e,21)+FLUXES(y_s:y_e,27)) &
                                / POOLS(y_s:y_e,4)) * wood_filter(y_s:y_e)) / dble(steps_per_year-sum(wood_hak(y_s:y_e)))
        ! Litter (fol+fine root)
        out_var5(i,4,y) = sum( ((FLUXES(y_s:y_e,13)+FLUXES(y_s:y_e,15)+FLUXES(y_s:y_e,22)+FLUXES(y_s:y_e,28)) &
                                / POOLS(y_s:y_e,5)) * lit_filter(y_s:y_e)) / dble(steps_per_year-sum(lit_hak(y_s:y_e)))
        ! Soil
        out_var5(i,5,y) = sum( ((FLUXES(y_s:y_e,14)+FLUXES(y_s:y_e,23)) &
                                / POOLS(y_s:y_e,6)) * som_filter(y_s:y_e)) / dble(steps_per_year-sum(som_hak(y_s:y_e)))
     end do

     ! Once the canopy has closed the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs. Similarly the litter pool input
     ! i.e. foliage and fine roots are likely to be in steady state.
     out_var4(i,1) = sum(FLUXES(:,4)+FLUXES(:,8)) ! Foliage
     out_var4(i,2) = sum(FLUXES(:,6)) ! Fine root
     out_var4(i,3) = sum(FLUXES(:,7)) ! Wood
     out_var4(i,4) = sum(FLUXES(1:nodays,10)+FLUXES(1:nodays,12) &
                        +FLUXES(1:nodays,24)+FLUXES(1:nodays,25)+FLUXES(1:nodays,26)) ! lit
     ! However soil C inputs are still changing as the wood pool is not in steady state.
     ! Therefore, at this point we can account for the fol, fine root, and disturbance inputs but NOT wood.
     ! The wood input is estimated later based on the steady state its steady state estimate
     out_var4(i,5) = sum(FLUXES(:,15)+FLUXES(1:nodays,27)+FLUXES(1:nodays,28)) ! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  ! NOTE: out_var3(:,1) is canopy residence time extracted directly in years from pars(5)
  out_var3(:,2:5) = (out_var3(:,2:5)*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)
  out_var5 = (out_var5*365.25d0)**(-1d0)               ! iter,(fol,root,wood,lit+litwood,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var4 = (out_var4 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the foliar, fine root, wood and litter steady states.
  out_var4(:,1:4) = out_var4(:,1:4) * out_var3(:,1:4) ! multiply by residence time in years
  ! Using the wood SS estimate (gC/m2) the steady state input to the som litter pool...
  out_var4(:,5) = (out_var4(:,5) + (out_var4(:,3) / out_var3(:,3))) * out_var3(:,5)

  ! deallocate harvested variable
  deallocate(extracted_C)

  ! return back to the subroutine then
  return

end subroutine rdaleccdealufires
