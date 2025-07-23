# How to edit models for thread safety
### src / model .f90 file 
Edit `src/<model>.f90` file , assited by `cardamom_model_type.py` .  See model 004 before and after example.  As written in the python script's top comment:
1. Insert the line `type model_working_variables` after the list of module parameter variables and before the list of non-parameter variables, before the line  `double precision :: minlwp = minlwp_default`
2. Insert the line `end type` at the close of list of module variables, before `contains` .
3. Get rid of all the variables in "public:" block except  CARBON_MODEL and read-only parameters USEd by likelihood file such as nos_soil_layers, top_soil_depth .  Add `mVs` to `public` list.
4. Close the file and run the script on it, e.g. 
		  `python cardamom_model_type.py LIBRARY/CARDAMOM_F/model/DALEC.A1.C1.D2.F2.H2.P1.004/src/DALEC.A1.C1.D2.F2.H2.P1.004.f90` .
	- This inserts `mv%` in front of all variables in the `type model_working_variables` block and adds an argument `mV` to subroutine definitions and subroutine calls.
	- This creates a new file `<...>_editted`
5. Open the new `_editted` file and check for successful `mV` insertions. Move the `_editted` file to the original filename
6. Try compiling cardamom with the model of interest.    It will not compile , but since the model file is in the compile order first it should compile past the model file and get errors in some other file such as (now not compatible) MODEL_LIKELIHOOD file.
7. Fix common problems
   - ``do mV%soil_layer = 1, nos_soil_layers`` - here a loop counter `soil_layers` happened to have the same name as a variable in `mV` and was wrongly editted .  Change to a different loop counter variable name.
   - Errors related to `find_gs_iWUE` function:
	- *Inside* subroutine `calculate_stomatal_conductance`, insert a function definition
				  ```fortran
				   subroutine calculate_stomatal_conductance (mV)
				  ...
				   contains
				    double precision function find_gs_iWUE_(x)
				      double precision, intent(in):: x
				      find_gs_iWUE_ = find_gs_iWUE(x, mV)
				    
				  end function
				  ```
				- and use this variant name when passing the function to `zbrent` *only* (all other references to `find_gs_iWUE` in `calculate_stomatal_conductance` unchanged)
					  ```fortran
					                  mV%stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
					                                                find_gs_iWUE_, mV%minimum_conductance, mV%potential_conductance, tol_gs*mV%lai, mV%iWUE_step*0.10d0)
					  ```
				- All other calls to `find_gs_iWUE` should have the second argument `mV`, may need to add this
			- Errors related to ``water_retention_saxton_eqns`` : Do the same for `water_retention_saxton_eqns` inside `calculate_field_capacity`
				  ```fortran
				    subroutine calculate_field_capacity (mV)
				  
				      use brent_zero, only: zbrent
				  
				      ! field capacity calculations for saxton eqns !
				  
				      implicit none
				  
				        type(model_working_variables):: mV
				  
				      ! local variables..
				      integer:: i
				      double precision:: x1, x2
				  
				      x1 = 0.1d0; x2 = 0.7d0  ! low/high guess
				      do i = 1, nos_soil_layers+1
				         mV%water_retention_pass = i
				         ! field capacity is water content at which SWP = -10 kPa
				         mV%field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
				                                     water_retention_saxton_eqns_, x1, x2, 0.001d0, 0d0 )
				      enddo
				  
				      contains
				        double precision function water_retention_saxton_eqns_(x)
				          double precision, intent(in):: x
				          water_retention_saxton_eqns_ = water_retention_saxton_eqns(x, mV)
				        end function
				  
				    end subroutine calculate_field_capacity
				  ```
			- Same for any other functions passed to `zbrent` : only single-argument functions can be passed to zbrent, so we have to define a single-argument function as a wrapper around the functions with `mV` argument.
8. After `end type`, devlare an array of `model_working_variables` structs :
		  ```fortran
		  type(model_working_variables), allocatable, dimension(:):: mVs
		  ```
9. Make subroutine `initialize_mv`
		  ```fortran
		   subroutine initialize_mv(mV, nodays, nomet, nopars)
		      !! For a single chain's model_working_varibles type object mV, allocate arrays
		      !! and calculate initial values.  
		      use cardamom_structures, only: DATAin
		      implicit none
		      type(model_working_variables), intent(out):: mV
		      integer, intent(in):: nodays, nomet, nopars
		  
		      ! copy these arrays from global, read-only DATAin struct :
		      double precision:: deltat(nodays)     ! time step in decimal days
		      double precision:: met(nomet, nodays)  ! met drivers
		      double precision:: lat
		  
		      integer:: n
		  
		      deltat = DATAin%deltat
		      met = DATAin%met
		      lat = DATAin%lat
		          mV%soil_frac_sand = DATAin%soil_frac_sand
		          mV%soil_frac_clay = DATAin%soil_frac_clay
		  
		          ! allocate variables dimension which are fixed per site only the once
		          allocate(mV%deltat_1(nodays), mV%wSWP_time(nodays), mV%rSWP_time(nodays), mV%gs_demand_supply_ratio(nodays), &
		                   mV%gs_total_canopy(nodays), mV%gb_total_canopy(nodays), mV%canopy_par_MJday_time(nodays), &
		                   mV%soil_par_MJday_time(nodays), &
		                   mV%daylength_hours(nodays), mV%daylength_seconds(nodays), mV%daylength_seconds_1(nodays), &
		                   mV%rainfall_time(nodays), mV%cica_time(nodays), mV%root_depth_time(nodays), mV%snow_storage_time(nodays))
		  
		        !
		          ! Timing variables which are needed first
		          !
		  
		          mV%deltat_1 = deltat**(-1d0)
		  
		          !
		          ! Iteration independent variables using functions and thus need to be in a loop
		          !
		  
		          ! first those linked to the time period of the analysis
		          do n = 1, nodays
		             ! check positive values only for rainfall input
		             mV%rainfall_time(n) = max(0d0, met(7, n))
		             ! calculate daylength in hours and seconds
		             call calculate_daylength((met(6, n)-(deltat(n)*0.5d0)), lat, mV)
		             mV%daylength_hours(n) = mV%dayl_hours; mV%daylength_seconds(n) = mV%dayl_seconds
		          end do
		  
		       ! calculate inverse for each time step in seconds
		          mV%daylength_seconds_1 = mV%daylength_seconds ** (-1d0)
		          ! fraction of temperture period above freezing
		          mV%airt_zero_fraction_time = (met(3, :)-0d0) / (met(3, :)-met(2, :))
		  
		          ! number of time steps per year
		          mV%steps_per_year = nint(dble(nodays)/(sum(deltat)*0.002737851d0))
		          ! mean days per step
		          mV%mean_days_per_step = sum(deltat) / dble(nodays)
		  
		          !
		          ! Initialise the water model
		          !
		  
		          ! zero variables not done elsewhere
		          mV%total_water_flux = 0d0; mV%water_flux_mmolH2Om2s = 0d0
		          ! initialise some time invarient parameters
		          call saxton_parameters(mV%soil_frac_clay, mV%soil_frac_sand, mV)
		          call initialise_soils(mV%soil_frac_clay, mV%soil_frac_sand, mV)
		          ! save the initial conditions for later
		          mV%soil_waterfrac_initial = mV%soil_waterfrac
		          mV%SWP_initial = mV%SWP
		          mV%field_capacity_initial = mV%field_capacity
		          mV%porosity_initial = mV%porosity
		    end subroutine
		  
		  ```
		  with type like this example.  It contains
		- setting `deltat`, `met`, and `lat` which are needed in code below from `DATAin`
		- setting `mv%soil_frac_sand` and `mv%soil_frac_clay` arrays from DATAin
		- AND setup code moved here from `CARBON_MODEL`'s `if (.not.allocated(mV%deltat_1)) then` block, including allocation of all other arrays in this model's `model_working_variables` type
				- problem with `call update_soil_initial_conditions(pars(24), mV)` : delete this line from `initialize_mv`.  It belongs in the `else` branch of `CARBON_MODEL` - `if (.not.allocated(mV%deltat_1))` only .
9. After moving it to `initialize_mv`, remove the initialization code from `CARBON_MODEL`'s `if (.not.allocated(mV%deltat_1)) then` branch and delete the if/else .
10. Done when we get compile error about `MODEL_LIKELIHOOD.f90` file instead of model file
### likelihood / MODEL_LIKELIHOOD.f90 file
1. remove ` type (EDCDIAGNOSTICS), save:: EDCD` at the top of module
2. delete subroutine `find_edc_initial_values` , this is now in `_main_utils` .  Delete from `public` list.
3. add `edc_model_likelihood` to `public::` list
4. change all use ``MCMCOPT, only: PI`` to  `use model_shared, only: PI`
5. change all `use carbon_model_mod, only: carbon_model` to `use carbon_model_mod, only: carbon_model,  mVs`
6. add `thread_id` last argument and declaration `integer, intent(in), optional:: thread_id` to `edc_model_likelihood`
7. add `thread_id` last argument and declaration `integer, intent(in), optional:: thread_id` to `model_sanity_check`
8. add `thread_id` last argument and declaration `integer, intent(in), optional:: thread_id` to `model_likelihood`, `scaled_model_likelihood`
9. add internal variable `type (EDCDIAGNOSTICS):: EDCD` to `edc_model_likelihood`
10. add last argument `EDCD` and declaration `type (EDCDIAGNOSTICS), intent(inout):: EDCD` to  `assess_EDC2`
11. add `EDCD` last argument to calls to `assess_EDC2`
12. change `PI%parini` first argument of `call model_sanity_check()` to `PARS`
13. add `mVs(thread_id)` as last argument to all `call carbon_model`
14. in `model_sanity_check`, `model_likelihood`, `scaled_model_likelihood`, `edc_model_likelihood` add local variable declarations
```fortran
	double precision,dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
	double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
	double precision,dimension(datain%nodays, datain%nodiags)::  M_DIAGS
```
and pass these to `carbon_model`, `assess_EDC2` instead of `DATAin%M_FLUXES`, `DATAin%M_POOLS`, `DATAin%M_DIAGS`
  	
   - for error
```
		1052 |                      ,DATAin%M_FLUXES, DATAin%M_POOLS, DATAin%M_DIAGS &
		     |                      1
			Error: Variable ‘datain’ is PROTECTED and cannot appear in a variable definition context (actual argument to INTENT = OUT/INOUT) at (1)
```

  - also reference `M_DIAGS` instead of `DATAin%M_DIAGS` etc in the body of `model_sanity_check`
