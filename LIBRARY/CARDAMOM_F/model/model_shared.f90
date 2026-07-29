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
! Module contains subroutines and variables common across the DALEC varients.
!
! This code was written by Jason Klebes (Jason.Klebes@ed.ac.uk, University of Edinburgh)
! Subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module model_shared
  use samplers_shared, only: PARINFO  
  use cardamom_structures, only: DATA_TYPE

  public

  type(PARINFO)  :: PI  ! should not be writted to except by pars_info() !  
                        ! type definition differs by model and lines in _PARS.f90 file
                        ! If we want parallel runs 
                        ! now Contains read-only description of the model parameters only

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine initialize_parinfo()
     use MODEL_PARAMETERS, only:  pars_info

     implicit none 

     ! split from read_pari_data

     ! Local variables
     integer:: i

     ! load parameter max/min information, npars 
     call pars_info(PI)

     ! Check allocation completed
     if (.not. allocated(PI%parmin) .or. .not. allocated(PI%parmax)) then
         write(*,*)"The pars_info() of the chosen model did not correctly allocate PI%parmax or PI%parmin"
         stop 
     end if

!     if (.not. allocated(PI%parmin)) then 
!         allocate(PI%parmin(PI%npars))
!         PI%parmin = 0d0
!     endif 
!     if (.not. allocated(PI%parmax)) then 
!         allocate(PI%parmax(PI%npars))
!         PI%parmax = 0d0 
!     endif
     ! Begin allocating parameter info
     if (.not. allocated(PI%parfix)) then 
         allocate(PI%parfix(PI%npars))
         PI%parfix = .false.
     endif
     if (.not. allocated(PI%paradj)) then 
         allocate(PI%paradj(PI%npars))
         PI%paradj = 0d0 
     endif

     ! force zero
     !PI%parini = 0d0  ! no longer exist-use DATAIN directly to -> pars in MCOUT
     !PI%parfix = 0d0; 
     ! PI%parvar = 0d0;  ! -> stats
     PI%paradj = 0d0  ! still exists, calcualted here
     !PI%covariance = 0d0; PI%iC = 0d0  ! -> stats TODO split some of these to initialize_stats
    
     ! For log-normalisation procedure, no parameter can be <= 0.
     ! To facilitate easy of setting parameter ranges to real values
     ! we here instead calculate the adjustment need to ensure positive only values
     where (PI%parmin <= 0d0) PI%paradj = abs(PI%parmin) + 1d0

     ! defining initial MHMCMC stepsize and standard deviation
     ! PI%parvar = 1d0; PI%Nparvar = 0d0
     ! Covariance matrix cannot be set to zero therefore set initial value to a
     ! small positive value along to variance access
     ! PI%covariance = 0d0; PI%meanpar = 0d0; PI%cov = .false. ; PI%use_multivariate = .false.
     !do i = 1, PI%npars
     !    PI%covariance(i, i) = 1d0
     !end do
     ! report back to user
     !write(*,*) "Created field for parameter and covariances"
   
  end subroutine initialize_parinfo
  !
  !------------------------------------------------------------------
  !
  subroutine initialize_carbon_model(n_chains)
     use cardamom_structures, only: DATAin
     use carbon_model_memory, only: mVs
     use carbon_model_mod, only: initialize_mv

     ! prepare N model_working_variables type objects to hold seperate sets of persistent values for
     ! each independent parallel chain.
     ! Must have DATAin filled first.

     ! Arguements
     integer, intent(in), optional:: n_chains

     ! Local variables
     integer:: n_chains_
     integer:: i

     ! Load number of chains as appropriate
     if (present(n_chains)) then
         n_chains_ = n_chains
     else
         n_chains_ = 1
     endif

     ! Allocate number of model memories for each chain
     allocate(mVs(n_chains))

     ! Initialise each chains memory
     do i = 1, n_chains_
         call initialize_mv(mVs(i), DATAin%nodays, DATAin%nomet, DATAin%nopars, DATAin%deltat, &
                            DATAin%soil_frac_sand, DATAin%soil_frac_clay, DATAin%met, DATAin%lat) 
     end do

  end subroutine initialize_carbon_model
  !
  !------------------------------------------------------------------
  !
  subroutine destroy_carbon_model()
     use carbon_model_memory, only: mVs

     ! deallocate members of model_working_variables struct(s) mVs
     if (allocated(mVs)) deallocate(mVs)

  end subroutine destroy_carbon_model
  !
  !------------------------------------------------------------------
  !  
end module model_shared
