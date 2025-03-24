module model_shared
  use samplers_shared, only: PARINFO  
  public:: PI 

  type(PARINFO), protected:: PI  ! should not be writted to except by pars_info() !  
                                ! If we want parallel runs 
                                ! now Contains read-only description of the model parameters only
  contains

  ! split from read_pari_data
  subroutine initialize_parinfo()
    use MODEL_PARAMETERS, only: pars_info
    implicit none 
    integer:: i

    ! load parameter max/min information, npars 
    call pars_info()

  ! Begin allocating parameter info
    if (.not. allocated(PI%parmin)) then 
         allocate(PI%parmin(PI%npars))
         PI%parmin = 0d0
    endif 
    if (.not. allocated(PI%parmax)) then 
        allocate(PI%parmax(PI%npars))
        PI%parmax = 0d0 
    endif
    allocate(&
            !PI%parini(PI%npars) &
            PI%parfix(PI%npars), & ! never used
            !PI%parvar(PI%npars), 
            PI%paradj(PI%npars) &
            !,PI%covariance(PI%npars, PI%npars), PI%mean_par(PI%npars) &
            !,PI%iC(PI%npars, PI%npars)&
            )


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
    ! PI%covariance = 0d0; PI%mean_par = 0d0; PI%cov = .false. ; PI%use_multivariate = .false.
    !do i = 1, PI%npars
    !    PI%covariance(i, i) = 1d0
    !end do
    ! report back to user
    !write(*,*) "Created field for parameter and covariances"
   

  end subroutine
  
end module
