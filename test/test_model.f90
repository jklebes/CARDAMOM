module test_model
  !!Testing the example model DALEC.A1.C1.D2.F2.H2.P1
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use test_functions
  use test_math, only: approx
  use random_uniform  
  use model_shared, only: initialize_carbon_model, destroy_carbon_model
  use cardamom_MHMCMC
  use model_shared, only: PI
  use cardamom_io, only: initialize
  use cardamom_main_utils
  use OMP_LIB
  implicit none
  private

  public:: collect_modeltests

  character(len=*), parameter:: infile = "/home/jklebes/CARDAMOM/test/data/UK_baseline_sites_AliceHolt.bin"
  ! TODO hard coded path

contains

!> Collect all exported unit tests
subroutine collect_modeltests(testsuite)
  implicit none
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("model_initialize", test_model_initialize), &
    new_unittest("model_repeat_evaluation", test_model_repeat_evaluation), &
    new_unittest("model_group_repeat_evaluation", test_model_group_repeat_evaluation_3), &
    new_unittest("model_group_repeat_evaluation_3", test_model_group_repeat_evaluation) &
    ]

end subroutine collect_modeltests


subroutine test_model_initialize(error)
  use CARBON_MODEL_MOD, only: mvs
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains 
  type(MCMC_OUTPUT):: MCOUT
  nchains = 1 
    call initialize(infile) 
    call initialize_carbon_model(nchains)
    call initialize_stats(MCOUT, PI%npars)
    ! check initialized?
    call check(error, allocated(mVs(1)%rainfall_time))
    call destroy_carbon_model(nchains)
    call check(error, .not. allocated(mVs))
end subroutine test_model_initialize

subroutine test_model_repeat_evaluation(error)
  use CARBON_MODEL_MOD, only: mvs, carbon_model
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
  !! Similar to original cardamom "Sanity check" function
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains 
    double precision, dimension(PI%npars):: PARS
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: pools1, pools2
    double precision, dimension(DATAin%nodays, DATAin%nofluxes):: fluxes1, fluxes2
    double precision, dimension(DATAin%nodays, DATAin%nodiags):: diags1, diags2
    double precision:: pool_error, flux_error, diag_error
    ! These are here because they're kept at this level in model_likelihood.f90 files, 
    ! and they are there in DALEC models because this text insertion was the simplest way to edit 37 models
integer:: seed
    type(UNIF_VECTOR):: random_uniform
  nchains = 1
    call initialize(infile) 
    call initialize_carbon_model(nchains)
    PARS = DATAin%parpriors(1:PI%npars)
    seed = irand()  ! TODO record later
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes1, pools1, diags1,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(1))
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes2, pools2, diags2,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(1))
    flux_error = sum(abs(fluxes1-fluxes2))
    pool_error = sum(abs(pools1-pools2))
    diag_error = sum(abs(diags1-diags2))
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    call check(error, diag_error < .0000000001)
    call destroy_carbon_model(nchains)
end subroutine test_model_repeat_evaluation

subroutine test_model_repeat_and_independent_evaluation(error)
  !! Similar to original cardamom "Sanity check" function
  implicit none
  type(error_type), allocatable, intent(out):: error
end subroutine test_model_repeat_and_independent_evaluation


subroutine test_model_group_repeat_evaluation(error)
  use CARBON_MODEL_MOD, only: mvs, carbon_model
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
  !! Similar to original cardamom "Sanity check" function
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains 
    double precision, dimension(PI%npars):: PARS
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: pools1, pools2
    double precision, dimension(DATAin%nodays, DATAin%nofluxes):: fluxes1, fluxes2
    double precision, dimension(DATAin%nodays, DATAin%nodiags):: diags1, diags2
    double precision:: pool_error, flux_error, diag_error
integer:: seed, i, clock
    type(UNIF_VECTOR):: random_uniform
    nchains = 4
    call initialize(infile) 
    PARS = DATAin%parpriors(1:PI%npars)
    seed = irand()  
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    call initialize_carbon_model(nchains)
    do i = 1, nchains
    !call initialize_carbon_model(nchains)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes1, pools1, diags1,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    !call destroy_carbon_model(nchains)
    !call initialize_carbon_model(nchains)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes2, pools2, diags2,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    !call destroy_carbon_model(nchains)
    flux_error = sum(abs(fluxes1-fluxes2))
    pool_error = sum(abs(pools1-pools2))
    diag_error = sum(abs(diags1-diags2))
    write(*,*) i, flux_error, pool_error
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    call check(error, diag_error < .0000000001)
    end do
    call destroy_carbon_model(nchains)
end subroutine test_model_group_repeat_evaluation

subroutine test_model_group_repeat_evaluation_3(error)
  !! Similar to original cardamom "Sanity check" function
  !! This variant with 4 x 3 evaluations instead of 4 x 2 passes
  !! when group_repeat_evaluation doesn't, mysterious
  use CARBON_MODEL_MOD, only: mvs, carbon_model
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains 
    double precision, dimension(PI%npars):: PARS
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: pools1, pools2, pools3
    double precision, dimension(DATAin%nodays, DATAin%nofluxes):: fluxes1, fluxes2, fluxes3
    double precision, dimension(DATAin%nodays, DATAin%nodiags):: diags1, diags2, diags3
    double precision:: pool_error, flux_error, diag_error
integer:: seed, i
    type(UNIF_VECTOR):: random_uniform
    nchains = 4
    call initialize(infile) 
    call initialize_carbon_model(nchains)
    PARS = DATAin%parpriors(1:PI%npars)
    seed = irand()  ! TODO record later
    write(*,*) "seed", seed
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    do i = 1, nchains
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes1, pools1, diags1,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes2, pools2, diags2,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    flux_error = sum(abs(fluxes1-fluxes2))
    pool_error = sum(abs(pools1-pools2))
    diag_error = sum(abs(diags1-diags2))
    write(*,*) i, flux_error, pool_error
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    call check(error, diag_error < .0000000001)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes3, pools3, diags3,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    flux_error = sum(abs(fluxes3-fluxes2))
    pool_error = sum(abs(pools3-pools2))
    diag_error = sum(abs(diags1-diags2))
    write(*,*) i, flux_error, pool_error
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    end do
    call destroy_carbon_model(nchains)
end subroutine test_model_group_repeat_evaluation_3

subroutine test_model_cross_repeat_evaluation(error)
  !! Multiple threads that change threadid and evaluate the state again with mV arrays
  !! previously used by a different thread
  implicit none
  type(error_type), allocatable, intent(out):: error
end subroutine test_model_cross_repeat_evaluation

end module test_model
