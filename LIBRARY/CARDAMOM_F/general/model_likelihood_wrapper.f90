module model_likelihood_wrapper

! jklebes 2024
! This module wraps the module_likelihood subroutine found in 
! model/<chosen model>/likelihood/model_likelihood.f90 
! (whichever is chosen for compilation) to a standard-shaped function
! model_likelihood : params vector -> loglikelihood value
! that can be used as input to cardamom_samplers 

use model_likelihood_module, only: model_likelihood, &
                            find_edc_initial_values, &
          sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood!, log_model_likelihood_dtm
use iso_c_binding
implicit none

contains

function model_likelihood_fct(params) result(loglikelihood)
  double precision, intent(inout), dimension(:):: params  ! TODO should be intent out, reform model_likelihood
  double precision:: loglikelihood

  double precision:: ML_obs_out, ML_prior_out

  ! call the function to write to ML_obs_out and ML_prior_out
  call model_likelihood(params, ML_obs_out, ML_prior_out)

  ! for the purpose of running samplers we are only interested in the sum
  loglikelihood = ML_obs_out+ML_prior_out
end function

! For R : subroutine -> C void function, npars is given additionally
subroutine model_likelihood_fct2(params, npars, loglikelihood) bind(c, name="C_modellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(in), dimension(npars):: params
  real(c_double), dimension(npars):: params2
  real(c_double), intent(out):: loglikelihood

  double precision:: ML_obs_out, ML_prior_out
  params2 = params !TODO fix intent at model
  !TODO can we use an expected_npars from model files?
  !if (npars .neq. expected_npars) then
    !write(*,*) "Error : Passed ", npars, "from R (as indicated by second argument to modellikelihood), but this model takes a vector
    !of " , expected_npars, "values."
    ! TODO force exit ?
!else:
  ! call the function to write to ML_obs_out and ML_prior_out
  call model_likelihood(params2, ML_obs_out, ML_prior_out)

  ! for the purpose of running samplers we are only interested in the sum
  loglikelihood = ML_obs_out+ML_prior_out
!endif 
end subroutine

! TODO helpers for use from R-can we get a function to return/print 
! model name/code, expected_npars ?

! return the compiled model's expected unmber of parameters
subroutine get_npars(npars) bind(c, name="C_getmodelnpars")
  use iso_c_binding
  use model_parameters , only: pars_info
  use MCMCOPT, only: PI
  implicit none
  integer(c_int), intent(out)  :: npars
  call pars_info() !ideally this would bt written in an object, not a function in model_parameters
  npars =  PI%npars
end subroutine


! get the compiled model's parmin list for R
subroutine get_parmin(npars, parmin) bind(c, name="C_getmodelparmin")
  use iso_c_binding
  use model_parameters , only: pars_info
  use MCMCOPT, only: PI
  implicit none
  integer(c_int) , intent(in) :: npars
  real(c_double), dimension(npars), intent(out)  :: parmin
  call pars_info()
  parmin =  PI%parmin
end subroutine

subroutine get_parmax(npars, parmax) bind(c, name="C_getmodelparmax")
    use iso_c_binding
    use model_parameters , only: pars_info
    use MCMCOPT, only: PI
    implicit none
    integer(c_int) , intent(in) :: npars
    real(c_double), dimension(npars), intent(out)  :: parmax
    call pars_info()
    parmax =  PI%parmax
    end subroutine

subroutine inintialize_example_FI_Hyy() bind(c, name = "C_TMP_initialize")
  ! trigger read_binary_data(hard coded filename) for testing, 
  ! later do this better
  use cardamom_io, only: initialize
  implicit none
  call initialize("/home/jklebes/cardamom_profiling/FI-Hyy_example/DATA/FI-Hyy_example_FI-Hyy.bin")
end subroutine

end module
