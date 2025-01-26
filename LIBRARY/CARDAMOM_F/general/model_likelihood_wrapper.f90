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

! For R : npars is given additionally
function model_likelihood_fct2(params, npars) result(loglikelihood) bind(c, name="modellikelihood")
  integer, intent(in)  :: npars
  double precision, intent(inout), dimension(npars):: params
  double precision:: loglikelihood

  double precision:: ML_obs_out, ML_prior_out
  !TODO can we use an expected_npars from model files?
  !if (npars .neq. expected_npars) then
    !write(*,*) "Error : Passed ", npars, "from R (as indicated by second argument to modellikelihood), but this model takes a vector
    !of " , expected_npars, "values."
    ! TODO force exit ?
!else:
  ! call the function to write to ML_obs_out and ML_prior_out
  call model_likelihood(params, ML_obs_out, ML_prior_out)

  ! for the purpose of running samplers we are only interested in the sum
  loglikelihood = ML_obs_out+ML_prior_out
!endif 
end function

! TODO helpers for use from R-can we get a function to return/print 
! model name/code, expected_npars ?

end module
