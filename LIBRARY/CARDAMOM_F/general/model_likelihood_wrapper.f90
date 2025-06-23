module model_likelihood_wrapper

! jklebes 2024
! This module wraps the module_likelihood subroutine found in 
! model/<chosen model>/likelihood/model_likelihood.f90 
! (whichever is chosen for compilation) to a standard-shaped function
! model_likelihood : params vector -> loglikelihood value
! that can be used as input to cardamom_samplers 

use model_likelihood_module, only: model_likelihood, &
                            edc_model_likelihood, & 
          sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood!, log_model_likelihood_dtm
use iso_c_binding
implicit none

contains


! model_likelihood, wrapped to be 
! shaped like the generic function to hand to both cardamom_samplers and (via C binding) R
! For R : has to be subroutine ( -> C void function), have to give npars 
! 
subroutine model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_modellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out
  integer(c_int), intent(in), optional:: id

  ! call the function to write to ML_obs_out and ML_prior_out
  call model_likelihood(params, ML_obs_out, ML_prior_out, id)

  ! for the purpose of running samplers we are only interested in the sum
  loglikelihood = ML_obs_out+ML_prior_out

end subroutine

! variants
subroutine log_model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_logmodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out
  integer(c_int), intent(in), optional:: id

  call log_model_likelihood(params, ML_obs_out, ML_prior_out, id)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine

subroutine sub_model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_submodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out
  integer(c_int), intent(in), optional:: id

  call sub_model_likelihood(params, ML_obs_out, ML_prior_out, id)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine

subroutine sqrt_model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_sqrtmodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out
  integer(c_int), intent(in), optional:: id

  call sqrt_model_likelihood(params, ML_obs_out, ML_prior_out, id)
  loglikelihood = ML_obs_out+ML_prior_out
end subroutine

! the function for starting loops, which finds a set of parameters fulfilling otherwise
! hard boundary conditions with a softer potential
subroutine edc_model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_edcmodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out
  integer(c_int), intent(in), optional:: id

  call edc_model_likelihood(params, ML_obs_out, ML_prior_out, id)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine

end module
