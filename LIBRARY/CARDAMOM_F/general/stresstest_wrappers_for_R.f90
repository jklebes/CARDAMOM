module cardamom_Rinterfaces
use iso_c_binding
use MHMCMC_StressTests
contains
! helpers for use from R-


subroutine initialize_stresstest_circle() bind(c, name = "C_initialize_stresstest_circle")
  use model_shared, only: PI, initialize_parinfo
  implicit none
  character(len = 350):: infile, outfile
  infile = ""
  outfile = "Circle"
  call prepare_for_stress_test(infile, outfile)
end subroutine

subroutine get_stresstest_parmax(npars, parmax) bind(c, name="C_getstresstestparmax")
  use model_shared, only: PI
  implicit none
  integer(c_int), intent(in):: npars
  real(c_double), intent(inout), dimension(npars):: parmax
  parmax = PI%parmax
end subroutine

subroutine get_stresstest_parmin(npars, parmin) bind(c, name="C_getstresstestparmin")
  use model_shared, only: PI
  implicit none
  integer(c_int), intent(in):: npars
  real(c_double), intent(inout), dimension(npars):: parmin
  parmin = PI%parmin
end subroutine

end module
