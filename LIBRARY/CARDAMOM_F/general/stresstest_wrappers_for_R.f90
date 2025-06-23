module cardamom_Rinterfaces
use iso_c_binding
use MHMCMC_StressTests
contains
! helpers for use from R-


subroutine initialize_stresstest_circlex() bind(c, name = "C_initialize_stresstest_circle")
  use model_shared, only: PI, initialize_parinfo
  implicit none
  character(len = 350):: infile, outfile
  infile = ""
  outfile = "Circle"
  call prepare_for_stress_test(infile, outfile)
end subroutine

end module
