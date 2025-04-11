module test_wrappers_for_R
  !! Helper functions for R scripts
  !! if they don't run from fortran they definitely won't run from R!
  !! Test here to see info on segfaults
  !! Behavior from R is hard to debug: fails intermittently and no info
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use cardamom_Rinterfaces
  implicit none
  private

  public:: collect_test_wrappers_for_R

contains

subroutine collect_test_wrappers_for_R(testsuite)
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("stresstest setup", test_initialize_stresstest_circle), &
    new_unittest("stresstest initialize twice", test_repeat_initialize) &
    ]

end subroutine collect_test_wrappers_for_R

subroutine test_initialize_stresstest_circle(error)
  !! this should be able to run first in a program, set up and allocate
  !! all objects needed in stresstest (PI, DATAin, MCOUT, ...)
  use model_shared, only: PI
  use cardamom_structures, only: DATAin
  type(error_type), allocatable, intent(out):: error
  call initialize_stresstest_circle()
  !! expect PI is initialized to expected values 
  call check(error, PI%parmin(1), -10d0)
  call check(error, PI%parmin(10), 1d0)
  call check(error, PI%parmax(3), 300d0)
  call check(error, PI%parmax(9), 300d0)
  ! datain was set, and "Circle" routine was selected
  call check(error, DATAin%noobs, 9)
end subroutine 


subroutine test_repeat_initialize(error)
  !! check we don't get crashes when called initialize_stresstest_circle twice in a R session
  type(error_type), allocatable, intent(out):: error
  call initialize_stresstest_circle()
  call initialize_stresstest_circle()
end subroutine 

end module test_wrappers_for_R
