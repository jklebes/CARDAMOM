module cardamom_Rinterfaces

contains
! helpers for use from R-
! TODO different module and file

! return the compiled model's expected unmber of parameters
subroutine get_npars(npars) bind(c, name="C_getmodelnpars")
  use iso_c_binding
  use model_parameters, only: pars_info
  use model_shared, only: PI
  implicit none
  integer(c_int), intent(out)  :: npars
  call pars_info(PI)
  npars =  PI%npars
end subroutine


! get the compiled model's parmin list for R
subroutine get_parmin(npars, parmin) bind(c, name="C_getmodelparmin")
  use iso_c_binding
  use model_parameters, only: pars_info
  use model_shared, only: PI
  implicit none
  integer(c_int), intent(in):: npars
  real(c_double), dimension(npars), intent(out)  :: parmin
  call pars_info(PI)
  parmin =  PI%parmin
end subroutine

subroutine get_parmax(npars, parmax) bind(c, name="C_getmodelparmax")
    use iso_c_binding
    use model_parameters, only: pars_info
    use model_shared, only: PI
    implicit none
    integer(c_int), intent(in):: npars
    real(c_double), dimension(npars), intent(out)  :: parmax
    call pars_info(PI)
    parmax =  PI%parmax
    end subroutine

subroutine inintialize_example_FI_Hyy() bind(c, name = "C_TMP_initialize")
  ! trigger read_binary_data(hard coded filename) for testing, 
  ! later do this better
  use cardamom_io, only: initialize
  implicit none
  character(len = 350)  :: filename
  filename = "/home/jklebes/cardamom_profiling/FI-Hyy_example/DATA/FI-Hyy_example_FI-Hyy.bin"
  ! TODO add results
  call initialize(filename)
end subroutine
end module
