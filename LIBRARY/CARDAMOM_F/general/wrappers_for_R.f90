module cardamom_Rinterfaces
   implicit none(external)
   public

contains
! helpers for use from R-
! TODO different module and file

! return the compiled model's expected unmber of parameters
   subroutine get_npars(npars) bind(c, name="C_getmodelnpars")
      use iso_c_binding
      !use model_parameters, only: pars_info
      use model_shared, only: PI
      implicit none(type, external)
      integer(c_int), intent(out)  :: npars
      !call pars_info(PI)
      npars = PI%npars
   end subroutine get_npars

! get the compiled model's parmin list for R
   subroutine get_parmin(npars, parmin) bind(c, name="C_getmodelparmin")
      use iso_c_binding
      !use model_parameters, only: pars_info
      use model_shared, only: PI
      implicit none(type, external)
      integer(c_int), intent(in):: npars
      real(c_double), dimension(npars), intent(out)  :: parmin
      !call pars_info(PI)
      parmin = PI%parmin
   end subroutine get_parmin

   subroutine get_parmax(npars, parmax) bind(c, name="C_getmodelparmax")
      use iso_c_binding
      !use model_parameters, only: pars_info
      use model_shared, only: PI
      implicit none(type, external)
      integer(c_int), intent(in):: npars
      real(c_double), dimension(npars), intent(out)  :: parmax
      !call pars_info(PI)
      parmax = PI%parmax
   end subroutine get_parmax

   subroutine initialize_cardamom(nchains) bind(c, name="C_initialize_model")
      use iso_c_binding
      use cardamom_io, only: initialize
      use model_shared, only: initialize_carbon_model
      implicit none(type, external)
      !character(kind = c_char, len = *), intent(in):: datain_filename
      character(kind=c_char, len=350):: filename
      integer(c_int), intent(in), optional:: nchains
      integer(c_int):: nchains_
      if (.not. present(nchains)) then
         nchains_ = 1_c_int
      else
         nchains_ = nchains
      end if
      !filename = trim(datain_filename)//c_null_char
      ! it turns out using string from R to C is hard so I'm hardcoding it here for now
      filename = "/home/jklebes/CARDAMOM/test/data/UK_baseline_sites_AliceHolt.bin"
      call initialize(filename)
      nchains_ = 4
      call initialize_carbon_model(nchains_)
   end subroutine initialize_cardamom

   subroutine initialize_stresstest_circle() bind(c, name="C_initialize_stresstest_circle")
      use model_shared, only: PI, initialize_parinfo
      use MHMCMC_StressTests
      implicit none(type, external)
      character(len=350):: infile, outfile
      infile = ""
      outfile = "Circle"
      call prepare_for_stress_test(infile, outfile)
   end subroutine initialize_stresstest_circle

end module cardamom_Rinterfaces
