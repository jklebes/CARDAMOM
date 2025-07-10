module cardamom_main_utils
   implicit none
   public

contains

   subroutine initialize_stats(MCOUT, npars)
      use MHMCMC, only: MCMC_OUTPUT
      integer, intent(in):: npars
      type(MCMC_OUTPUT), intent(inout):: MCOUT
      allocate (MCOUT%covariance(npars, npars), MCOUT%parvar(npars), MCOUT%meanpar(npars))
      call reset_stats(MCOUT, npars)
   end subroutine
   subroutine reset_stats(MCOUT, npars)
      use MHMCMC, only: MCMC_OUTPUT
      type(MCMC_OUTPUT), intent(inout):: MCOUT
      integer, intent(in):: npars
      integer:: n
      MCOUT%parvar = 1d0; MCOUT%Nparvar = 0d0
      ! Covariance matrix cannot be set to zero therefore set initial
      ! value to a small positive value along to variance access
      MCOUT%covariance = 0d0; MCOUT%meanpar = 0d0; MCOUT%cov = .false.
      MCOUT%use_multivariate = .false.
      do n = 1, npars
         MCOUT%covariance(n, n) = 1d0
      end do
   end subroutine
end module
