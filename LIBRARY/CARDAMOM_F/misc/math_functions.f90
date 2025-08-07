!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
! assimilate observations and ecological theory to retrieve parameters for the 
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently. 
! Copyright (C) 2024  University of Edinburgh, 
!                     Mathew Williams (mat.williams@ed.ac.uk), 
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see < https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains functions needed to mathematical calculations in CARDAMOM
! 
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module math_functions

  !!!!!!!!!!!
  ! Module contains functions needed to mathematical calculations in CARDAMOM
  !!!!!!!!!!!

  implicit none

  ! assume default private
  private

  ! make explicit bits we want others to see
  public:: randn, &
            random_normal, random_uniform, &
            random_multivariate, &
            matrix_vector_func, &
            linear_model_gradient

  !!!!!!!!!!!
  ! Subroutines rand(), narray() and rnstrt() are from:
  !!!!!!!!!!!

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2000-09-10  Time: 16:37:48
  ! Latest revision-16 January 2003

  ! FORTRAN 77 version of "ran_array"
  ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
  !       including the MODIFICATIONS made in the 9th printing (2002)
  ! ********* see the book for explanations and caveats! *********
  ! Author: Steve Kifowit
  ! http://ourworld.compuserve.com/homepages/steve_kifowit
  ! with modifications by Alan Miller to rnarry and rnstrt based upon
  ! Knuth's code.

  ! For Donald Knuth's Fortran 77 versions, go to:
  ! http://www-cs-faculty.stanford.edu/~knuth/programs
  ! Look for frng.f and frngdb.f

  ! NOTE: that minimum number of values to be returned is 100
  !!!!!!!!!!!

  !!!!!!!!!!!
  ! Function randn()
  !!!!!!!!!!!

  ! Code is a modified version of that found in the uniform distribution generator from Numerical receipes
  ! Modified to give 0-1 unform on input of value 0 and normal distribution (mean = 0 sd = 1) on input of 1
  ! Modified by TLS

  !!!!!!!!!!!

  ! randn() related seed value
  double precision:: idum

  ! rand(), narray(), rnstrt() related
  integer, parameter  :: kk = 100, ll = 37, mm = 2**30, tt = 70, kkk = kk+kk-1
  integer, save       :: ranx(kk)

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine matrix_vector_func(uplo, n, alpha, A, lda, X, incx, beta, Y, incy)

    ! Performs the matrix-vector operation
    ! y := alpha*A*x+beta*y, 
    ! where alpha and beta are scalars, x and y are n element vectors and
    ! A is an n by n symmetric matrix.

    !  Arguments:
    !  ==========
    !
    ! intent(in):: UPLO
    !          UPLO is CHARACTER*1
    !          On entry, UPLO specifies whether the upper or lower
    !          triangular part of the array A is to be referenced as
    !          follows:
    !
    !              UPLO = 'U' or 'u'   Only the upper triangular part of A
    !                                  is to be referenced.
    !
    !              UPLO = 'L' or 'l'   Only the lower triangular part of A
    !                                  is to be referenced.!
    !
    ! intent(in):: N
    !           N is INTEGER
    !           On entry, N specifies the order of the matrix A.
    !           N must be at least zero.
    !
    ! intent(in):: ALPHA
    !           ALPHA is DOUBLE PRECISION.
    !           On entry, ALPHA specifies the scalar alpha.
    !
    ! intent(in):: A
    !           A is DOUBLE PRECISION array, dimension ( LDA, N )
    !           Before entry with  UPLO = 'U' or 'u', the leading n by n
    !           upper triangular part of the array A must contain the upper
    !           triangular part of the symmetric matrix and the strictly
    !           lower triangular part of A is not referenced.
    !           Before entry with UPLO = 'L' or 'l', the leading n by n
    !           lower triangular part of the array A must contain the lower
    !           triangular part of the symmetric matrix and the strictly
    !           upper triangular part of A is not referenced.
    !
    ! intent(in):: LDA
    !           LDA is INTEGER
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. LDA must be at least
    !           max( 1, n ).
    !
    ! intent(in):: X
    !           X is DOUBLE PRECISION array, dimension at least
    !           ( 1 + ( n-1 )*abs( INCX ) ).
    !           Before entry, the incremented array X must contain the n
    !           element vector x.
    !
    ! intent(in):: INCX
    !           INCX is INTEGER
    !           On entry, INCX specifies the increment for the elements of
    !           X. INCX must not be zero.
    !
    ! intent(in):: BETA
    !           BETA is DOUBLE PRECISION.
    !           On entry, BETA specifies the scalar beta. When BETA is
    !           supplied as zero then Y need not be set on input.
    !
    ! intent(inout):: Y
    !           Y is DOUBLE PRECISION array, dimension at least
    !           ( 1 + ( n-1 )*abs( INCY ) ).
    !           Before entry, the incremented array Y must contain the n
    !           element vector y. On exit, Y is overwritten by the updated
    !           vector y.
    !
    ! intent(in):: INCY
    !           INCY is INTEGER
    !           On entry, INCY specifies the increment for the elements of
    !           Y. INCY must not be zero.
    !
    !  Authors:
    !  ========
    !
    ! author Univ. of Tennessee
    ! author Univ. of California Berkeley
    ! author Univ. of Colorado Denver
    ! author NAG Ltd.
    !
    ! Date: December 2016
    !
    ! Further Details:
    ! =====================
    !
    !  Level 2 Basic Linear Algebra Subprograms (BLAS) routines.
    !  The vector and matrix arguments are not referenced when N = 0, or M = 0
    !
    !  -- Written on 22-October-1986.
    !     Jack Dongarra, Argonne National Lab.
    !     Jeremy Du Croz, Nag Central Office.
    !     Sven Hammarling, Nag Central Office.
    !     Richard Hanson, Sandia National Labs.
    !
    !  -- Code modifed for inclusion into CARDAMOM 27th June 2019
    !     T. Luke Smallman (UoE; t.l.smallman@ed.ac.uk)
    !
    !  =====================================================================
    !
    !  -- Reference BLAS level2 routine (version 3.7.0) --
    !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     December 2016

    ! Arguments
    double precision, intent(in):: alpha, beta
    integer, intent(in):: incx, incy, lda, n
    character, intent(in):: uplo
    double precision, intent(in):: X(n)
    double precision, intent (in):: A(lda, n)
    double precision, intent (inout):: Y(n)

    !
    !  =====================================================================
    !

    ! local variables
    double precision, parameter:: one = 1d0, zero = 0d0
    double precision:: tmp1, tmp2
    integer:: i, info, ix, iy, j, jx, jy, kx, ky

    ! Test the input parameters.

    info = 0
    if (trim(uplo) /= 'U' .and. trim(uplo) /= 'L') THEN
        info = 1
    else if (n <= 0) then
        info = 2
    else if (lda < max(1, n)) then
        info = 5
    else if (incx == 0) then
        info = 7
    else if (incy == 0) then
        info = 10
    end if
    if (info /= 0) then
        print*,"Inputs to matrix_vector_func not correct-error code = ",info
        stop
    end if

    !     Quick return if possible.
    if (alpha == zero .and. beta == one) return

    !
    !     Set up the start points in  X  and  Y.
    !
    if (incx > 0) then
        kx = 1
    else
        kx = 1 - (n-1)*incx
    end if
    if (incy > 0) then
        ky = 1
    else
        ky = 1 - (n-1)*incy
    end if

    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through the triangular part
    !     of A.
    !
    !     First form  y := beta*y.

    if (beta /= one) then
        if (incy == 1) then
            if (beta == zero) then
                do i = 1, n
                    y(i) = zero
                end do
            else
                DO i = 1, n
                    y(i) = beta*y(i)
                end do
            end if
        else
            iy = ky
            if (beta == zero) then
                do i = 1, n
                    y(iy) = zero
                    iy = iy+incy
                end do
            else
                do i = 1, n
                    y(iy) = beta*y(iy)
                    iy = iy+incy
                end do
            end if
        end if
    end if
    ! have done beta effect, no more change will occur if alpha is zero therefore return
    if (alpha == zero) return

    ! Now apply alpha component on calculation

    if (trim(uplo) == 'U') then

       ! Form  y  when A is stored in upper triangle.

       if ((incx .eq. 1) .and. (incy .eq. 1)) then

           do j = 1, n
              tmp1 = alpha*x(j)
              tmp2 = zero
              do i = 1, j-1
                 y(i) = y(i) + tmp1*A(i, j)
                 tmp2 = tmp2+a(i, j)*x(i)
              end do
              y(j) = y(j) + tmp1*A(j, j) + alpha*tmp2
           end do

       else
           jx = kx
           jy = ky
           do j = 1, n
              tmp1 = alpha*x(jx)
              tmp2 = zero
              ix = kx
              iy = ky
              do i = 1, j-1
                 y(iy) = y(iy) + tmp1*a(i, j)
                 tmp2 = tmp2+A(i, j)*x(ix)
                 ix = ix+incx
                 iy = iy+incy
              end do
              y(jy) = y(jy) + tmp1*A(j, j) + alpha*tmp2
              jx = jx+incx
              jy = jy+incy
           end do
       end if

    else  ! assume we must be in lower triangle

    ! Form  y  when A is stored in lower triangle.

       if (incx .eq. 1 .and. incy .eq. 1) then
           do j = 1, n
              tmp1 = alpha*x(j)
              tmp2 = zero
              y(j) = y(j) + tmp1*A(j, j)
              do i = j+1, n
                 y(i) = y(i) + tmp1*A(i, j)
                 tmp2 = tmp2+A(i, j)*x(i)
              end do
              y(j) = y(j) + alpha*tmp2
           end do
       else  ! incx .eq. 1 .and. incy .eq. 1
           jx = kx
           jy = ky
           do j = 1, n
              tmp1 = alpha*x(jx)
              tmp2 = zero
              y(jy) = y(jy) + tmp1*A(j, j)
              ix = jx
              iy = jy
              do i = j+1, n
                 ix = ix+incx
                 iy = iy+incy
                 y(iy) = y(iy) + tmp1*A(i, j)
                 tmp2 = tmp2+A(i, j)*x(ix)
              end do
              y(jy) = y(jy) + alpha*tmp2
              jx = jx+incx
              jy = jy+incy
           end do

       end if  ! incx .eq. 1 .and. incy .eq. 1

    end if ! (trim(uplo) == 'U')

    ! now return finally to the user
    return

  end subroutine matrix_vector_func
  !
  !--------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x, y, interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer:: interval  ! the total number of variables being regressed over
    double precision, dimension(interval):: x, y 

    ! declare local variables
    double precision:: sum_x, sum_y, sumsq_x, sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
    !                      / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )
    ! Linear regression done as single line to reduce assignment requirements
    linear_model_gradient = ( (dble(interval)*sum(x*y)) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sum(x*x)) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !------------------------------------------------------------------
  !
  double precision function randn(option)

    ! From Numerical Recipes p271 Press et al., 1986 2nd Edition Chapter 7, 
    ! Random Numbers function returns real random number between 0-1 based on an initial start
    ! point (ran1). The start point (default = -1) is reinitialised every time the model runs
    ! providing the same distribution each run. To ensure random numbers each time use the
    ! sum of the current system time.

    ! This code was modified based on CARDAMOM routines provided by A. A. Bloom, 
    ! available at github.com/CARDAMOM-framework/CARDAMOM_2.1.6c
    ! (contact abloom@jpl.nasa.gov for access)

    implicit none
    integer:: IA, IM, IQ, IR, NTAB, NDIV, option
    double precision:: AM, EPS, RNMX, const, r1, r2, pi
    parameter(IA = 16807, IM = 2147483647, AM = 1d0/dble(IM), IQ = 127773, &
      & IR = 2836, NTAB = 32, NDIV = 1+(IM-1)/NTAB, EPS = 1.2d-30, RNMX = 1d0-EPS)
    integer:: j, k, iv(NTAB), iy
    SAVE iv, iy
    DATA iv/NTAB*0/
    DATA iy/0/

    const = 1d0
    pi = 3.141592653589793d0

    if (option == 0) then
        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum, const)
            do j = (NTAB+8), 1, -1
               k = nint(idum/dble(IQ))
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum/dble(IQ))
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)

        ! output now
        randn = min(AM*dble(iy), RNMX)
        return

    else

        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum, const)
            do j = (NTAB+8), 1, -1
               k = nint(idum/dble(IQ))
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum)/IQ
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)
        r1 = max(min(AM*dble(iy), RNMX), 1d-30)

        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum, const)
            do j = (NTAB+8), 1, -1
               k = nint(idum)/IQ
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum)/IQ
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)
        r2 = max(min(AM*dble(iy), RNMX), 1d-30)

        ! output now
        randn = sqrt(-2d0*log(r1)) * cos(2d0*pi*r2)
        return

    endif

  end function randn
  !
  !--------------------------------------------------------------------
  !
  subroutine random_normal(uniform, random_length, uniform_random_vector, fn_val)

    ! Generate a random normal deviate using the polar method.
    ! Reference: Marsaglia, G. & Bray, T.A. 'A convenient method for
    ! generating normal variables',
    ! Siam Rev., vol.6, 260-264, 1964.
    ! Code modified from that created by Alan Miller
    ! (https://jblevins.org/mirror/amiller/rnorm.f90, last updated February 2004)

    implicit none

    ! arguments
    integer, intent(in):: random_length
    integer, intent(inout):: uniform
    double precision, intent(out):: fn_val
    double precision, dimension(random_length), intent(inout):: uniform_random_vector

    ! Local variables
    double precision:: u, sumsq
    double precision, save:: v, sln
    logical, save:: second = .false.
    double precision, parameter:: one = 1d0, vsmall = tiny( one )
    double precision, parameter:: sample_mean = 0d0, sample_std = 1d0

    if (second) then

        ! If second, use the second random number generated on last call
        second = .false.
        fn_val = sample_mean+sample_std*v * sln

    else

        ! First call; generate a pair of random normals
        second = .true. ; sumsq = 0d0
        do while (sumsq > one .or. sumsq == 0d0)
           u = uniform_random_vector(uniform); uniform = uniform+1
           v = uniform_random_vector(uniform); uniform = uniform+1
           u = u*2d0-one
           v = v*2d0-one
           sumsq = u*u + v*v
           if (uniform >= random_length) then
               call random_uniform(uniform_random_vector, random_length)
               uniform = 1
           endif
        end do
        sln = sqrt(-2d0*log(sumsq) / sumsq)
        fn_val = sample_mean+sample_std*u * sln
    end if

    return

  end subroutine random_normal
  !
  !--------------------------------------------------------------------
  !
  subroutine random_uniform(u, n)

    ! Generate an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)  :: n  ! number of random values wanted
    double precision, intent(out):: u(n)  ! output vector

    ! Local array
    integer, allocatable, dimension(:)  :: aa

    ! allocate memory
    allocate(aa(n))

    call rnarry(aa, n)
    u(1:n) = scale( dble(aa), -30)

    ! tidy
    deallocate(aa)

    return

  end subroutine random_uniform
  !
  !--------------------------------------------------------------------
  !
  subroutine random_multivariate ( m, n, u, u_size, uniform_random_vector, a, mu, x )

    !
    !  Discussion:
    !
    !    The multivariate normal distribution for the M dimensional vector X
    !    has the form:
    !
    !      pdf(X) = (2*pi*det(A))**(-M/2) * exp(-0.5*(X-MU)'*inverse(A)*(X-MU))
    !
    !    where MU is the mean vector, and A is a positive definite symmetric
    !    matrix called the variance-covariance matrix.
    !
    !  Licensing: This code is distributed under the GNU LGPL license.
    !
    !  Last Modified: Thu 07 Aug 2025 15:44:49 BST
    !
    !  Original Author: John Burkardt (07 December 2009)
    !
    !  Modified for coupling into CARDAMOM (03 May 2019):
    !    T. L. Smallman (t.l.smallman@ed.ac.uk)
    !
    !  Parameters:
    !
    !    Input, M, the dimension of the space.
    !    Input, N, the number of points.
    !    Input, A(M, M), the variance-covariance
    !    Matrix. A must be positive definite symmetric.
    !    Input,  MU(M), the mean vector.
    !
    !    Output, X(M), the points.
    !

    implicit none

    ! arguments
    integer, intent(in):: m, & ! number of parameters
                           n, & ! number of multivariate samples wanted per parameter
                      u_size    ! vector length for the uniform_random_vector
    integer, intent(inout):: u  ! current position within the random uniform vector
    double precision, intent(in):: a(m, m), mu(m)
    double precision, intent(out):: x(m, n)
    double precision, dimension(u_size), intent(inout):: uniform_random_vector

    ! local variables
    integer:: info, i, j
    double precision:: r(m, m)

    !
    !  Compute the upper triangular Cholesky factor R of the variance-covariance
    !  matrix.
    !

    ! Requires variance-covariance matrix, however as the matrix is over written, 
    ! make a duplicate
    r(1:m, 1:m) = a(1:m, 1:m)
    call cholesky_factor ( m, r, info )

    ! error checking from Cholesky Factor calculation
    if ( info /= 0 ) then

        ! Normal_multivariate-Fatal error!
        ! The variance-covariance matrix is not positive definite symmetric.
        ! Non-multivariate normal sample returned instead as default.'

        do i = 1, m  ! number of parameters searching for
           do j = 1, n  ! number of samples needed per parameter
              ! Draw random number from normal distribution
              call random_normal(u, u_size, uniform_random_vector, x(i, j))
              ! Updated variance based on variance from covariance matrix
              x(i, j) = x(i, j) * a(i, i)
           end do  ! j
        end do  ! i

        ! return early
        return

    end if

    !
    !  Get an MxN matrix of samples of the 1D normal distribution with mean 0
    !  and variance 1.
    !

    do i = 1, m
       do j = 1, n
          ! draw random number from normal distribution
          call random_normal(u, u_size, uniform_random_vector, x(i, j))
       end do  ! j
    end do  ! i

    !
    !  Compute R' * X.
    !  We actually carry out this computation in the equivalent form X' * R.
    !

    ! Whole multivariate estimated via matrix multplication
    ! Each desires set of combinations wanted it iterated
    do j = 1, n
       x(1:m, j) = mu(1:m) + matmul ( x(1:m, j), r(1:m, 1:m) )
    end do

    return

  end subroutine random_multivariate
  !
  !--------------------------------------------------------------------
  !
end module math_functions
