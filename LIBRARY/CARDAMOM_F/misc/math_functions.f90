
module math_functions

  !!!!!!!!!!!
  ! Module contains functions needed to mathematical calculations in CARDAMOM
  !!!!!!!!!!!

  implicit none

  ! assume default private
  private

  ! make explicit bits we want others to see
  public :: randn, std, idum, covariance_matrix, &
            random_normal, random_uniform, rnstrt, &
            random_multivariate, increment_covariance_matrix, &
            par2nor, nor2par

  !!!!!!!!!!!
  ! Subroutines rand(), narray() and rnstrt() are from:
  !!!!!!!!!!!

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2000-09-10  Time: 16:37:48
  ! Latest revision - 16 January 2003

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
  double precision :: idum

  ! rand(), narray(), rnstrt() related
  integer, parameter  :: kk = 100, ll = 37, mm = 2**30, tt = 70, kkk = kk+kk-1
  integer, save       :: ranx(kk)

  contains

  !
  !--------------------------------------------------------------------
  !
  subroutine covariance_matrix(PARSALL,mean_par,npars,naccepted,covariance)

    ! Subroutine to estimate the covariance matrix
    ! Cov(X, Y) = Σ ( Xi - X ) ( Yi - Y ) / (N-1)
    ! X = mean for parameter 1
    ! Y = mean parameter 2
    ! Xi = ith member of the vector
    ! N = number of parameters accepted so far

    implicit none

    ! Arguments
    integer, intent(in) :: npars, naccepted
    double precision, intent(in) :: PARSALL(npars,naccepted)
    double precision, intent(out) :: mean_par(npars), covariance(npars,npars)

    ! local variables
    integer :: i
    double precision :: deviances(npars,naccepted)

    ! calculate components needed for covariance
    do i = 1, npars
       mean_par(i) = sum(PARSALL(i,:)) / dble(naccepted)
       ! estimate deviance
       deviances(i,:) = PARSALL(i,:) - mean_par(i)
    end do

    ! use matrix multiplication to estimate covariance
    ! NOTE: that naccepted-1 makes this the sample covariance
    covariance = matmul(deviances,transpose(deviances)) * dble(naccepted-1)**(-1)

    ! return to user
    return

  end subroutine covariance_matrix
  !
  !--------------------------------------------------------------------
  !
  subroutine increment_covariance_matrix(PARSALL,mean_par,npars,cur,new,covariance)

    ! Subroutine for incremental update of a the covariance matrix
    ! CMOUT = CM*(N-1)/(N-1+ar) + (N*M'*M-(N+ar)*Mi'*Mi+x'*x*ar)/(N-1+ar)
    ! M  = mean vector for parameters
    ! Mi = new mean vector for updated covariance_matrix
    ! ar = number of new parameters to be added
    ! N = number of parameters accepted so far

    implicit none

    ! Arguments
    integer, intent(in) :: npars, new
    double precision, intent(in) :: PARSALL(npars,new)
    double precision, intent(inout) :: cur, mean_par(npars), covariance(npars,npars)

    ! local variables
    integer :: n, i, j
    double precision :: new_mean_par(npars), nnew

    nnew = 1d0
    ! loop through each accepted parameter set...
    do n = 1, new
       ! ...estimate the new mean value for each parameter...
       new_mean_par = ((mean_par * cur) + (PARSALL(:,n) * nnew)) &
                    / (cur + nnew)
       ! ...update the covariance matrix with each new parameter vector in turn
       do i = 1, npars
          do j = 1, npars
             covariance(i,j) = covariance(i,j)*(cur-1d0)/(cur-1d0+nnew) &
                             + (cur*mean_par(i)*mean_par(j)- &
                               (cur+nnew)*new_mean_par(i)*new_mean_par(j) + &
                                nnew*PARSALL(i,n)*PARSALL(j,n))/(cur-1d0+nnew)
          end do ! j = 1, npars
       end do ! i = 1, npars
       ! update running totals and mean for the next iteration
       cur = cur + 1 ; mean_par = new_mean_par
    end do ! new_accepted

    ! return to user
    return

  end subroutine increment_covariance_matrix
  !
  !--------------------------------------------------------------------
  !
  double precision function std(a,n)

    ! Function to determine the standard deviation
    ! inputs are the vector of values and number of values included

    implicit none

    ! declare inputs
    integer, intent(in) :: n ! number of values in vector
    double precision, dimension(n), intent(in) :: a

    ! declare local variables
    integer :: i
    double precision :: mean, sq_diff_sum, diff, variance, sample

    ! if no length has been returned then provide value which ensures crash (i.e.
    ! infinity)
!    if (n == 0) then
!        std = 0d0
!        write(*,*) "no sample size has been provided to std function"
!        return
!    endif

    ! multiple use variable
    sample = dble(n)
    ! first calculate the mean
    mean = sum(a) / sample
    ! ensure zero values
    sq_diff_sum = 0d0
    ! calculate cumulative square difference
    do i = 1, n
       diff = a(i)-mean
       sq_diff_sum = sq_diff_sum + diff**2d0
    end do
    ! calculate the variance
    variance = sq_diff_sum / (sample-1d0)
    ! return the standard deviation
    std = sqrt(variance)

    return

  end function std
  !
  !------------------------------------------------------------------
  !
  subroutine par2nor(niter,initial_par,min_par,max_par,out_par)

    ! functions to normalised log parameter values and return them back to
    ! un-normalised value.

    ! converting parameters on log scale between 0-1 for min/max values
    implicit none
    integer, intent(in) :: niter     ! number of iterations in current vector
    double precision, intent(in) :: min_par, max_par
    double precision, dimension(niter), intent(in) :: initial_par
    double precision, dimension(niter), intent(out) :: out_par

    ! then normalise
    out_par = (initial_par-min_par)/(max_par-min_par)

    ! explicit return
    return

  end subroutine par2nor
  !
  !---------------------and vise versa ------------------------------
  !
  subroutine nor2par(niter,initial_par,min_par,max_par,out_par)

    ! Converting values back from normalised (0-1) to 'real' numbers

    implicit none
    integer, intent(in) :: niter     ! number of iterations in current vector
    double precision, intent(in) :: min_par, max_par
    double precision, dimension(niter), intent(in) :: initial_par
    double precision, dimension(niter), intent(out) :: out_par

    ! ...then un-normalise without logs as we cross zero and logs wont work
    out_par = min_par+(max_par-min_par)*initial_par

    ! explicit return
    return

  end subroutine nor2par
  !
  !------------------------------------------------------------------
  !
  double precision function randn(option)

    ! from Numerical Receipes p271 Press et al., 1986 2nd Edition Chapter 7,
    ! Random Numbers function returns real random number between 0-1 based on an initial start
    ! point (ran1).
    ! The start point (default = -1) is reinitialised every time the model runs
    ! providing the same distribution each run
    ! To ensure random numbers each time use the sum of the current system time

    ! modified based on blooms C code to alter range of random numbers

    implicit none
    integer :: IA,IM,IQ,IR,NTAB,NDIV,option
    double precision :: AM,EPS,RNMX,const,r1,r2,pi
    parameter(IA = 16807,IM = 2147483647,AM=1d0/dble(IM),IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-30,RNMX=1d0-EPS)
    integer :: j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/
    DATA iy /0/

    const = 1d0
    pi = 3.141592653589793d0

    if (option == 0) then
        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum,const)
            do j = (NTAB+8), 1, -1
               k = nint(idum/dble(IQ))
               idum=dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum + dble(IM)
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
        randn=min(AM*dble(iy), RNMX)
        return

    else

        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum,const)
            do j = (NTAB+8),1,-1
               k = nint(idum/dble(IQ))
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k=nint(idum)/IQ
        idum=dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum+dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)
        r1 = max(min(AM*dble(iy), RNMX),1d-30)

        if (idum < 0d0 .or. iy == 0) then
            idum = max(-idum,const)
            do j = (NTAB+8),1,-1
               k = nint(idum)/IQ
               idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
               if (idum < 0d0) idum = idum+dble(IM)
               if (j < NTAB) iv(j) = nint(idum)
            enddo
            iy = iv(1)
        endif
        k = nint(idum)/IQ
        idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
        if (idum < 0d0) idum = idum + dble(IM)
        j = 1+iy/NDIV
        iy = iv(j)
        iv(j) = nint(idum)
        r2 = max(min(AM*dble(iy), RNMX),1d-30)

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
    ! Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for
    ! generating normal variables',
    ! Siam Rev., vol.6, 260-264, 1964.

    implicit none

    ! arguments
    integer, intent(in) :: random_length
    integer, intent(inout) :: uniform
    double precision, intent(out) :: fn_val
    double precision, dimension(random_length), intent(inout) :: uniform_random_vector

    ! Local variables
    double precision :: u, sumsq
    double precision, save :: v, sln
    logical, save   :: second = .false.
    double precision, parameter :: one = 1d0, vsmall = tiny( one )
    double precision, parameter :: sample_mean = 0d0, sample_std = 1d0

    if (second) then

        ! If second, use the second random number generated on last call
        second = .false.
        fn_val = sample_mean + sample_std * v * sln

    else

        ! First call; generate a pair of random normals
        second = .true. ; sumsq = 0d0
        do while (sumsq > one .or. sumsq == 0d0)
           u = uniform_random_vector(uniform) ; uniform = uniform + 1
           v = uniform_random_vector(uniform) ; uniform = uniform + 1
           u = u * 2d0 - one
           v = v * 2d0 - one
           sumsq = u*u + v*v
           if (uniform >= random_length) then
               call random_uniform(uniform_random_vector,random_length)
               uniform = 1
           endif
        end do
        sln = sqrt(-2d0 * log(sumsq) / sumsq)
        fn_val = sample_mean + sample_std * u * sln
    end if

    return

  end subroutine random_normal
  !
  !--------------------------------------------------------------------
  !
  subroutine random_uniform(u, n)

    ! Generate an array of n double precision values between 0 and 1.

    integer, intent(in)  :: n ! number of random values wanted
    double precision, intent(out) :: u(n) ! output vector

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
    !  Last Modified: 03/05/2019
    !
    !  Original Author: John Burkardt
    !
    !  Modified for coupling into CARDAMOM:
    !    T. L. Smallman (t.l.smallman@ed.ac.uk)
    !
    !  Parameters:
    !
    !    Input, M, the dimension of the space.
    !    Input, N, the number of points.
    !    Input, A(M,M), the variance-covariance
    !    matrix. A must be positive definite symmetric.
    !    Input,  MU(M), the mean vector.
    !
    !    Output, X(M), the points.
    !

    implicit none

    ! arguments
    integer, intent(in) :: m, & ! number of parameters
                           n, & ! number of multivariate samples wanted per parameter
                      u_size    ! vector length for the uniform_random_vector
    integer, intent(inout) :: u ! current position within the random uniform vector
    double precision, intent(in) :: a(m,m), mu(m)
    double precision, intent(out) :: x(m,n)
    double precision, dimension(u_size), intent(inout) :: uniform_random_vector

    ! local variables
    integer :: info, i, j
    double precision :: r(m,m)

    !
    !  Compute the upper triangular Cholesky factor R of the variance-covariance
    !  matrix.
    !

    ! Requires variance-covariance matrix, however as the matrix is over written,
    ! make a duplicate
    r(1:m,1:m) = a(1:m,1:m)
    call cholesky_factor ( m, r, info )

    ! error checking from Cholesky Factor calculation
    if ( info /= 0) then

        ! Normal_multivariate - Fatal error!
        ! The variance-covariance matrix is not positive definite symmetric.
        ! Non-multivariate normal sample returned instead as default.'

        do i = 1, m ! number of parameters searching for
           do j = 1, n ! number of samples needed per parameter
              ! Draw random number from normal distribution
              call random_normal(u, u_size, uniform_random_vector, x(i,j))
              ! Updated variance based on variance from covariance matrix
              x(i,j) = x(i,j) * a(i,i)
           end do ! j
        end do ! i

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
          call random_normal(u, u_size, uniform_random_vector, x(i,j))
       end do ! j
    end do ! i

    !
    !  Compute R' * X.
    !  We actually carry out this computation in the equivalent form X' * R.
    !

    ! Whole multivariate estimated via matrix multplication
    ! Each desires set of combinations wanted it iterated
    do j = 1, n
       x(1:m,j) = mu(1:m) + matmul ( x(1:m,j), r(1:m,1:m) )
    end do

    return

  end subroutine random_multivariate
  !
  !--------------------------------------------------------------------
  !
  subroutine cholesky_factor ( n, a, info )

    !
    !  Discussion:
    !
    !    Subroutine calculates the symmetric positive definite
    !    matrix and its inverse. The Cholesky factor is an
    !    upper triangular matrix.
    !
    !    Only the diagonal and upper triangle of the square array are used.
    !    For clarity, the lower triangle is set to zero.
    !
    !    The positive definite symmetric matrix A has a Cholesky factorization
    !    of the form:
    !
    !      A = R' * R
    !
    !    where R is an upper triangular matrix with positive elements on
    !    its diagonal.  This routine overwrites the matrix A with its
    !    factor R.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Last Modified:
    !
    !    03/05/2019
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    !    FORTRAN90 version by John Burkardt.
    !    Modified for coupling to CARDAMOM by T. L. Smallman (t.l.smallman@ed.ac.uk)
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !  Parameters:
    !
    !    Input, integer N, the order of the matrix.
    !
    !    Input/output, double precision A(N,N).
    !    On output, the Cholesky factor R.
    !
    !    Output, integer INFO, error flag.
    !    0, normal return.
    !    K, error condition. The principal minor of order K is not
    !    positive definite, and the factorization was not completed.
    !

    implicit none

    ! arguments
    integer, intent(in) :: n
    integer, intent(out) :: info
    double precision, intent(inout) :: a(n,n)

    ! local arguments
    integer :: i, j, k
    double precision :: s

    ! Loop through the matrix along one dimension
    do j = 1, n

       ! doing the upper triangle only
       do k = 1, j - 1
          a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
       end do

       s = a(j,j) - sum ( a(1:j-1,j)**2 )

       ! error checking
       if ( s <= 0.0D+00 ) then
           info = j
           return
       end if

       ! final calculation of Cholesky
       a(j,j) = sqrt ( s )

    end do ! j = 1, n

    info = 0

    !
    !  Since the Cholesky factor is upper right corner only, be sure to
    !  zero out the lower triangle.
    !

    do i = 1, n
       do j = 1, i-1
          a(i,j) = 0.0D+00
       end do
    end do

    return

  end subroutine cholesky_factor
  !
  !--------------------------------------------------------------------
  !
  subroutine rnarry(aa, n)

    ! Generate an array of n integers between 0 and 2^30-1.

    integer, intent(in)   :: n
    integer, intent(out)  :: aa(n)

    ! Local variables
    integer  :: j

    aa(1:kk) = ranx(1:kk)
    do j = kk + 1, n
       aa(j) = aa(j-kk) - aa(j-ll)
       if (aa(j) < 0) aa(j) = aa(j) + mm
    end do
    do j = 1, ll
       ranx(j) = aa(n+j-kk) - aa(n+j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do
    do j = ll+1, kk
       ranx(j) = aa(n+j-kk) - ranx(j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do

    return

  end subroutine rnarry
  !
  !--------------------------------------------------------------------
  !
  subroutine rnstrt(seed)

    ! Initialize integer array ranx using the input seed.

    integer, intent(in)  :: seed

    ! Local variables
    integer  :: x(kkk), j, ss, sseed, t

    if (seed < 0) then
        sseed = mm - 1 - mod(-1-seed, mm)
    else
        sseed = mod(seed, mm)
    end if
    ss = sseed - mod(sseed,2) + 2
    do j = 1, kk
       x(j) = ss
       ss = ishft(ss, 1)
       if (ss >= mm) ss = ss - mm + 2
    end do
    x(kk+1:kkk) = 0
    x(2) = x(2)+1
    ss = sseed
    t = tt - 1
10  do j = kk, 2, -1
       x(j+j-1) = x(j)
    end do
    do j = kkk, kk + 1, -1
       x(j-(kk-ll)) = x(j-(kk-ll)) - x(j)
       if (x(j-(kk-ll)) < 0) x(j-(kk-ll)) = x(j-(kk-ll)) + mm
       x(j-kk) = x(j-kk) - x(j)
       if (x(j-kk) < 0) x(j-kk) = x(j-kk) + mm
    end do
    if (mod(ss,2) == 1) then
        do j = kk, 1, -1
           x(j+1) = x(j)
        end do
        x(1) = x(kk+1)
        x(ll+1) = x(ll+1) - x(kk+1)
        if (x(ll+1) < 0) x(ll+1) = x(ll+1) + mm
    end if
    if (ss /= 0) THEN
        ss = ishft(ss, -1)
    else
        t = t - 1
    end if
    if (t > 0) GO TO 10

    do j = 1, ll
       ranx(j+kk-ll) = x(j)
    end do
    do j = ll+1, kk
       ranx(j-ll) = x(j)
    end do

    do j = 1, 10
       call rnarry(x,kkk)
    end do

    return
  end subroutine rnstrt
  !
  !--------------------------------------------------------------------
  !
end module math_functions
