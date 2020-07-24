!!
!> @brief      modern versions of numerical recipes routines
!!
!! Collection of updated versions of subroutines and functions from the numerical recipes book
!!
!! @author     Rodrigo Navarro Perez
!!
module num_recipes

use precisions, only : dp
use constants, only : pi, i_

implicit none

private

public :: dfridr, context, func, sphbes, bessjy, cmplx_log_gamma, legendre_poly, spherical_harmonic

!!
!> @brief      generic type for data in function callbacks
!!
!! The idea is for this type to be generic enough to be used by
!! different functions
!!
!! @author     Rodrigo Navarro Perez
!!
type context
    real(dp) :: a !< a real number
    real(dp) :: b !< a real number
    real(dp) :: c !< a real number
    real(dp) :: d !< a real number
    real(dp) :: e !< a real number
    integer :: i !< an integer
    integer :: j !< an integer
    integer :: k !< an integer
    integer :: l !< an integer
    real(dp), allocatable :: x(:) !< an array of reals
    integer, allocatable :: i_x(:) !< an array of integers
    character(len=1024) :: string !< a string
    logical :: log !< a logical variable
    procedure(), pointer, nopass :: proc
end type

!!
!> @brief      interface for function callbacks
!!
!! @author     Rodrigo Navarro Perez
!!
interface
    real(dp) function func(x, data)
        use precisions, only: dp
        import context
        implicit none
        real(dp), intent(in) :: x !< point at which a function will be evaluated
        type(context), intent(in) :: data !< contains data for function evaluation
    end function
end interface



contains

!!
!> @brief      Logarithm of gamma function
!!
!! Series approximation to the Logarithm of a gamma function for a complex argument.
!!
!! \f$ \log \left[ \Gamma (z) \right] \f$, where \f$ z \f$ is complex
!!
!! @author     Rodrigo Navarro Perez
!!
complex(dp) function cmplx_log_gamma(z) result(l_gamma)
    implicit none
    complex(dp), intent(in) :: z !< a complex number

    integer, parameter :: n_coeff = 6
    real(dp), parameter :: lanczos_coefficients(0:n_coeff) = [1.000000000190015_dp, &
        76.18009172947146_dp, -86.50532032941677_dp, 24.01409824083091_dp, -1.231739572450155_dp, &
        0.1208650973866179e-2_dp, -0.5395239384953e-5_dp]
    integer, parameter :: gamma = 5
    complex (dp) :: sum
    integer :: i

    if (real(z) < 0) then
        stop 'cmplx_log_gamma is only valid for Re(z) >= 0'
    endif

    sum = lanczos_coefficients(0)
    do i = 1, n_coeff
        sum = sum + lanczos_coefficients(i)/(z + i)
    enddo

    l_gamma = (z + 0.5_dp)*log(z + gamma + 0.5_dp) - (z + gamma + 0.5_dp) + log(sqrt(2*pi)*sum/z)
end function cmplx_log_gamma

!!
!> @brief      Spherical harmonic function
!!
!! Calculates the spherical harmonic
!! \f$Y_l^m(\theta, \varphi) = (-1)^m \sqrt{\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}} P_{lm}(\cos \theta) e^{im\varphi}\f$,
!! where \f$ P_{lm}(x) \f$ is the usual Legendre polynomial.
!!
!! @author     Rodrigo Navarro Perez
!!
complex(dp) function spherical_harmonic(l, m, theta, phi) result(Ylm)
    implicit none
    integer, intent(in) :: l !< orbital momentum quantum number
    integer, intent(in) :: m !< magnetic quantum number
    real(dp), intent(in) :: theta !< polar angle
    real(dp), intent(in) :: phi !< azimuthal angle

    real(dp) :: ratio, factorial, x, Plm
    integer :: lmm, i

    if (m == 0) then
        ratio = 1._dp
    else
        lmm = l - abs(m)
        factorial = lmm + 1
        do i = 2, 2*abs(m)
            factorial = factorial*(lmm + i)
        enddo
        ratio = 1/factorial
    endif
    x = cos(theta)
    Plm = legendre_poly(l, abs(m), x)
    Ylm = sqrt((2*l + 1)*ratio/(4*pi))*Plm*exp(i_*m*phi)
    if (m < 0) Ylm = (-1)**m*Ylm
end function spherical_harmonic


!!
!> @brief      Legendre polynomial
!!
!! Uses a recursive relation to calculate the usual Legendre polynomial
!!
!! The algorithm has been modified to save
!! the values from previous calls and build up from them. This eliminates the
!! need to start from the beginning when using same x and m values 
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
real(dp) function legendre_poly(l, m, x) result(r)
    implicit none
    integer, intent(in) :: l !< \f$l\f$ index. Has to be smaller than \f$m\f$ and 2000 
    integer, intent(in) :: m !< \fm\f$ index. Has to bee positive
    real(dp), intent(in) :: x !< a real number between \f$-1\f$ and \f$1\f$
    !------- place holder variables--------------------------
    real(dp), save :: x_pre = -1.1
    integer, save :: m_pre = -1, l_largest = -1
    integer, parameter :: L_MAX = 2000 ! max number of polynomial iterations
    real(dp), save, dimension(0:L_MAX) :: pl_array = 0 ! initialize to 0
    !-------------------------------------------------------
    integer ::  i
    real(dp) :: p_mm, factor, odd_factorial, p_mmp1, p_ml

    ! set
    if (m < 0) stop 'legendre_poly is only valid for m >= 0'
    if (l < m) stop 'legendre_poly is only valid for l >= m'
    if (l > L_MAX) stop 'legendre_poly limit of 2000 reached'
    if (abs(x) > 1) stop 'legendre_poly is only valid for -1 <= x <= 1'
    p_ml=0
    if (x==x_pre .and. m==m_pre) then
        !Same arguments were used last time, we can use previous calculations
        if (l <= l_largest) then
            !The polynomial for this l has been calculated, use it
            r = pl_array(l)
        else
            ! The The polynomial for this l has NOT been calculated
            if (l == m+1) then
                ! There's only one previous calculation
                p_mm = pl_array(l-1)
                p_ml = x*(2*m + 1)*p_mm
                pl_array(l) = p_ml
            else
                ! build up from previous values
                do i = l_largest + 1, l
                    p_mm = pl_array(i-2)
                    p_mmp1 = pl_array(i-1)
                    p_ml = (x*(2*i - 1)*p_mmp1 - (i + m - 1)*p_mm)/(i-m) ! build l+1 using l-2 and l-1
                    pl_array(i) = p_ml ! save value for future use
                enddo
            endif
            r = p_ml ! set final l
            l_largest = l
        endif
    else
        !Arguments are different, we need to start from scracth
        pl_array = 0
        p_mm = 1._dp
        if (m > 0) then
            factor = sqrt(1 - x**2)
            odd_factorial = 1
            do i = 1, m
                p_mm = -p_mm*odd_factorial*factor
                odd_factorial = odd_factorial + 2
            enddo
        endif
        if (l == m) then
            r = p_mm
            ! save result
            pl_array(l) = p_mm
        else
            p_mmp1 = x*(2*m + 1)*p_mm
            if (l == m + 1) then
                r = p_mmp1
                ! save result
                pl_array(l) = p_mmp1
            else
                p_ml = p_mmp1
                do i = m + 2, l
                    p_ml = (x*(2*i - 1)*p_mmp1 - (i + m - 1)*p_mm)/(i-m) ! build l+1 using l-2 and l-1
                    pl_array(i) = p_ml ! save value for future use
                    p_mm = p_mmp1 ! set l-2 to l-1
                    p_mmp1 = p_ml ! set l-1 to l+1
                    ! repeat until you reach l
                enddo
                r = p_ml ! set final l
            endif
        endif
        l_largest = l
        x_pre = x
        m_pre = m
    endif
    ! THIS IS THE ORIGINAL UNMODIFIED ALGORITHM
    ! p_mm = 1._dp
    ! if (m > 0) then
    !     factor = sqrt(1 - x**2)
    !     odd_factorial = 1
    !     do i = 1, m
    !         p_mm = -p_mm*odd_factorial*factor
    !         odd_factorial = odd_factorial + 2
    !     enddo
    ! endif
    ! if (l == m) then
    !     r = p_mm ! set l-2
    ! else
    !     p_mmp1 = x*(2*m + 1)*p_mm
    !     if (l == m + 1) then
    !         r = p_mmp1 ! set l-1
    !     else
    !         p_ml = p_mmp1 ! set current to l-1
    !         do i = m + 2, l
    !             p_ml = (x*(2*i - 1)*p_mmp1 - (i + m - 1)*p_mm)/(i-m) ! build l+1 using l-2 and l-1
    !             p_mm = p_mmp1 ! set l-2 to l-1
    !             p_mmp1 = p_ml ! set l-1 to l+1
    !             ! repeat until you reach l
    !         enddo
    !         r = p_ml ! set final l
    !     endif
    ! endif
end function legendre_poly

!!
!> @brief      numerical derivative of a function
!!
!! Returns the derivative of a function f at a point x by Ridders’ method of polynomial
!! extrapolation. The value h is input as an estimated initial step-size; it need not be small,
!! but rather should be an increment in x over which f changes substantially. An estimate
!! of the error in the derivative is returned as err.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine dfridr(f, data, x, h, df_dx, err)
    implicit none
    procedure(func) :: f !< function whose derivative will be calculated
    type(context), intent(in) :: data !< context type containing the data necessary to evaluate f
    real(dp), intent(in) :: x !< point at which the derivative will be evaluated
    real(dp), intent(in) :: h !< estimated initial step-size. should be an increment in x over which f changes substantially
    real(dp), intent(out) :: df_dx !< derivative of f at x
    real(dp), intent(out) :: err !< estimate of the error in the derivative

    integer, parameter  :: ntab = 10
    real(dp), parameter :: con = 1.4_dp, con2 = con**2, big = 1.e+30_dp, safe = 2._dp
    integer :: i, j
    real(dp) :: errt,fac
    real(dp) :: hh
    real(dp) :: a(1:ntab, 1:ntab)

    if (h == 0._dp) then
        stop 'h must be nonzero in dfridr'
    endif

    hh = h
    a(1,1) = (f(x + hh, data) - f(x - hh, data))/(2*hh)
    err = big

    do i = 2, ntab
        hh = hh/con
        a(1,i) = (f(x + hh, data) - f(x - hh, data))/(2*hh)
        fac = con2
        do j = 2, i
            a(j, i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac - 1)
            fac = con2*fac
            errt = max(abs(a(j,i) - a(j-1,i)), abs(a(j,i) - a(j-1,i-1)))
            if (errt <= err) then
                err = errt
                df_dx = a(j,i)
            endif
        enddo
        if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) return
    enddo
end subroutine dfridr


!!
!> @brief      spherical Bessel functions
!!
!! Modified from Numerical Recipes
!!
!! Calculates the spherical Bessel functions \f$ j_n(x) \f$, \f$ y_n(x) \f$, and their derivatives
!! \f$ j_n'(x) \f$, \f$ y_n'(x) \f$ for an integer \f$ n \f$.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine sphbes(n, x, sj, sy, sjp, syp)
    implicit none
    integer, intent(in) :: n !< order of the spherical Bessel functions.
    real(dp), intent(in) :: x !< point at which the Bessel functions are evaluated.
    real(dp), intent(out) :: sj !< regular spherical Bessel function \f$ j_n(x) \f$.
    real(dp), intent(out) :: sy !< irregular spherical Bessel function \f$ y_n(x) \f$.
    real(dp), intent(out) :: sjp !< derivative of the regular spherical Bessel function \f$ j_n'(x) \f$.
    real(dp), intent(out) :: syp !< derivative of the irregular spherical Bessel function \f$ y_n'(x) \f$.

    real(dp) :: factor, order, rj, rjp, ry, ryp

    if ( n < 0 .or. x <= 0 ) stop 'bad arguments in sphbes'
    order = n + 1/2._dp
    call bessjy(x, order, rj, ry, rjp, ryp)
    factor = sqrt(pi/(2*x))
    sj = factor*rj
    sy = factor*ry
    sjp = factor*rjp - sj/(2*x)
    syp = factor*ryp - sy/(2*x)
end subroutine sphbes


!!
!> @brief      Bessel functions of fractional order
!!
!! Modified from Numerical Recipes
!!
!! Returns the Bessel functions \f$ J_\nu \f$ , \f$ Y_\nu \f$ and their derivatives \f$ J_\nu' \f$,
!! \f$ Y_\nu' \f$ for positive x and for \f$ \nu \geq 0 \f$.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine bessjy(x, nu, rj, ry, rjp, ryp)
    implicit none
    real(dp), intent(in) :: x !< point at which the Bessel functions are evaluated
    real(dp), intent(in) :: nu !< Order of the Bessel functions
    real(dp), intent(out) :: rj !< regular spherical Bessel function \f$ J_\nu(x) \f$.
    real(dp), intent(out) :: ry !< irregular spherical Bessel function \f$ Y_\nu(x) \f$.
    real(dp), intent(out) :: rjp !< derivative of regular spherical Bessel function \f$ J_\nu'(x) \f$.
    real(dp), intent(out) :: ryp !< derivative of irregular spherical Bessel function \f$ Y_\nu'(x) \f$.

    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 1.e-16_dp, fpmin = 10*tiny(1._dp), xmin = 2._dp

    integer ::  i, isign, l, nl
    real(dp) :: a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2, &
        fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl, rjl1, rjmu, rjp1, &
        rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1, temp, w, x2, xi, xi2, xmu, xmu2

    if (x <=0 .or. nu <= 0) stop 'bad arguments in bessjy'
    if (x <= xmin ) then
        nl = int(nu + 0.5_dp)
    else
        nl = max(0, int(nu - x + 1.5_dp))
    endif
    xmu = nu - nl
    xmu2 = xmu*xmu
    xi = 1/x
    xi2 = 2*xi
    w = xi2/pi
    isign = 1
    h = nu*xi
    if (h < fpmin) h = fpmin
    b = xi2*nu
    d = 0
    c = h
    do i = 1, maxit
        b = b + xi2
        d = b - d
        if (abs(d) < fpmin) d = fpmin
        c = b - 1/c
        if (abs(c) < fpmin) c = fpmin
        d = 1/d
        del = c*d
        h = del*h
        if (d < 0) isign = -isign
        if (abs(del - 1) < eps) exit
    enddo
    if (i >= maxit) stop 'x too large in bessjy; try asymptotic expansion'
    rjl = isign*fpmin
    rjpl = h*rjl
    rjl1 = rjl
    rjp1 = rjpl
    fact = nu*xi
    do l = nl, 1, -1
        rjtemp = fact*rjl + rjpl
        fact = fact - xi
        rjpl = fact*rjtemp - rjl
        rjl = rjtemp
    enddo
    if (rjl == 0._dp) rjl = eps
    f = rjpl/rjl
    if (x < xmin) then
        x2 = x/2
        pimu = pi*xmu
        if (abs(pimu) < eps)then
            fact = 1
        else
            fact = pimu/sin(pimu)
        endif
        d = -log(x2)
        e = xmu*d
        if (abs(e) < eps)then
            fact2 = 1
        else
            fact2 = sinh(e)/e
        endif
        call beschb(xmu, gam1, gam2, gampl, gammi)
        ff = 2/pi*fact*(gam1*cosh(e) + gam2*fact2*d)
        e = exp(e)
        p = e/(gampl*pi)
        q = 1/(e*pi*gammi)
        pimu2 = pimu/2
        if (abs(pimu2) < eps)then
            fact3 = 1
        else
            fact3 = sin(pimu2)/pimu2
        endif
        r = pi*pimu2*fact3*fact3
        c = 1
        d = -x2*x2
        sum = ff + r*q
        sum1 = p
        do i = 1, maxit
            ff = (i*ff + p + q)/(i*i - xmu2)
            c = c*d/i
            p = p/(i - xmu)
            q = q/(i + xmu)
            del = c*(ff + r*q)
            sum = sum + del
            del1 = c*p - i*del
            sum1 = sum1 + del1
            if (abs(del) < (1 + abs(sum))*eps) exit
        enddo
        if (i >= maxit) stop 'bessy series failed to converge'
        rymu = -sum
        ry1 = -sum1*xi2
        rymup = xmu*xi*rymu - ry1
        rjmu = w/(rymup - f*rymu)
    else
        a = 0.25_dp - xmu2
        p = -0.5_dp*xi
        q = 1
        br = 2*x
        bi = 2
        fact = a*xi/(p*p + q*q)
        cr = br + q*fact
        ci = bi + p*fact
        den = br*br + bi*bi
        dr = br/den
        di = -bi/den
        dlr = cr*dr - ci*di
        dli = cr*di + ci*dr
        temp = p*dlr - q*dli
        q = p*dli + q*dlr
        p = temp
        do i = 2, maxit
            a = a + 2*(i - 1)
            bi = bi + 2
            dr = a*dr + br
            di = a*di + bi
            if (abs(dr) + abs(di) < fpmin) dr = fpmin
            fact = a/(cr*cr + ci*ci)
            cr = br + cr*fact
            ci = bi - ci*fact
            if (abs(cr) + abs(ci) < fpmin) cr = fpmin
            den = dr*dr + di*di
            dr = dr/den
            di = -di/den
            dlr = cr*dr - ci*di
            dli = cr*di + ci*dr
            temp = p*dlr - q*dli
            q = p*dli + q*dlr
            p = temp
            if (abs(dlr-1) + abs(dli) < eps) exit
        enddo
        if (i >= maxit) stop 'cf2 failed in bessjy'
        gam = (p - f)/q
        rjmu = sqrt(w/((p - f)*gam + q))
        rjmu = sign(rjmu,rjl)
        rymu = rjmu*gam
        rymup = rymu*(p + q/gam)
        ry1 = xmu*xi*rymu - rymup
    endif
    fact = rjmu/rjl
    rj = rjl1*fact
    rjp = rjp1*fact
    do i = 1,nl
        rytemp = (xmu + i)*xi2*ry1 - rymu
        rymu = ry1
        ry1 = rytemp
    enddo
    ry = rymu
    ryp = nu*xi*rymu - ry1
end subroutine bessjy

!!
!> @brief      Chebyshev expansion for gamma terms in bessjy
!!
!! Modified from Numerical Recipes
!!
!! Evaluates the \f$ \Gamma_1(x) \f$ and \f$ \Gamma_2(x) \f$ terms for the bessjy calculation by
!! Chebyshev expansion for \f$ |x| \leq 1/2 \f$. Also returns \f$ 1/\Gamma(1 + x ) \f$ and
!! \f$ 1/\Gamma(1 − x ) \f$.
!!
!! \f$ \Gamma_1(x) = \frac{1}{2x} \left[\frac{1}{\Gamma(1-x)} - \frac{1}{\Gamma(1+x)} \right] \f$
!!
!! \f$ \Gamma_2(x) = \frac{1}{2} \left[\frac{1}{\Gamma(1-x)} - \frac{1}{\Gamma(1+x)} \right] \f$
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine beschb(x, gam1, gam2, gampl, gammi)
    implicit none
    real(dp), intent(in) :: x !< Point at which the \f$ \Gamma \f$ terms will be evaluated
    real(dp), intent(out) :: gam1 !< \f$ \Gamma_1(x) \f$
    real(dp), intent(out) :: gam2 !< \f$ \Gamma_2(x) \f$
    real(dp), intent(out) :: gampl !< \f$ 1/\Gamma_(1+x) \f$
    real(dp), intent(out) :: gammi !< \f$ 1/\Gamma_(1-x) \f$
    real(dp) :: xx
    real(dp), parameter, dimension(1:7) :: c1 = [-1.142022680371172_dp, 6.516511267076e-3_dp, &
        3.08709017308e-4_dp, -3.470626964e-6_dp, 6.943764e-9_dp, 3.6780e-11_dp, -1.36e-13_dp]
    real(dp), parameter, dimension(1:8) :: c2 = [1.843740587300906_dp, -.076852840844786_dp, &
        1.271927136655e-3_dp, -4.971736704e-6_dp, -3.3126120e-8_dp, 2.42310e-10_dp, -1.70e-13_dp, &
        -1.e-15_dp]

    xx = 8.d0*x*x - 1
    gam1 = chebev(-1._dp, 1._dp, c1, xx)
    gam2 = chebev(-1._dp, 1._dp, c2, xx)
    gampl = gam2 - x*gam1
    gammi = gam2 + x*gam1
end subroutine beschb

!!
!> @brief      Chebyshev evaluation
!!
!! Modified from Numerical Recipes
!!
!! c(:) is an array of Chebyshev coefficients, usually the output from chebft (which must have been
!! called with the same a and b ). The Chebyshev polynomial
!!
!! \f$  \sum_{k=1}^m = c_k T_{k-1}(y) - c_1/2 \f$
!!
!! is evaluated at a point \f$ y = [ x − ( b + a )/2]/[( b − a )/2] \f$.
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function chebev(a, b, c, x) result(r)
    implicit none
    real(dp), intent(in) :: a !< interval lower limit
    real(dp), intent(in) :: b !< interval upper limit
    real(dp), intent(in) :: c(:) !< Chebyshev coefficients
    real(dp), intent(in) :: x !< point where the polynomial es evaluated
    integer :: m, j
    real(dp) :: d, dd, sv, y, y2

    m = size(c)
    if ((x - a)*(x - b) > 0) stop 'x not in range in chebev'
    d = 0
    dd = 0
    y = (2*x - a - b)/(b - a)
    y2 = 2*y
    do j = m, 2, -1
        sv = d
        d = y2*d - dd + c(j)
        dd = sv
    enddo
    r = y*d - dd + c(1)/2
end function chebev

end module num_recipes
