!!
!> @brief      modern versions of numerical recipes routines
!!
!! Collection of updated versions of subroutines and functions from the numerical recipes book
!!
!! @author     Rodrigo Navarro Perez
!!
module num_recipes

use precisions, only : dp
use constants, only : pi

implicit none

private

public :: dfridr, context, func, sphbes, bessjy

!!
!> @brief      generic type for data in function callbacks
!!
!! @author     Rodrigo Navarro Perez
!!
type context
    real(dp) :: a, b, c, d
    integer :: i, j, k, l
    real(dp), allocatable :: x(:)
    integer, allocatable :: i_x(:)
    character(len=1024) :: string
    logical :: log
end type

!!
!> @brief      interface for function callbacks
!!
!! 
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
!! Calculates the spherical Bessel functions \f$ j_n(x) \f$, \f$ y_n(x) \f$, and their derivatives 
!! \f$ j_n'(x) \f$, \f$ y_n'(x) \f$ for an integer \$f n \f$.
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


subroutine bessjy(x, xnu, rj, ry, rjp, ryp)
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(in) :: xnu
    real(dp), intent(out) :: rj
    real(dp), intent(out) :: ry
    real(dp), intent(out) :: rjp
    real(dp), intent(out) :: ryp

    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 1.e-16_dp, fpmin = 10*tiny(1._dp), xmin = 2._dp

    integer ::  i, isign, l, nl
    real(dp) :: a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2, &
        fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl, rjl1, rjmu, rjp1, &
        rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1, temp, w, x2, xi, xi2, xmu, xmu2
    
    if (x <=0 .or. xnu <= 0) stop 'bad arguments in bessjy'
    if (x <= xmin ) then
        nl = int(xnu + 0.5_dp)
    else
        nl = max(0, int(xnu - x + 1.5_dp))
    endif
    xmu = xnu - nl
    xmu2 = xmu*xmu
    xi = 1/x
    xi2 = 2*xi
    w = xi2/pi
    isign = 1
    h = xnu*xi
    if (h < fpmin) h = fpmin
    b = xi2*xnu
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
    fact = xnu*xi
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
    ryp = xnu*xi*rymu - ry1
end subroutine bessjy

subroutine beschb(x, gam1, gam2, gampl, gammi)
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(out) :: gam1
    real(dp), intent(out) :: gam2
    real(dp), intent(out) :: gampl
    real(dp), intent(out) :: gammi
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

real(dp) function chebev(a, b, c, x) result(r)
    implicit none
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c(:)
    real(dp), intent(in) :: x
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