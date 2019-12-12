!!
!> @brief      modern versions of numerical recipes routines
!!
!! Collection of updated versions of subroutines and functions from the numerical recipes book
!!
!! @author     Rodrigo Navarro Perez
!!
module num_recipes

use precisions, only : dp

implicit none

private

public :: dfridr, context

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
   
end module num_recipes