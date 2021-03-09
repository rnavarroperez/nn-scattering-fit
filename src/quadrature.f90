!!
!> @brief      Tools for quadrature
!!
!! Subroutine and functions to perform integrals
!! 
!! @author     Rodrigo Navarro Perez
!!
module quadrature
use precisions, only: dp
implicit none

private
public :: booles_quadrature
contains
    
!!
!> @brief      Booles quadrature
!!
!! Given an array containing equidistant evaluations of a function and
!! the distance between evaluation points, uses the Booles rule to
!! calculate the integral of the evaluated function.
!!
!! Given the number of points necessary for the Booles rule, the 
!! size of the array has to be a multiple of 4 minus 1
!!
!! @return     Integral of a function
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in), dimension(1:) :: fx !< Array of equidistant evaluations of a function
    real(dp), intent(in) :: delta_x !< distance between evaluation points

    integer :: fx_size, i

    fx_size = size(fx)

    if (mod(fx_size-1,4) /= 0) then
        print *, 'fx array size in booles_quadrature has to be a multiple of 4 plus 1'
        stop
    endif

    s = 0._dp

    do i=1,fx_size-1,4
        s = s + booles_rule(fx(i:i+4), delta_x)
    enddo
end function booles_quadrature

!!
!> @brief      Booles rule
!!
!! Given an array containing 5 equidistant evaluations of a function and
!! the distance between evaluation points, calculates Booles rule using those 5 points.
!!
!! @return     Booles rule for 5 points
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in), dimension(1:) :: fx !< Array of 5 equidistant evaluations of a function
    real(dp), intent(in) :: delta_x !< distance between evaluation points

    real(dp) :: fx0, fx1, fx2, fx3, fx4

    if( size(fx) /= 5 ) then
        print *, 'fx array size in booles_rule has to be 5'
        stop
    endif

    fx0 = fx(1)
    fx1 = fx(2)
    fx2 = fx(3)
    fx3 = fx(4)
    fx4 = fx(5)
    
    s = delta_x*2/45._dp*(7*fx4 + 32*fx3 + 12*fx2 + 32*fx1 + 7*fx0)
end function booles_rule

end module quadrature