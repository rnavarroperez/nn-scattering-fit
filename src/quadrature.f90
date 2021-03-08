module quadrature
use precisions, only: dp
implicit none

private
public :: booles_quadrature
contains
    
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in), dimension(1:) :: fx
    real(dp), intent(in) :: delta_x

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

real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in), dimension(1:) :: fx
    real(dp), intent(in) :: delta_x

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