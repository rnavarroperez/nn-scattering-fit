module special_functions

use precisions, only : dp

implicit none

private
public :: bessel_i0, bessel_i1, bessel_k0, bessel_k1

contains

real(dp) function bessel_k0(x) result(k0)
    implicit none
    real(dp), intent(in) :: x
    real(dp), parameter, dimension(1:7) :: p = [-0.57721566_dp, 0.42278420_dp, 0.23069756_dp, &
        0.3488590e-1_dp, 0.262698e-2_dp, 0.10750e-3_dp, 0.74e-5_dp]
    real(dp), parameter, dimension(1:7) :: q = [1.25331414_dp, -0.7832358e-1_dp, 0.2189568e-1_dp, &
        -0.1062446e-1_dp, 0.587872e-2_dp, -0.251540e-2_dp, 0.53208e-3_dp]
    real(dp), parameter :: cut_off = 2._dp
    real(dp) :: y

    if (x <= 2._dp) then
        y = (x/cut_off)**2
        k0 = (-log(x/cut_off)*bessel_i0(x)) + (p(1) + y*(p(2) + y*(p(3) + y*(p(4) + y*(p(5) + y*(p(6) + y*p(7)))))))
    else
        y = cut_off/x
        k0 = (exp(-x)/sqrt(x))*(q(1) + y*(q(2) + y*(q(3) + y*(q(4) + y*(q(5) + y*(q(6) + y*q(7)))))))
    endif
    
end function bessel_k0

real(dp) function bessel_k1(x) result(k1)
    implicit none
    real(dp) :: x
    real(dp), parameter, dimension(1:7) :: p = [1.0_dp, 0.15443144_dp, -0.67278579_dp, &
         -0.18156897_dp, -0.1919402e-1_dp, -0.110404e-2_dp, -0.4686e-4_dp]
    real(dp), parameter, dimension(1:7) :: q = [1.25331414_dp, 0.23498619_dp, &
         -0.3655620e-1_dp, 0.1504268e-1_dp, -0.780353e-2_dp, 0.325614e-2_dp, &
         -0.68245e-3_dp]
    real(dp), parameter :: cut_off = 2._dp
    real(dp) :: y

    if (x <= 2._dp) then
        y = (x/cut_off)**2
        k1 = (log(x/cut_off)*bessel_i1(x)) + 1/x*(p(1) + y*(p(2) + y*(p(3) + y*(p(4) + y*(p(5) + y*(p(6) + y*p(7)))))))
    else
        y = cut_off/x
        k1 = (exp(-x)/sqrt(x))*(q(1) + y*(q(2) + y*(q(3) + y*(q(4) + y*(q(5) + y*(q(6) + y*q(7)))))))
    endif    
end function bessel_k1

real(dp) function bessel_i0(x) result(i0)
    implicit none
    real(dp), intent(in) :: x
    real(dp), parameter, dimension(1:7) :: p = [1._dp, 3.5156229_dp, 3.0899424_dp, 1.2067492_dp, 0.2659732_dp, &
        0.360768e-1_dp, 0.45813e-2_dp]
    real(dp), parameter, dimension(1:9) :: q = [0.39894228_dp, 0.1328592e-1_dp, 0.225319e-2_dp, -0.157565e-2_dp, &
        0.916281e-2_dp, -0.2057706e-1_dp, 0.2635537e-1_dp, -0.1647633e-1_dp, 0.392377e-2_dp]
    real(dp), parameter :: cut_off = 15._dp/4._dp
    real(dp) :: y

    if (abs(x) < cut_off) then
        y = (x/cut_off)**2
        i0 = p(1) + y*(p(2) + y*(p(3) + y*(p(4) + y*(p(5) + y*(p(6) + y*p(7))))))
    else
        y = cut_off/abs(x)
        i0 = (exp(abs(x))/sqrt(abs(x)))*(q(1) + y*(q(2) + y*(q(3) + y*(q(4) + y*(q(5) + y*(q(6) + y*(q(7) + y*(q(8) &
            + y*q(9)))))))))
    endif

end function bessel_i0

real(dp) function bessel_i1(x) result(i1)
    implicit none
    real(dp), intent(in) :: x
    real(dp), parameter, dimension(1:7) :: p = [0.5_dp, 0.87890594_dp, 0.51498869_dp, 0.15084934_dp, 0.2658733e-1_dp, &
        0.301532e-2_dp, 0.32411e-3_dp]
    real(dp), parameter, dimension(1:9) :: q = [0.39894228_dp, -0.3988024e-1_dp, -0.362018e-2_dp, 0.163801e-2_dp, &
        -0.1031555e-1_dp, 0.2282967e-1_dp, -0.2895312e-1_dp, 0.1787654e-1_dp, -0.420059e-2_dp]
    real(dp), parameter :: cut_off = 15._dp/4._dp
    real(dp) :: y

    if (abs(x) < cut_off) then
        y = (x/cut_off)**2
        i1 = x*(p(1) + y*(p(2) + y*(p(3) + y*(p(4) + y*(p(5) + y*(p(6) + y*p(7)))))))
    else
        y = cut_off/abs(x)
        i1 = (exp(abs(x))/sqrt(abs(x)))*(q(1) + y*(q(2) + y*(q(3) + y*(q(4) + y*(q(5) + y*(q(6) + y*(q(7) + y*(q(8) &
            + y*q(9)))))))))
        if(x < 0) i1 = -i1
    endif

end function bessel_i1
    
end module special_functions