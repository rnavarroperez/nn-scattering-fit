module deuteron
use precisions, only : dp
use nn_phaseshifts, only : nn_local_model
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi, alpha
implicit none



contains


subroutine deuteron_binding(k, model, parameters, be, dbe)
    implicit none
    real(dp), intent(in) :: k
    type(nn_local_model), intent(in) :: model
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(out) :: be
    real(dp), intent(out), allocatable, dimension(:) :: dbe

    real(dp) ::  r
    real(dp), parameter :: mu = 2*m_p*m_n/(m_p + m_n)
    real(dp), dimension(1:5, 1:2) :: v_pw
    real(dp), allocatable, dimension(:, :, :) :: dv_pw
    real(dp) :: a_alpha, b_alpha, c_alpha, d_alpha
    real(dp) :: a_beta, b_beta, c_beta, d_beta
    real(dp) :: beta
    integer :: counter

    a_alpha = -1._dp
    b_alpha = -1._dp
    c_alpha = 0._dp
    d_alpha = 0._dp

    a_beta = 0._dp
    b_beta = 0._dp
    c_beta = -1._dp
    d_beta = -1._dp

    counter = 0
    r = model%r_max + model%dr/2._dp
    do
       if(r < 0._dp) exit
       call model%potential(parameters, r, 'np', v_pw, dv_pw)
       v_pw = v_pw*mu*model%dr/(hbar_c**2)
       counter = counter + 1
       call reverse_variable_wave(k, r, v_pw(3:5, 2), a_alpha, b_alpha, c_alpha, d_alpha)
       call reverse_variable_wave(k, r, v_pw(3:5, 2), a_beta, b_beta, c_beta, d_beta)
       r = r - model%dr 
    enddo

    beta = -d_alpha/d_beta
    be = (b_alpha + beta*b_beta)
    ! print*, k, d_alpha, d_beta, beta, be**2, counter, int(model%r_max/model%dr)+1  
    allocate(dbe, mold=parameters)
    dbe = 0
end subroutine deuteron_binding

subroutine reverse_variable_wave(k, r, lambdas, a, b, c, d)
    implicit none
    real(dp), intent(in) :: k
    real(dp), intent(in) :: r
    real(dp), intent(in), dimension(1:3) :: lambdas
    real(dp), intent(inout) :: a
    real(dp), intent(inout) :: b
    real(dp), intent(inout) :: c
    real(dp), intent(inout) :: d

    real(dp) :: j0, y0, j2, y2
    real(dp) :: l_3s1, l_ep1, l_3d1
    real(dp) :: diff_b, diff_d
    j0 = j0_bound(k*r)
    y0 = y0_bound(k*r)
    j2 = j2_bound(k*r)
    y2 = y2_bound(k*r)
    l_3s1 = lambdas(1)
    l_ep1 = lambdas(2)
    l_3d1 = lambdas(3)
    
    diff_b = -j0*(l_3s1*(a*j0 + b*y0) + l_ep1*(c*j2 + d*y2))/k
    diff_d = -j2*(l_3d1*(c*j2 + d*y2) + l_ep1*(a*j0 + b*y0))/k

    a = a - diff_b*y0/j0
    b = b + diff_b
    c = c - diff_d*y2/j2
    d = d + diff_d
end subroutine reverse_variable_wave

real(dp) function j0_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    r = sinh(x)
end function j0_bound

real(dp) function y0_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    r = -cosh(x)
end function y0_bound

real(dp) function j2_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    r = sinh(x) - 3*cosh(x)/x + 3*sinh(x)/(x**2)
end function j2_bound

real(dp) function y2_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    r = -cosh(x) + 3*sinh(x)/x -3*cosh(x)/(x**2)
end function y2_bound

end module deuteron