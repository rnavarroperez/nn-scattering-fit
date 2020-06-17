module deuteron
use precisions, only : dp
use nn_phaseshifts, only : nn_local_model
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi, alpha
implicit none

private

public :: binding_energy

type :: ds_wave_function
    integer :: n_radii
    real(dp) :: gamma
    real(dp), allocatable, dimension(:) :: a_s
    real(dp), allocatable, dimension(:) :: b_s
    real(dp), allocatable, dimension(:) :: a_d
    real(dp), allocatable, dimension(:) :: b_d
    real(dp), allocatable, dimension(:) :: radii
    real(dp) :: dr
end type ds_wave_function

contains

subroutine binding_energy(model, parameters, be, dbe)
    implicit none
    type(nn_local_model), intent(in) :: model
    real(dp), intent(in), dimension(:) ::  parameters
    real(dp), intent(out) :: be
    real(dp), intent(out), allocatable, dimension(:) :: dbe

    real(dp) :: k, r, v, w
    type(ds_wave_function) :: wavefunc
    integer :: i
    real(dp), parameter :: mu = 2*m_p*m_n/(m_p + m_n)
    real(dp), dimension(1:5, 1:2) :: v_pw
    real(dp), allocatable, dimension(:, :, :) :: dv_pw

    k = wave_number(model, parameters)
    be = m_p + m_n - sqrt(m_p**2 - (k*hbar_c)**2) - sqrt(m_n**2 - (k*hbar_c)**2)
    wavefunc = wave_function(k, model, parameters)
    call normalize(wavefunc)
    allocate(dbe, mold=parameters)
    dbe = 0._dp
    do i = 1, wavefunc%n_radii
        r = wavefunc%radii(i)
        call model%potential(parameters, r, 'np', v_pw, dv_pw)
        dv_pw = dv_pw*mu*model%dr/(hbar_c**2)
        v = s_wave_function(wavefunc, i)
        w = d_wave_function(wavefunc, i)
        dbe = dbe + dv_pw(:, 3, 2)*v**2 + 2*dv_pw(:, 4, 2)*v*w + dv_pw(:, 5, 2)*w**2
    enddo
    dbe = -0.5_dp*(1/sqrt(m_p**2 - (k*hbar_c)**2) + 1/sqrt(m_p**2 - (k*hbar_c)**2))*hbar_c**2*dbe
end subroutine binding_energy

subroutine normalize(wavefunc)
    implicit none
    type(ds_wave_function), intent(inout) :: wavefunc

    real(dp) :: r, k, s
    real(dp) :: a_s, a_d, a_norm
    integer :: i, n_points
    n_points = wavefunc%n_radii
    ! one half from the trapezoid rule plus one quarter from the triangle from r = 0 to r = r_1 = 0.5*dr
    s = 0.75_dp*probability_density(wavefunc, 1)
    ! trapezoid rule for all intermediate points
    do i = 2, n_points - 1
        s = s + probability_density(wavefunc, i)
    enddo
    ! one half of the last point from trapezoid rule
    s = s + 0.5_dp*probability_density(wavefunc, n_points)
    ! overall dr factor from the trapezoid rule
    s = s*wavefunc%dr
    a_s = wavefunc%a_s(n_points)
    a_d = wavefunc%a_d(n_points)
    r = wavefunc%radii(n_points)
    k = wavefunc%gamma
    ! Integral from last grid point to infinity where the wave function has exponential decay
    s = s + exp(-2*k*r)*(a_s**2/(2*k) + a_d**2*(6 + k*r*(12 + k*r*(6 + k*r)))/(2*k**4*r**3))
    a_norm = sqrt(1._dp/s)
    wavefunc%a_s = a_norm*wavefunc%a_s
    wavefunc%b_s = a_norm*wavefunc%b_s
    wavefunc%a_d = a_norm*wavefunc%a_d
    wavefunc%b_d = a_norm*wavefunc%b_d
end subroutine normalize

real(dp) function probability_density(wavefunc, i) result(p)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc
    integer, intent(in) :: i

    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in probability_density is incompatible with wavefunc grid'
    endif
    p = s_wave_function(wavefunc, i)**2 + d_wave_function(wavefunc, i)**2
end function probability_density

real(dp) function s_wave_function(wavefunc, i) result(v)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc
    integer, intent(in) :: i

    real(dp) :: k ,r
    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in probability_density is incompatible with wavefunc grid'
    endif
    k = wavefunc%gamma
    r = wavefunc%radii(i)
    v = wavefunc%a_s(i)*j0_bound(k*r) + wavefunc%b_s(i)*y0_bound(k*r)
end function s_wave_function

real(dp) function d_wave_function(wavefunc, i) result(w)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc
    integer, intent(in) :: i

    real(dp) :: k ,r
    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in probability_density is incompatible with wavefunc grid'
    endif
    k = wavefunc%gamma
    r = wavefunc%radii(i)
    w = wavefunc%a_d(i)*j2_bound(k*r) + wavefunc%b_d(i)*y2_bound(k*r)
end function d_wave_function

real(dp) function wave_number(model, parameters, tolerance) result(r)
    implicit none
    type(nn_local_model), intent(in) :: model
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in), optional :: tolerance

    real(dp), parameter :: default_tol = 1.e-10_dp
    real(dp), parameter :: k_init_l = 0.18_dp
    real(dp), parameter :: k_init_h = 0.30_dp
    real(dp) :: k_low
    real(dp) :: k_high
    real(dp) :: bs_low, bs_high, k_0, bs_0, slope, dk
    real(dp) :: local_tolerance

    if (present(tolerance)) then
        local_tolerance = tolerance
    else
        local_tolerance = default_tol
    end if
    k_low = k_init_l
    k_high = k_init_h
    bs_low = irregularity(k_low, model, parameters)
    bs_high = irregularity(k_high, model, parameters)

    if (sign(1._dp, bs_low) == sign(1._dp, bs_high)) then
        print*, k_low, bs_low,  k_high, bs_high     
        stop 'deuteron wave function irregularity is not bracketed in wave_number function'
    endif
    do
        slope = (bs_low - bs_high)/(k_low - k_high)
        dk = -bs_low/slope
        k_0 = -(bs_low - slope*k_low)/slope
        bs_0 = irregularity(k_0, model, parameters)
        if ((abs(k_low-k_0) < local_tolerance) .or. (abs(k_high-k_0) < local_tolerance)) exit
        if (sign(1._dp, bs_low) == sign(1._dp, bs_0)) then
            k_low = k_0
            bs_low = bs_0
        else
            k_high = k_0
            bs_high = bs_0
        endif
    enddo
    r = k_0
        
end function wave_number

real(dp) function irregularity(k, model, parameters) result(r)
    implicit none
    real(dp), intent(in) :: k
    type(nn_local_model), intent(in) :: model
    real(dp), intent(in), dimension(:) :: parameters

    type(ds_wave_function) :: deuteron_wf

    deuteron_wf = wave_function(k, model, parameters)
    r = deuteron_wf%b_s(0)
end function irregularity

type(ds_wave_function) function wave_function(k, model, parameters) result(r)
    implicit none
    real(dp), intent(in) :: k
    type(nn_local_model), intent(in) :: model
    real(dp), intent(in), dimension(:) :: parameters
    
    real(dp) ::  r_i
    real(dp), parameter :: mu = 2*m_p*m_n/(m_p + m_n)
    real(dp), dimension(1:5, 1:2) :: v_pw
    real(dp), allocatable, dimension(:, :, :) :: dv_pw
    real(dp), allocatable, dimension(:) :: a_alpha, b_alpha, c_alpha, d_alpha
    real(dp), allocatable, dimension(:) :: a_beta, b_beta, c_beta, d_beta
    integer :: i
    real(dp) :: beta

    r%n_radii = int(model%r_max/model%dr)
    r%gamma = k
    allocate(r%radii(1:r%n_radii))
    allocate(r%a_s(0:r%n_radii))
    allocate(r%b_s, r%a_d, r%b_d, a_alpha, b_alpha, c_alpha, d_alpha, a_beta, b_beta, c_beta, d_beta, mold=r%a_s)
    r%dr = model%dr
    a_alpha(r%n_radii) = -1._dp
    b_alpha(r%n_radii) = -1._dp
    c_alpha(r%n_radii) = 0._dp
    d_alpha(r%n_radii) = 0._dp

    a_beta(r%n_radii) = 0._dp
    b_beta(r%n_radii) = 0._dp
    c_beta(r%n_radii) = -1._dp
    d_beta(r%n_radii) = -1._dp

    do i = r%n_radii , 1, -1
        a_alpha(i-1) = a_alpha(i)
        b_alpha(i-1) = b_alpha(i)
        c_alpha(i-1) = c_alpha(i)
        d_alpha(i-1) = d_alpha(i)

        a_beta(i-1) = a_beta(i)
        b_beta(i-1) = b_beta(i)
        c_beta(i-1) = c_beta(i)
        d_beta(i-1) = d_beta(i)

        r_i = model%dr*(i - 0.5_dp)
        call model%potential(parameters, r_i, 'np', v_pw, dv_pw)
        v_pw = v_pw*mu*model%dr/(hbar_c**2)
        call reverse_variable_wave(k, r_i, v_pw(3:5, 2), a_alpha(i-1), b_alpha(i-1), c_alpha(i-1), d_alpha(i-1))
        call reverse_variable_wave(k, r_i, v_pw(3:5, 2),  a_beta(i-1),  b_beta(i-1),  c_beta(i-1),  d_beta(i-1))
        r%radii(i) = r_i
    enddo

    beta = -d_alpha(0)/d_beta(0)
    r%a_s = (a_alpha + beta*a_beta)
    r%b_s = (b_alpha + beta*b_beta)
    r%a_d = (c_alpha + beta*c_beta)
    r%b_d = (d_alpha + beta*d_beta)
    
end function wave_function

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