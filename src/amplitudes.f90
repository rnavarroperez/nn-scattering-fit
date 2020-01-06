module amplitudes
use precisions, only : dp
use nn_phaseshifts, only : eta_prime
use constants, only : i_, m_p => proton_mass, hbar_c, alpha, pi, m_e => electron_mass, &
    mu_p => mu_proton
use num_recipes, only : cmplx_log_gamma, spherical_harmonic, kronecker_delta
implicit none

private

public :: saclay_amplitudes

contains

subroutine saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, &
    d_c, d_d, d_e)
    implicit none
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: theta
    character(len=2), intent(in) :: reaction
    real(dp), intent(in) :: phases(:, :)
    real(dp), intent(in) :: d_phases(:, :, :)
    complex(dp), intent(out) :: a
    complex(dp), intent(out) :: b
    complex(dp), intent(out) :: c
    complex(dp), intent(out) :: d
    complex(dp), intent(out) :: e
    complex(dp), intent(out), allocatable :: d_a(:)
    complex(dp), intent(out), allocatable :: d_b(:)
    complex(dp), intent(out), allocatable :: d_c(:)
    complex(dp), intent(out), allocatable :: d_d(:)
    complex(dp), intent(out), allocatable :: d_e(:)

    complex(dp) :: m_000, m_100, m_110, m_101, m_111, m_11m1
    complex(dp), allocatable, dimension(:) :: d_m_000, d_m_100, d_m_110, d_m_101, d_m_111, d_m_11m1
    integer :: n_params

    n_params = size(d_phases, 1)
    allocate(d_a(1: n_params))
    d_a = (0, 0)
    allocate(d_b, d_c, d_d, d_e, source = d_a)
    call partial_wave_amplitude_sum(0, 0, 0, k_cm, theta, reaction, phases, d_phases, m_000,  d_m_000)
    call partial_wave_amplitude_sum(1, 0, 0, k_cm, theta, reaction, phases, d_phases, m_100,  d_m_100)
    call partial_wave_amplitude_sum(1, 1, 0, k_cm, theta, reaction, phases, d_phases, m_110,  d_m_110)
    call partial_wave_amplitude_sum(1, 0, 1, k_cm, theta, reaction, phases, d_phases, m_101,  d_m_101)
    call partial_wave_amplitude_sum(1, 1, 1, k_cm, theta, reaction, phases, d_phases, m_111,  d_m_111)
    call partial_wave_amplitude_sum(1, 1,-1, k_cm, theta, reaction, phases, d_phases, m_11m1, d_m_11m1)


    a = 0.5_dp*(m_111 + m_100 - m_11m1)
    b = 0.5_dp*(m_111 + m_000 + m_11m1)
    c = 0.5_dp*(m_111 - m_000 + m_11m1)
    if (theta == 0 .or. theta == pi) then
        d = (-m_111 + m_100 + m_11m1)/(2*cos(theta)) 
    else
        d = -(m_110 + m_101)/(sqrt(2._dp)*sin(Theta))
    endif
    e = i_*(m_110 - m_101)/sqrt(2._dp)

    d_a = 0.5_dp*(d_m_111 + d_m_100 - d_m_11m1)
    d_b = 0.5_dp*(d_m_111 + d_m_000 + d_m_11m1)
    d_c = 0.5_dp*(d_m_111 - d_m_000 + d_m_11m1)
    if (theta == 0 .or. theta == pi) then
        d_d = (-d_m_111 + d_m_100 + d_m_11m1)/(2*cos(theta))
    else
        d_d = -(d_m_110 + d_m_101)/(sqrt(2._dp)*sin(Theta)) 
    endif
    d_e = i_*(d_m_110 - d_m_101)/sqrt(2._dp)
end subroutine saclay_amplitudes


subroutine partial_wave_amplitude_sum(s, ms, mj, k_cm, theta, reaction, phases, d_phases, m, d_m)
    implicit none
    integer, intent(in) :: s
    integer, intent(in) :: ms
    integer, intent(in) :: mj
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: theta
    character(len=2), intent(in) :: reaction
    real(dp), intent(in) :: phases(:, :)
    real(dp), intent(in) :: d_phases(:, :, :)
    complex(dp), intent(out) :: m
    complex(dp), intent(out), allocatable  :: d_m(:)

    integer :: n_params, j_max, j, l, lp
    real(dp) :: etap, sigma_0, sigma_l, sigma_lp, Ylp, cg_1, cg_2
    complex(dp) :: sm
    complex(dp), allocatable :: d_sm(:)

    m = (0, 0)
    n_params = size(d_phases, 1)
    allocate(d_m(1: n_params))
    d_m = (0, 0)

    if (reaction == 'pp') then
        etap = eta_prime(k_cm)
        sigma_0 = coulomb_sigma_l(0._dp, etap)
    endif

    j_max = size(phases, 2) - 1

    do j = 0, j_max
        do l = j-1, j+1
            if (reaction == 'pp') sigma_l = coulomb_sigma_l(real(l, kind=dp), etap)
            do lp = j-1, j+1
                if (reaction == 'pp') sigma_lp = coulomb_sigma_l(real(lp, kind=dp), etap)
                if (lp >= abs(mj - ms) .and. lp >= 0 .and. l >= 0) then
                    Ylp = real(spherical_harmonic(lp, mj-ms, theta, 0._dp))
                    call s_matrix(lp, l, j, s, phases, d_phases, k_cm, etap, reaction, sm, d_sm)
                    cg_1 = clebsch_gordan(lp, s , j, ms, mj)
                    cg_2 = clebsch_gordan(l, s , j, mj, mj)
                    m = m + 2*Ylp*cg_1*(i_**(l - lp))*exp(i_*(sigma_lp - sigma_0))*sm &
                        *exp(i_*(sigma_l-sigma_0))*cg_2*sqrt(2*l + 1._dp)*sqrt(4*pi)/(2*i_*k_cm)
                    d_m = d_m + 2*Ylp*cg_1*(i_**(l - lp))*exp(i_*(sigma_lp - sigma_0))*d_sm &
                        *exp(i_*(sigma_l-sigma_0))*cg_2*sqrt(2*l + 1._dp)*sqrt(4*pi)/(2*i_*k_cm)
                endif
            enddo
        enddo
    enddo
end subroutine partial_wave_amplitude_sum

real(dp) function clebsch_gordan(l_, s_ , j_, ms_, mj_) result(cgc)
    implicit none
    integer :: l_, s_ , j_, ms_, mj_
    real(dp) :: l, s, j, ms, mj
    l = real(l_, kind = dp)
    s = real(s_, kind = dp)
    j = real(j_, kind = dp)
    ms = real(ms_, kind = dp)
    mj = real(mj_, kind = dp)

    if (abs(ms_) > 1 .or. abs(mj_) > j_ .or. abs(mj_ - ms_) > l_) then
        cgc = 0._dp
        return
    endif

    select case(s_)
    case (0)
        if (l_ == j_) then
            cgc = 1._dp
        else
            cgc = 0._dp
        endif
    case (1)
        select case(ms_)
        case (1)
            if (j_ == l_ + 1) then
                cgc = sqrt((l + mj)*(l + mj + 1)/(2*l + 1)/(2*l + 2))
            elseif (j_ == l_) then
                cgc = -sqrt((l + mj)*(l - mj + 1)/(2*l)/(l + 1))
            elseif (j_ == l_ - 1) then
                cgc = sqrt((l - mj)*(l - mj + 1)/(2*l)/(2*l + 1))
            else
                cgc = 0._dp
            endif
        case (0)
            if (j_ == l_ + 1) then
                cgc = sqrt((l - mj + 1)*(l + mj + 1)/(2*l + 1)/(l + 1))
            elseif (j_ == l_ .and. l_ > 0) then
                cgc = sqrt(mj**2/l/(l+1))
            elseif (j_ == l_ - 1) then
                cgc = -sqrt((l - mj)*(l + mj)/l/(2*l + 1))
            else
                cgc = 0._dp
            endif
        case (-1)
            if (j_ == l_ + 1) then
                cgc = sqrt((l - mj)*(l - mj + 1)/(2*l + 1)/(2*l + 2))
            elseif (j_ == l_) then
                cgc = sqrt((l - mj)*(l + mj + 1)/(2*l)/(l + 1))
            elseif (j_ == l_ - 1) then
                cgc = sqrt((l + mj + 1)*(l + mj)/(2*l)/(2*l + 1))
            else
                cgc = 0._dp
            endif
        case default
            stop 'ms_ has to be -1, 0, or 1 in clebsch_gordan'
        end select
    case default
        stop 's_ has to be 0 or 1 in clebsch_gordan'
    end select
    
end function clebsch_gordan

subroutine s_matrix(lp, l, j, s, phases, d_phases, k_cm, eta, reaction, sm, d_sm)
    implicit none
    integer, intent(in) :: lp
    integer, intent(in) :: l
    integer, intent(in) :: j
    integer, intent(in) :: s
    real(dp), intent(in) :: phases(:,:)
    real(dp), intent(in) :: d_phases(:, :, :)
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: eta
    character(len=2), intent(in) :: reaction
    complex(dp), intent(out) :: sm
    complex(dp), intent(out), allocatable :: d_sm(:)

    integer :: n_params
    real(dp) :: alphap, lambda, sigma_l, sigma_lambda, rho_l, t_lab, nu, tau00, tau0, lambdap, &
        sigma_lp, sigma_lambdap, rho_lp
    real(dp), allocatable :: mm_phases(:, :)

    n_params = size(d_phases, 1)
    allocate(d_sm(1: n_params))
    d_sm = (0, 0)

    allocate(mm_phases, mold = phases)
    mm_phases = 0

    select case(reaction)
    case('pp')
        alphap = 2*k_cm*eta/m_p*hbar_c
        lambda = (-1 + sqrt(1 + 4*l*(l+1) - 4*alpha*alphap))/2._dp
        sigma_l = coulomb_sigma_l(real(l, kind=dp), eta)
        sigma_lambda = coulomb_sigma_l(lambda, eta)
        rho_l = sigma_lambda - sigma_l + (l - lambda)*pi/2._dp

        lambdap = (-1 + sqrt(1 + 4*lp*(lp+1) - 4*alpha*alphap))/2._dp
        sigma_lp = coulomb_sigma_l(real(lp, kind = dp), eta)
        sigma_lambdap = coulomb_sigma_l(lambdap, eta)
        rho_lp = sigma_lambdap - sigma_lp + (lp - lambdap)*pi/2._dp

        t_lab = 2/m_p*(k_cm*hbar_c)**2
        nu = 4*m_e**2/(m_p*t_lab)
        tau00 = -alpha*eta/(6*pi)*(0.5_dp*(log(2._dp/nu))**2 + 1.7615_dp - 0.2804_dp*log(2._dp/nu))
        tau0 = alpha*eta/(3*Pi)*(5.2_dp*eta - 3.83_dp*eta*sqrt(nu/2._dp) - 0.29_dp*eta**2 &
             + (1.202_dp*eta**2 - 2.36_dp*eta*sqrt(nu/2._dp))*log(2._dp/nu)) + tau00
        if (s == 1) then
            call mm_phaseshifts(k_cm, eta, mm_phases)
        endif
    case('np')
        rho_l = 0._dp
        rho_lp = 0._dp
        tau0 = 0._dp
    case default
        stop 's_matrix only takes pp and np as reaction channel'
    end select

    select case(s)
    case (0)
        if (lp == l .and. l == j) then
            sm = (exp(2*i_*phases(1, j+1)) - 1)*exp(2*i_*(rho_l + tau0))
            d_sm = 2*i_*exp(2*i_*phases(1, j+1))*d_phases(:, 1, j+1)*exp(2*i_*(rho_l + tau0))
        else
            sm = 1 - kronecker_delta(lp, l)
            d_sm = (0, 0)
        endif        
    case (1)
        if (l == lp) then
            if (l == j) then
                sm = (exp(2*i_*phases(2, j+1)) - 1)*exp(2*i_*rho_l)*exp(2*i_*mm_phases(2, j+1))
                d_sm = 2*i_*exp(2*i_*phases(2, j+1))*d_phases(:, 2, j+1)*exp(2*i_*rho_l)&
                    *exp(2*i_*mm_phases(2, j+1))
            else if (l == j-1) then
                sm = (cos(2*phases(4, j+1))*exp(2*i_*phases(3, j+1)) - 1)*exp(2*i_*rho_l)&
                    *exp(2*i_*mm_phases(3, j+1))
                d_sm = (-2*sin(2*phases(4, j+1))*d_phases(:, 4, j+1)*exp(2*i_*phases(3, j+1)) &
                    + cos(2*phases(4, j+1))*2*i_*exp(2*i_*phases(3, j+1))*d_phases(:, 3, j+1))&
                    *exp(2*i_*rho_l)*exp(2*i_*mm_phases(3, j+1))
            else if (l == j+1) then
                sm = (cos(2*phases(4, j+1))*exp(2*i_*phases(5, j+1)) - 1)*exp(2*i_*rho_l)&
                    *exp(2*i_*mm_phases(5, j+1))
                d_sm = (-2*sin(2*phases(4, j+1))*d_phases(:, 4, j+1)*exp(2*i_*phases(5, j+1)) &
                    + cos(2*phases(4, j+1))*2*i_*exp(2*i_*phases(5, j+1))*d_phases(:, 5, j+1))&
                    *exp(2*i_*rho_l)*exp(2*i_*mm_phases(5, j+1))
            else
                sm = -1
                d_sm = (0, 0)
            endif
        else if (abs(l-lp) == 2) then
            sm = i_*sin(2*phases(4, j+1))*exp(i_*(phases(3, j+1) + phases(5, j+1)))&
                *exp(i_*(rho_l + rho_lp))*exp(i_*(mm_phases(3, j+1) + mm_phases(5, j+1)))
            d_sm = i_*(2*cos(2*phases(4, j+1))*d_phases(:, 4, j+1)&
                *exp(i_*(phases(3, j+1) + phases(5, j+1))) + sin(2*phases(4, j+1))*i_&
                *exp(i_*(phases(3, j+1) + phases(5, j+1)))*(d_phases(:, 3, j+1) &
                + d_phases(:, 5, j+1)))*exp(i_*(rho_l+rho_lp)) &
                *exp(i_*(mm_phases(3, j+1) + mm_phases(5, j+1)))
        else
            sm = -kronecker_delta(lp, l)
            d_sm = (0, 0)
        endif
    case default
        stop 's has to be zero or one in s_matrix'
    end select
end subroutine s_matrix

subroutine mm_phaseshifts(k_cm, eta, mm_phases)
    implicit none
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: eta
    real(dp), intent(out) :: mm_phases(:, :)

    real(dp) :: f_T, f_ls, I_ll, I_lp2lp2, I_llp2
    integer :: l, j_max, i

    j_max = size(mm_phases, 2)

    f_T = -alpha*mu_p**2/(4*m_p**2)
    f_ls = -alpha*(8*mu_p - 2)/(4*m_p**2)
    mm_phases = 0._dp
    l = -1
    I_lp2lp2 = mm_coulomb_Ill(l+2, eta)

    mm_phases(5, 1) = ((2*l + 6)/(2*l + 3._dp)*f_T + (l + 3)*f_ls)*I_lp2lp2

    do i = 2, j_max - 1, 2
        l = i - 1
        I_ll = I_lp2lp2
        mm_phases(2, i) = -(2*f_T - f_ls)*I_ll
        mm_phases(3, i+1) = -(-2*l/(2*l+3._dp)*f_T + l*f_ls)*I_ll
        I_llp2 = mm_coulomb_Illp2(l, eta)
        mm_phases(4, i+1) = -(6*sqrt((l + 1._dp)*(l + 2))/(2*l + 3._dp)*f_T)*I_llp2
        I_lp2lp2 = mm_coulomb_Ill(l+2, eta)
        mm_phases(5, i+1) = ((2*l + 6)/(2*l + 3._dp)*f_T + (l + 3)*f_ls)*I_lp2lp2
    enddo
    mm_phases = m_p*k_cm*mm_phases*hbar_c    
end subroutine mm_phaseshifts

real(dp) function mm_coulomb_Ill(l, eta) result(I_ll)
    implicit none
    integer, intent(in) :: l
    real(dp), intent(in) :: eta
    integer :: i
    real(dp) ::  sum
    sum = 0._dp
    do i = 0, l
       sum = sum + 1._dp/(i**2 + eta**2)
    endDo
    I_ll = 1._dp/(2._dp*l*(l+1)) + (1 - pi*eta + pi*eta/tanh(pi*eta) - 2*eta**2*sum)/&
        (2._dp*l*(l+1)*(2*l+1))
end function mm_coulomb_Ill

real(dp) function mm_coulomb_Illp2(l, eta) result(I_llp2)
    implicit none
    integer, intent(in) :: l
    real(dp), intent(in) :: eta
    I_llp2 = 1._dp/(6*abs(l + 1 +i_*eta)*abs(l + 2 + i_*eta))
end function mm_coulomb_Illp2

real(dp) function coulomb_sigma_l(l, eta) result(sig_l)
    implicit none
    real(dp), intent(in) :: l
    real(dp), intent(in) :: eta
    complex(dp) :: z
    z = 1 + l + i_*eta
    sig_l = aimag(cmplx_log_gamma(z))
end function coulomb_sigma_l

    
end module amplitudes