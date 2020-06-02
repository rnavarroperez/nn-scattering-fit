!!
!> @brief      NN scattering amplitudes
!!
!! Module to calculate the nucleon-nucleon scattering amplitude given a laboratory
!! energy, scattering angle and the corresponding NN phase-shifts
!!
!! @author     Rodrigo Navarro Perez
!! @author     Raul L Bernal-Gonzalez
!!
module amplitudes

use precisions, only : dp
use nn_phaseshifts, only : eta_prime, momentum_cm
use constants, only : i_, m_p => proton_mass, hbar_c, alpha, pi, m_e => electron_mass, &
    mu_p => mu_proton, m_n => neutron_mass, mu_n => mu_neutron
use num_recipes, only : cmplx_log_gamma, spherical_harmonic, kronecker_delta, legendre_poly
implicit none

real(dp), parameter :: f_T = -alpha*mu_p**2/(4*m_p**2) !< Tensor factor in magnetic moment amplitude
real(dp), parameter :: f_ls = -alpha*(8*mu_p - 2)/(4*m_p**2) !< spin-orbit factor in magnetic moment amplitude

private

public :: saclay_amplitudes, em_amplitudes!, f_amplitudes, df_amplitudes

!!
!> @brief      interface for S matrix subroutine
!!
!! @author     Rodrigo Navarro Perez
!!
interface

    subroutine s_matrix_elements(l_prime, l, j, s, k_cm, eta, reaction, phases, d_phases, sm, d_sm)
        use precisions, only : dp
        implicit none
        integer, intent(in) :: l_prime !< orbital angular momentum quantum number \f$ l' \f$
        integer, intent(in) :: l !< orbital angular momentum quantum number \f$ l\f$
        integer, intent(in) :: j !< total angular momentum quantum number \f$ j \f$
        integer, intent(in) :: s !< spin quantum number
        real(dp), intent(in) :: k_cm !< center of mass momentum, in fm\f$^{-1}\f$
        real(dp), intent(in) :: eta !< Sommerfeld parameter
        character(len=2), intent(in) :: reaction !< reaction channel 'pp' or 'np'
        real(dp), intent(in) :: phases(:,:) !< NN scattering phase-shifts corresponding to the given center of mass momentum
        real(dp), optional, intent(in) :: d_phases(:, :, :) !< derivatives of phase-shifts with respect to the fitting parameters
        complex(dp), intent(out) :: sm !< S matrix element for the given quantum numbers
        complex(dp), optional, intent(out), allocatable :: d_sm(:) !< derivatives of the S matrix element with respect to the fitting parameters
    end subroutine s_matrix_elements
end interface

contains

!!
!> @brief      Calculate electromagnetic amplitudes
!!
!! Wrapper function to calculate either pp or np electromagnetic
!! amplitudes for a given laboratory energy and scattering angle
!!
!! @return     electromagnetic amplitude
!!
!! @author     Rodrigo Navarro Perez
!!
function em_amplitudes(t_lab, angle, channel) result(r)
    implicit none
    real(dp), intent(in) :: t_lab !< Laboratory energy in MeV
    real(dp), intent(in) :: angle !< Scattering angle in degrees
    character(len=*), intent(in) :: channel !< reaction channel ('pp' or 'np')
    complex(dp), dimension(1:5) :: r
    real(dp) :: k_cm, theta
    complex(dp) :: a, b, c, d, e
    k_cm = momentum_cm(t_lab, channel)
    theta = angle*pi/180
    select case(trim(channel))
    case ('pp')
        call em_pp_amplitudes(k_cm, theta, a, b, c, d, e)
    case ('np')
        call em_np_amplitudes(k_cm, theta, a, b, c, d, e)
    case default
        stop 'invalid reacttion channel in em_amplitudes'
    end select
    r = [a, b, c, d, e]
end function em_amplitudes

!!
!> @brief      calculate Saclay amplitudes
!!
!! Calculates the NN scattering amplitude in the Saclay parametrization for a given c.m. momentum, scattering angle,
!! reaction channel and the corresponding phase-shifts.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, &
    d_c, d_d, d_e)
    implicit none
    real(dp), intent(in) :: k_cm !< C.M. momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: theta !< scattering angle in radians
    character(len=2), intent(in) :: reaction !< reaction channel (pp or np)
    real(dp), intent(in) :: phases(:, :) !< corresponding phase-shifts
    real(dp), intent(in) :: d_phases(:, :, :) !< derivatives of the phase-shifts
    complex(dp), intent(out) :: a !< Saclay parameter \f$a\f$
    complex(dp), intent(out) :: b !< Saclay parameter \f$b\f$
    complex(dp), intent(out) :: c !< Saclay parameter \f$c\f$
    complex(dp), intent(out) :: d !< Saclay parameter \f$d\f$
    complex(dp), intent(out) :: e !< Saclay parameter \f$e\f$
    complex(dp), intent(out), allocatable :: d_a(:) !< derivatives of the Saclay parameter \f$a\f$
    complex(dp), intent(out), allocatable :: d_b(:) !< derivatives of the Saclay parameter \f$b\f$
    complex(dp), intent(out), allocatable :: d_c(:) !< derivatives of the Saclay parameter \f$c\f$
    complex(dp), intent(out), allocatable :: d_d(:) !< derivatives of the Saclay parameter \f$d\f$
    complex(dp), intent(out), allocatable :: d_e(:) !< derivatives of the Saclay parameter \f$e\f$

    complex(dp) :: m_000, m_100, m_110, m_101, m_111, m_11m1
    complex(dp), allocatable, dimension(:) :: d_m_000, d_m_100, d_m_110, d_m_101, d_m_111, d_m_11m1
    integer :: n_params

    n_params = size(d_phases, 1)
    allocate(d_a(1: n_params))
    d_a = 0
    allocate(d_b, d_c, d_d, d_e, source = d_a)
    call partial_wave_amplitude_sum(s_matrix, 0, 0, 0, k_cm, theta, reaction, phases, d_phases, m_000,  d_m_000)
    call partial_wave_amplitude_sum(s_matrix, 1, 0, 0, k_cm, theta, reaction, phases, d_phases, m_100,  d_m_100)
    call partial_wave_amplitude_sum(s_matrix, 1, 1, 0, k_cm, theta, reaction, phases, d_phases, m_110,  d_m_110)
    call partial_wave_amplitude_sum(s_matrix, 1, 0, 1, k_cm, theta, reaction, phases, d_phases, m_101,  d_m_101)
    call partial_wave_amplitude_sum(s_matrix, 1, 1, 1, k_cm, theta, reaction, phases, d_phases, m_111,  d_m_111)
    call partial_wave_amplitude_sum(s_matrix, 1, 1,-1, k_cm, theta, reaction, phases, d_phases, m_11m1, d_m_11m1)

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

!!
!> @brief      Calculates the pp electromagnetic amplitude
!! 
!! Given the C.M. momentum and scattering angle, calculates the corresponding pp electromagnetic amplitude.
!! See section III A of PhysRevC.88.064002 for more details.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine em_pp_amplitudes(k_cm, theta, a, b, c, d, e)
    implicit none
    real(dp), intent(in) :: k_cm !< C.M. momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: theta !< scattering angle in radians
    complex(dp), intent(out) :: a !< Saclay parameter \f$a\f$
    complex(dp), intent(out) :: b !< Saclay parameter \f$b\f$
    complex(dp), intent(out) :: c !< Saclay parameter \f$c\f$
    complex(dp), intent(out) :: d !< Saclay parameter \f$d\f$
    complex(dp), intent(out) :: e !< Saclay parameter \f$e\f$

    real(dp) :: etap
    complex(dp) :: f_c, f_cpi, f_c2, f_c2pi, f_vp, f_vppi, m_c0, m_c1, a_coulomb, b_coulomb, c_coulomb, &
        m_00, m_11, m_10, m_01, m_1m1, a_mm, b_mm, c_mm, d_mm, e_mm

    etap = eta_prime(k_cm)
    f_c   = f_coulomb(k_cm, etap, theta)
    f_cpi = f_coulomb(k_cm, etap, pi - theta)
    f_c2   = f_coulomb2(k_cm, etap, theta)
    f_c2pi = f_coulomb2(k_cm, etap, pi - theta)
    f_vp   = f_vacuum_polarization(k_cm, etap, theta)
    f_vppi = f_vacuum_polarization(k_cm, etap, pi - theta)
    m_c0 = (f_c + f_cpi) + (f_c2 + f_c2pi) + (f_vp + f_vppi)
    m_c1 = (f_c - f_cpi) + (f_c2 - f_c2pi) + (f_vp - f_vppi)
    a_coulomb = m_c1
    b_coulomb = 0.5_dp*(m_c1 + m_c0)
    c_coulomb = 0.5_dp*(m_c1 - m_c0)
    call mm_amplitudes(k_cm, etap, theta, m_00, m_11, m_10, m_01, m_1m1)
    a_mm = 0.5_dp*(m_11 + m_00 - m_1m1)
    b_mm = 0.5_dp*(m_11 + m_1m1)
    c_mm = 0.5_dp*(m_11 + m_1m1)
    d_mm = -(m_10 + m_01)/(sqrt(2._dp)*sin(theta))
    if(theta == 0._dp .or. theta == pi) then
        d_mm = (-m_11 + m_00 + m_1m1)/(2._dp*cos(theta))
    endIf
    e_mm = i_*(m_10 - m_01)/sqrt(2._dp)
    a = a_coulomb + a_mm
    b = b_coulomb + b_mm
    c = c_coulomb + c_mm
    d = d_mm
    e = e_mm
end subroutine em_pp_amplitudes

!!
!> @brief      Calculates the np electromagnetic amplitude
!! 
!! Given the C.M. momentum and scattering angle, calculates the corresponding np electromagnetic amplitude.
!! See section 34 of PhysRevC.88.064002 for more details.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine em_np_amplitudes(k_cm, theta, a, b, c, d, e)
    implicit none
    real(dp), intent(in) :: k_cm !< C.M. momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: theta !< scattering angle in radians
    complex(dp), intent(out) :: a !< Saclay parameter \f$a\f$
    complex(dp), intent(out) :: b !< Saclay parameter \f$b\f$
    complex(dp), intent(out) :: c !< Saclay parameter \f$c\f$
    complex(dp), intent(out) :: d !< Saclay parameter \f$d\f$
    complex(dp), intent(out) :: e !< Saclay parameter \f$e\f$

    real(dp) :: etap, f1p, f2p, f1n, f2n, sm, tm, mp_fm, mn_fm
    f1p = 1._dp
    f1n = 0._dp
    mp_fm = m_p/hbar_c !in fm-1
    mn_fm = m_n/hbar_c !in fm-1
    f2p = (mu_p - 1)/(2*mp_fm) !in fm
    f2n = mu_n/(2*mn_fm) !in fm

    etap = eta_prime(k_cm)

    tm = -2*k_cm**2*(1 - cos(theta)) !in fm-2
    sm = 2*sqrt((k_cm**2 + mn_fm**2)*(k_cm**2 + mp_fm**2)) + 2*k_cm**2 + mn_fm**2 + mp_fm**2 !in fm-2
    a = alpha/(tm*sqrt(sm))*((f1n*f1p + tm*f2n*f2p)*(sm - mn_fm**2 - mp_fm**2 &
        + tm/(8*sm*k_cm**2)*((sm - (mn_fm + mp_fm)**2)*(3*sm - (mn_fm - mp_fm)**2) &
        + 2*(sm - (mn_fm - mp_fm)**2)*(sqrt(sm) - mn_fm - mp_fm)**2) + tm**2/(16*sm*k_cm**4)&
        *(sm - (mn_fm - mp_fm)**2)*(sqrt(sm) - mn_fm - mp_fm)**2) &
        + (f1n*f2p + f2n*f1p)*tm*(2*sqrt(sm) - mn_fm - mp_fm + tm/(2*k_cm**2)*(sqrt(sm)-mn_fm-mp_fm)))
    b = alpha/(tm*sqrt(sm))*((f1n*f1p - tm*f2n*f2p)*(sm - mn_fm**2 - mp_fm**2 &
        + tm/(8*sm*k_cm**2)*(sm + (mn_fm - mp_fm)**2)*(sm - (mn_fm + mp_fm)**2)) + (f1n*f2p - f2n*f1p)*tm*(mn_fm - mp_fm))
    c = alpha/(2*sqrt(sm))*(f1n + 2*mn_fm*f2n)*(f1p + 2*mp_fm*f2p)
    d = -c
    e = -i_*alpha*sin(theta)/(tm*sqrt(sm))*((f1n*f1p + tm*f2n*f2p)*(sm - mn_fm**2 - mp_fm**2 &
        - (mn_fm + mp_fm)/(2*sqrt(sm))*(sm - (mn_fm - mp_fm)**2) &
        + (tm*(sqrt(sm) - mn_fm - mp_fm))/(2*(sqrt(sm) + mn_fm + mp_fm))) &
        + (f1n*f2p + f2n*f1p)*(2*k_cm**2*sqrt(sm) + tm*(sqrt(sm) - mn_fm - mp_fm)))

end subroutine em_np_amplitudes

!!
!> @brief      Coulomb scattering amplitude
!!
!! Calculates the scattering amplitude corresponding to the energy dependent Coulomb potential
!! \f$ f_{C1, k}(\theta) = - \frac{\eta}{k} \frac{e^{-i\eta\ln[(1-\cos\theta)/2]}}{1-\cos\theta} \f$.
!!
!! See equation 20 in Phys. Rev. C 91, 029901 (2015) for more details
!!
!! @return     Coulomb scattering amplitude
!!
!! @author     Rodrigo Navarro Perez
!!
complex(dp) function f_coulomb(k_cm, eta, theta) result(f_c)
    implicit none
    real(dp), intent(in) :: k_cm !< Center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: eta !< Energy dependent Sommerfeld parameter (dimensionless)
    real(dp), intent(in) :: theta !< Scattering angle in radians
    f_c = -eta/k_cm*exp(-i_*eta*log((1 - cos(theta))/2._dp))/(1 - cos(theta))
end function f_coulomb


!!
!> @brief      Two photon exchange scattering amplitude
!!
!! Calculates the scattering amplitude corresponding to the energy dependent two photon exchange potential
!! \f$ f_{C2, k}(\theta) = \frac{1}{2ik} \sum_l (2_l+1) e^{2i(\sigma_l - \sigma_0)} (e^{2i\rho_l}-1) P_l(\theta) \f$,
!! where \f$\sigma_l\f$ is the Coulomb phase-shift, \f$\rho_l\f$ is the two photon exchange phase-shift, and 
!! \f$ P_l(\theta) \f$ are the usual Legendre polynomials.
!!
!! See equations 22 in Phys. Rev. C 91, 029901 (2015) for more details
!! 
!! @return     Two photon exchange scattering amplitude
!!
!! @author     Rodrigo Navarro Perez
!!
complex(dp) function f_coulomb2(k_cm, eta, theta) result(f_c)
    implicit none
    real(dp), intent(in) :: k_cm !< Center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: eta !< Energy dependent Sommerfeld parameter (dimensionless)
    real(dp), intent(in) :: theta !< Scattering angle in radians

    integer, parameter :: l_max = 1000
    integer :: l
    real(dp) :: lambda, sigma_0, sigma_l, sigma_lambda, rho, x, alphap, p_l_0, p_lm2_0, p_lm1_0

    f_c = (0._dp, 0._dp)
    sigma_0 = coulomb_sigma_l(0._dp, eta)
    x = cos(theta)
    alphap = 2*k_cm*eta/m_p*hbar_c

    ! l == 0
    l = 0
    lambda = (-1 + sqrt(1 + 4*l*(l+1) - 4*alpha*alphap))/2._dp
    sigma_l = sigma_0
    sigma_lambda = coulomb_sigma_l(lambda, eta)
    rho = sigma_lambda - sigma_l + (l - lambda)*pi/2._dp
    p_l_0 = legendre_poly(l, 0, x)
    p_lm2_0 = p_l_0
    f_c = (exp(2*i_*rho) - 1)*p_l_0/(2*i_*k_cm)

    ! l == 1
    l = 1
    lambda = (-1 + sqrt(1 + 4*l*(l+1) - 4*alpha*alphap))/2._dp
    sigma_l = coulomb_sigma_l(real(l, kind = dp), eta)
    sigma_lambda = coulomb_sigma_l(lambda, eta)
    rho = sigma_lambda - sigma_l + (l - lambda)*pi/2._dp
    p_l_0 = legendre_poly(l, 0, x)
    p_lm1_0 = p_l_0
    f_c = f_c + (2*l + 1)*exp(2*i_*(sigma_l - sigma_0))*(exp(2*i_*rho) - 1)*p_l_0/(2*i_*k_cm)

    do l = 2, l_max
        lambda = (-1 + sqrt(1 + 4*l*(l+1) - 4*alpha*alphap))/2._dp
        sigma_l = coulomb_sigma_l(real(l, kind = dp), eta)
        sigma_lambda = coulomb_sigma_l(lambda, eta)
        rho = sigma_lambda - sigma_l + (l - lambda)*pi/2._dp
        p_l_0 = (x*(2*l - 1)*p_lm1_0 - (l - 1)*p_lm2_0)/l
        f_c = f_c + (2*l + 1)*exp(2*i_*(sigma_l - sigma_0))*(exp(2*i_*rho) - 1)*p_l_0/(2*i_*k_cm)
        p_lm2_0 = p_lm1_0
        p_lm1_0 = p_l_0
    enddo
end function f_coulomb2

!!
!> @brief      Vacuum polarization scattering amplitude
!!
!! Calculates a series expansion to second leading order of scattering amplitude corresponding to the
!! energy dependent vacuum polarization potential.
!!
!! See equation 22 to 27 in Phys. Rev. C 91, 029901 (2015) for more details
!! 
!! @return     Vacuum polarization scattering amplitude
!!
!! @author     Rodrigo Navarro Perez
!!
complex(dp) function f_vacuum_polarization(k_cm, eta, theta) result(f_vp)
    implicit none
    real(dp), intent(in) :: k_cm !< Center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: eta !< Energy dependent Sommerfeld parameter (dimensionless)
    real(dp), intent(in) :: theta !< Scattering angle in radians

    real(dp) :: f_vp_0, re_f_vp_1, im_f_vp_1, x, nu, y, F, k_MeV

    k_MeV = k_cm*hbar_c
    x = cos(theta)
    nu = 2*m_e**2/(k_MeV**2)
    y = nu/(1 - x)
    F = -5/3._dp + y + sqrt(1 + y)*(1 - 0.5_dp*y)*log((sqrt(1 + y) + 1)/(sqrt(1 + y) -1))
    f_vp_0 = -alpha*eta*F/(3*pi*(1 - x))
    re_f_vp_1 = 4*eta**2*alpha/(3*pi*(1 - x))*sqrt((1 - x)/(1 + x))*(atan(sqrt((1 + x)/(1 - x))) &
        - atan(sqrt((nu*(1 + x))/(2*(1 - x)))))
    im_f_vp_1 = alpha*eta**2/(3*pi*(1 - x))*log(1/y)*(log(k_MeV/m_e) - 3._dp/2._dp*log(2/(1 - x)))
    f_vp = f_vp_0  + re_f_vp_1 + i_*im_f_vp_1
    f_vp = f_vp/k_cm
end function f_vacuum_polarization

!!
!> @brief      pp magnetic moment scattering amplitude
!!
!! Calculates the pp magnetic moment scattering amplitude for a given center of mass momentum, and scattering angle
!!
!! See equations 15 and 29 to 33 in Phys. Rev. C 91, 029901 (2015) for more details
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine mm_amplitudes(k_cm, eta, theta, m_00, m_11, m_10, m_01, m_1m1)
    implicit none
    real(dp), intent(in) :: k_cm !< Center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: eta !< Energy dependent Sommerfeld parameter (dimensionless)
    real(dp), intent(in) :: theta !< Scattering angle in radians
    complex(dp), intent(out) :: m_00 !< \f$ M^1_{00} \f$ pp magnetic moment amplitude
    complex(dp), intent(out) :: m_11 !< \f$ M^1_{11} \f$ pp magnetic moment amplitude
    complex(dp), intent(out) :: m_10 !< \f$ M^1_{10} \f$ pp magnetic moment amplitude
    complex(dp), intent(out) :: m_01 !< \f$ M^1_{01} \f$ pp magnetic moment amplitude
    complex(dp), intent(out) :: m_1m1 !< \f$ M^1_{1-1} \f$ pp magnetic moment amplitude

    integer, parameter :: j_max = 1001
    real(dp) :: mm_phases(1:5, 1:j_max), mm_phases_ls(1:5, 1:j_max)
    character(len=2), parameter :: reaction = 'pp'
    complex(dp) :: z_ls

    call mm_phaseshifts(k_cm, eta, mm_phases)
    call mm_phases_ls_term(k_cm, mm_phases_ls)
    call partial_wave_amplitude_sum(s_matrix_mm, 1, 0, 0, k_cm, theta, reaction, mm_phases, m = m_00)
    call partial_wave_amplitude_sum(s_matrix_mm, 1, 1, 1, k_cm, theta, reaction, mm_phases, m = m_11)
    call partial_wave_amplitude_sum(s_matrix_mm, 1, 1,-1, k_cm, theta, reaction, mm_phases, m = m_1m1)
    mm_phases = mm_phases - mm_phases_ls
    call partial_wave_amplitude_sum(s_matrix_mm, 1, 1, 0, k_cm, theta, reaction, mm_phases, m = m_10)
    call partial_wave_amplitude_sum(s_matrix_mm, 1, 0, 1, k_cm, theta, reaction, mm_phases, m = m_01)
    call m01_ls_contribution(k_cm, eta, theta, mm_phases_ls, m_01)
    z_ls = -m_p*f_ls/(sin(theta)*sqrt(2._dp))*(exp(-i_*eta*log((1 - cos(theta))/2._dp)) &
        + exp(-i_*eta*log((1 + cos(theta))/2._dp)) - 1)*hbar_c
    m_10 = m_10 + z_ls
    m_01 = m_01 - z_ls

end subroutine mm_amplitudes

!!
!> @brief      Spin-Orbit contribution to the pp magnetic moment phase-shifts
!!
!! Calculates only the spin orbit contribution to the pp magnetic moment phase-shifts
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine mm_phases_ls_term(k_cm, mm_phases)
    implicit none
    real(dp), intent(in) :: k_cm !< Center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(out) :: mm_phases(:, :) !< pp magnetic moment phase-shifts

    integer :: l, i
    real(dp) :: k_MeV, I_lp2lp2, I_ll

    mm_phases = 0._dp

    k_MeV = k_cm*hbar_c
    l = -1
    I_lp2lp2 = 1._dp/(2*(l+2)*((l+2)+1._dp))

    mm_phases(5, 1) = m_p*k_MeV*(l + 3)*f_ls*I_lp2lp2
    I_ll = I_lp2lp2

    do i = 2, size(mm_phases, 2) - 1, 2
        l = i - 1
        mm_phases(3, i+1) = -m_p*k_MeV*l*f_ls*I_ll
        I_lp2lp2 = 1/(2._dp*(l + 2)*(l + 3))
        mm_phases(5, i+1) = m_p*k_MeV*(l+3)*f_ls*I_lp2lp2
        I_ll = I_lp2lp2
    enddo
end subroutine mm_phases_ls_term

!!
!> @brief     Spin-Orbit contribution to the \f$ M^1_{01} \f$ pp magnetic moment amplitude
!!
!! Calculate only the spin orbit contribution to the \f$ M^1_{01} \f$ pp magnetic moment amplitude
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine m01_ls_contribution(k_cm, eta, theta, mm_phases, m_01)
    implicit none
    real(dp), intent(in) :: k_cm !< Center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: eta !< Energy dependent Sommerfeld parameter (dimensionless)
    real(dp), intent(in) :: theta !< Scattering angle in radians
    real(dp), intent(in) :: mm_phases(:, :) !< pp magnetic moment phase-shifts without LS terms
    complex(dp), intent(inout) :: m_01 !< \f$ M^1_{01} \f$ pp magnetic moment amplitude

    integer :: l
    real(dp) :: x, sigma_l, sigma_0
    complex(dp) :: sm_lp1, sm_lm1
    character(len=2), parameter :: reaction = 'pp'

    x = cos(theta)

    sigma_0 = coulomb_sigma_l(0._dp, eta)
    do l=1, size(mm_phases, 2) - 2, 2
        sigma_l = coulomb_sigma_l(real(l, kind = dp), eta)
        call s_matrix_mm(l, l, l+1, 1, k_cm, eta, reaction, mm_phases, sm = sm_lp1)
        call s_matrix_mm(l, l, l-1, 1, k_cm, eta, reaction, mm_phases, sm = sm_lm1)
        m_01 = m_01 + legendre_poly(l, 1, x)*(1/(l+1._dp)*sm_lp1 + 1._dp/l*sm_lm1)*exp(2*i_*(sigma_l - sigma_0))&
            /(i_*k_cm*sqrt(2._dp))
    enddo

end subroutine m01_ls_contribution

!!
!> @brief     S matrix approximation for pp magnetic moment amplitude 
!!
!! Gives the \f$ S_{\rm MM, pp} - 1 \approx 2 i K_{\rm MM, pp} \f$
!! 
!! The k_cm, eta, reaction, d_phases, and d_sm variables are only used to have the same
!! interface as the s_matrix subroutine so that both can be given as arguments to 
!! partial_wave_amplitude_sum
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine s_matrix_mm(l_prime, l, j, s, k_cm, eta, reaction, phases, d_phases, sm, d_sm)
    implicit none
    integer, intent(in) :: l_prime !< \f$l'\f$ quantum number
    integer, intent(in) :: l !< \f$l\f$ quantum number
    integer, intent(in) :: j !< \f$j\f$ quantum number
    integer, intent(in) :: s !< \f$s\f$ quantum number
    real(dp), intent(in) :: k_cm !< C.M. momentum in \f$^{-1}\f$
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta'\f$
    character(len=2), intent(in) :: reaction !< reaction channel (pp or np)
    real(dp), intent(in) :: phases(:,:) !< corresponding phase-shifts
    real(dp), optional, intent(in) :: d_phases(:, :, :) !< derivatives of the phase-shifts
    complex(dp), intent(out) :: sm !< matrix element
    complex(dp), optional, intent(out), allocatable :: d_sm(:) !< derivatives of matrix element

    if (present(d_sm) .and. .not. present(d_phases)) then
        stop 'if derivatives of m are requested in s_matrix_mm, derivatives of phases need to be provided'
        ! the lines below are to avoid -Wunused-dummy-argument warning during compilation,
        ! they're never really executed
        if (reaction == 'pp') then
            sm = k_cm + eta
        endif
    endif

    select case(s)
    case (0)
        if (l_prime == l .and. l == j) then
            sm = 2*i_*phases(1, j+1)
        else
            sm = 1 - kronecker_delta(l_prime, l)
        endif
    case (1)
        if (l == l_prime) then
            if (l == j) then
                sm = 2*i_*phases(2, j+1)
            else if (l == j-1) then
                sm = 2*i_*phases(3, j+1)
            else if (l == j+1) then
                sm = 2*i_*phases(5, j+1)
            else
                sm = -1
            endif
        else if (abs(l-l_prime) == 2) then
            sm = 2*i_*phases(4, j+1)
        else
            sm = -kronecker_delta(l_prime, l)
        endif
    case default
        stop 's has to be zero or one in s_matrix_mm'
    end select
end subroutine s_matrix_mm

!!
!> @brief      sums over partial waves to calculate the amplitude
!!
!! Calculates a particular Wolfenstein parameter of the NN scattering amplitude determined by the set of quantum numbers
!! \f$s\f$, \f$m_s\f$ and \f$m_s'\f$
!!
!! The parameter is given by
!! \f[M_{m_s',m_s}^s(\theta) = \frac{1}{2ik} \sum_{J,l',l} \sqrt{4\pi(2l+1)} Y_{m_s'-m_s}^{l'}(\theta,0)
!!    C_{m_s-m_s',m_s',m_s}^{l',s,J} i^{l-l'} (S_{l,l'}^{J,s}-\delta_{l',l})C_{0,m_s,m_s}^{l,s,J},\f]
!!
!! where \f$Y_{m}^{l}(\theta,\varphi)\f$ are the spherical harmonics, \f$C_{m_l,m_s,m_j}^{l,s,j}\f$ are the
!! Clebsch-Gordan coefficients and \f$S_{l,l'}^{j,s}\f$ are the scattering matrix elements (i.e. phase-shifts)
!!
!! See equation 15 in Phys. Rev. C 91, 029901 (2015) for more details
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine partial_wave_amplitude_sum(s_mat, s, ms, mj, k_cm, theta, reaction, phases, d_phases, m, d_m)
    implicit none
    procedure(s_matrix_elements) :: s_mat !< S matrix elements (subroutine)
    integer, intent(in) :: s !< \f$s\f$ quantum number
    integer, intent(in) :: ms !< \f$m_s\f$ quantum number
    integer, intent(in) :: mj !< \f$m_j\f$ quantum number
    real(dp), intent(in) :: k_cm !< C.M. momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: theta !< scattering angle in radians
    character(len=2), intent(in) :: reaction !< reaction channel (pp or np)
    real(dp), intent(in) :: phases(:, :) !< corresponding phase-shifts
    real(dp), optional, intent(in) :: d_phases(:, :, :) !< derivatives of the phase-shifts
    complex(dp), intent(out) :: m !< Wolfenstein parameter
    complex(dp), optional, intent(out), allocatable  :: d_m(:) !< derivatives of the Wolfenstein parameter

    integer :: n_params, j_max, j, l, l_prime
    real(dp) :: etap, sigma_0, sigma_l, sigma_lp, Ylp, cg_1, cg_2
    complex(dp) :: sm
    complex(dp), allocatable :: d_sm(:)

    if (present(d_m) .and. .not. present(d_phases)) then
        stop 'if derivatives of m are requested in partial_wave_amplitude_sum, derivatives of phases need to be provided'
    endif

    m = (0, 0)
    n_params = size(d_phases, 1)
    if (present(d_m)) then
        allocate(d_m(1: n_params))
        d_m = (0, 0)
    endif

    sigma_l = 0._dp
    sigma_lp = 0._dp

    if (reaction == 'pp') then
        etap = eta_prime(k_cm)
        sigma_0 = coulomb_sigma_l(0._dp, etap)
    else
        etap = 0._dp
        sigma_0 = 0._dp
    endif

    j_max = size(phases, 2) - 1

    do j = 0, j_max
        do l = j-1, j+1
            if (reaction == 'pp') sigma_l = coulomb_sigma_l(real(l, kind=dp), etap)
            do l_prime = j-1, j+1
                if (reaction == 'pp') sigma_lp = coulomb_sigma_l(real(l_prime, kind=dp), etap)
                if (l_prime >= abs(mj - ms) .and. l_prime >= 0 .and. l >= 0) then
                    Ylp = real(spherical_harmonic(l_prime, mj-ms, theta, 0._dp))
                    if (present(d_m)) then
                        call s_mat(l_prime, l, j, s, k_cm, etap, reaction, phases, d_phases, sm, d_sm)
                    else
                        call s_mat(l_prime, l, j, s, k_cm, etap, reaction, phases, sm = sm)
                    endif

                    cg_1 = clebsch_gordan(l_prime, s , j, ms, mj)
                    cg_2 = clebsch_gordan(l, s , j, mj, mj)
                    m = m + 2*Ylp*cg_1*(i_**(l - l_prime))*exp(i_*(sigma_lp - sigma_0))*sm &
                        *exp(i_*(sigma_l-sigma_0))*cg_2*sqrt(2*l + 1._dp)*sqrt(4*pi)/(2*i_*k_cm)
                    if (present(d_m)) then
                        d_m = d_m + 2*Ylp*cg_1*(i_**(l - l_prime))*exp(i_*(sigma_lp - sigma_0))*d_sm &
                            *exp(i_*(sigma_l-sigma_0))*cg_2*sqrt(2*l + 1._dp)*sqrt(4*pi)/(2*i_*k_cm)
                    endif
                endif
            enddo
        enddo
    enddo
    if (reaction == 'np') then
        m = 0.5_dp*m
        if (present(d_m)) then
            d_m = 0.5_dp*d_m
        endif
    endif
end subroutine partial_wave_amplitude_sum

!!
!> @brief      Clebsch-Gordan coefficient
!!
!! Calculates the Clebsch-Gordan coefficient \f$C_{ms,mj}^{l,s,j}\f$
!!
!! @return     Clebsch-Gordan coefficient \f$C_{ms,mj}^{l,s,j}\f$
!!
!! @author     Rodrigo Navarro
!!
real(dp) function clebsch_gordan(l_, s_ , j_, ms_, mj_) result(cgc)
    implicit none
    integer, intent(in) :: l_ !< \f$l\f$ quantum number
    integer, intent(in) :: s_ !< \f$s\f$ quantum number
    integer, intent(in) :: j_ !< \f$j\f$ quantum number
    integer, intent(in) :: ms_ !< \f$m_s\f$ quantum number
    integer, intent(in) :: mj_ !< \f$m_j\f$ quantum number
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

!!
!> @brief      S mastrix elements
!!
!! Given a set of quantum number \f$l'\f$, \f$l\f$, \f$j\f$, \f$s\f$ extracts the corresponding element from the NN
!! scattering matrix \f$S - \delta_{l, l'}\f$
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine s_matrix(l_prime, l, j, s, k_cm, eta, reaction, phases, d_phases, sm, d_sm)
    implicit none
    integer, intent(in) :: l_prime !< \f$l'\f$ quantum number
    integer, intent(in) :: l !< \f$l\f$ quantum number
    integer, intent(in) :: j !< \f$j\f$ quantum number
    integer, intent(in) :: s !< \f$s\f$ quantum number
    real(dp), intent(in) :: k_cm !< C.M. momentum in \f$^{-1}\f$
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta'\f$
    character(len=2), intent(in) :: reaction !< reaction channel (pp or np)
    real(dp), intent(in) :: phases(:,:) !< corresponding phase-shifts
    real(dp), optional, intent(in) :: d_phases(:, :, :) !< derivatives of the phase-shifts
    complex(dp), intent(out) :: sm !< matrix element
    complex(dp), optional, intent(out), allocatable :: d_sm(:) !< derivatives of matrix element

    integer :: n_params
    real(dp) :: alphap, lambda, sigma_l, sigma_lambda, rho_l, t_lab, nu, tau00, tau0, lambdap, &
        sigma_lp, sigma_lambdap, rho_lp
    real(dp), allocatable :: mm_phases(:, :)

    if (present(d_sm) .and. .not. present(d_phases)) then
        stop 'if derivatives of sm are requested in s_matrix, derivatives of phases need to be provided'
    endif

    n_params = size(d_phases, 1)
    if (present(d_sm)) then
        allocate(d_sm(1: n_params))
        d_sm = (0, 0)
    endif

    allocate(mm_phases, mold = phases)
    mm_phases = 0

    select case(reaction)
    case('pp')
        alphap = 2*k_cm*eta/m_p*hbar_c
        lambda = (-1 + sqrt(1 + 4*l*(l+1) - 4*alpha*alphap))/2._dp
        sigma_l = coulomb_sigma_l(real(l, kind=dp), eta)
        sigma_lambda = coulomb_sigma_l(lambda, eta)
        rho_l = sigma_lambda - sigma_l + (l - lambda)*pi/2._dp

        lambdap = (-1 + sqrt(1 + 4*l_prime*(l_prime+1) - 4*alpha*alphap))/2._dp
        sigma_lp = coulomb_sigma_l(real(l_prime, kind = dp), eta)
        sigma_lambdap = coulomb_sigma_l(lambdap, eta)
        rho_lp = sigma_lambdap - sigma_lp + (l_prime - lambdap)*pi/2._dp

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
        if (l_prime == l .and. l == j) then
            sm = (exp(2*i_*phases(1, j+1)) - 1)*exp(2*i_*(rho_l + tau0))
            if (present(d_sm)) then
                d_sm = 2*i_*exp(2*i_*phases(1, j+1))*d_phases(:, 1, j+1)*exp(2*i_*(rho_l + tau0))
            endif
        else
            sm = 1 - kronecker_delta(l_prime, l)
            if (present(d_sm)) then
                d_sm = (0, 0)
            endif
        endif
    case (1)
        if (l == l_prime) then
            if (l == j) then
                sm = (exp(2*i_*phases(2, j+1)) - 1)*exp(2*i_*rho_l)*exp(2*i_*mm_phases(2, j+1))
                if (present(d_sm)) then
                    d_sm = 2*i_*exp(2*i_*phases(2, j+1))*d_phases(:, 2, j+1)*exp(2*i_*rho_l)&
                        *exp(2*i_*mm_phases(2, j+1))
                endif
            else if (l == j-1) then
                sm = (cos(2*phases(4, j+1))*exp(2*i_*phases(3, j+1)) - 1)*exp(2*i_*rho_l)&
                    *exp(2*i_*mm_phases(3, j+1))
                if (present(d_sm)) then
                    d_sm = (-2*sin(2*phases(4, j+1))*d_phases(:, 4, j+1)*exp(2*i_*phases(3, j+1)) &
                        + cos(2*phases(4, j+1))*2*i_*exp(2*i_*phases(3, j+1))*d_phases(:, 3, j+1))&
                        *exp(2*i_*rho_l)*exp(2*i_*mm_phases(3, j+1))
                endif
            else if (l == j+1) then
                sm = (cos(2*phases(4, j+1))*exp(2*i_*phases(5, j+1)) - 1)*exp(2*i_*rho_l)&
                    *exp(2*i_*mm_phases(5, j+1))
                if (present(d_sm)) then
                    d_sm = (-2*sin(2*phases(4, j+1))*d_phases(:, 4, j+1)*exp(2*i_*phases(5, j+1)) &
                        + cos(2*phases(4, j+1))*2*i_*exp(2*i_*phases(5, j+1))*d_phases(:, 5, j+1))&
                        *exp(2*i_*rho_l)*exp(2*i_*mm_phases(5, j+1))
                endif
            else
                sm = -1
                if (present(d_sm)) then
                    d_sm = (0, 0)
                endif
            endif
        else if (abs(l-l_prime) == 2) then
            sm = i_*sin(2*phases(4, j+1))*exp(i_*(phases(3, j+1) + phases(5, j+1)))&
                *exp(i_*(rho_l + rho_lp))*exp(i_*(mm_phases(3, j+1) + mm_phases(5, j+1)))
            if (present(d_sm)) then
                d_sm = i_*(2*cos(2*phases(4, j+1))*d_phases(:, 4, j+1)&
                    *exp(i_*(phases(3, j+1) + phases(5, j+1))) + sin(2*phases(4, j+1))*i_&
                    *exp(i_*(phases(3, j+1) + phases(5, j+1)))*(d_phases(:, 3, j+1) &
                    + d_phases(:, 5, j+1)))*exp(i_*(rho_l+rho_lp)) &
                    *exp(i_*(mm_phases(3, j+1) + mm_phases(5, j+1)))
            endif
        else
            sm = -kronecker_delta(l_prime, l)
            if (present(d_sm)) then
                d_sm = (0, 0)
            endif
        endif
    case default
        stop 's has to be zero or one in s_matrix'
    end select
end subroutine s_matrix

!!
!> @brief      magnetic moment pp phase-shifts
!!
!! Calculates the partial-wave \f$K\f$ matrix elements for the magnetic moment pp amplitude
!!
!! The \f$K\f$ matrix is defined as \f$S-1=2iK(1-iK)^{-1} \f$ where \f$S\f$ is the the scattering matrix determined by
!! the pp magnetic moment potential.
!!
!! See Eq. 29 in PhysRevC.88.064002 for more details.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine mm_phaseshifts(k_cm, eta, mm_phases)
    implicit none
    real(dp), intent(in) :: k_cm !< C.M. momentum in \f$^{-1}\f$
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta\f$
    real(dp), intent(out) :: mm_phases(:, :) !< magnetic momentum phase-shifts

    real(dp) :: I_ll, I_lp2lp2, I_llp2 !f_T, f_ls,
    integer :: l, j_max, i

    j_max = size(mm_phases, 2)


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


!!
!> @brief      \f$I_{l,l}\f$ term for mm phase-shifts
!!
!! Calculates the \f$I_{l,l}\f$ term necessary for the calculation of the magnetic moment phases
!!
!! The \f$I_{l,l}\f$ term results as an integral of the \f$1/r^3\f$ dependence on the magnetic moment potential and is
!! given by
!!
!! \f[I_{l,l} = \frac{1}{2l(l+1)} + \frac{1-\pi \eta + \pi \eta \coth(\pi \eta)-2\eta^2\sum_{n=0}^l (n^2 + \eta^2)^{-1}}
!!   {2l(l+1)(2l+1)}\f]
!!
!! See Eq. 31 of PhysRevC.88.064002 for more details.
!!
!! @return     \f$I_{l,l}\f$ term for mm phase-shifts
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function mm_coulomb_Ill(l, eta) result(I_ll)
    implicit none
    integer, intent(in) :: l !< \f$l\f$ quantum number
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta\f$
    integer :: i
    real(dp) ::  sum
    sum = 0._dp
    do i = 0, l
       sum = sum + 1._dp/(i**2 + eta**2)
    endDo
    I_ll = 1._dp/(2._dp*l*(l+1)) + (1 - pi*eta + pi*eta/tanh(pi*eta) - 2*eta**2*sum)/&
        (2._dp*l*(l+1)*(2*l+1))
end function mm_coulomb_Ill

!!
!> @brief      \f$I_{l,l+2}\f$ term for mm phase-shifts
!!
!! Calculates the \f$I_{l,l+2}\f$ term necessary for the calculation of the magnetic moment phases
!!
!! The \f$I_{l,l+2}\f$ term results as an integral of the \f$1/r^3\f$ dependence on the magnetic moment potential and
!! is given by
!!
!! \f[I_{l,l+2} = \frac{1}{6} |l+1+i\eta|^{-1} |l+2+i\eta|^{-1}\f]
!!
!! @return     \f$I_{l,l+2}\f$ term for mm phase-shifts
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function mm_coulomb_Illp2(l, eta) result(I_llp2)
    implicit none
    integer, intent(in) :: l !< \f$l\f$ quantum number
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta\f$
    I_llp2 = 1._dp/(6*abs(l + 1 +i_*eta)*abs(l + 2 + i_*eta))
end function mm_coulomb_Illp2

!!
!> @brief      Coulomb phase-shift
!!
!! Calculates the Coulomb phase-shift \f$\sigma_l \f$
!!
!! The phase-shift is given by
!! \f[\sigma_l = {\rm arg} \Gamma(l+1+i\eta)\f]
!!
!! @return     Coulomb phase-shift \f$\sigma_l\f$
!!
real(dp) function coulomb_sigma_l(l, eta) result(sig_l)
    implicit none
    real(dp), intent(in) :: l !< \f$l\f$ quantum number
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta\f$
    complex(dp) :: z
    z = 1 + l + i_*eta
    sig_l = aimag(cmplx_log_gamma(z))
end function coulomb_sigma_l

! The functions below were written to test the derivatives of saclay_parameters against numerical calculations
! and require the av18 module. All analytic calculations of the derivatives match the numerical ones.
! since the amplitudes module should not depend on a specific module for the NN interaction the functions
! are commented and left here for reference.




end module amplitudes
