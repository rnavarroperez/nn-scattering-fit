module amplitudes
use precisions, only : dp
use nn_phaseshifts, only : eta_prime
use constants, only : i_, m_p => proton_mass, hbar_c, alpha, pi, m_e => electron_mass, &
    mu_p => mu_proton
use num_recipes, only : cmplx_log_gamma, spherical_harmonic, kronecker_delta
implicit none

private

public :: saclay_amplitudes!, f_amplitudes, df_amplitudes

contains

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
    complex(dp), intent(out) :: a !< Sacalay parameter \f$a\f$
    complex(dp), intent(out) :: b !< Sacalay parameter \f$b\f$
    complex(dp), intent(out) :: c !< Sacalay parameter \f$c\f$
    complex(dp), intent(out) :: d !< Sacalay parameter \f$d\f$
    complex(dp), intent(out) :: e !< Sacalay parameter \f$e\f$
    complex(dp), intent(out), allocatable :: d_a(:) !< derivatives of the Sacalay parameter \f$a\f$
    complex(dp), intent(out), allocatable :: d_b(:) !< derivatives of the Sacalay parameter \f$b\f$
    complex(dp), intent(out), allocatable :: d_c(:) !< derivatives of the Sacalay parameter \f$c\f$
    complex(dp), intent(out), allocatable :: d_d(:) !< derivatives of the Sacalay parameter \f$d\f$
    complex(dp), intent(out), allocatable :: d_e(:) !< derivatives of the Sacalay parameter \f$e\f$

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
subroutine partial_wave_amplitude_sum(s, ms, mj, k_cm, theta, reaction, phases, d_phases, m, d_m)
    implicit none
    integer, intent(in) :: s !< \f$s\f$ quantum number
    integer, intent(in) :: ms !< \f$m_s\f$ quantum number
    integer, intent(in) :: mj !< \f$m_j\f$ quantum number
    real(dp), intent(in) :: k_cm !< C.M. momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: theta !< scattering angle in radians
    character(len=2), intent(in) :: reaction !< reaction channel (pp or np)
    real(dp), intent(in) :: phases(:, :) !< corresponding phase-shifts
    real(dp), intent(in) :: d_phases(:, :, :) !< derivatives of the phase-shifts
    complex(dp), intent(out) :: m !< Wolfenstein parameter
    complex(dp), intent(out), allocatable  :: d_m(:) !< derivatives of the Wolfenstein parameter

    integer :: n_params, j_max, j, l, lp
    real(dp) :: etap, sigma_0, sigma_l, sigma_lp, Ylp, cg_1, cg_2
    complex(dp) :: sm
    complex(dp), allocatable :: d_sm(:)

    m = (0, 0)
    n_params = size(d_phases, 1)
    allocate(d_m(1: n_params))
    d_m = (0, 0)

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
    if (reaction == 'np') then
        m = 0.5_dp*m
        d_m = 0.5_dp*d_m
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
!! scattering matrix \f$S - \detla_{l, l'}\f$
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine s_matrix(lp, l, j, s, phases, d_phases, k_cm, eta, reaction, sm, d_sm)
    implicit none
    integer, intent(in) :: lp !< \f$l'\f$ quantum number
    integer, intent(in) :: l !< \f$l\f$ quantum number
    integer, intent(in) :: j !< \f$j\f$ quantum number
    integer, intent(in) :: s !< \f$s\f$ quantum number
    real(dp), intent(in) :: phases(:,:) !< corresponding phase-shifts
    real(dp), intent(in) :: d_phases(:, :, :) !< derivatives of the phase-shifts
    real(dp), intent(in) :: k_cm !< C.M. momentum in \f$^{-1}\f$
    real(dp), intent(in) :: eta !< Sommerfeld parameter \f$\eta'\f$
    character(len=2), intent(in) :: reaction !< reaction channel (pp or np)
    complex(dp), intent(out) :: sm !< matrix element
    complex(dp), intent(out), allocatable :: d_sm(:) !< derivatives of matrix element

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


! !!
! !> @brief      wrapper function for saclay_amplitudes
! !!
! !! This wrapper function is used to test the derivatives of the saclay_amplitudes subroutine.
! !! The generic data of type context is used to receive all the arguments necessary to call 
! !! saclay_amplitudes. The same data of type context is used to receive which parameter will
! !! be varied by the dfridr subroutine and which partial wave will be returned.
! !!
! !! @returns    NN phase-shift at an specific lab energy and partial wave
! !!
! !! @author     Rodrigo Navarro Perez
! !!
! real(dp) function f_amplitudes(x, data) result(r)
!     use num_recipes, only : context
!     use av18, only : av18_all_partial_waves
!     use nn_phaseshifts, only : all_phaseshifts
!     implicit none
!     real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
!     type(context), intent(in) :: data !< data structure with all the arguments for saclay_amplitudes

!     real(dp), allocatable :: ap(:)
!     real(dp) :: t_lab, r_max, dr, theta, k_cm, phases(1:5,1:20)
!     real(dp), allocatable :: d_phases(:, :, :)
!     integer :: i_target, i_parameter
!     character(len=2) :: reaction
!     complex(dp) :: a, b, c, d, e
!     complex(dp), allocatable, dimension(:) :: d_a, d_b, d_c, d_d, d_e

!     allocate(ap, source = data%x)
!     t_lab = data%a
!     r_max = data%b 
!     dr = data%c
!     theta = data%d
!     reaction = trim(data%string)
!     i_parameter = data%i
!     i_target = data%j

!     ap(i_parameter) = x
!     call all_phaseshifts(av18_all_partial_waves, ap, t_lab, reaction, r_max, dr, k_cm, phases, d_phases)
!     call saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, d_c, d_d, d_e)

!     select case (i_target)
!     case (1)
!         r = real(a)
!     case (2)
!         r = aimag(a)
!     case (3)
!         r = real(b)
!     case (4)
!         r = aimag(b)
!     case (5)
!         r = real(c)
!     case (6)
!         r = aimag(c)
!     case (7)
!         r = real(d)
!     case (8)
!         r = aimag(d)
!     case (9)
!         r = real(e)
!     case (10)
!         r = aimag(e)
!     case default
!         r = 0
!     end select
! end function f_amplitudes

! !!
! !> @brief      wrapper function for the derivatives of saclay_amplitudes
! !!
! !! This wrapper function is used to test the derivatives of the saclay_amplitudes subroutine.
! !! The generic data of type context is used to receive all the arguments necessary to call 
! !! saclay_amplitudes. The same data of type context is used to receive which parameter will
! !! be varied by the dfridr subroutine and which partial wave will be returned.
! !!
! !! @returns    derivatives of a NN phase-shift at an specific lab energy and partial wave
! !!
! !! @author     Rodrigo Navarro Perez
! !!
! function df_amplitudes(data) result(r)
!     use num_recipes, only : context
!     use av18, only : av18_all_partial_waves
!     use nn_phaseshifts, only : all_phaseshifts
!     implicit none
!     type(context), intent(in) :: data !< data structure with all the arguments for saclay_amplitudes
!     real(dp), allocatable :: r(:)

!     real(dp), allocatable :: ap(:)
!     real(dp) :: t_lab, r_max, dr, theta, k_cm, phases(1:5, 1:20)
!     real(dp), allocatable :: d_phases(:, :, :)
!     integer :: i_target, i_parameter
!     character(len=2) :: reaction 
!     complex(dp) :: a, b, c, d, e
!     complex(dp), allocatable, dimension(:) :: d_a, d_b, d_c, d_d, d_e

!     allocate(ap, source = data%x)
!     t_lab = data%a
!     r_max = data%b 
!     dr = data%c
!     theta = data%d
!     reaction = trim(data%string)
!     i_parameter = data%i
!     i_target = data%j

!     call all_phaseshifts(av18_all_partial_waves, ap, t_lab, reaction, r_max, dr, k_cm, phases, d_phases)
!     call saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, d_c, d_d, d_e)

!     allocate(r, mold = ap)

!     select case (i_target)
!     case (1)
!         r = real(d_a)
!     case (2)
!         r = aimag(d_a)
!     case (3)
!         r = real(d_b)
!     case (4)
!         r = aimag(d_b)
!     case (5)
!         r = real(d_c)
!     case (6)
!         r = aimag(d_c)
!     case (7)
!         r = real(d_d)
!     case (8)
!         r = aimag(d_d)
!     case (9)
!         r = real(d_e)
!     case (10)
!         r = aimag(d_e)
!     case default
!         r = 0
!     end select
! end function df_amplitudes

    
end module amplitudes