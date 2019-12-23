module nn_scattering

use precisions, only : dp
use num_recipes, only : sphbes
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi

implicit none

private

public :: all_phaseshifts

interface
    subroutine nn_potential(ap, r, reaction, v_pw, dv_pw)
        use precisions, only : dp
        implicit none
        real(dp), intent(in) :: ap(:)
        real(dp), intent(in) :: r
        character(len=2), intent(in) :: reaction
        real(dp), intent(out) :: v_pw(:, :)
        real(dp), allocatable, intent(out) :: dv_pw(:, :, :)
    end subroutine nn_potential
end interface

contains

subroutine all_phaseshifts(model, params, t_lab, reaction, r_max, dr, phases, d_phases)
    implicit none
    procedure(nn_potential) :: model
    real(dp), intent(in) :: params(:)
    real(dp), intent(in) :: t_lab
    character(len=2), intent(in) :: reaction
    real(dp), intent(in) :: r_max
    real(dp), intent(in) :: dr
    real(dp), intent(out) :: phases(:, :)
    real(dp), allocatable, intent(out) :: d_phases(:, :, :)

    integer :: n_params, n_waves, j_max
    integer :: ij, j
    real(dp) :: r, k, mu, alfa_1, alfa_2, ps_eigen(1:3), ps_bar(1:3)
    real(dp), allocatable :: v_pw(:, :), dv_pw(:, :, :)
    real(dp), allocatable :: singlets(:), triplets(:)
    real(dp), allocatable :: a1(:), a2(:), b1(:), b2(:), c1(:), c2(:), d1(:), d2(:)

    n_params = size(params)
    n_waves = size(phases, 1)
    j_max = size(phases, 2)

    if (n_waves /= 5) stop 'incorrect number of waves for v_pw in all_phaseshifts'

    allocate(d_phases(1:n_params, 1:n_waves, 1:j_max))
    d_phases = 0._dp
    allocate(v_pw(1:n_waves, 1:j_max))
    allocate(singlets(1:j_max))
    singlets = 0
    allocate(triplets, source = singlets)
    allocate(a1(1:j_max - 1))
    a1 = 0
    allocate(a2, b1, b2, c1, c2, d1, d2, source = a1)
    a1 = 1
    c2 = 1

    select case (reaction)
    case ('pp')
        mu = m_p
    case ('np')
        mu = 2*m_p*m_n/(m_p + m_n)
    case ('nn')
        mu = m_n
    case default
        stop 'incorrect reaction channel in all_phaseshifts'
    end select

    k = momentum_cm(t_lab, reaction)
    r = dr/2
    do
        if( r > r_max) exit
        call model(params, r, reaction, v_pw, dv_pw)
        v_pw = v_pw*mu*dr/(hbar_c**2)
        dv_pw = dv_pw*mu*dr/(hbar_c**2)
        call uncoupled_variable_phase(0, k, r, v_pw(1, 1), singlets(1))
        call uncoupled_variable_phase(1, k, r, v_pw(5, 1), triplets(1))
        do ij = 2, j_max
            j = ij - 1
            call uncoupled_variable_phase(j, k, r, v_pw(1, ij), singlets(ij))
            call uncoupled_variable_phase(j, k, r, v_pw(2, ij), triplets(ij))
            call coupled_variable_phase(j, k, r, v_pw(3:5, ij), a1(j), b1(j), c1(j), d1(j))
            call coupled_variable_phase(j, k, r, v_pw(3:5, ij), a2(j), b2(j), c2(j), d2(j))
        enddo
        r = r + dr
    enddo
    phases(1, :) = atan(singlets(1))
    phases(5, 1) = atan(triplets(1))
    do ij = 2, j_max
        j = ij - 1
        phases(1, ij) = atan(singlets(ij))
        phases(2, ij) = atan(triplets(ij))
        call solve_alfas(a1(j), b1(j), c1(j), d1(j), a2(j), b2(j), c2(j), d2(j), alfa_1, alfa_2)
        ps_eigen = eigen_phases(a1(j), b1(j), c1(j), d1(j), a2(j), b2(j), c2(j), d2(j), alfa_1, alfa_2)
        ps_bar = eigen_2_bar(ps_eigen)
        phases(3:5, ij) = ps_bar
    enddo

end subroutine all_phaseshifts


subroutine coupled_variable_phase(j, k, r, lambdas, a, b, c, d)
    implicit none
    integer, intent(in) :: j
    real(dp), intent(in) :: k
    real(dp), intent(in) :: r
    real(dp), intent(in) :: lambdas(1:3)
    real(dp), intent(inout) :: a
    real(dp), intent(inout) :: b
    real(dp), intent(inout) :: c
    real(dp), intent(inout) :: d

    real(dp) :: j_hat_m1, j_hat_p1, y_hat_m1, y_hat_p1, lambda_jm1, lambda_j, lambda_jp1, &
        lin_comb_ab, lin_comb_cd, diff_b, diff_d, sj, sy, sjp, syp

    lambda_jm1 = lambdas(1)
    lambda_j   = lambdas(2)
    lambda_jp1 = lambdas(3)

    call sphbes(j-1, r*k, sj, sy, sjp, syp)
    j_hat_m1 = sj*r*k
    y_hat_m1 = sy*r*k
    call sphbes(j+1, r*k, sj, sy, sjp, syp)
    j_hat_p1 = sj*r*k
    y_hat_p1 = sy*r*k

    lin_comb_ab = a*j_hat_m1 + b*y_hat_m1
    lin_comb_cd = c*j_hat_p1 + d*y_hat_p1
    diff_b = (lambda_jm1*lin_comb_ab + lambda_j*lin_comb_cd)/k
    diff_d = (lambda_jp1*lin_comb_cd + lambda_j*lin_comb_ab)/k

    a = a - diff_b*y_hat_m1
    b = b + diff_b*j_hat_m1
    c = c - diff_d*y_hat_p1
    d = d + diff_d*j_hat_p1

end subroutine coupled_variable_phase

!!
!> @brief      variable phase equation for uncoupled waves
!!
!! Applies the variable phase equation to continue the integration of the uncoupled partial waves
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine uncoupled_variable_phase(l, k, r, lambda, tan_delta)
    implicit none
    integer, intent(in) :: l !< orbital angular momentum quantum number
    real(dp), intent(in) ::  k !< center of mass momentum (in units of fm\f$^{-1}\f$)
    real(dp), intent(in) :: r !< radius at which the matching is being done
    real(dp), intent(in) :: lambda !< lambda from the potential
    real(dp), intent(inout) :: tan_delta !< tangent of the wave function

    real(dp) :: j_hat, y_hat, phi, denominator, sj, sy, sjp, syp

    call sphbes(l, r*k, sj, sy, sjp, syp)
    j_hat = sj*r*k
    y_hat = sy*r*k
    phi = j_hat - tan_delta*y_hat
    denominator = 1 - lambda*y_hat*phi/k
    tan_delta = (tan_delta - lambda*j_hat*phi/k)/denominator
end subroutine uncoupled_variable_phase

!!
!> @brief      center of mass momentum
!!
!! Given the energy in the LAB frame and the type of reaction, calculates the center of mass
!! momentum in units of fm\f$^{-1}\f$
!!
!! @returns    center of mass in fm\f$^{-1}\f$ 
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function momentum_cm(t_lab, reaction) result(k)
    implicit none
    real(dp), intent(in) :: t_lab !< Energy in LAB frame in MeV
    character(len=2), intent(in) :: reaction !< Type of reaction. pp, np or nn
    real(dp) :: mu
    select case (reaction)
    case ('pp')
        mu = m_p/2
        k = sqrt(mu*t_lab)/hbar_c
    case ('np')
        k = sqrt((m_p**2*t_lab*(t_lab + 2*m_n))/((m_p + m_n)**2 + 2*t_lab*m_p))/hbar_c 
    case ('nn')
        mu = m_n/2
        k = sqrt(mu*t_lab)/hbar_c        
    case default
        stop 'incorrect reaction channel in av18_all_partial_waves'
    end select
end function momentum_cm

subroutine solve_alfas(a1, b1, c1, d1, a2, c2, b2, d2, alfa_1, alfa_2)
    implicit none
    real(dp), intent(in) :: a1, b1, c1, d1, a2, c2, b2, d2
    real(dp), intent(out) :: alfa_1, alfa_2
    real(dp) :: ar, br, cr, radical 
    ar = a2*d2 - c2*b2
    br = a1*d2 + a2*d1 - c1*b2 - c2*b1
    cr = a1*d1 - c1*b1
    radical = br**2 - 4*ar*cr
    if (ar == 0._dp) then
        alfa_1 = -cr/br
        alfa_2 = alfa_1
    elseif (radical > 0) then
        alfa_1 = (-br + sqrt(radical))/(2*ar)
        alfa_2 = (-br - sqrt(radical))/(2*ar)
    else
        print*, 'WARNING: No real solutions in solve_alfas'
        alfa_1 = 0
        alfa_2 = 0
    endif
end subroutine solve_alfas

function eigen_phases(a1, b1, c1, d1, a2, c2, b2, d2, alfa_1, alfa_2) result(ps_eigen)
    implicit none
    real(dp), intent(in) :: a1, b1, c1, d1, a2, c2, b2, d2, alfa_1, alfa_2
    real(dp) :: ps_eigen(1:3)
    ps_eigen(1) = -atan((b1 + alfa_1*b2)/(a1 + alfa_1*a2))
    ps_eigen(2) =  atan((d1 + alfa_1*d2)/(b1 + alfa_1*b2))
    ps_eigen(3) = -atan((d1 + alfa_2*d2)/(c1 + alfa_2*c2))
end function eigen_phases

function eigen_2_bar(ps_eigen) result(ps_bar)
    implicit none
    real(dp), intent(in) :: ps_eigen(1:3)
    real(dp) :: ps_bar(1:3)

    real(dp) :: sin_2eps, sin_diff, arg_arcsin, numerator, denominator, fraction

    sin_2eps = sin(2*ps_eigen(2))
    sin_diff = sin(ps_eigen(1) - ps_eigen(3))
    arg_arcsin = sin_2eps*sin_diff
    ps_bar(2) = asin(arg_arcsin)/2

    numerator = tan(2*ps_bar(2))
    denominator = tan(2*ps_eigen(2))
    fraction = numerator/denominator

    if (ps_eigen(1) - ps_eigen(3) > pi/2) then
        ps_bar(1) = (ps_eigen(1) + ps_eigen(3) + pi - asin(fraction))/2
        ps_bar(3) = (ps_eigen(1) + ps_eigen(3) - pi + asin(fraction))/2
    else
        ps_bar(1) = (ps_eigen(1) + ps_eigen(3) + asin(fraction))/2
        ps_bar(3) = (ps_eigen(1) + ps_eigen(3) - asin(fraction))/2
    endif

end function eigen_2_bar

end module nn_scattering