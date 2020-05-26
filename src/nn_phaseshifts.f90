module nn_phaseshifts

use precisions, only : dp
use num_recipes, only : sphbes
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi, alpha
use coulombwf, only : coul90
implicit none

private

public :: all_phaseshifts, eta_prime, momentum_cm, nn_potential!, f_all_phaseshifts, df_all_phaseshifts

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

!!
!> @brief      phaseshifts for all partial waves
!!
!! given a model for a local NN potential with its corresponding phenomenological parameters,
!! uses the variable phase equation to integrate the reduced Schrödinger equation and calculates
!! the phaseshifts (i.e. S-matrix) at all partial waves up to a maximum total angular momentum
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine all_phaseshifts(model, params, t_lab, reaction, r_max, dr, phases, d_phases)
    implicit none
    procedure(nn_potential) :: model !< local NN potential (subroutine)
    real(dp), intent(in) :: params(:) !< phenomenological parameters for the local potential
    real(dp), intent(in) :: t_lab !< laboratory energy of the scattering in MeV
    character(len=2), intent(in) :: reaction !< reaction channel (pp, np or nn)
    real(dp), intent(in) :: r_max !< maximum radius of integration
    real(dp), intent(in) :: dr !< integration step
    real(dp), intent(out) :: phases(:, :) !< phaseshifts in all partial waves
    real(dp), allocatable, intent(out) :: d_phases(:, :, :) !< derivatives of the partial waves with respect to the potential parameters

    real(dp) :: k_cm !< exchange momentum in center of mass frame. In units of fm \f$^{-1}\f$
    integer :: n_params, n_waves, j_max
    integer :: ij, j
    real(dp) :: r, mu, alfa_1, alfa_2, ps_eigen(1:3)
    real(dp), allocatable :: v_pw(:, :), dv_pw(:, :, :)
    real(dp), allocatable :: singlets(:), triplets(:), d_singlets(:, :), d_triplets(:, :)
    real(dp), allocatable :: a1(:), a2(:), b1(:), b2(:), c1(:), c2(:), d1(:), d2(:)
    real(dp), allocatable, dimension(:, :) :: d_a1, d_a2, d_b1, d_b2, d_c1, d_c2, d_d1, d_d2, d_ps_eigen
    real(dp), allocatable, dimension(:) :: d_alfa_1, d_alfa_2

    n_params = size(params)
    n_waves = size(phases, 1)
    j_max = size(phases, 2)

    phases = 0

    if (n_waves /= 5) stop 'incorrect number of waves for v_pw in all_phaseshifts'

    allocate(d_phases(1:n_params, 1:n_waves, 1:j_max))
    d_phases = 0
    allocate(v_pw(1:n_waves, 1:j_max))
    allocate(singlets(1:j_max))
    singlets = 0
    allocate(triplets, source = singlets)
    allocate(d_singlets(1:n_params, 1:j_max))
    d_singlets = 0
    allocate(d_triplets, source = d_singlets)
    allocate(a1(1:j_max - 1))
    a1 = 0
    allocate(a2, b1, b2, c1, c2, d1, d2, source = a1)
    allocate(d_a1(1:n_params, 1:j_max - 1))
    d_a1 = 0
    allocate(d_a2, d_b1, d_b2, d_c1, d_c2, d_d1, d_d2, source = d_a1)

    a1 = 1
    c2 = 1

    allocate(d_alfa_1(1:n_params))
    d_alfa_1 = 0
    allocate(d_alfa_2, source = d_alfa_1)

    allocate(d_ps_eigen(1:n_params, 1:3))
    d_ps_eigen = 0

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
    k_cm = momentum_cm(t_lab, reaction)
    r = dr/2
    do
        if( r > r_max) exit
        call model(params, r, reaction, v_pw, dv_pw)
        if (reaction == 'pp') then
            call add_coulomb(r, k_cm, v_pw)
        endif
        v_pw = v_pw*mu*dr/(hbar_c**2)
        dv_pw = dv_pw*mu*dr/(hbar_c**2)
        call uncoupled_variable_phase(0, k_cm, r, v_pw(1, 1), dv_pw(:, 1, 1), singlets(1), d_singlets(:,1))
        call uncoupled_variable_phase(1, k_cm, r, v_pw(5, 1), dv_pw(:, 5, 1), triplets(1), d_triplets(:,1))
        do ij = 2, j_max
            j = ij - 1
            if (reaction == 'np') then
                call uncoupled_variable_phase(j, k_cm, r, v_pw(1, ij), dv_pw(:, 1, ij), singlets(ij), d_singlets(:,ij))
                call uncoupled_variable_phase(j, k_cm, r, v_pw(2, ij), dv_pw(:, 2, ij), triplets(ij), d_triplets(:,ij))
                call coupled_variable_phase(j, k_cm, r, v_pw(3:5, ij), dv_pw(:, 3:5, ij), a1(j), b1(j), &
                    c1(j), d1(j), d_a1(:, j), d_b1(:, j), d_c1(:, j), d_d1(:, j))
                call coupled_variable_phase(j, k_cm, r, v_pw(3:5, ij), dv_pw(:, 3:5, ij), a2(j), b2(j), &
                    c2(j), d2(j), d_a2(:, j), d_b2(:, j), d_c2(:, j), d_d2(:, j))
            elseif (mod(j, 2) == 1) then
                call uncoupled_variable_phase(j, k_cm, r, v_pw(2, ij), dv_pw(:, 2, ij), triplets(ij), d_triplets(:,ij))
            else
                call uncoupled_variable_phase(j, k_cm, r, v_pw(1, ij), dv_pw(:, 1, ij), singlets(ij), d_singlets(:,ij))
                call coupled_variable_phase(j, k_cm, r, v_pw(3:5, ij), dv_pw(:, 3:5, ij), a1(j), b1(j), &
                    c1(j), d1(j), d_a1(:, j), d_b1(:, j), d_c1(:, j), d_d1(:, j))
                call coupled_variable_phase(j, k_cm, r, v_pw(3:5, ij), dv_pw(:, 3:5, ij), a2(j), b2(j), &
                    c2(j), d2(j), d_a2(:, j), d_b2(:, j), d_c2(:, j), d_d2(:, j))
            endif
        enddo
        r = r + dr
    enddo
    if (reaction == 'pp') then
        call model(params, r, reaction, v_pw, dv_pw)
        call add_coulomb(r, k_cm, v_pw)
        v_pw = v_pw*mu*dr/(hbar_c**2)
        dv_pw = dv_pw*mu*dr/(hbar_c**2)
        v_pw(2, 1) = v_pw(5, 1)
        dv_pw(:, 2 , 1) = dv_pw(:, 5, 1)
        call match_uncoupled_waves(0, k_cm, r, v_pw(1, :), dv_pw(:, 1, :), singlets, d_singlets)
        call match_uncoupled_waves(1, k_cm, r, v_pw(2, :), dv_pw(:, 2, :), triplets, d_triplets)
        call match_coupled_waves(k_cm, r, v_pw(3:5, :), dv_pw(:, 3:5, :), a1, b1, c1, d1, d_a1, d_b1, d_c1, d_d1)
        call match_coupled_waves(k_cm, r, v_pw(3:5, :), dv_pw(:, 3:5, :), a2, b2, c2, d2, d_a2, d_b2, d_c2, d_d2)
    endif
    phases(1, 1) = atan(singlets(1))
    phases(5, 1) = atan(triplets(1))
    d_phases(:, 1, 1) = 1/(1 + singlets(1)**2)*d_singlets(:, 1)
    d_phases(:, 5, 1) = 1/(1 + triplets(1)**2)*d_triplets(:, 1)
    do ij = 2, j_max
        j = ij - 1
        phases(1, ij) = atan(singlets(ij))
        phases(2, ij) = atan(triplets(ij))
        d_phases(:, 1, ij) = 1/(1 + singlets(ij)**2)*d_singlets(:, ij)
        d_phases(:, 2, ij) = 1/(1 + triplets(ij)**2)*d_triplets(:, ij)
        if (reaction == 'np' .or. mod(j,2) == 0) then
            call solve_alfas(a1(j), b1(j), c1(j), d1(j), a2(j), b2(j), c2(j), d2(j), d_a1(:, j), &
                d_b1(:, j), d_c1(:, j), d_d1(:, j), d_a2(:, j), d_b2(:, j), d_c2(:, j), d_d2(:, j), &
                alfa_1, alfa_2, d_alfa_1, d_alfa_2)
            call eigen_phases(a1(j), b1(j), c1(j), d1(j), a2(j), b2(j), c2(j), d2(j), d_a1(:, j), &
                d_b1(:, j), d_c1(:, j), d_d1(:, j), d_a2(:, j), d_b2(:, j), d_c2(:, j), d_d2(:, j), &
                alfa_1, alfa_2, d_alfa_1, d_alfa_2, ps_eigen, d_ps_eigen)
            if (ps_eigen(1) < 0 .and. t_lab < 30 .and. j == 1) ps_eigen(1) = ps_eigen(1) + pi
            call eigen_2_bar(ps_eigen, d_ps_eigen, phases(3:5, ij), d_phases(:, 3:5, ij))
        endif
    enddo

end subroutine all_phaseshifts


!!
!> @brief      variable phase equation integration in coupled channels
!!
!! Performs a single step in the integration of the variable phase equation for coupled channels
!! with a total angular momentum j
!!
!! See equation B27 in Phys. Rev. C 88 (2013) 064002
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine coupled_variable_phase(j, k, r, lambdas, d_lambdas, a, b, c, d, d_a, d_b, d_c, d_d)
    implicit none
    integer, intent(in) :: j !< total angular momentum quantum number
    real(dp), intent(in) :: k !< center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< integration radius in fm
    real(dp), intent(in) :: lambdas(1:3) !< lambda strength coefficients in fm\f$^{-2}\f$
    real(dp), intent(in) :: d_lambdas(:, :) !< derivatives of the lambda strength coefficients
    real(dp), intent(inout) :: a !< \f$A\f$ parameter
    real(dp), intent(inout) :: b !< \f$B\f$ parameter
    real(dp), intent(inout) :: c !< \f$C\f$ parameter
    real(dp), intent(inout) :: d !< \f$D\f$ parameter
    real(dp), intent(inout) :: d_a(:) !< derivatives of the \f$A\f$ parameter
    real(dp), intent(inout) :: d_b(:) !< derivatives of the \f$B\f$ parameter
    real(dp), intent(inout) :: d_c(:) !< derivatives of the \f$C\f$ parameter
    real(dp), intent(inout) :: d_d(:) !< derivatives of the \f$D\f$ parameter

    real(dp) :: j_hat_m1, j_hat_p1, y_hat_m1, y_hat_p1, lambda_jm1, lambda_j, lambda_jp1, &
        lin_comb_ab, lin_comb_cd, diff_b, diff_d, sj, sy, sjp, syp
    real(dp), allocatable, dimension(:) :: d_lambda_jm1, d_lambda_j, d_lambda_jp1, &
        d_lin_comb_ab, d_lin_comb_cd, d_diff_b, d_diff_d
    integer n_params

    n_params = size(d_lambdas, 1)
    allocate(d_lambda_jm1(1:n_params))
    d_lambda_jm1 = 0
    allocate(d_lambda_j, d_lambda_jp1, d_lin_comb_ab, d_lin_comb_cd, d_diff_b, d_diff_d, &
        source = d_lambda_jm1)

    lambda_jm1 = lambdas(1)
    lambda_j   = lambdas(2)
    lambda_jp1 = lambdas(3)

    d_lambda_jm1 = d_lambdas(:, 1)
    d_lambda_j   = d_lambdas(:, 2)
    d_lambda_jp1 = d_lambdas(:, 3)

    call sphbes(j - 1, r*k, sj, sy, sjp, syp)
    j_hat_m1 = sj*r*k
    y_hat_m1 = sy*r*k
    call sphbes(j + 1, r*k, sj, sy, sjp, syp)
    j_hat_p1 = sj*r*k
    y_hat_p1 = sy*r*k

    lin_comb_ab = a*j_hat_m1 + b*y_hat_m1
    lin_comb_cd = c*j_hat_p1 + d*y_hat_p1
    diff_b = (lambda_jm1*lin_comb_ab + lambda_j*lin_comb_cd)/k
    diff_d = (lambda_jp1*lin_comb_cd + lambda_j*lin_comb_ab)/k

    d_lin_comb_ab = d_a*j_hat_m1 + d_b*y_hat_m1
    d_lin_comb_cd = d_c*j_hat_p1 + d_d*y_hat_p1
    d_diff_b = (d_lambda_jm1*lin_comb_ab + lambda_jm1*d_lin_comb_ab + &
        d_lambda_j*lin_comb_cd + lambda_j*d_lin_comb_cd)/k
    d_diff_d = (d_lambda_jp1*lin_comb_cd + lambda_jp1*d_lin_comb_cd + &
        d_lambda_j*lin_comb_ab + lambda_j*d_lin_comb_ab)/k

    a = a - diff_b*y_hat_m1
    b = b + diff_b*j_hat_m1
    c = c - diff_d*y_hat_p1
    d = d + diff_d*j_hat_p1

    d_a = d_a - d_diff_b*y_hat_m1
    d_b = d_b + d_diff_b*j_hat_m1
    d_c = d_c - d_diff_d*y_hat_p1
    d_d = d_d + d_diff_d*j_hat_p1

end subroutine coupled_variable_phase

!!
!> @brief      variable phase equation for uncoupled waves
!!
!! Performs a single step in the integration of the variable phase equation for uncoupled channels
!! with a linear angular momentum l
!!
!! See equation B14 in Phys. Rev. C 88 (2013) 064002
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine uncoupled_variable_phase(l, k, r, lambda, d_lambda, tan_delta, d_tan_delta)
    implicit none
    integer, intent(in) :: l !< orbital angular momentum quantum number
    real(dp), intent(in) ::  k !< center of mass momentum (in units of fm\f$^{-1}\f$)
    real(dp), intent(in) :: r !< integration radius in fm
    real(dp), intent(in) :: lambda !< lambda strength coefficient in fm\f$^{-2}\f$
    real(dp), intent(in) :: d_lambda(:) !< derivatives of lambda with respect of model parameters
    real(dp), intent(inout) :: tan_delta !< tangent of the phase shift
    real(dp), intent(inout) :: d_tan_delta(:) !< derivatives of tan_delta with respect of model parameters

    real(dp) :: j_hat, y_hat, phi, numerator, denominator, sj, sy, sjp, syp
    real(dp), allocatable :: d_phi(:), d_numerator(:), d_denominator(:)

    allocate(d_phi, d_numerator, d_denominator, mold = d_lambda)

    call sphbes(l, r*k, sj, sy, sjp, syp)
    j_hat = sj*r*k
    y_hat = sy*r*k
    phi = j_hat - tan_delta*y_hat
    d_phi = -d_tan_delta*y_hat
    numerator = tan_delta - lambda*j_hat*phi/k
    d_numerator = d_tan_delta - (d_lambda*phi + lambda*d_phi)*j_hat/k
    denominator = 1 - lambda*y_hat*phi/k
    d_denominator = -(d_lambda*phi + lambda*d_phi)*y_hat/k
    tan_delta = numerator/denominator
    d_tan_delta = (d_numerator*denominator - numerator*d_denominator)/denominator**2
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
    character(len=2), intent(in) :: reaction !< Type of reaction. (pp, np or nn)
    select case (reaction)
    case ('pp')
        k = sqrt(m_p/2*t_lab)/hbar_c
    case ('np')
        k = sqrt((m_p**2*t_lab*(t_lab + 2*m_n))/((m_p + m_n)**2 + 2*t_lab*m_p))/hbar_c
    case ('nn')
        k = sqrt(m_n/2*t_lab)/hbar_c
    case default
        stop 'incorrect reaction channel in av18_all_partial_waves'
    end select
end function momentum_cm

!!
!> @brief      find alfa parameters in coupled channels
!!
!! After integrating the variable phase equation for coupled channels, solves a quadratic equation
!! for the alpha parameters that determine the eigen phaseshifts.
!!
!! See equation B25 in Phys. Rev. C 88 (2013) 064002
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine solve_alfas(a1, b1, c1, d1, a2, b2, c2, d2, d_a1, d_b1, d_c1, d_d1, d_a2, d_b2, d_c2, &
    d_d2, alfa_1, alfa_2, d_alfa_1, d_alfa_2)
    implicit none
    real(dp), intent(in) :: a1 !< the \f$A_1\f$ parameter
    real(dp), intent(in) :: b1 !< the \f$B_1\f$ parameter
    real(dp), intent(in) :: c1 !< the \f$C_1\f$ parameter
    real(dp), intent(in) :: d1 !< the \f$D_1\f$ parameter
    real(dp), intent(in) :: a2 !< the \f$A_2\f$ parameter
    real(dp), intent(in) :: b2 !< the \f$B_2\f$ parameter
    real(dp), intent(in) :: c2 !< the \f$C_2\f$ parameter
    real(dp), intent(in) :: d2 !< the \f$D_2\f$ parameter
    real(dp), intent(in) :: d_a1(:) !< derivatives of the \f$A_1\f$ parameter
    real(dp), intent(in) :: d_b1(:) !< derivatives of the \f$B_1\f$ parameter
    real(dp), intent(in) :: d_c1(:) !< derivatives of the \f$C_1\f$ parameter
    real(dp), intent(in) :: d_d1(:) !< derivatives of the \f$D_1\f$ parameter
    real(dp), intent(in) :: d_a2(:) !< derivatives of the \f$A_2\f$ parameter
    real(dp), intent(in) :: d_b2(:) !< derivatives of the \f$B_2\f$ parameter
    real(dp), intent(in) :: d_c2(:) !< derivatives of the \f$C_2\f$ parameter
    real(dp), intent(in) :: d_d2(:) !< derivatives of the \f$D_2\f$ parameter
    real(dp), intent(out) :: alfa_1 !< the \f$\alpha_1\f$ solution
    real(dp), intent(out) :: alfa_2 !< the \f$\alpha_2\f$ solution
    real(dp), intent(out) :: d_alfa_1(:) !< derivatives of the \f$\alpha_1\f$ solution
    real(dp), intent(out) :: d_alfa_2(:) !< derivatives of the \f$\alpha_2\f$ solution
    real(dp) :: ar, br, cr, radical, numerator, denominator
    real(dp), allocatable, dimension(:) :: d_ar, d_br, d_cr, d_radical, d_numerator, d_denominator
    allocate(d_ar, mold = d_a1)
    ar = 0._dp
    allocate(d_br, d_cr, d_radical, d_numerator, d_denominator, source = d_a1)
    ar = a2*d2 - c2*b2
    br = a1*d2 + a2*d1 - c1*b2 - c2*b1
    cr = a1*d1 - c1*b1

    d_ar = d_a2*d2 + a2*d_d2 - d_c2*b2 - c2*d_b2
    d_br = d_a1*d2 + a1*d_d2 + d_a2*d1 + a2*d_d1 - d_c1*b2 - c1*d_b2 - d_c2*b1 - c2*d_b1
    d_cr = d_a1*d1 + a1*d_d1 - d_c1*b1 - c1*d_b1

    radical = br**2 - 4*ar*cr
    d_radical = 2*br*d_br - 4*(d_ar*cr + ar*d_cr)
    if (ar == 0._dp) then
        numerator = -cr
        denominator = br
        d_numerator = -d_cr
        d_denominator = d_br
        alfa_1 = numerator/denominator
        alfa_2 = alfa_1
        d_alfa_1 = (d_numerator*denominator - numerator*d_denominator)/denominator**2
        d_alfa_2 = d_alfa_1
    elseif (radical > 0) then
        numerator = -br + sqrt(radical)
        d_numerator = -d_br + d_radical/(2*sqrt(radical))
        denominator = 2*ar
        d_denominator = 2*d_ar
        alfa_1 = numerator/denominator
        d_alfa_1 = (d_numerator*denominator - numerator*d_denominator)/denominator**2

        numerator = -br - sqrt(radical)
        d_numerator = -d_br - d_radical/(2*sqrt(radical))
        alfa_2 = numerator/denominator
        d_alfa_2 = (d_numerator*denominator - numerator*d_denominator)/denominator**2
    else
        print*, 'WARNING: No real solutions in solve_alfas'
        alfa_1 = 0
        alfa_2 = 0
        d_alfa_1 = 0
        d_alfa_2 = 0
    endif
end subroutine solve_alfas

!!
!> @brief      Eigen phaseshifts for coupled channel
!!
!! Calculates the Eigen phaseshifts of a coupled partial waves channel after the variable phases
!! equation has been integrated.
!!
!! See equation B25 in Phys. Rev. C 88 (2013) 064002
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine eigen_phases(a1, b1, c1, d1, a2, b2, c2, d2, d_a1, d_b1, d_c1, d_d1, d_a2, d_b2, d_c2, &
    d_d2, alfa_1, alfa_2, d_alfa_1, d_alfa_2, ps_eigen, d_ps_eigen)
    implicit none
    real(dp), intent(in) :: a1 !< the \f$A_1\f$ parameter
    real(dp), intent(in) :: b1 !< the \f$B_1\f$ parameter
    real(dp), intent(in) :: c1 !< the \f$C_1\f$ parameter
    real(dp), intent(in) :: d1 !< the \f$D_1\f$ parameter
    real(dp), intent(in) :: a2 !< the \f$A_2\f$ parameter
    real(dp), intent(in) :: b2 !< the \f$B_2\f$ parameter
    real(dp), intent(in) :: c2 !< the \f$C_2\f$ parameter
    real(dp), intent(in) :: d2 !< the \f$D_2\f$ parameter
    real(dp), intent(in) :: d_a1(:) !< derivatives of the \f$A_1\f$ parameter
    real(dp), intent(in) :: d_b1(:) !< derivatives of the \f$B_1\f$ parameter
    real(dp), intent(in) :: d_c1(:) !< derivatives of the \f$C_1\f$ parameter
    real(dp), intent(in) :: d_d1(:) !< derivatives of the \f$D_1\f$ parameter
    real(dp), intent(in) :: d_a2(:) !< derivatives of the \f$A_2\f$ parameter
    real(dp), intent(in) :: d_b2(:) !< derivatives of the \f$B_2\f$ parameter
    real(dp), intent(in) :: d_c2(:) !< derivatives of the \f$C_2\f$ parameter
    real(dp), intent(in) :: d_d2(:) !< derivatives of the \f$D_2\f$ parameter
    real(dp), intent(in) :: alfa_1 !< the \f$\alpha_1\f$ solution
    real(dp), intent(in) :: alfa_2 !< the \f$\alpha_2\f$ solution
    real(dp), intent(out) :: d_alfa_1(:) !< derivatives of the \f$\alpha_1\f$ solution
    real(dp), intent(out) :: d_alfa_2(:) !< derivatives of the \f$\alpha_2\f$ solution
    real(dp), intent(out) :: ps_eigen(1:3) !< eigen phaseshifts in a coupled channel
    real(dp), intent(out) :: d_ps_eigen(:, :) !< derivatives of the eigen phaseshifts in a coupled channel
    real(dp) :: numerator, denominator, argument
    real(dp), allocatable, dimension(:) :: d_numerator, d_denominator, d_argument
    allocate(d_numerator, mold = d_a1)
    d_numerator = 0
    allocate(d_denominator, d_argument, source = d_numerator)
    numerator = b1 + alfa_1*b2
    denominator = a1 + alfa_1*a2
    argument = numerator/denominator
    d_numerator = d_b1 + d_alfa_1*b2 + alfa_1*d_b2
    d_denominator = d_a1 + d_alfa_1*a2 + alfa_1*d_a2
    d_argument = (d_numerator*denominator - numerator*d_denominator)/denominator**2
    ps_eigen(1) = -atan(argument)
    d_ps_eigen(:, 1) = -1/(1 + argument**2)*d_argument

    numerator = d1 + alfa_1*d2
    denominator = b1 + alfa_1*b2
    argument = numerator/denominator
    d_numerator = d_d1 + d_alfa_1*d2 + alfa_1*d_d2
    d_denominator = d_b1 + d_alfa_1*b2 + alfa_1*d_b2
    d_argument = (d_numerator*denominator - numerator*d_denominator)/denominator**2
    ps_eigen(2) =  atan(argument)
    d_ps_eigen(:, 2) = 1/(1 + argument**2)*d_argument

    numerator = d1 + alfa_2*d2
    denominator = c1 + alfa_2*c2
    argument = numerator/denominator
    d_numerator = d_d1 + d_alfa_2*d2 + alfa_2*d_d2
    d_denominator = d_c1 + d_alfa_2*c2 + alfa_2*d_c2
    d_argument = (d_numerator*denominator - numerator*d_denominator)/denominator**2
    ps_eigen(3) = -atan(argument)
    d_ps_eigen(:, 3) = -1/(1 + argument**2)*d_argument
end subroutine eigen_phases

!!
!> @brief      eigen to nuclear bar conversion
!!
!! Converts the set of phaseshifts in a coupled channel from the eigen to the nuclear bar
!! parametrization
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine eigen_2_bar(ps_eigen, d_ps_eigen, ps_bar, d_ps_bar)
    implicit none
    real(dp), intent(in) :: ps_eigen(1:3) !< eigen phaseshifts in a coupled channel
    real(dp), intent(in) :: d_ps_eigen(:, :) !< derivatives of the eigen phaseshifts in a coupled channel
    real(dp), intent(out) :: ps_bar(1:3) !< nuclear bar phaseshifts in a coupled channel
    real(dp), intent(out) :: d_ps_bar(:, :) !< derivatives of the nuclear bar phaseshifts in a coupled channel

    real(dp) :: sin_2eps, sin_diff, arg_arcsin, numerator, denominator, fraction
    real(dp), allocatable, dimension(:) :: d_sin_2eps, d_sin_diff, d_arg_arcsin, d_numerator, &
        d_denominator, d_fraction
    integer :: sign

    allocate(d_sin_2eps, mold = d_ps_eigen(:, 1))
    d_sin_2eps = 0
    allocate(d_sin_diff, d_arg_arcsin, d_numerator, d_denominator, d_fraction, source = d_sin_2eps)

    sin_2eps = sin(2*ps_eigen(2))
    d_sin_2eps = cos(2*ps_eigen(2))*2*d_ps_eigen(:, 2)
    sin_diff = sin(ps_eigen(1) - ps_eigen(3))
    d_sin_diff = cos(ps_eigen(1) - ps_eigen(3))*(d_ps_eigen(:, 1) - d_ps_eigen(:, 3))
    arg_arcsin = sin_2eps*sin_diff
    d_arg_arcsin = d_sin_2eps*sin_diff + sin_2eps*d_sin_diff
    ps_bar(2) = asin(arg_arcsin)/2
    d_ps_bar(:, 2) = 1/(2*sqrt(1 - arg_arcsin**2))*d_arg_arcsin

    numerator = tan(2*ps_bar(2))
    d_numerator = (1/cos(2*ps_bar(2))**2)*2*d_ps_bar(:, 2)
    denominator = tan(2*ps_eigen(2))
    d_denominator = (1/cos(2*ps_eigen(2))**2)*2*d_ps_eigen(:, 2)
    fraction = numerator/denominator
    d_fraction = (d_numerator*denominator - numerator*d_denominator)/denominator**2

    if (ps_eigen(1) - ps_eigen(3) > pi/2) then
        ps_bar(1) = (ps_eigen(1) + ps_eigen(3) + pi - asin(fraction))/2
        ps_bar(3) = (ps_eigen(1) + ps_eigen(3) - pi + asin(fraction))/2
        sign = -1
    else
        ps_bar(1) = (ps_eigen(1) + ps_eigen(3) + asin(fraction))/2
        ps_bar(3) = (ps_eigen(1) + ps_eigen(3) - asin(fraction))/2
        sign = +1
    endif
    d_ps_bar(: ,1) = (d_ps_eigen(:, 1) + d_ps_eigen(:, 3) + sign/sqrt(1 - fraction**2)*d_fraction)/2
    d_ps_bar(: ,3) = (d_ps_eigen(:, 1) + d_ps_eigen(:, 3) - sign/sqrt(1 - fraction**2)*d_fraction)/2

end subroutine eigen_2_bar

!!
!> @brief      Matches the uncoupled asymptotic solution the Coulomb wave function
!!
!! When integrating the uncoupled variable phase equation in the pp channel, the wave function
!! in the last integration step is matched form the free particle wave functions (reduced spherical
!! Bessel functions) to the Coulomb wave functions.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine match_uncoupled_waves(s, k, r, lambdas, d_lambdas, tan_deltas, d_tan_deltas)
    implicit none
    integer, intent(in) :: s !< spin quantum number
    real(dp), intent(in) :: k !< center of mass momentum (in units of fm\f$^{-1}\f$)
    real(dp), intent(in) :: r !< integration radius in fm
    real(dp), intent(in) :: lambdas(:) !< lambda strength coefficients for all uncoupled waves in fm\f$^{-2}\f$
    real(dp), intent(in) :: d_lambdas(:, :)  !< derivatives of lambda strength coefficients for all uncoupled waves
    real(dp), intent(inout) :: tan_deltas(:) !< tangent of the phase shift for all uncoupled waves
    real(dp), intent(inout) :: d_tan_deltas(:, :) !< derivatives of the tangent of the phase shift for all uncoupled waves
    real(dp) :: etap, lambda, eta0, numerator, denominator, diff
    real(dp), allocatable :: d_lambda(:), d_numerator(:), d_denominator(:), d_diff(:)
    integer :: l_max, l, ifail, i, lqm
    real(dp), dimension(:), allocatable :: FC, GC, FCP, GCP, jc, yc, jcp, ycp
    real(dp) :: F, G, jh, yh, jhp, yhp, Fp, Gp
    ifail = 0
    l_max = size(lambdas)

    allocate(d_lambda, mold = d_lambdas(:, 1))
    d_lambda = 0
    allocate(d_numerator, d_denominator, d_diff, source = d_lambda)

    if (l_max /= size(tan_deltas)) then
        stop 'lambdas and tan_deltas must have the same size in match_uncoupled_waves'
    endif

    allocate(FC(0:l_max))
    FC = 0
    allocate(GC, FCP, GCP, jc, yc, jcp, ycp, source = FC)
    eta0 = 0._dp
    etap = eta_prime(k)
    call COUL90(r*k, eta0, 0._dp, l_max, jc, yc, jcp, ycp, 1, ifail)
    if (ifail /= 0) stop 'coul90 fail'
    call COUL90(r*k, etap, 0._dp, l_max, fc, gc, fcp, gcp, 0, ifail)
    if (ifail /= 0) stop 'coul90 fail'

    do l = 0, l_max - 1
        if (mod(s+l,2) == 0 .or. (s == 1 .and. l == 0)) then
            if (s == 1 .and. l==0) then
                lqm = 1
            else
                lqm = l
            endif
            i = l + 1
            lambda = lambdas(i)
            d_lambda = d_lambdas(:, i)
            jh = jc(lqm)*r*k
            yh = yc(lqm)*r*k
            jhp = jcp(lqm)*r*k + jc(lqm)
            yhp = ycp(lqm)*r*k + yc(lqm)
            F = fc(lqm)
            G = gc(lqm)
            Fp = fcp(lqm)
            Gp = gcp(lqm)
            diff = jh - tan_deltas(i)*yh
            d_diff = -d_tan_deltas(:, i)*yh
            numerator   =  lambda*diff*F + k*(F*jhp - jh*Fp + tan_deltas(i)*(yh*Fp - F*yhp))
            denominator = -lambda*diff*G + k*(jh*Gp - g*jhp + tan_deltas(i)*(G*yhp - yh*Gp))
            d_numerator   =  (d_lambda*diff + lambda*d_diff)*F + k*d_tan_deltas(:, i)*(yh*Fp - F*yhp)
            d_denominator = -(d_lambda*diff + lambda*d_diff)*G + k*d_tan_deltas(:, i)*(G*yhp - yh*Gp)
            tan_deltas(i) = numerator/denominator
            d_tan_deltas(:, i) = (d_numerator*denominator - numerator*d_denominator)/denominator**2
        endif
    enddo

end subroutine match_uncoupled_waves

!!
!! @brief      energy dependent Sommerfeld parameter
!!
!! Calculates the energy dependent Sommerfeld parameter
!! \f$\eta' = \alpha M_p / (2 k) \left(\frac{1 + 2 k^2/M_p^2}{\sqrt{1 + k^2/M_p^2}} \right) \f$,
!! where \f$ \alpha \f$ is the fine structure constant, \f$ M_p \f$ is the mass of the proton
!! and \f$ k \f$ is the center of mass momentum
!!
!! @return     energy dependent Sommerfeld parameter
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function eta_prime(k) result(etap)
    implicit none
    real(dp), intent(in) :: k !< Center of mass momentum in fm\f$^{-1}\f$
    etap = m_p*alpha/(2*k*hbar_c)*(1 + 2*(k*hbar_c)**2/m_p**2)/sqrt(1 + (k*hbar_c)**2/m_p**2)
end function eta_prime

!!
!> @brief      Matches the coupled asymptotic solution to the Coulomb wave function
!!
!! When integrating the coupled variable phase equation in the pp channel, the wave function
!! in the last integration step is matched form the free particle wave functions (reduced spherical
!! Bessel functions) to the Coulomb wave functions.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine match_coupled_waves(k, r, lambdas, d_lambdas, a, b, c, d, d_a, d_b, d_c, d_d)
    implicit none
    real(dp), intent(in) :: k !< center of mass momentum (in units of fm\f$^{-1}\f$)
    real(dp), intent(in) :: r !< integration radius in fm
    real(dp), intent(in) :: lambdas(:, :) !< lambda strength coefficients for all coupled waves in fm\f$^{-2}\f$
    real(dp), intent(in) :: d_lambdas(:, :, :) !< derivatives lambda strength coefficients for all coupled waves
    real(dp), intent(inout) :: a(:) !< the \f$A\f$ parameter for all coupled waves
    real(dp), intent(inout) :: b(:) !< the \f$B\f$ parameter for all coupled waves
    real(dp), intent(inout) :: c(:) !< the \f$C\f$ parameter for all coupled waves
    real(dp), intent(inout) :: d(:) !< the \f$D\f$ parameter for all coupled waves
    real(dp), intent(inout) :: d_a(:, :) !< derivatives of the \f$A\f$ parameter for all coupled waves
    real(dp), intent(inout) :: d_b(:, :) !< derivatives of the \f$B\f$ parameter for all coupled waves
    real(dp), intent(inout) :: d_c(:, :) !< derivatives of the \f$C\f$ parameter for all coupled waves
    real(dp), intent(inout) :: d_d(:, :) !< derivatives of the \f$D\f$ parameter for all coupled waves
    integer :: j_max, ifail, j, ij
    real(dp) :: eta0, etap, ljm1, lj, ljp1, jhjm1, yhjm1, jhpjm1, yhpjm1, jhjp1, yhjp1, jhpjp1, &
        yhpjp1, Fjm1, Gjm1, Fpjm1, Gpjm1, Fjp1, Gjp1, Fpjp1, Gpjp1, Bf, Df
    real(dp), dimension(:), allocatable :: FC, GC, FCP, GCP, jc, yc, jcp, ycp
    real(dp), dimension(:), allocatable :: d_ljm1, d_lj, d_ljp1, d_bf, d_df

    j_max = size(lambdas, 2)

    if (size(a) /= size(b) .or. size(a) /= size(c) .or. size(a) /= size(d)) then
        stop 'incorrect array size in match_coupled_waves'
    endif

    if (j_max /= size(a) + 1 ) stop 'incorrect array size in match_coupled_waves'

    if (size(lambdas,1) /= 3) stop 'incorrect size for lambdas in match_uncoupled_waves'

    allocate(FC(0:j_max))
    FC = 0
    allocate(GC, FCP, GCP, jc, yc, jcp, ycp, source = FC)
    eta0 = 0._dp
    etap = eta_prime(k)
    call COUL90(r*k, eta0, 0._dp, j_max, jc, yc, jcp, ycp, 1, ifail)
    if (ifail /= 0) stop 'coul90 fail'
    call COUL90(r*k, etap, 0._dp, j_max, fc, gc, fcp, gcp, 0, ifail)
    if (ifail /= 0) stop 'coul90 fail'

    allocate(d_ljm1, mold = d_a(:, 1))
    d_ljm1 = 0
    allocate(d_lj, d_ljp1, d_bf, d_df, source = d_ljm1)

    do ij = 2, j_max
        j = ij - 1
        ljm1 = lambdas(1, ij)
        lj   = lambdas(2, ij)
        ljp1 = lambdas(3, ij)

        d_ljm1 = d_lambdas(:, 1, ij)
        d_lj   = d_lambdas(:, 2, ij)
        d_ljp1 = d_lambdas(:, 3, ij)


        jhjm1 = jc(j-1)*r*k
        yhjm1 = yc(j-1)*r*k
        jhpjm1 = jcp(j-1)*r*k + jc(j-1)
        yhpjm1 = ycp(j-1)*r*k + yc(j-1)
        jhjp1 = jc(j+1)*r*k
        yhjp1 = yc(j+1)*r*k
        jhpjp1 = jcp(j+1)*r*k + jc(j+1)
        yhpjp1 = ycp(j+1)*r*k + yc(j+1)

        Fjm1 = fc(j-1)
        Gjm1 = gc(j-1)
        Fpjm1 = fcp(j-1)
        Gpjm1 = gcp(j-1)
        Fjp1 = fc(j+1)
        Gjp1 = gc(j+1)
        Fpjp1 = fcp(j+1)
        Gpjp1 = gcp(j+1)

        Bf = Fjm1*(ljm1*(A(j)*jhjm1 + B(j)*yhjm1) + lj*(C(j)*jhjp1 + D(j)*yhjp1))/k &
            -B(j)*(yhjm1*Fpjm1 - yhpjm1*Fjm1) - A(j)*(jhjm1*Fpjm1 - jhpjm1*Fjm1)
        Df = Fjp1*(ljp1*(C(j)*jhjp1 + D(j)*yhjp1) + lj*(A(j)*jhjm1 + B(j)*yhjm1))/k &
            -D(j)*(yhjp1*Fpjp1 - yhpjp1*Fjp1) - C(j)*(jhjp1*Fpjp1 - jhpjp1*Fjp1)

        d_Bf = Fjm1*(d_ljm1*(A(j)*jhjm1 + B(j)*yhjm1) + ljm1*(d_A(:, j)*jhjm1 + d_B(:, j)*yhjm1) &
                     + d_lj*(C(j)*jhjp1 + D(j)*yhjp1) +   lj*(d_C(:, j)*jhjp1 + d_D(:, j)*yhjp1))/k &
              -d_B(:, j)*(yhjm1*Fpjm1 - yhpjm1*Fjm1) - d_A(:, j)*(jhjm1*Fpjm1 - jhpjm1*Fjm1)
        d_Df = Fjp1*(d_ljp1*(C(j)*jhjp1 + D(j)*yhjp1) + ljp1*(d_C(:, j)*jhjp1 + d_D(:, j)*yhjp1) &
                     + d_lj*(A(j)*jhjm1 + B(j)*yhjm1) +   lj*(d_A(:, j)*jhjm1 + d_B(:, j)*yhjm1))/k &
              -d_D(:, j)*(yhjp1*Fpjp1 - yhpjp1*Fjp1) - d_C(:, j)*(jhjp1*Fpjp1 - jhpjp1*Fjp1)

        A(j) = (A(j)*jhjm1 + B(j)*yhjm1 + Bf*Gjm1)/Fjm1
        C(j) = (C(j)*jhjp1 + D(j)*yhjp1 + Df*Gjp1)/Fjp1
        B(j) = Bf
        D(j) = Df

        d_A(:, j) = (d_A(:, j)*jhjm1 + d_B(:, j)*yhjm1 + d_Bf*Gjm1)/Fjm1
        d_C(:, j) = (d_C(:, j)*jhjp1 + d_D(:, j)*yhjp1 + d_Df*Gjp1)/Fjp1
        d_B(:, j) = d_Bf
        d_D(:, j) = d_Df
    enddo

end subroutine match_coupled_waves

!!
!> @brief      Adds energy dependent Coulomb term to pp potential
!!
!! Adds an energy dependent Coulomb term to the pp potential in all partial waves
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine add_coulomb(r, k, v_pw)
    implicit none
    real(dp), intent(in) :: r !< potential radius in fm
    real(dp), intent(in) :: k !< center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(out) :: v_pw(:, :) !< pp potential for all partial waves in MeV
    integer :: i
    real(dp) :: v_coul, fcoulr, br, kmev, alphap
    real(dp), parameter :: b=4.27_dp, small=0.e-5_dp
    kmev = k*hbar_c
    alphap = alpha*(1+2*kmev**2/m_p**2)/sqrt(1+kmev**2/m_p**2)
    br = b*r
    if (r < small) then
       fcoulr = 5*b/16
    else
        fcoulr = (1 - (1 +   11*br/16   + 3*br**2/16 + br**3/48)*exp(-br))/r
    end if

    v_coul = alphap*hbar_c*fcoulr

    do i = 1, size(v_pw, 2)
        v_pw(1:3, i) = v_pw(1:3, i) + v_coul
        v_pw(5, i) = v_pw(5, i) + v_coul
    enddo
end subroutine add_coulomb

! The functions below were written to test the derivatives of all_phaseshifts against numerical calculations
! and require the av18 module. All analytic calculations of the derivatives match the numerical ones.
! since the nn_phaseshifts module should not depend on a specific module for the NN interaction the functions
! are commented and left here for reference.

! !!
! !> @brief      wrapper function for all_phaseshifts
! !!
! !! This wrapper function is used to test the derivatives of the all_phaseshifts subroutine.
! !! The generic data of type context is used to receive all the arguments necessary to call
! !! all_phaseshifts. The same data of type context is used to receive which parameter will
! !! be varied by the dfridr subroutine and which partial wave will be returned.
! !!
! !! @returns    NN phase-shift at an specific lab energy and partial wave
! !!
! !! @author     Rodrigo Navarro Perez
! !!
! real(dp) function f_all_phaseshifts(x, data) result(r)
!     use num_recipes, only : context
!     use av18, only : av18_all_partial_waves
!     implicit none
!     real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
!     type(context), intent(in) :: data !< data structure with all the arguments for av18_operator

!     real(dp), allocatable :: ap(:)
!     real(dp) :: t_lab, r_max, dr
!     real(dp), allocatable :: phases(:, :), d_phases(:, :, :)
!     integer :: i_target, i_parameter, ic, ij
!     character(len=2) :: reaction

!     allocate(ap, source = data%x)
!     t_lab = data%a
!     r_max = data%b
!     dr = data%c
!     reaction = trim(data%string)
!     i_parameter = data%i
!     i_target = data%j

!     ap(i_parameter) = x

!     ic = mod(i_target - 1, 5) + 1
!     ij = 1 + (i_target - 1)/5

!     allocate(phases(1:5, ij))

!     call all_phaseshifts(av18_all_partial_waves, ap, t_lab, reaction, r_max, dr, phases, d_phases)

!     r = phases(ic, ij)
! end function f_all_phaseshifts

! !!
! !> @brief      wrapper function for the derivatives of all_phaseshifts
! !!
! !! This wrapper function is used to test the derivatives of the all_phaseshifts subroutine.
! !! The generic data of type context is used to receive all the arguments necessary to call
! !! all_phaseshifts. The same data of type context is used to receive which parameter will
! !! be varied by the dfridr subroutine and which partial wave will be returned.
! !!
! !! @returns    derivatives of a NN phase-shift at an specific lab energy and partial wave
! !!
! !! @author     Rodrigo Navarro Perez
! !!
! function df_all_phaseshifts(data) result(r)
!     use num_recipes, only : context
!     use av18, only : av18_all_partial_waves
!     implicit none
!     type(context), intent(in) :: data !< data structure with all the arguments for av18_operator
!     real(dp), allocatable :: r(:)

!     real(dp), allocatable :: ap(:)
!     real(dp) :: t_lab, r_max, dr
!     real(dp), allocatable :: phases(:, :), d_phases(:, :, :)
!     integer :: i_target, i_parameter, ic, ij
!     character(len=2) :: reaction

!     allocate(ap, source = data%x)
!     t_lab = data%a
!     r_max = data%b
!     dr = data%c
!     reaction = trim(data%string)
!     i_parameter = data%i
!     i_target = data%j

!     ic = mod(i_target - 1, 5) + 1
!     ij = 1 + (i_target - 1)/5

!     allocate(phases(1:5, ij))

!     call all_phaseshifts(av18_all_partial_waves, ap, t_lab, reaction, r_max, dr, phases, d_phases)

!     r = d_phases(:, ic, ij)
! end function df_all_phaseshifts

end module nn_phaseshifts