module nn_scattering

use precisions, only : dp
use num_recipes, only : sphbes

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
    real(dp) :: r
    real(dp), allocatable :: v_pw(:, :), dv_pw(:, :, :)

    phases = 0

    n_params = size(params)
    n_waves = size(phases, 1)
    j_max = size(phases, 2)

    if (n_waves /= 5) stop 'incorrect number of waves for v_pw in all_phaseshifts'

    allocate(d_phases(1:n_params, 1:n_waves, 1:j_max))
    d_phases = 0._dp
    allocate(v_pw(1:n_waves, 1:j_max))

    r = dr/2 + t_lab*0
    do
        if( r > r_max) exit
        call model(params, r, reaction, v_pw, dv_pw)
        call uncoupled_variable_phase(0, 1.0_dp, r, v_pw(1, 1), phases(1, 1))
        do ij = 2, j_max
            j = ij - 1
        enddo
        r = r + dr
    enddo
    
end subroutine all_phaseshifts

subroutine uncoupled_variable_phase(l, k, r, lambda, tan_delta)
    implicit none
    integer, intent(in) :: l !< orbital angular momentum quantum number
    real(dp), intent(in) ::  k !< center of mass momentum (in units of fm\f$^{-1}\f$)
    real(dp), intent(in) :: r !< radius at which the matching is being done
    real(dp), intent(in) :: lambda !< lambda from the potential
    real(dp), intent(inout) :: tan_delta !< tangent of the wave function

    real(dp) :: j_hat, y_hat, phi, denominator, sj, sy, sjp, syp


    !call RedSphBes(l,r*k,F,G)
    call sphbes(l, r*k, sj, sy, sjp, syp)
    j_hat = sj*r*k
    y_hat = sy*r*k
    phi = j_hat - tan_delta*y_hat
    denominator = 1 - lambda*y_hat*phi/k
    tan_delta = (tan_delta - lambda*j_hat*phi/k)/denominator


end subroutine uncoupled_variable_phase
    
end module nn_scattering