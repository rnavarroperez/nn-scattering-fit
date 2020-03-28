module observables
use av18, only: av18_all_partial_waves, default_params, n_parameters
use nn_phaseshifts, only: all_phaseshifts, momentum_cm
use amplitudes, only: saclay_amplitudes
use precisions, only: dp
use constants
implicit none
public observable
private

integer, parameter :: n_obs = 26
integer, parameter :: max_j = 20
character(len=4), dimension(1:n_obs), parameter :: &
       obs_types = ['DSG ','DT  ','AYY ','D   ','P   ','AZZ ','R   '&
         ,'RT  ','RPT ','AT  ','D0SK','NSKN','NSSN','NNKK','A   '&
         ,'AXX ','CKP ','RP  ','MSSN','MSKN','AZX ','AP  ','DTRT'&
         ,'SGT ','SGTT','SGTL']

contains


!> @brief Calculates a NN scattering observable for a given laboratory energy
!! in MeV and scattering angle in degrees.
!!
!! The observable to be calculated is determined by the type integer
!! and corresponds to the list of observable labels in the strObs
!! array.
!!
!! To avoid recalculating phase-shifts (the more time consuming part
!! of the calculation) the phases are only calculated if t and tprev
!! are different. At the end of the subroutine tprev is updated to
!! the value of t
!!
!! The derivative of the observable with respect of the parameters is stored
!! on the ap array
!!
!! @author Raul L Bernal-Gonzalez
subroutine observable(t_lab, pre_t_lab, angle, type, reac, obs, d_obs)
    implicit none
    real(dp), intent(in) :: t_lab ! laboratory energy
    real(dp), intent(inout) :: pre_t_lab ! previous value of t_lab
    real(dp), intent(in) :: angle ! scattering angle in degrees
    integer, intent(in) :: type ! index to indicate the type of observable
    character(len=2), intent(in) :: reac ! reaction channel
    real(dp), intent(out) :: obs ! NN scattering observable
    real(dp), intent(out), dimension(1:n_parameters) :: d_obs ! derivitive of observbles
    real(dp) :: k, eta, theta
    real(dp), parameter :: r_max = 2 , dr = .002 ! integration radius and step in fm
    save k, eta

    if(t_lab /= pre_t_lab) then
        k = momentum_cm(t_lab, reac)
        !call all_phaseshifts(av18_all_partial_waves, default_params, t_lab, reac, r_max, dr, k, phases, d_phases)
    end if
    theta = angle*pi/180.0_dp


end subroutine observable
end module observables
