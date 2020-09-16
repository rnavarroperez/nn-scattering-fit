!!
!> @brief      Chi Square
!!
!! This module contains necessary functions
!! to evaluate and minimize the chi square
!!
!! @author Rodrigo Navarro Perez
!! @author Raul L Bernal-Gonzalez
!!
module chi_square
use precisions, only: dp
use exp_data, only : nn_experiment
use observables, only: observable, kinematics
use delta_shell, only : nn_model
use amplitudes, only : em_amplitudes
implicit none

private
public :: simple_chi_square
contains



    !!
    !> @brief      chi_square
    !!
    !! Calculates the chi square merit function for all experiments
    !! given the provided parameters and model
    !!
    !! @author Rodrigo Navarro Perez
    !! @author Raul L Bernal-Gonzalez
    !!
subroutine chi_square(experiments, potential_parameters, model, chi2)
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: experiments !< input experiment data
    real(dp), intent(in) :: potential_parameters(:) !< potential model parameters
    type(nn_model), intent(in) :: model !< potential model
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model

    type(kinematics) :: kine
    real(dp) :: obs
    real(dp), allocatable :: d_obs(:)
    integer :: n_points, i, j
    real(dp) :: exp_val, sigma, residual

    chi2 = 0
    n_points = 0
    do i = 1, size(experiments)
        if (experiments(i)%rejected) cycle
        kine%channel = experiments(i)%channel
        kine%type = experiments(i)%obs_type
        do j = 1, experiments(i)%n_data
            kine%t_lab = experiments(i)%data_points(j)%t_lab
            kine%angle = experiments(i)%data_points(j)%theta
            kine%em_amplitude = experiments(i)%data_points(j)%em_amplitude
            call observable(kine, potential_parameters, model, obs, d_obs)
            exp_val = experiments(i)%data_points(j)%value
            sigma = experiments(i)%data_points(j)%stat_error
            residual = (exp_val - obs)/sigma
            chi2 = chi2 + residual**2
            n_points = n_points + 1
        enddo
    enddo
end subroutine simple_chi_square

end module chi_square
