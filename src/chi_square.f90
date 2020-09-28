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
public :: calc_chi_square
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
subroutine calc_chi_square(experiments, potential_parameters, model, chi2)
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: experiments !< input experiment data
    real(dp), intent(in) :: potential_parameters(:) !< potential model parameters
    type(nn_model), intent(in) :: model !< potential model
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model

    type(kinematics) :: kine
    real(dp) :: exp_val, sigma, obs, z_scale, chi_sys_error_cont
    real(dp) :: znum, zden, sys_error, nth_chi2
    real(dp), allocatable :: d_obs(:)
    logical :: float
    integer :: i, j

    ! initialize chi-squares at 0
    chi2 = 0
    nth_chi2 = 0
    !n_points = 0
    ! initialize znum and zden
    znum = 0
    zden = 0

    do i = 1, size(experiments) ! goes through all the experiments
        if (experiments(i)%rejected) cycle
        kine%channel = experiments(i)%channel
        kine%type = experiments(i)%obs_type
        sys_error = experiments(i)%sys_error
        do j = 1, experiments(i)%n_data ! goes through all points of nth experiment
            kine%t_lab = experiments(i)%data_points(j)%t_lab
            kine%angle = experiments(i)%data_points(j)%theta
            kine%em_amplitude = experiments(i)%data_points(j)%em_amplitude
            exp_val = experiments(i)%data_points(j)%value
            sigma = experiments(i)%data_points(j)%stat_error
            call observable(kine, potential_parameters, model, obs, d_obs)
            znum = znum + (exp_val*obs)/sigma**2
            zden = zden + (obs/sigma)**2
        end do
        call calc_z_scale(sys_error, znum, zden, z_scale, chi_sys_error_cont, float)
        ! calculate chi-square for nth experiment
        nth_chi2 = ((exp_val - obs*z_scale)/sigma)**2
        ! check if data is floated
        if(float) then ! add the systematic error contribution to the chi-square
             nth_chi2 = nth_chi2 + chi_sys_error_cont
        end if
        ! total chi-square is the sum of all nth experiments
        chi2 = chi2 + nth_chi2
        ! reset valued for znum and zden
        znum = 0
        zden = 0
        ! n_points = n_points + 1
    end do

end subroutine calc_chi_square

! calculate Z scale
subroutine calc_z_scale(sys_error, znum, zden, z_scale, chi_sys_error_cont, float)
    implicit none
    real(dp), intent(in) :: sys_error
    real(dp), intent(in) :: znum
    real(dp), intent(in) :: zden
    real(dp), intent(out) :: z_scale
    real(dp), intent(out) :: chi_sys_error_cont
    logical, intent(out) :: float
    real(dp) :: num, den

    if(sys_error == 0) then ! check if data is absolute
        z_scale = 1._dp
        float = .false.
    else if (sys_error > 0.25_dp) then ! check for large systematic error
        z_scale = znum/zden
        float = .false.
    else ! if systematic error is greater than 0 but less than 0.25
        num = znum + 1._dp/(sys_error)**2
        den = zden + 1._dp/(sys_error)**2
        z_scale = num/den
        float = .true.
        chi_sys_error_cont = ((z_scale-1._dp)/sys_error)**2
        if (chi_sys_error_cont > 9._dp) then
            float = .false.
            z_scale = znum/zden
        end if
    end if
end subroutine calc_z_scale

end module chi_square
