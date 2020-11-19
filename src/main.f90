!!
!> @brief fits a NN potential to scattering data
!!
!! Uses the Levenberg Marquardt algorithm to adjust the parameters of a NN interaction and
!! reproduce experimental data collected since 1954.
!!
!! @author Rodrigo Navarro Perez
!!
program nn_fit
use precisions, only : dp
use av18, only : default_params, av18_all_partial_waves
use delta_shell, only : nn_model
use exp_data, only : nn_experiment, read_database, init_ex_em_amplitudes
!use chi_square, only: calc_chi_square
use chi_optimization, only: lavenberg_marquardt

implicit none

real(dp), parameter :: r_max = 12.5_dp
real(dp), parameter :: dr = 0.01_dp
type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: experiments
real(dp), allocatable :: covariance(:,:), new_parameters(:)
real(dp), allocatable :: alpha(:,:), beta(:)
integer :: n_points
real(dp) :: chi2

model%potential => av18_all_partial_waves
model%r_max = r_max
model%dr = dr
model%potential_type = 'local'

allocate(experiments(1:2))
call read_database('database/granada_database.dat', experiments)
call init_ex_em_amplitudes(experiments)

!call calc_chi_square(experiments, default_params, model, n_points, chi2, alpha, beta)
!print*, 'Before minimization: ', chi2, n_points, chi2/n_points
call lavenberg_marquardt(experiments, default_params, model, n_points, chi2, covariance, new_parameters)
print*, 'after minimization: ', chi2, n_points, chi2/n_points
end program nn_fit
