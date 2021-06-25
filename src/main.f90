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
use delta_shell, only : nn_model
use exp_data, only : nn_experiment!, read_database, init_ex_em_amplitudes
use optimization, only: lavenberg_marquardt, setup_optimization
use string_functions, only : mask_to_string
implicit none

type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: database
real(dp), allocatable :: covariance(:,:)
real(dp), allocatable, dimension(:) :: parameters
logical, allocatable, dimension(:) :: mask
real(dp) :: chi2
integer :: n_points

call setup_optimization(model, parameters, mask, database)

call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
print*, 'after minimization: ', chi2, n_points, chi2/n_points
end program nn_fit
