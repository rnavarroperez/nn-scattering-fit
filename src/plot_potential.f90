program plot_potential
use precisions, only : dp
use delta_shell, only : nn_model
use exp_data, only : nn_experiment
use optimization, only : setup_optimization, covariance_matrix
use chi_square, only : total_chi_square
use read_write, only : plot_potential_components
implicit none
type(nn_model) :: model
real(dp), allocatable, dimension(:) :: parameters
logical, allocatable, dimension(:) :: mask
type(nn_experiment), allocatable, dimension(:) :: database
logical :: save_results
character(len=1024) :: output_file
integer :: n_points
real(dp) :: chi2
real(dp), allocatable, dimension(:, :) :: alpha, covariance
real(dp), allocatable, dimension(:) :: beta
real(dp), parameter :: r_min = 0.0_dp, r_max = 2.0_dp, r_step = 0.0078125_dp

call setup_optimization(model, parameters, mask, database, save_results, output_file)
call total_chi_square(database, parameters, mask, model, n_points, chi2, alpha, beta)
covariance = covariance_matrix(alpha, mask)
call plot_potential_components(model, parameters, covariance, r_min, r_max, r_step, output_file)

end program plot_potential