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
use av18, only : default_params, av18_all_partial_waves, display_parameters
use av18_compatibility, only : read_marias_format
use delta_shell, only : nn_model
use exp_data, only : nn_experiment, read_database, init_ex_em_amplitudes
use optimization, only: lavenberg_marquardt
use randomize_exp
implicit none

real(dp), parameter :: r_max = 12.5_dp
real(dp), parameter :: dr = 0.01_dp
type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: experiments
real(dp), allocatable :: covariance(:,:)
real(dp), allocatable, dimension(:) :: parameters
real(dp), allocatable, dimension(:) :: new_parameters
logical, allocatable, dimension(:) :: mask
real(dp) :: chi2
integer :: n_points

model%potential => av18_all_partial_waves
model%display_subroutine => display_parameters
model%r_max = r_max
model%dr = dr
model%potential_type = 'local'

allocate(parameters, source=default_params)

allocate(experiments(1:2))
call read_database('database/granada_database.dat', experiments)
call init_ex_em_amplitudes(experiments)

allocate(mask(1: size(parameters)))
mask = .true.
call bootstrap(experiments, mask, model, parameters, new_parameters, chi2, n_points) 
!call lavenberg_marquardt(experiments, mask, model, parameters, n_points, chi2, covariance)
print*, 'after minimization: ', chi2, n_points, chi2/n_points
end program nn_fit
