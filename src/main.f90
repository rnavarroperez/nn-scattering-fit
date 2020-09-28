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
use av18, only : default_params, av18_all_partial_waves, n_parameters
use delta_shell, only : nn_model
use exp_data, only : nn_experiment, read_database, init_ex_em_amplitudes
use observables, only : kinematics, observable
use amplitudes, only : em_amplitudes
use chi_square, only: calc_chi_square
implicit none

real(dp), parameter :: r_max = 12.5_dp
real(dp), parameter :: dr = 0.01_dp
type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: experiments
type(kinematics) :: kine
integer :: i, j
real(dp) :: obs,residual, exp_val, sigma, chi2
integer :: n_points
real(dp), allocatable, dimension(:) :: d_obs

model%potential => av18_all_partial_waves
model%r_max = r_max
model%dr = dr
model%potential_type = 'local'

allocate(experiments(1:2))
call read_database('database/granada_database.dat', experiments)
call init_ex_em_amplitudes(experiments)
call calc_chi_square(experiments, default_params, model, chi2)
print*, chi2

end program nn_fit
