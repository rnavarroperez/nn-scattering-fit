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
use delta_shell, only : nn_local_model
use exp_data, only : nn_experiment, read_database
use observables, only : kinematics, observable
use amplitudes, only : em_amplitudes
implicit none

real(dp), parameter :: r_max = 12.5_dp
real(dp), parameter :: dr = 0.01_dp
type(nn_local_model) :: model
type(nn_experiment), allocatable, dimension(:) :: experiments
type(kinematics) :: kine
integer :: i, j
real(dp) :: obs,residual, chi2, exp_val, sigma
integer :: n_points
real(dp), allocatable, dimension(:) :: d_obs

model%potential => av18_all_partial_waves
model%r_max = r_max
model%dr = dr

allocate(experiments(1:2))
call read_database('database/granada_database.dat', experiments)
chi2 = 0
n_points = 0
do i = 1, size(experiments)
    if (experiments(i)%rejected) cycle
    kine%channel = experiments(i)%channel
    kine%type = experiments(i)%obs_type
    do j = 1, experiments(i)%n_data
        kine%t_lab = experiments(i)%data_points(j)%t_lab
        kine%angle = experiments(i)%data_points(j)%theta
        kine%em_amplitude = em_amplitudes(kine%t_lab, kine%angle, kine%channel)
        call observable(kine, default_params, model, obs, d_obs)
        exp_val = experiments(i)%data_points(j)%value
        sigma = experiments(i)%data_points(j)%stat_error
        residual = (exp_val - obs)/sigma
        chi2 = chi2 + residual**2
        n_points = n_points + 1
        print*, kine%channel, kine%type, kine%t_lab, kine%angle, exp_val, obs, residual, chi2/n_points
    enddo
enddo

end program nn_fit
