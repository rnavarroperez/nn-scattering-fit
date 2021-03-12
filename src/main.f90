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
use random_num, only: box_muller_num
use randomize_exp, only: randomize_experiment, calc_mean, calc_stan_dev
implicit none

real(dp), parameter :: r_max = 12.5_dp
real(dp), parameter :: dr = 0.01_dp
type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: experiments
type(nn_experiment) :: new_exp, old_exp
real(dp), allocatable :: covariance(:,:)
real(dp), allocatable, dimension(:) :: parameters
real(dp), dimension(14) :: check
logical, allocatable, dimension(:) :: mask
integer :: n_points
real(dp) :: chi2
real(dp) :: x, y, mean, stan_dev
integer :: i
integer, parameter :: n_samples = 10000

model%potential => av18_all_partial_waves
model%display_subroutine => display_parameters
model%r_max = r_max
model%dr = dr
model%potential_type = 'local'


call box_muller_num(x)
!PRINT*, x
!DO i = 1,1000
 !       call box_muller_num(x)
 !       PRINT*, x
!END DO
!call generator_100_num
!call verify_box_muller_num(n_samples)

!allocate(parameters, mold=default_params)
!call read_marias_format('av18.bob.in', parameters)

!allocate(experiments(1:2))
call read_database('database/granada_database.dat', experiments)
old_exp = experiments(1)
!new_exp = old_exp
new_exp = randomize_experiment(old_exp)
open(8, file = 'check_experiment.dat', status = 'unknown')
do i = 1, new_exp%n_data
        y = (old_exp%data_points(i)%value - new_exp%data_points(i)%value) / &
                (old_exp%data_points(i)%stat_error + &
                new_exp%data_points(i)%stat_error)
        !y = old_exp%data_points(i)%value
        check(i) = y
        write(8,*) y
end do
mean = calc_mean(check)
stan_dev = calc_stan_dev(check)
write(8,*) mean
write(8,*) stan_dev
close(8)
!call init_ex_em_amplitudes(experiments)

!allocate(mask(1: size(parameters)))
!mask = .true.
!call lavenberg_marquardt(experiments, mask, model, parameters, n_points, chi2, covariance)
!print*, 'after minimization: ', chi2, n_points, chi2/n_points
end program nn_fit
