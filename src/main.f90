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
use read_write, only : write_optimization_results
implicit none

type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: database
real(dp), allocatable :: covariance(:,:)
real(dp), allocatable, dimension(:) :: parameters
real(dp), allocatable, dimension(:) :: initial_parameters
logical, allocatable, dimension(:) :: mask
real(dp), parameter, dimension(1:3) :: target_shape = [1.9_dp, 0.525_dp, 0.21_dp]
real(dp), dimension(1:3) :: delta_p
integer, parameter :: n_steps = 256
logical :: save_results
character(len=1024) :: output_name
real(dp) :: chi2
integer :: n_points, i, unit

call setup_optimization(model, parameters, mask, database, save_results, output_name)
delta_p = (target_shape - parameters(52:54))/n_steps
allocate(initial_parameters, source=parameters) !make a copy of the initial parameters to later save them
open(newunit=unit, file='adiabatic_fits_17_19.dat')
do i = 1, n_steps
    parameters(52:54) = parameters(52:54) + delta_p
    call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
    write(unit, *) parameters(52:54), chi2, n_points, chi2/n_points
enddo
close(unit)
output_name = 'av19_cutoff19'
if (save_results) then
    call write_optimization_results(model, initial_parameters, parameters, mask, chi2, n_points, &
        covariance, output_name)
endif
end program nn_fit
