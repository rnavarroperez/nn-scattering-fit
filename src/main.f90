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
use read_write, only : write_chiral_kernals !, write_optimization_results
implicit none

! type(nn_model) :: model
! type(nn_experiment), allocatable, dimension(:) :: database
! real(dp), allocatable :: covariance(:,:)
! real(dp), allocatable, dimension(:) :: parameters
! real(dp), allocatable, dimension(:) :: initial_parameters
! logical, allocatable, dimension(:) :: mask
! logical :: save_results
! character(len=1024) :: output_name
! integer :: n_points

character(len=*) :: file_name
real(dp) :: r

! call setup_optimization(model, parameters, mask, database, save_results, output_name)
! allocate(initial_parameters, source=parameters) !make a copy of the initial parameters to later save them
! call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
! if (save_results) then
!     call write_optimization_results(model, initial_parameters, parameters, mask, chi2, n_points, &
!         covariance, output_name)
! endif

r = 0.1_dp

do
    if (r >= 12.1_dp) exit

    call write_chiral_kernals(r, file_name)
    r = r + 2.0_dp
    file_name = 'Chiral Integrands, r is '
end do

end program nn_fit
