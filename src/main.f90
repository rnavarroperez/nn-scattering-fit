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
use read_write, only : write_chiral_kernels, write_chiral_integrals !, write_optimization_results
use chiral_potential, only : vf_integral, vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9 
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

character(len=50) :: file_name
real(dp) :: r!, r_max

! call setup_optimization(model, parameters, mask, database, save_results, output_name)
! allocate(initial_parameters, source=parameters) !make a copy of the initial parameters to later save them
! call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
! if (save_results) then
!     call write_optimization_results(model, initial_parameters, parameters, mask, chi2, n_points, &
!         covariance, output_name)
! endif

r = 0.1_dp
! r_max = 12.1_dp

! FOR write_chiral_integrals
write(file_name, *) 'chiral_integrals'
call write_chiral_integrals(r, file_name)

! do
!     if (r > r_max) exit

!     ! FOR write_chiral_kernels
!     write(file_name, '(a,i0.3,a)') 'chiral_integrands_r_', int(10*r), '.dat'
!     call write_chiral_kernels(r, file_name)

!     !FOR printing chiral integrals in terminal
!     print'(f11.1,9es20.9)', r, vf_integral(vf_1, r), vf_integral(vf_2, r), vf_integral(vf_3, r), &
!             vf_integral(vf_4, r), vf_integral(vf_5, r), vf_integral(vf_6, r), &
!             vf_integral(vf_7, r), vf_integral(vf_8, r), vf_integral(vf_9, r)

!     r = r + 2.0_dp

! end do

end program nn_fit
