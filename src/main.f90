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
use read_write, only : write_chiral_kernels, write_chiral_integrals, write_long_range_chiral_potentials, write_optimization_results
use long_range_chiral_potentials, only : vf_integral, vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9, chiral_integrals
! use short_range_chiral_potential, only : 
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

! call setup_optimization(model, parameters, mask, database, save_results, output_name)
! allocate(initial_parameters, source=parameters) !make a copy of the initial parameters to later save them
! call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
! if (save_results) then
!     call write_optimization_results(model, initial_parameters, parameters, mask, chi2, n_points, &
!         covariance, output_name)
! endif


!! THE FOLLOWING CODE IS FOR SHORT_RANGE_CHIRAL.f90 (-Ky Putnam):

    real(dp), dimension(1:28) :: short_lecs
    ! real(dp) :: r

    short_lecs(1) = 2.936041 !C_s
    short_lecs(2) = -0.4933897 !C_T

    short_lecs(3) = -0.1013462 !C_1
    short_lecs(4) = -0.1444844 !C_2
    short_lecs(5) = -0.03647634 !C_3
    short_lecs(6) = -0.01630825 !C_4
    short_lecs(7) = -0.006658100 !C_5
    short_lecs(8) = -0.06176835 !C_6
    short_lecs(9) = -0.9578191 !C_7

    short_lecs(10) = -0.03102824 !D_1
    short_lecs(11) = -0.004438695 !D_2
    short_lecs(12) = -0.01351171 !D_3
    short_lecs(13) = -0.0007084459 !D_4
    short_lecs(14) = 0.01110108 !D_5
    short_lecs(15) = -0.008598857 !D_6
    short_lecs(16) = -0.05367908 !D_7
    short_lecs(17) = 0.03119241 !D_8
    short_lecs(18) = 0.03281636 !D_9
    short_lecs(19) = -0.08647128 !D_10
    short_lecs(20) = -0.01167788 !D_11

    short_lecs(21) = 0.009575695 !C_0_IV
    short_lecs(22) = 0.02194758 !C_0_IT

    short_lecs(23) = -0.001550501 !C_1_IT
    short_lecs(24) = -0.008354679 !C_2_IT
    short_lecs(25) = -0.006682746 !C_3_IT
    short_lecs(26) = 0.01276971 !C_4_IT

    short_lecs(27) = 0.6 !R_s
    short_lecs(28) = 0.8 !R_L

    ! r = 0.01_dp

    ! call short_range_chiral_potential(r, short_lecs, v_short)

    

!! THE FOLLOWING CODE IS FOR CHIRAL_POTENTIAL.f90 (-Ky Putnam):
!! write_chiral_integrals, write_chiral_kernels,
!! and write_all_potentials_functions are located in read_write.f90

    ! r = 0.1_dp
    ! r_max = 12.1_dp

    ! for write_chiral_kernels
    ! character(len=50) :: file_name
    ! real(dp) :: r, r_max
    ! real(dp), dimension(9) :: mu_integrals

    call write_long_range_chiral_potentials()

    ! ! FOR write_chiral_integrals
    ! write(file_name, *) 'chiral_integrals'
    ! call write_chiral_potentials(r, file_name)

    ! do
    !     if (r > r_max) exit

    ! !     call chiral_integrals(r, mu_integrals)
    ! ! 
    ! !     ! FOR write_chiral_kernels
    !     write(file_name, '(a,i0.3,a)') 'chiral_integrands_r_', int(10*r), '.dat'
    !     call write_chiral_kernels(r, file_name)

    !     r = r + 2.0_dp

    ! enddo

end program nn_fit
