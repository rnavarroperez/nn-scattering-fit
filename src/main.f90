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
use read_write, only : write_short_range_chiral_potentials, write_chiral_kernels, write_chiral_integrals, &
                       write_long_range_chiral_potentials, write_optimization_results
use long_range_chiral_potentials, only : vf_integral, vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9, chiral_integrals, &
                                         long_range_potentials
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

    short_lecs(1) = 2.936041_dp !C_s
    short_lecs(2) = -0.4933897_dp !C_T

    short_lecs(3) = -0.1013462_dp !C_1
    short_lecs(4) = -0.1444844_dp !C_2
    short_lecs(5) = -0.03647634_dp !C_3
    short_lecs(6) = -0.01630825_dp !C_4
    short_lecs(7) = -0.006658100_dp !C_5
    short_lecs(8) = -0.06176835_dp !C_6
    short_lecs(9) = -0.9578191_dp !C_7

    short_lecs(10) = -0.03102824_dp !D_1
    short_lecs(11) = -0.004438695_dp !D_2
    short_lecs(12) = -0.01351171_dp !D_3
    short_lecs(13) = -0.0007084459_dp !D_4
    short_lecs(14) = 0.01110108_dp !D_5
    short_lecs(15) = -0.008598857_dp !D_6
    short_lecs(16) = -0.05367908_dp !D_7
    short_lecs(17) = 0.03119241_dp !D_8
    short_lecs(18) = 0.03281636_dp !D_9
    short_lecs(19) = -0.08647128_dp !D_10
    short_lecs(20) = -0.01167788_dp !D_11

    short_lecs(21) = 0.009575695_dp !C_0_IV
    short_lecs(22) = 0.02194758_dp !C_0_IT

    short_lecs(23) = -0.001550501_dp !C_1_IT
    short_lecs(24) = -0.008354679_dp !C_2_IT
    short_lecs(25) = -0.006682746_dp !C_3_IT
    short_lecs(26) = 0.01276971_dp !C_4_IT

    short_lecs(27) = 0.6_dp !R_s
    short_lecs(28) = 0.8_dp !R_L

!! THE FOLLOWING CODE IS FOR CHIRAL_POTENTIAL.f90 (-Ky Putnam):
!! write_chiral_integrals, write_chiral_kernels,
!! and write_all_potentials_functions are located in read_write.f90

! call test_potentials()

    ! r = 0.1_dp
    ! r_max = 12.1_dp

    ! for write_chiral_kernels
    ! character(len=50) :: file_name
    ! real(dp) :: r, r_max
    ! real(dp), dimension(9) :: mu_integrals

    ! call write_long_range_chiral_potentials()
    call write_short_range_chiral_potentials(short_lecs)

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
