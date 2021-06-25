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
use exp_data, only : nn_experiment
use optimization, only : lavenberg_marquardt, setup_optimization
use randomize_exp, only : full_bootstrap
implicit none

type(nn_model) :: model
type(nn_experiment), allocatable, dimension(:) :: database
real(dp), allocatable :: covariance(:,:)
real(dp), allocatable, dimension(:) :: parameters
real(dp), allocatable, dimension(:) :: new_parameters
logical, allocatable, dimension(:) :: mask
real(dp) :: chi2
integer :: n_points
integer, parameter :: n_runs = 5
real(dp), allocatable, dimension(:) :: all_chi2
integer, allocatable, dimension(:) :: all_npoints
real(dp), allocatable :: all_parameters(:,:)
real(dp), allocatable ::alpha(:,:)
real(dp), allocatable :: beta(:)

call setup_optimization(model, parameters, mask, database)

call full_bootstrap(database, mask, model, parameters, n_runs,&
       all_chi2, all_npoints, all_parameters)
!call bootstrap(database, mask, model, parameters, new_parameters, chi2, n_points) 
!call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
!print*, 'after minimization: ', chi2, n_points, chi2/n_points

end program nn_fit
