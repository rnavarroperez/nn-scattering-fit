module bootstrap
use exp_data
use precisions, only: dp
use optimization
use delta_shell, only: nn_model
use chi_square, only : total_chi_square
use randomize_exp, only: randomize_database
IMPLICIT NONE
private
public full_bootstrap
contains

!!
!> @brief       Readjusts parameters with a randomized database
!!
!! Subroutine to readjust parameters using the lavenberg_marquardt 
!! subroutine. Parameters are readjusted with a randomized database 
!! from the randomize_database function.  new_parameters and new_exp
!! are allocated based on the size of parameters and experiments. 
!!
!! @author      Marielle Elizabeth Duran
!!
subroutine one_bootstrap(experiments, mask, model, parameters, new_parameters, chi2, n_points)
        IMPLICIT NONE
        type(nn_experiment), intent(in), dimension(:) :: experiments
        type(nn_experiment), allocatable,  dimension(:) :: new_exp
        LOGICAL, intent(in), dimension(:) :: mask

        type(nn_model), intent(in) :: model
        REAL(dp), intent(in) :: parameters(:)
        REAL(dp), allocatable, intent(out) :: new_parameters(:)
        REAL(dp), intent(out) :: chi2
        INTEGER, intent(out) :: n_points
        REAL(dp), allocatable :: covariance(:,:)

        allocate(new_parameters(1:SIZE(parameters)))
        allocate(new_exp(1:SIZE(experiments)))

        new_exp = randomize_database(experiments)              
        new_parameters = parameters
        call lavenberg_marquardt(new_exp, mask, model, new_parameters, n_points, chi2, covariance)
END subroutine one_bootstrap

!!
!> @brief       Implements full bootstrap method 
!!
!! Subroutine to implement the entire bootstrap method using the 
!! lavenberg-marquardt method and the total_chi_square subroutine. This
!! follows the same steps as the botstrap subroutine listed above,
!! except this subroutine uses the total_chi_square subroutine to 
!! find the chi square values. It also puts values of new parameters,
!! chi square, and number of points into three separate arrays. These
!! three arrays are written into the file, all_arrays.
!!
!! @author      Marielle Duran
!! 
subroutine full_bootstrap(old_exp, mask, model, parameters, n_runs, all_chi2, , all_npoints, all_parameters)
        IMPLICIT NONE
        type(nn_experiment), intent(in), dimension(:) :: old_exp
        type(nn_experiment), allocatable, dimension(:) :: new_exp
        LOGICAL, intent(in), dimension(:) :: mask
        type(nn_model), intent(in) :: model
        REAL(dp), intent(in) :: parameters(:)
        REAL(dp), allocatable :: new_parameters(:)
        REAL(dp) :: chi2
        INTEGER :: n_points
        REAL(dp), allocatable :: covariance(:,:)

        REAL(dp), allocatable, intent(out), dimension(:) :: all_chi2
        INTEGER, allocatable, intent(out), dimension(:) :: all_npoints
        REAL(dp), allocatable, intent(out) :: all_parameters(:,:)
        REAL(dp), allocatable ::alpha(:,:)
        REAL(dp), allocatable :: beta(:)
        INTEGER :: i, unit
        INTEGER, intent(in) :: n_runs

        allocate(new_parameters(1:SIZE(parameters)))
        allocate(new_exp(1:SIZE(old_exp)))
        allocate(all_chi2(1:n_runs))
        allocate(all_npoints(1:n_runs))
        allocate(all_parameters(1:SIZE(parameters),1:n_runs))

        !open(newunit=unit, file=output_name//'_bootstrap.dat', status='unknown')
        !        write(unit, *), '# '
        !        write(unit, *), '# potential parameters, chi square, number of points '
        !close(unit)
        DO i = 1, n_runs
                new_exp = randomize_database(old_exp)
                new_parameters = parameters
                call lavenberg_marquardt(new_exp, mask, model, new_parameters, n_points, chi2, covariance)
                call total_chi_square(old_exp, new_parameters, mask, model, n_points, chi2, alpha, beta)

                all_parameters(:,i) = new_parameters
                all_chi2(i) = chi2
                all_npoints(i) = n_points
                !open(newunit, )
                !write(unit,*) all_parameters(:,i), all_chi2(i), all_npoints(i)
                !close(unit)
        END DO
END subroutine full_bootstrap

END module bootstrap