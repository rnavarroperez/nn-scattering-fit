!!
!> @brief   chi_minimization
!!
!! Subroutines and functions to implement the Lavenberg-Marquardt method for
!! minimizing the chi-square
!!
!! @author Raul L Bernal-Gonzalez
!!
module chi_optimization
use precisions, only: dp
use exp_data, only: nn_experiment
use delta_shell, only: nn_model
use chi_square, only: calc_chi_square
implicit none

private
public :: lavenberg_marquardt
contains

!!
!> @brief   lavenberg-marquardt
!!
!! This subroutine contains the main loop that implements the
!! lavenberg-marquardt method. It will return the minimize chi-square, covariance,
!! and the optimal parameters for the input model
!!
!! @author Raul L Bernal-Gonzalez
!!
subroutine lavenberg_marquardt(experiments, model_parameters, model, n_points, chi2, covariance, new_parameters)
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: experiments !< input experiment data
    real(dp), intent(in) :: model_parameters(:) !< potential model parameters
    type(nn_model), intent(in) :: model !< potential model
    integer, intent(out) :: n_points !< number of data points in the total chi-square
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model
    real(dp), intent(out), allocatable :: covariance(:,:) !< estimated covariance
    real(dp), intent(out), allocatable :: new_parameters(:) !< optimize parameters

    real(dp) :: lambda, delta, prev_chi2
    real(dp), allocatable :: alpha_prime(:,:), alpha(:,:), beta(:), old_parameters(:)
    integer :: counter, factor, i, j, n_param, limit

    ! allocate new_parameters
    n_param = size(model_parameters)
    allocate(new_parameters, old_parameters, mold=model_parameters)
    allocate(covariance(n_param, n_param))
    allocate(alpha_prime, mold=covariance)
    ! initialize lambda
    lambda = 0.001_dp
    ! set term for negligible differnce
    delta = 0.01_dp
    ! factor to increase lambda by
    factor = 10
    ! condition counter
    limit = 0 ! stop after 2 consecutive negligible differences
    counter = 0
    ! first call outside the loop to initiate values using original parameters
    call calc_chi_square(experiments, model_parameters, model, n_points, chi2, alpha, beta)
    prev_chi2 = chi2
    old_parameters = model_parameters
    ! begin optimization loop
    do
        ! check exit condition
        if(limit >= 2) exit
        ! limit iterations for testing
        if(counter >= 5) exit
        ! set alpha-prime
        alpha_prime = set_alpha_prime(alpha, lambda)
        ! calculate new parameters
        call calc_new_parameters(alpha_prime, beta, old_parameters, new_parameters)
        call exit(0)
        ! save new parameters
        old_parameters = new_parameters
        ! calculate chi-square with new parameters
        call calc_chi_square(experiments, new_parameters, model, n_points, chi2, alpha, beta)
        ! compare chi-squares
        if((prev_chi2 - chi2) <= delta) then
            ! increase limit by 1 if there is a negligible difference
            limit = limit + 1
        else ! we want 2 consecutive negligible differences
            limit = 0
        end if
        ! determined whether to raise or lower lambda
        if(chi2 >= prev_chi2) then
            lambda = lambda*factor
        else if(chi2 < prev_chi2) then
            lambda = lambda/factor
        end if
        ! update prev chi
        prev_chi2 = chi2
        counter = counter + 1
    end do
    ! set covariance to last alpha calculated
    covariance = alpha
end subroutine lavenberg_marquardt

!!
!> @brief   calc_new_parameters
!!
!! The function uses lapack library to inverse alpha and solve the linear equation
!! alpha*delta = beta. It then adds delta to the parameters to create new parameters
!!
!! @author Raul L Bernal-Gonzalez
!!
subroutine calc_new_parameters(alpha, beta, parameters, new_params)
    implicit none
    real(dp), intent(in) :: alpha(:,:) !< alpha matrix
    real(dp), intent(in) :: beta(:) !< beta vector
    real(dp), intent(in) :: parameters(:) !< current parameters
    real(dp), intent(out), allocatable :: new_params(:) !< new parameters

    real(dp), allocatable :: delta_params(:), work(:,:)
    integer :: info, n_param, i, j

    allocate(new_params, mold=parameters)
    allocate(delta_params, mold=beta)
    allocate(work, mold=alpha)

    ! set work equal to alpha to prevent changing alpha
    work = alpha
    ! number of parameters
    n_param = size(parameters)
    ! use lapack to inverse work
    call dpotrf('U', n_param, work, n_param, info)
    ! check if call was successful
    if(info /= 0) then
        print*, 'Error calling dportf: ', info
        call exit(0)
    end if
    call dpotri('U', n_param, work, n_param, info)
    ! check if call was successful
    if(info /= 0) then
        print*, 'Error calling dportf: ', info
        call exit(0)
    end if

    ! rebuild the matrix from upper triangular half
    do i = 1, n_param
        do j = 1, n_param
                work(j,i) = work(i,j)
        end do
    end do
    ! multiply alpha inverse by beta to get deltas
    delta_params = matmul(work, beta)

    ! add deltas to get new parameters
    new_params = parameters + delta_params
end subroutine calc_new_parameters

!!
!> @brief set_alpha_prime
!!
!! Function to multiply the lambda fudge factor to the diagonal of alpha to
!! create alpha-prime
!!
!! @author Raul L Bernal-Gonzalez
function set_alpha_prime(alpha, lambda) result(alpha_prime)
    implicit none
    real(dp), intent(in) :: alpha(:,:)
    real(dp), intent(in) :: lambda
    real(dp), allocatable :: alpha_prime(:,:)

    integer :: i, j

    allocate(alpha_prime, mold=alpha)

    alpha_prime = alpha
    do i = 1, size(alpha, 1)
        do j = 1, size(alpha,2)
            if(i == j) then
                alpha_prime(i,j) = alpha(i,j)*(1 + lambda) ! diagonal
            end if
        end do
    end do
end function set_alpha_prime
end module
