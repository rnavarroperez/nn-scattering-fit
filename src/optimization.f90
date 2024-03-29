!!
!> @brief      optimization
!!
!! Subroutines and functions to implement the Lavenberg-Marquardt method for
!! minimizing the chi-square
!!
!! @author Raul L Bernal-Gonzalez
!!
module optimization
use precisions, only: dp
use exp_data, only: nn_experiment, read_database, init_ex_em_amplitudes
use delta_shell, only: nn_model
use av18, only: set_av18_potential, default_av18_mask=>default_mask
use chi_square, only: total_chi_square
use read_write, only: write_potential_setup, setup_from_namelist
implicit none

private
public :: lavenberg_marquardt, invert_alpha, setup_optimization, covariance_matrix, adiabatic_fit
contains

subroutine setup_optimization(model, parameters, mask, database, save_results, output_name)
    use iso_fortran_env, only : output_unit
    implicit none
    type(nn_model), intent(out) :: model
    real(dp), intent(out), allocatable, dimension(:) ::  parameters
    logical, intent(out), allocatable, dimension(:) :: mask
    type(nn_experiment), intent(out), allocatable, dimension(:) :: database
    logical, intent(out) :: save_results
    character(len=*), intent(out) :: output_name
    
    integer :: n_arguments
    character(len=1024) :: namelist_file
    character(len=1024) :: database_file = 'database/granada_database.dat'

    n_arguments = command_argument_count()

    if (n_arguments == 0) then
        print*, 'No namelist file was given as an argument, using default setup'
        print*, ''
        save_results = .true.
        output_name = 'av18'
        call set_av18_potential(model, parameters)
        allocate(mask(1:size(parameters)))
        mask = default_av18_mask
    else if(n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        call setup_from_namelist(namelist_file, model, parameters, mask, database_file, &
            save_results, output_name)
    else
        print*, 'The program takes either zero or one argument.'
        print*, 'See documentation for details.'
        print*, 'Stopping the program'
        stop
    endif
    call write_potential_setup(model, parameters, mask, output_unit)
    print*, 'Reading database from: ', trim(database_file)
    print*, ''
    call read_database(database_file, database)
    print*, 'Calculating EM amplitudes for each experimental point ...'
    call init_ex_em_amplitudes(database)
    print*, '... done'
    print*, 'Optimization setup finished'
    print*, ''
end subroutine setup_optimization

!!
!> @brief   lavenberg-marquardt
!!
!! This subroutine contains the main loop that implements the
!! lavenberg-marquardt method. It will return the minimize chi-square, covariance,
!! and the optimal parameters for the input model
!!
!! @author Raul L Bernal-Gonzalez
!!
subroutine lavenberg_marquardt(experiments, mask, model, parameters, n_points, chi2, covariance)
    use iso_fortran_env, only : output_unit
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: experiments !< input experiment data
    logical, intent(in), dimension(:) :: mask
    type(nn_model), intent(in) :: model !< potential model
    real(dp), intent(inout) :: parameters(:) !< potential model parameters
    integer, intent(out) :: n_points !< number of data points in the total chi-square
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model
    real(dp), intent(out), allocatable :: covariance(:,:) !< estimated covariance

    real(dp), allocatable :: alpha(:, :), beta(:)
    real(dp), allocatable :: prev_alpha(:, :), prev_beta(:), prev_parameters(:), alpha_prime(:, :)
    real(dp) :: lambda, prev_chi_ratio, chi_ratio
    integer :: limit, counter
    real(dp), parameter :: delta = 1.e-2_dp
    real(dp), parameter :: factor = 10._dp

    limit = 0
    counter = 0
    lambda = 1.e-3_dp
    call total_chi_square(experiments, parameters, mask, model, n_points, chi2, alpha, beta)
    chi_ratio = chi2/n_points
    
    call model%display_subroutine(parameters, mask, output_unit)
    print 1, 'chi^2:', chi2, 'N_data:', n_points, 'chi^2/N_data:', chi_ratio, 'counter:', &
        counter, 'lambda:', lambda, 'limit:', limit
   ! print*,
    
    allocate(prev_parameters, source=parameters)
    allocate(prev_alpha, source=alpha)
    allocate(prev_beta, source=beta)
    allocate(alpha_prime, mold=alpha)
    prev_chi_ratio = chi2/n_points
    do
        if(limit == 5 .or. lambda > 1.e+12_dp) exit
        alpha_prime = set_alpha_prime(prev_alpha, lambda)
        parameters = get_new_parameters(alpha_prime, prev_beta, prev_parameters, mask)
        call total_chi_square(experiments, parameters, mask, model, n_points, chi2, alpha, beta)
        chi_ratio = chi2/n_points
        ! determined whether to raise or lower lambda
        if(chi_ratio >= prev_chi_ratio) then
            lambda = lambda*factor
        else if(chi_ratio < prev_chi_ratio) then
            lambda = lambda/factor
            if((prev_chi_ratio - chi_ratio)*n_points <= delta) then
                ! increase limit by 1 if there is a negligible difference
                limit = limit + 1
            else ! we want consecutive negligible differences
                limit = 0
            end if
            prev_alpha = alpha
            prev_beta = beta
            prev_parameters = parameters
            prev_chi_ratio = chi_ratio
        end if
        counter = counter + 1
        call model%display_subroutine(parameters, mask, output_unit)
        print 1, 'chi^2:', chi2, 'N_data:', n_points, 'chi^2/N_data:', chi_ratio, 'counter:', &
            counter, 'lambda:', lambda, 'limit:', limit
        !print*,
    enddo
    alpha = prev_alpha
    parameters = prev_parameters
    covariance = covariance_matrix(alpha, mask)
    call model%display_subroutine(parameters, mask, output_unit, covariance)

1 format(1x,a,f20.8,2x,a,i5,2x,a,f20.8,2x,a,i5,2x,a,e11.4,2x,a,i2)
end subroutine lavenberg_marquardt

function covariance_matrix(alpha, mask) result(covariance)
    implicit none
    real(dp), intent(in), dimension(:, :) :: alpha
    logical, intent(in), dimension(:) :: mask
    real(dp), allocatable, dimension(:, :) :: covariance

    real(dp), allocatable, dimension(:, :) :: alpha_inverse
    integer :: i, j, k, l, n_parameters
    n_parameters = size(mask)
    allocate(alpha_inverse, mold = alpha)
    allocate(covariance(1:n_parameters, 1:n_parameters))
    
    alpha_inverse = invert_alpha(alpha)
    covariance = 0._dp
    l = 0
    do j = 1, n_parameters
        if (.not.mask(j)) cycle
        l = l + 1
        k = 0
        do i = 1, n_parameters
            if (.not.mask(i)) cycle
            k = k + 1
            covariance(i, j) = alpha_inverse(k, l)
        enddo
    enddo
end function covariance_matrix

!!
!> @brief   calc_new_parameters
!!
!! The function uses lapack library to inverse alpha and solve the linear equation
!! alpha*delta = beta. It then adds delta to the parameters to create new parameters
!!
!! @author Raul L Bernal-Gonzalez
!!
function get_new_parameters(alpha, beta, parameters, mask) result(new_parameters)
    implicit none
    real(dp), intent(in) :: alpha(:,:) !< alpha matrix
    real(dp), intent(in) :: beta(:) !< beta vector
    real(dp), intent(in) :: parameters(:) !< current parameters
    logical, intent(in) :: mask(:)
    real(dp), allocatable :: new_parameters(:) !< new parameters

    real(dp), allocatable :: delta_params(:), work(:,:)
    integer :: i, j

    allocate(new_parameters, source=parameters)
    allocate(delta_params, mold=beta)
    allocate(work, mold=alpha)

    
    ! multiply alpha inverse by beta to get deltas
    work = invert_alpha(alpha)
    delta_params = matmul(work, beta)
    ! add deltas to get new parameters
    j = 0
    do i = 1, size(parameters)
        if (.not.mask(i)) cycle
        j = j + 1
        new_parameters(i) = new_parameters(i) + delta_params(j)
    enddo
    
end function get_new_parameters

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

    integer :: i

    allocate(alpha_prime, source=alpha)

    do i = 1, size(alpha, 1)
        alpha_prime(i,i) = alpha_prime(i,i)*(1 + lambda) ! diagonal
    end do
end function set_alpha_prime

function invert_alpha(alpha) result(alpha_inv)
    implicit none
    real(dp), intent(in) :: alpha(:,:)
    real(dp), allocatable :: alpha_inv(:,:)
    integer :: i, j, info, n_param

    allocate(alpha_inv, mold=alpha)

    ! set alpha_inv equal to alpha to prevent changing alpha
    alpha_inv = alpha
    ! alpha is n_param x n_param matrix
    n_param = size(alpha,dim=1)

    ! invert alpha
    call dpotrf('U', n_param, alpha_inv, n_param, info)
    ! check if call was successful
    if(info /= 0) then
        print*, 'Error calling dportf: ', info
        stop
    end if
    call dpotri('U', n_param, alpha_inv, n_param, info)
    ! check if call was successful
    if(info /= 0) then
        print*, 'Error calling dpotri: ', info
        stop
    end if

    ! rebuild the matrix from upper triangular half
    do i = 1, n_param
        do j = i+1, n_param
                alpha_inv(j,i) = alpha_inv(i,j)
        end do
    end do
end function invert_alpha


subroutine adiabatic_fit(parameters, target_shape, n_steps, database, mask, model, &
    log_filename, initial_parameters, chi2, n_points, covariance)
    implicit none
    real(dp), intent(inout), dimension(:) :: parameters
    real(dp), intent(in), dimension(1:3) :: target_shape
    integer, intent(in) :: n_steps
    type(nn_experiment), intent(in), dimension(:) :: database
    logical, intent(in), dimension(:) :: mask
    type(nn_model), intent(in) :: model
    character(len=*), intent(in) :: log_filename
    real(dp), intent(out), allocatable, dimension(:) :: initial_parameters
    real(dp), intent(out) :: chi2
    integer, intent(out) :: n_points
    real(dp), intent(out), allocatable, dimension(:, :) :: covariance

    real(dp), dimension(1:3) :: delta_p
    integer :: unit, i
    

    delta_p = (target_shape - parameters(58:60))/n_steps
    allocate(initial_parameters, source=parameters) !make a copy of the initial parameters to later save them
    open(newunit=unit, file=log_filename)
    do i = 1, n_steps
        parameters(58:60) = parameters(58:60) + delta_p
        call lavenberg_marquardt(database, mask, model, parameters, n_points, chi2, covariance)
        write(unit, '(4f15.8,i8,f15.8)') parameters(58:60), chi2, n_points, chi2/n_points
    enddo
    close(unit)

end subroutine adiabatic_fit

end module optimization
