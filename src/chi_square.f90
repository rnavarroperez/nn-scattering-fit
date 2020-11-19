!!
!> @brief      Chi Square
!!
!! This module contains necessary functions
!! to evaluate and minimize the chi square
!!
!! @author Rodrigo Navarro Perez
!! @author Raul L Bernal-Gonzalez
!!
module chi_square
use precisions, only: dp
use exp_data, only : nn_experiment
use observables, only: observable, kinematics
use delta_shell, only : nn_model
use omp_lib
implicit none

private
public :: calc_chi_square
contains

!!
!> @brief   calc_chi_square
!!
!! This will calculate the chi-square along with the alpha, and beta
!! matrices for the current model parameters
!!
!! @author Rodrigo Navarro-Perez
!! @author Raul L Bernal-Gonzalez
!!
subroutine calc_chi_square(experiments, model_parameters, model, n_points, chi2, alpha, beta)
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: experiments !< input experiment data
    real(dp), intent(in) :: model_parameters(:) !< potential model parameters
    type(nn_model), intent(in) :: model !< potential model
    integer, intent(out) :: n_points !< number of data points in the total chi-square
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model
    real(dp), intent(out), allocatable :: alpha(:,:), beta(:) !< matrices to minimize chi-square

    real(dp), allocatable :: all_chi(:)
    integer, allocatable :: all_n_points(:)
    real(dp), allocatable :: all_alpha(:,:,:), all_beta(:,:)
    integer :: i, n_exps, n_param, j

    ! get the number of experiments
    n_exps = size(experiments)
    ! print*, 'exp ', n_exps
    ! get the number of parameters
    n_param = size(model_parameters)
    ! allocate memory
    allocate(all_chi(n_exps))
    allocate(all_n_points(n_exps))
    allocate(all_beta(n_exps, n_param), all_alpha(n_exps,n_param, n_param))
    allocate(alpha(n_param, n_param), beta(n_param))
    ! initialize all_alpha and all_beta
    all_alpha = 0
    all_beta = 0
    ! initialze arrays
    all_chi = 0
    all_n_points = 0
!---Parallel section------------------------------------------------------------
    !$omp parallel default(none) private(i, chi2, n_points, beta, alpha) &
    !$omp & shared(model_parameters, model, all_chi, all_n_points, all_alpha, all_beta, experiments)
    !$omp do schedule(dynamic)
    do i = 1, size(experiments) ! calculate chi-square, alpha, and beta for each experiment
        if (experiments(i)%rejected) cycle
        call sum_chi_square(experiments(i), model_parameters, model, n_points, chi2, alpha, beta)
        all_chi(i) = chi2
        all_n_points(i) = n_points
        all_alpha(i,:,:) = alpha
        all_beta(i,:) = beta
    end do
    !$omp end do
    !$omp end parallel
!---End parallel section--------------------------------------------------------
    ! reset alpha and beta
    alpha = 0
    beta = 0
    ! sum all experiments
    chi2 = sum(all_chi)
    n_points = sum(all_n_points)
    alpha = sum(all_alpha, dim=1)!, mask=.not.isnan(all_alpha))
    beta = sum(all_beta, dim=1)!, mask=.not.isnan(all_beta))

    ! do i = 1, size(alpha, 1)/3
    !     do j = 1, size(alpha, 1)/3
    !         write(*,"(f15.10)", advance='no') alpha(i,j)
    !     end do
    !     write(*,*)
    ! end do

    ! beta = sum(all_beta, 1) ! sum along the columns of all_beta
    ! do i = 1, size(all_beta, 1)
    !     do j = 1, size(all_beta, 2)
    !         print*,'b ', beta(j)
    !         beta(j) = beta(j) + all_beta(i,j)
    !         print*,'a ' ,beta(j)
    !     end do
    ! end do
end subroutine calc_chi_square

    !!
    !> @brief   sum_chi_square
    !!
    !! Sums the chi-square for all data in the given experiment,
    !! parameters, and model
    !!
    !! @author Rodrigo Navarro Perez
    !! @author Raul L Bernal-Gonzalez
    !!
subroutine sum_chi_square(experiment, model_parameters, model, n_points, chi2, alpha, beta)
    implicit none
    ! type(nn_experiment), intent(in), dimension(:) :: experiment !< input experiment data
    type(nn_experiment), intent(in) :: experiment !< input single experiment data
    real(dp), intent(in) :: model_parameters(:) !< potential model parameters
    type(nn_model), intent(in) :: model !< potential model
    integer, intent(out) :: n_points !< number of data points in the total chi square
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model
    real(dp), intent(out) :: alpha(:,:), beta(:)

    type(kinematics) :: kine
    real(dp), allocatable, dimension (:) :: exp_val, sigma, obs
    real(dp) :: z_scale, chi_sys_error_cont
    real(dp) :: znum, zden, sys_error
    real(dp), allocatable :: d_obs(:), all_d_obs(:,:)
    logical :: first_nan, first_nan_b

    logical :: float
    integer :: i, k, l, n_size, n_param

    first_nan = .true.
    first_nan_b = .true.
    ! initialize chi-squares at 0
    chi2 = 0
    n_points = 0
    ! initialize znum and zden
    znum = 0
    zden = 0
    ! initialize kine data per experiment
    kine%channel = experiment%channel
    kine%type = experiment%obs_type
    !if (kine%type == 'dbe' .or. kine%type == 'asl') return ! skipping deuteron binding energy during development process
    ! get systematic error per experiment
    sys_error = experiment%sys_error
    ! get number of experiments
    n_size = experiment%n_data
    ! get number of parameters
    n_param = size(model_parameters)
    ! allocate arrays
    allocate(exp_val(1: n_size))
    allocate(sigma, obs, mold=exp_val)
    allocate(all_d_obs(n_size, n_param))
    ! initialize alpha and beta
    alpha = 0
    beta = 0
    do i = 1, n_size ! goes through all points of nth experiment
        n_points = n_points + 1
        kine%t_lab = experiment%data_points(i)%t_lab
        kine%angle = experiment%data_points(i)%theta
        kine%em_amplitude = experiment%data_points(i)%em_amplitude
        exp_val(i) = experiment%data_points(i)%value
        sigma(i) = experiment%data_points(i)%stat_error
        call observable(kine, model_parameters, model, obs(i), d_obs)
        ! save the derivative of the corresponding observable
        all_d_obs(i, :) = d_obs
        ! calculate the numerator and denominator of the z-scale term
        znum = znum + (exp_val(i)*obs(i))/sigma(i)**2
        zden = zden + (obs(i)/sigma(i))**2
    end do
    if (n_size > 1) then
        call calc_z_scale(sys_error, znum, zden, z_scale, chi_sys_error_cont, float)
    else
        z_scale = 1
        chi_sys_error_cont = 0._dp
        float = .false.
    endif
    ! calculate chi-square for sinlge experiment
    chi2 = sum(((exp_val - obs*z_scale)/sigma)**2)
    ! check if data is floated
    if(float) then ! add the systematic error contribution to the chi-square
        chi2 = chi2 + chi_sys_error_cont
        n_points = n_points + 1
    end if
    ! calculate alpha and beta for single experiment
    do i = 1, n_size
        do k = 1, n_param
            do l = 1, n_param
                ! build alpha for single experiment
                ! derivative of one, times all the others
                alpha(k,l) = alpha(k,l) + ((all_d_obs(i,k)*all_d_obs(i,l))*z_scale**2)/sigma(i)**2
            end do
            ! build beta
            ! the derivative of chi-square with respect to the parameters
            beta(k) = beta(k) + ((exp_val(i) - obs(i)*z_scale)/sigma(i)**2)*(z_scale*all_d_obs(i,k))
        end do
    end do
end subroutine sum_chi_square

!!
!> @brief   calc_z_scale
!!
!! Caculates the Z scaling factor and the
!! contribution of the systematic error to
!! the chi-square
!!
!! @author Rodrigo Navarro Perez
!! @author Raul L Bernal-Gonzalez
!!
subroutine calc_z_scale(sys_error, znum, zden, z_scale, chi_sys_error_cont, float)
    implicit none
    real(dp), intent(in) :: sys_error !< experiment's systematic error
    real(dp), intent(in) :: znum !< denominator of Z scaling factor
    real(dp), intent(in) :: zden !< numerator of the Z scaling factor
    real(dp), intent(out) :: z_scale !< Z scaling factor
    real(dp), intent(out) :: chi_sys_error_cont !< contribution of the systematic error to the chi-square
    logical, intent(out) :: float !< if true, add the chi_sys_error_cont to the chi-squre point
    real(dp) :: num, den

    if(sys_error == 0) then ! check if data is absolute
        z_scale = 1._dp
        float = .false.
    else if (sys_error > 0.25_dp) then ! check for large systematic error
        z_scale = znum/zden
        float = .false.
    else ! if systematic error is greater than 0 but less than 0.25
        num = znum + 1._dp/(sys_error)**2
        den = zden + 1._dp/(sys_error)**2
        z_scale = num/den
        float = .true.
        chi_sys_error_cont = ((z_scale-1._dp)/sys_error)**2
        if (chi_sys_error_cont > 9._dp) then
            float = .false.
            z_scale = znum/zden
        end if
    end if
end subroutine calc_z_scale
end module chi_square
