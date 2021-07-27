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
public :: total_chi_square
contains

!!
!> @brief   total_chi_square
!!
!! This will calculate the chi-square along with the alpha, and beta
!! matrices for the current model parameters
!!
!! @author Rodrigo Navarro-Perez
!! @author Raul L Bernal-Gonzalez
!!
subroutine total_chi_square(experiments, parameters, mask, model, n_points, chi2, alpha, beta)
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: experiments !< input experiment data
    real(dp), intent(in), dimension(:) :: parameters !< potential model parameters
    logical, intent(in), dimension(:) :: mask
    type(nn_model), intent(in) :: model !< potential model
    integer, intent(out), dimension(:) :: n_points !< number of data points in the total chi-square
    real(dp), intent(out), dimension(:) :: chi2 !< chi square for given parameters and model
    real(dp), intent(out), allocatable :: alpha(:,:), beta(:) !< matrices to minimize chi-square

    real(dp), allocatable :: all_chi(:, :)
    integer, allocatable :: all_n_points(:, :)
    real(dp), allocatable :: all_alpha(:,:,:), all_beta(:,:)
    integer :: i, n_exps, n_active_parameters
    

    ! get the number of experiments
    n_exps = size(experiments)
    ! get the number of parameters
    n_active_parameters = count(mask)
    ! allocate memory
    allocate(all_chi(1:size(chi2), n_exps))
    allocate(all_n_points(1:size(chi2), n_exps))
    allocate(all_beta(n_active_parameters, n_exps), all_alpha(n_active_parameters, n_active_parameters, n_exps))
    ! initialize all_alpha and all_beta
    all_alpha = 0
    all_beta = 0
    ! initialze arrays
    all_chi = 0
    all_n_points = 0
!---Parallel section------------------------------------------------------------
    !$omp parallel default(none) private(i, chi2, n_points, beta, alpha) &
    !$omp & shared(parameters, mask, model, all_chi, all_n_points, all_alpha, all_beta, experiments)
    !$omp do schedule(dynamic)
    do i = 1, size(experiments) ! calculate chi-square, alpha, and beta for each experiment
        if (experiments(i)%rejected) cycle
        if (experiments(i)%channel == 'nn' .and. model%potential_type == 'delta_shell') cycle ! delta-shell potentials don't have a nn channel (yet)
        call experiment_chi_square(experiments(i), parameters, mask, model, chi2, alpha, beta, n_points)
        all_chi(:, i) = chi2
        all_n_points(:, i) = n_points
        all_alpha(:, :, i) = alpha
        all_beta(:, i) = beta
    end do
    !$omp end do
    !$omp end parallel
!---End parallel section--------------------------------------------------------
    if(allocated(alpha)) deallocate(alpha)
    if(allocated(beta)) deallocate(beta)
    allocate(alpha(1:n_active_parameters, 1:n_active_parameters))
    allocate(beta(1:n_active_parameters))
    ! sum all experiments
    chi2 = sum(all_chi, dim=2)
    n_points = sum(all_n_points, dim=2)
    alpha = sum(all_alpha, dim=3)!, mask=.not.isnan(all_alpha))
    beta = sum(all_beta, dim=2)!, mask=.not.isnan(all_beta))
end subroutine total_chi_square

subroutine experiment_chi_square(experiment, parameters, mask, model, chi2, alpha, beta, n_points)
    implicit none
    type(nn_experiment), intent(in) :: experiment
    real(dp), intent(in), dimension(:) :: parameters
    logical, intent(in), dimension(:) :: mask
    type(nn_model), intent(in) :: model
    real(dp), intent(out), dimension(:) :: chi2
    real(dp), intent(out), allocatable, dimension(:, :) :: alpha
    real(dp), intent(out), allocatable, dimension(:) :: beta
    integer, intent(out), dimension(:) :: n_points

    type(kinematics) :: kine
    real(dp), allocatable, dimension (:) :: exp_values, sigmas, theory_values
    real(dp) :: z_scale, chi2_sys_error
    real(dp) :: znum, zden, sys_error
    real(dp), allocatable :: d_theory(:), all_derivatives(:,:)
    logical :: float
    integer :: i!, j, k, l, m, n_active_parameters
    integer :: my_n_points
    real(dp) :: my_chi2

    my_n_points = experiment%n_data

    allocate(exp_values(1:my_n_points))
    allocate(sigmas, theory_values, mold=exp_values)
    allocate(all_derivatives(1:size(parameters), my_n_points))

    kine%channel = experiment%channel
    kine%type = experiment%obs_type
    ! get systematic error per experiment
    sys_error = experiment%sys_error
    
    znum = 0._dp
    zden = 0._dp
    do i = 1, my_n_points
        kine%t_lab = experiment%data_points(i)%t_lab
        kine%angle = experiment%data_points(i)%theta
        kine%em_amplitude = experiment%data_points(i)%em_amplitude
        exp_values(i) = experiment%data_points(i)%value
        sigmas(i) = experiment%data_points(i)%stat_error
        call observable(kine, parameters, model, theory_values(i), d_theory)
        ! save the derivative of the corresponding observable
        all_derivatives(:, i) = d_theory
        ! calculate the numerator and denominator of the z-scale term
        znum = znum + (exp_values(i)*theory_values(i))/sigmas(i)**2
        zden = zden + (theory_values(i)/sigmas(i))**2
    enddo
    call calc_z_scale(sys_error, znum, zden, z_scale, chi2_sys_error, float)

    call calculate_alpha_beta(exp_values, theory_values, sigmas, all_derivatives, z_scale, mask, alpha, beta)

    ! calculate chi-square for sinlge experiment
    my_chi2 = sum(((exp_values - theory_values*z_scale)/sigmas)**2)
    ! check if data is floated
    if(float) then ! add the systematic error contribution to the chi-square
        my_chi2 = my_chi2 + chi2_sys_error
        my_n_points = my_n_points + 1
    end if

    chi2 = 0._dp
    n_points = 0
    select case(trim(kine%channel))
    case('pp')
        chi2(1) = my_chi2
        n_points(1) = my_n_points
    case('np')
        select case(trim(kine%type))
        case('dbe')
            chi2(3) = my_chi2
            n_points(3) = my_n_points
        case('asl')
            chi2(4) = my_chi2
            n_points(4) = my_n_points
        case default
            chi2(2) = my_chi2
            n_points(2) = my_n_points
        end select
    case('nn')
        chi2(5) = my_chi2
        n_points(5) = my_n_points
    case default
        print*, 'Unrecognized reaction channel in experiment_chi_square'
        stop
    end select

end subroutine experiment_chi_square

subroutine calculate_alpha_beta(exp_values, theory_values, sigmas, all_derivatives, z_scale, mask, alpha, beta)
    implicit none
    real(dp), intent(in), dimension(:) :: exp_values
    real(dp), intent(in), dimension(:) :: theory_values
    real(dp), intent(in), dimension(:) :: sigmas
    real(dp), intent(in), dimension(:, :) :: all_derivatives
    real(dp), intent(in) :: z_scale
    logical, intent(in), dimension(:) :: mask
    real(dp), intent(out), allocatable, dimension(:, :) :: alpha
    real(dp), intent(out), allocatable, dimension(:) :: beta
    
    integer :: n_active_parameters, n_parameters
    integer :: i, j, k, l, m

    n_active_parameters = count(mask)
    n_parameters = size(all_derivatives, 1)

    allocate(alpha(1:n_active_parameters, 1:n_active_parameters))
    allocate(beta(1:n_active_parameters))
    alpha = 0._dp
    beta = 0._dp
    do i = 1, size(exp_values)
        j = 0
        do k = 1, n_parameters
            if (.not.mask(k)) cycle
            j = j + 1
            l = 0
            do m = 1, n_parameters
                if (.not.mask(m)) cycle
                l = l + 1
                alpha(l, j) = alpha(l, j) + ((all_derivatives(m, i)*all_derivatives(k, i))*z_scale**2)/sigmas(i)**2
            enddo
            beta(j) = beta(j) + ((exp_values(i) - theory_values(i)*z_scale)/sigmas(i)**2)*(z_scale*all_derivatives(k, i))
        enddo
    enddo  

end subroutine calculate_alpha_beta

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
        chi_sys_error_cont = 0._dp
    else if (sys_error > 0.25_dp) then ! check for large systematic error
        z_scale = znum/zden
        float = .false.
        chi_sys_error_cont = 0._dp
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
