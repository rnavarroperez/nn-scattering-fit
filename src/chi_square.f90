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
    integer, intent(out) :: n_points !< number of data points in the total chi-square
    real(dp), intent(out) :: chi2 !< chi square for given parameters and model
    real(dp), intent(out), allocatable :: alpha(:,:), beta(:) !< matrices to minimize chi-square

    real(dp), allocatable :: all_chi(:)
    integer, allocatable :: all_n_points(:)
    real(dp), allocatable :: all_alpha(:,:,:), all_beta(:,:)
    integer :: i, n_exps, n_active_parameters
    type(nn_experiment) :: my_experiment
    

    ! get the number of experiments
    n_exps = size(experiments)
    ! get the number of parameters
    n_active_parameters = count(mask)
    ! allocate memory
    allocate(all_chi(n_exps))
    allocate(all_n_points(n_exps))
    allocate(all_beta(n_active_parameters, n_exps), all_alpha(n_active_parameters, n_active_parameters, n_exps))
    ! initialize all_alpha and all_beta
    all_alpha = 0
    all_beta = 0
    ! initialze arrays
    all_chi = 0
    all_n_points = 0
!---Parallel section------------------------------------------------------------
    !$omp parallel default(none) private(i, chi2, n_points, beta, alpha, my_experiment) &
    !$omp & shared(parameters, mask, model, all_chi, all_n_points, all_alpha, all_beta, experiments)
    !$omp do schedule(dynamic)
    do i = 1, size(experiments) ! calculate chi-square, alpha, and beta for each experiment
        if (experiments(i)%rejected) cycle
        if (experiments(i)%channel == 'nn' .and. model%potential_type == 'delta_shell') cycle ! delta-shell potentials don't have a nn channel (yet)
        if (experiments(i)%obs_type == 'dbe' .and. (.not. model%fit_deuteron) ) cycle !skip deuteron if the model says so
        call filter_by_energy(experiments(i), model%t_lab_limit, my_experiment)
        call experiment_chi_square(my_experiment, parameters, mask, model, chi2, alpha, beta, n_points)
        all_chi(i) = chi2
        all_n_points(i) = n_points
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
    chi2 = sum(all_chi)
    n_points = sum(all_n_points)
    alpha = sum(all_alpha, dim=3)!, mask=.not.isnan(all_alpha))
    beta = sum(all_beta, dim=2)!, mask=.not.isnan(all_beta))
end subroutine total_chi_square

subroutine filter_by_energy(full_experiment, t_lab_limit, filtered_experiment)
    implicit none
    type(nn_experiment), intent(in) :: full_experiment
    real(dp), intent(in) :: t_lab_limit
    type(nn_experiment), intent(out) :: filtered_experiment
    
    integer :: filtered_n_data, i, counter

    filtered_n_data = count(full_experiment%data_points(:)%t_lab <= t_lab_limit)
    filtered_experiment%n_data = filtered_n_data
    filtered_experiment%sys_error = full_experiment%sys_error
    filtered_experiment%obs_type = full_experiment%obs_type
    filtered_experiment%channel = full_experiment%channel
    filtered_experiment%rejected = full_experiment%rejected
    filtered_experiment%year = full_experiment%year
    filtered_experiment%reference = full_experiment%reference
    allocate(filtered_experiment%data_points(1:filtered_n_data))
    counter = 0
    do i = 1, full_experiment%n_data
        if (full_experiment%data_points(i)%t_lab <= t_lab_limit) then
            counter = counter + 1
            filtered_experiment%data_points(counter) = full_experiment%data_points(i)
        endif
    enddo
end subroutine filter_by_energy

subroutine experiment_chi_square(experiment, parameters, mask, model, chi2, alpha, beta, n_points)
    implicit none
    type(nn_experiment), intent(in) :: experiment
    real(dp), intent(in), dimension(:) :: parameters
    logical, intent(in), dimension(:) :: mask
    type(nn_model), intent(in) :: model
    real(dp), intent(out) :: chi2
    real(dp), intent(out), allocatable, dimension(:, :) :: alpha
    real(dp), intent(out), allocatable, dimension(:) :: beta
    integer, intent(out) :: n_points

    type(kinematics) :: kine
    real(dp), allocatable, dimension (:) :: exp_values, sigmas, theory_values
    real(dp) :: z_scale, chi2_sys_error
    real(dp) :: znum, zden, sys_error
    real(dp), allocatable :: d_theory(:), all_derivatives(:,:)
    logical :: float
    integer :: i!, j, k, l, m, n_active_parameters

    n_points = experiment%n_data

    allocate(exp_values(1:n_points))
    allocate(sigmas, theory_values, mold=exp_values)
    allocate(all_derivatives(1:size(parameters), n_points))

    kine%channel = experiment%channel
    kine%type = experiment%obs_type
    ! get systematic error per experiment
    sys_error = experiment%sys_error
    
    znum = 0._dp
    zden = 0._dp
    do i = 1, n_points
        kine%t_lab = experiment%data_points(i)%t_lab
        kine%angle = experiment%data_points(i)%theta
        kine%wave = experiment%data_points(i)%wave
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
    chi2 = sum(((exp_values - theory_values*z_scale)/sigmas)**2)
    ! check if data is floated
    if(float) then ! add the systematic error contribution to the chi-square
        chi2 = chi2 + chi2_sys_error
        n_points = n_points + 1
    end if

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
