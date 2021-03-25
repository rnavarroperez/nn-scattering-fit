module ode_solver
use precisions, only : dp
use delta_shell, only : nn_model
implicit none
private

public :: solve_runge_kutta_4, solve_runge_kutta_2

interface
    function func(r, t, work, model) result(f)
        use precisions, only : dp
        use delta_shell, only : nn_model
        implicit none
        real(dp), intent(in) :: r(:), t, work(:)
        type(nn_model), intent(in) :: model
        real(dp), allocatable :: f(:)
    end function func
end interface

contains

subroutine solve_runge_kutta_4(f, work, model, t_i, t_f, n, r_i, t, r)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: work(:)
    type(nn_model), intent(in) :: model
    real(dp), intent(in) :: t_i
    real(dp), intent(in) :: t_f
    integer, intent(in) :: n
    real(dp), intent(in) :: r_i(:)
    real(dp), intent(out), allocatable :: t(:), r(:,:)

    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), r_sol(:)
    integer :: n_variables, i
    real(dp) :: h, t_sol

    n_variables = size(r_i)

    if(allocated(t)) deallocate(t)
    if(allocated(t)) deallocate(t)

    allocate(k1(1:n_variables), k2(1:n_variables), k3(1:n_variables), k4(1:n_variables))
    allocate(t(1:n))
    allocate(r(1:n_variables, 1:n))
    allocate(r_sol(1:n_variables))

    h = (t_f - t_i)/n
    r_sol = r_i
    t_sol = t_i

    do i=1,n
        k1 = h*f(r_sol, t_sol, work, model)
        k2 = h*f(r_sol + 0.5_dp*k1, t_sol + 0.5_dp*h, work, model)
        k3 = h*f(r_sol + 0.5_dp*k2, t_sol + 0.5_dp*h, work, model)
        k4 = h*f(r_sol + k3, t_sol + h, work, model)
        r_sol = r_sol + (k1 + 2*k2 + 2*k3 + k4)/6._dp
        t_sol = t_sol + h
        t(i) = t_sol
        r(:, i) = r_sol
    enddo

end subroutine solve_runge_kutta_4

subroutine solve_runge_kutta_2(f, work, model, t_i, t_f, n, r_i, t, r)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: work(:)
    type(nn_model), intent(in) :: model
    real(dp), intent(in) :: t_i
    real(dp), intent(in) :: t_f
    integer, intent(in) :: n
    real(dp), intent(in) :: r_i(:)
    real(dp), intent(out), allocatable :: t(:), r(:,:)

    real(dp), allocatable :: k1(:), k2(:), r_sol(:)
    integer :: n_variables, i
    real(dp) :: h, t_sol

    n_variables = size(r_i)

    if(allocated(t)) deallocate(t)
    if(allocated(t)) deallocate(t)

    allocate(k1(1:n_variables), k2(1:n_variables))
    allocate(t(1:n))
    allocate(r(1:n_variables, 1:n))
    allocate(r_sol(1:n_variables))

    h = (t_f - t_i)/n
    r_sol = r_i
    t_sol = t_i

    do i=1,n
        k1 = h*f(r_sol, t_sol, work, model)
        k2 = h*f(r_sol + 0.5_dp*k1, t_sol + 0.5_dp*h, work, model)
        r_sol = r_sol + k2
        t_sol = t_sol + h
        t(i) = t_sol
        r(:, i) = r_sol
    enddo

end subroutine solve_runge_kutta_2

    
end module ode_solver