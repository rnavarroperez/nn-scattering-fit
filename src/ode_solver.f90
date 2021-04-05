!!
!> @brief      Solving Ordinary Differential Equations
!!
!! Subroutines to solve ODES with Runge Kutta
!!
!! @author     Rodrigo Navarro Perez
!!
module ode_solver
use precisions, only : dp
use delta_shell, only : nn_model
implicit none
private

public :: solve_runge_kutta_4, solve_runge_kutta_2

!!
!> @brief      Interface for receiving the ODE to solve
!!
!! The ODE to solve is receiving by the Runge Kutta algorithms
!! as a function that returns the derivative of the n_variables
!! that are being solved numerically
!!
!! A small modification had to be made from the more general
!! RK implementation to allow for a NN potential to be part 
!! of the arguments
!!
!> @author     Rodrigo Navarro Perez
!!
interface
    function func(r, t, work, model) result(f)
        use precisions, only : dp
        use delta_shell, only : nn_model
        implicit none
        real(dp), intent(in) :: r(:) !< Dependent variable(s) that we are solving for
        real(dp), intent(in) :: t !< Independent variable
        real(dp), intent(in) :: work(:) !< Array with everything necessary to calculate the derivatives of the dependent variable
        type(nn_model), intent(in) :: model !< NN potential (added to work with the NN scattering code)
        real(dp), allocatable :: f(:) !< Derivative(s) of the dependent variable(s)
    end function func
end interface

contains

!!
!> @brief      4th order Runge Kutta
!!
!! Implementation of 4h order Runge Kutta to solve a set of Ordinary Differential Equations
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine solve_runge_kutta_4(f, work, model, t_i, t_f, n, r_i, t, r)
    implicit none
    procedure(func) :: f !< Function that calculates the derivatives of the dependent variable(s)
    real(dp), intent(in) :: work(:) !< Array with everything necessary to calculate the derivatives of the dependent variable
    type(nn_model), intent(in) :: model !< NN potential (added to work with the NN scattering code)
    real(dp), intent(in) :: t_i !< Initial value of the independent variable
    real(dp), intent(in) :: t_f !< Final value of the independent variable
    integer, intent(in) :: n !< Number of points where the solution will be evaluated
    real(dp), intent(in) :: r_i(:) !< Initial value(s) of the dependent variable(s)
    real(dp), intent(out), allocatable :: t(:) !< Sampled values of the independent variable
    real(dp), intent(out), allocatable :: r(:,:) !< Sampled values of the dependent variable(s). The numerical solution

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

!!
!> @brief      2nd order Runge Kutta
!!
!! Implementation of 2nd order Runge Kutta to solve a set of Ordinary Differential Equations
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine solve_runge_kutta_2(f, work, model, t_i, t_f, n, r_i, t, r)
    implicit none
    procedure(func) :: f !< Function that calculates the derivatives of the dependent variable(s)
    real(dp), intent(in) :: work(:) !< Array with everything necessary to calculate the derivatives of the dependent variable
    type(nn_model), intent(in) :: model !< NN potential (added to work with the NN scattering code)
    real(dp), intent(in) :: t_i !< Initial value of the independent variable
    real(dp), intent(in) :: t_f !< Final value of the independent variable
    integer, intent(in) :: n !< Number of points where the solution will be evaluated
    real(dp), intent(in) :: r_i(:) !< Initial value(s) of the dependent variable(s)
    real(dp), intent(out), allocatable :: t(:) !< Sampled values of the independent variable
    real(dp), intent(out), allocatable :: r(:,:) !< Sampled values of the dependent variable(s). The numerical solution

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