!!
!> @brief      compare analytical derivatives vs numerical ones
!!
!! @author     Rodrigo Navarro Perez
!!
module test_derivatives

use precisions, only : dp
use num_recipes, only : dfridr, context, func

implicit none

private

public :: test_all_derivatives, f_av18, df_av18, context

!!
!> @brief      interface for calling back functions that return arrays
!!
!! @author     Rodrigo Navarro Perez
!!
interface
    function func_array(data) result(r)
        use precisions, only: dp
        import context
        implicit none
        type(context), intent(in) :: data !< contains data for function evaluation
        real(dp), allocatable :: r(:) !< array to be returned by the function
    end function
end interface

contains

!!
!> @brief      tests all the derivatives of a given function
!!
!! Given a function f that depends on a set of parameters and the corresponding function df that
!! calculates the derivatives of such f with respect of all the parameters, compares such
!! derivatives with a numerical approximation.
!!
!! The function might have more than one target (i.e. the procedure being tested results in an 
!! array). 
!!
!! First the df function is used to calculate all the analytic derivatives of an specific target
!! with respect of all the parameters. Then all the numerical derivatives are calculated and 
!! compared with analytic results
!!
!! If the absolute difference between analytic and numerical, or the error in the numerical 
!! calculation, is greater than \f$ 10^{-10}\f$, the target, parameter, difference and error are 
!! printed in screen. Otherwise, a message indicates that all calculations are consistent with each
!! other.
!!
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine test_all_derivatives(f, df, data, n_targets, n_parameters)
    implicit none
    procedure(func) :: f !< function whose derivatives will be tested
    procedure(func_array) :: df !< function that returns analytical derivatives
    type(context), intent(inout) :: data !< data used to evaluate f and df
    integer, intent(in) :: n_targets !< number of targets in f
    integer, intent(in) :: n_parameters !< number of parameters that f depends on
    
    real(dp), parameter :: small = 1.e-010_dp
    real(dp), allocatable :: d_analytic(:)
    real(dp) :: x, h, df_dx, err, diff
    integer :: i, j, counter

    h = 0.1_dp
    counter = 0
    do j = 1, n_targets
        data%j = j
        d_analytic = df(data)
        do i = 1, n_parameters
            data%i = i
            x = data%x(i)
            call dfridr(f, data, x, h, df_dx, err)
            diff = abs(d_analytic(i) - df_dx)
            if (diff > small .or. err > small) then
                print*, j, i, d_analytic(i) - df_dx, err
                counter = counter + 1
            endif
        enddo
    enddo
    if (counter == 0) then
        print*, 'all analytic derivatives match numerical calculation'
    endif

end subroutine test_all_derivatives

!!
!> @brief      wrapper function for the av18 potential
!!
!! This wrapper function is used to test the derivatives of the av18_operator subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call the 
!! av18_operator subroutine. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which operator will be returned.
!!
!! @returns    one of the operators of the av18 potential at an specific radius
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_av18(x, data) result(r)
    use av18, only : av18_operator, n_parameters, n_operators
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for av18_operator

    real(dp) :: ap(1:n_parameters)
    real(dp) :: radius
    real(dp) :: v_nn(1:n_operators)
    real(dp) :: dv_nn(1:n_parameters, 1:n_operators)
    integer :: i_target, i_parameter

    ap = data%x
    radius = data%a
    i_parameter = data%i
    i_target = data%j

    ap(i_parameter) = x
    call av18_operator(ap, radius, v_nn, dv_nn)

    r = v_nn(i_target)
       
end function f_av18

!!
!> @brief      wrapper function for the derivatives of the av18 potential
!!
!! This wrapper function is used to test the derivatives of the av18_operator subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call the 
!! av18_operator subroutine. The same data of type context is used to receive which operator
!! will be returned
!!
!! @returns    the derivatives of one of the operators of the av18 potential at an specific radius
!!
!! @author     Rodrigo Navarro Perez
!!
function df_av18(data) result(r)
    use av18, only : av18_operator, n_parameters, n_operators
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for av18_operator
    real(dp), allocatable :: r(:)

    real(dp) :: ap(1:n_parameters)
    real(dp) :: radius
    real(dp) :: v_nn(1:n_operators)
    real(dp) :: dv_nn(1:n_parameters, 1:n_operators)
    integer :: i_target

    allocate(r(1:n_parameters))

    ap = data%x
    radius = data%a
    i_target = data%j

    call av18_operator(ap, radius, v_nn, dv_nn)

    r = dv_nn(:, i_target)
    
end function df_av18
    
end module test_derivatives