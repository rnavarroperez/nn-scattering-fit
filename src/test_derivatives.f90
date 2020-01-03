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

public :: test_all_derivatives, context

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
                print'(2i5, 4f18.11)', j, i, d_analytic(i), df_dx, diff, err
                counter = counter + 1
            endif
        enddo
        if (counter == 0) then
            print*, 'all analytic derivatives match numerical calculation in target', j
        endif
        counter = 0
    enddo
end subroutine test_all_derivatives

end module test_derivatives