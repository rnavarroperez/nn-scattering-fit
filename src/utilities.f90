!!
!> @brief      abstract operations
!!
!! Set of functions and subroutines to perform abstract operations that
!! is not particular to any program. Things like reallocating arrays, 
!! or a kronecker delta, etc
!!
!! @author     Rodrigo Navarro Perez
!!
module utilities
use precisions, only : dp
implicit none
private
public :: double_2darray_allocation, trim_2d_array, int_to_logical, kronecker_delta
contains

!!
!> @brief      converts an integer to a logical
!!
!! Converts an integer to a logical. If the integer is zero,
!! false is returned; true is returned for any other value
!!
!! @return     a logical representation of an integer
!!
logical function int_to_logical(i) result(r)
    implicit none
    integer, intent(in) :: i !< an integer
    if (i == 0) then
        r = .false.
    else
        r = .true.
    endif    
end function int_to_logical

!!
!> @brief      Kronecker_delta delta
!!
!! The usual Kronecker delta \f$\delta_{i,j} = 1 \; {\rm if } \; i=j, \;  0\; {\rm if} \; i\neq j \f$
!!
!! @author     Rodrigo Navarro Perez
!!
integer function kronecker_delta(i, j) result(delta)
    implicit none
    integer, intent(in) :: i !< a integer
    integer, intent(in) :: j !< another integer
    if (i == j) then
        delta = 1
    else
        delta = 0
    endif
end function kronecker_delta

!!
!> @brief      doubles the size of a rank 2 array
!!
!! Given a 2 rank array of reals, returns an array with it's second dimension doubled
!!
!! The elements of the original array are kept in the first half of the new array 
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine double_2darray_allocation(array)
    implicit none
    real(dp), intent(inout), allocatable, dimension(:, :) :: array !< 2d array to extend on the second index

    real(dp), allocatable, dimension(:, :) :: temp
    integer, dimension(1:2) :: array_shape
    array_shape = shape(array)
    allocate(temp(1: array_shape(1), 1:2*array_shape(2)))
    temp(:, 1:array_shape(2)) = array
    call move_alloc(temp, array)
end subroutine double_2darray_allocation


!!
!> @brief      trims the second dimension of rank 2 array
!!
!! Given a rank 2 array of reals and a cut off index, trims along the second
!! dimension of the array all the elements above the cut off
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine trim_2d_array(cut_off, array)
    implicit none
    integer, intent(in) :: cut_off !< Size of the 2d array to keep (second index)
    real(dp), intent(inout), allocatable, dimension(:, :) :: array !< 2d array to trim

    real(dp), allocatable, dimension(:, :) :: temp
    if (cut_off > size(array, 2)) then
        stop 'cut_off has to be smaller than array second index size in trim_2d_array'
    endif
    allocate(temp(1:size(array,1), 1:cut_off))
    temp = array(:, 1:cut_off)
    call move_alloc(temp, array)
end subroutine trim_2d_array

end module utilities