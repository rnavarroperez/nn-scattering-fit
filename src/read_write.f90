!!
!> @brief Output and input handling
!!
!! These subroutines will deal with reading in data and outputting results
!!
!! @author Raul L Bernal-Gonzalez
!!
module read_write

use precisions, only: dp
use amplitudes, only: em_np_amplitudes, em_pp_amplitudes
use nn_phaseshifts, only: momentum_cm
use constants, only: pi
implicit none

private

public :: print_em_np_amplitudes, print_em_pp_amplitudes

contains

    !!
    !> @brief      Subroutine to test the output of em_np_amplitudes
    !!
    !! Makes a grid of t_lab from 1-350 in steps of 50 and theta from
    !! 10pi/180 - 180pi/180 in steps of 10
    !!
    !! @author     Raul L Bernal-Gonzalez
    !!
subroutine print_em_np_amplitudes()
    implicit none
    real(dp) :: momentum, angle, r
    complex(dp) :: a_s, b_s, c_s, d_s, e_s
    integer :: i, j, unit
    character(len=128) :: fmt ! format specifier

    ! open file for output
    open(newunit=unit, file='new_out_np.dat', status='unknown')
    ! format
    fmt = '(I3, 4x, F21.15, 5(5x, E21.15,SP, E21.15, "i"))'
    ! file header
    write(unit, '(2(A,4x), 5(40x, A))') 't_lab', 'theta', 'a', 'b', 'c', 'd', 'e'

    do i = 50, 350, 50
        r = real(i, kind=dp)
        momentum = momentum_cm(r, 'np')
        do j = 10, 180, 10
            angle =  j*pi/180.0_dp
            call em_np_amplitudes(momentum, angle, a_s, b_s, c_s, d_s, e_s)
            write(unit, fmt) i, angle, a_s, b_s, c_s, d_s, e_s
        end do
    end do
end subroutine print_em_np_amplitudes

subroutine print_em_pp_amplitudes()
        implicit none
        real(dp) :: momentum, angle, r
        complex(dp) :: a_s, b_s, c_s, d_s, e_s
        integer :: i, j, unit
        character(len=128) :: fmt ! format specifier

        ! open file for output
        open(newunit=unit, file='new_out_pp.dat', status='unknown')
        ! format
        fmt = '(I3, 4x, F21.15, 5(5x, E21.15,SP, E21.15, "i"))'
        ! file header
        write(unit, '(2(A,4x), 5(40x, A))') 't_lab', 'theta', 'a', 'b', 'c', 'd', 'e'

        do i = 50, 350, 50
            r = real(i, kind=dp)
            momentum = momentum_cm(r, 'pp')
            do j = 10, 180, 10
                angle =  j*pi/180.0_dp
                call em_pp_amplitudes(momentum, angle, a_s, b_s, c_s, d_s, e_s)
                write(unit, fmt) i, angle, a_s, b_s, c_s, d_s, e_s
            end do
        end do
end subroutine print_em_pp_amplitudes


end module read_write
