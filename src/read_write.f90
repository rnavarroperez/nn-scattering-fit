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
! use observables
use String_Functions, only: lower
use constants, only: pi
implicit none

private

public :: print_em_np_amplitudes, print_em_pp_amplitudes!, print_observables

contains

    !!
    !> @brief      Subroutines to test the output of em_np_amplitudes
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
    fmt = '(I3, 4x, F21.15, 5(5x, E21.15,SP, E21.15, "j"))'
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
        fmt = '(I3, 4x, F21.15, 5(5x, E21.15,SP, E21.15, "j"))'
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

!!
!> @brief      Subroutine to test the output of observables
!!
!! Makes a grid of t_lab from 1-350 in steps of 50 and theta from
!! 10pi/180 - 180pi/180 in steps of 10
!!
!! @author     Raul L Bernal-Gonzalez
!!
! subroutine print_observables()
!     implicit none
!     character(len=4), dimension(1:26), parameter :: &
!            obs_types = ['DSG ','DT  ','AYY ','D   ','P   ','AZZ ','R   '&
!              ,'RT  ','RPT ','AT  ','D0SK','NSKN','NSSN','NNKK','A   '&
!              ,'AXX ','CKP ','RP  ','MSSN','MSKN','AZX ','AP  ','DTRT'&
!              ,'SGT ','SGTT','SGTL'] !< array of possible observable characters
!     !---------------------------------------------------------------------------
!     real(dp) :: momentum, angle, obs
!     real(dp) :: pre = -1
!     real(dp), allocatable :: d_obs(:)
!     character(len=128) type
!     integer :: i, j, k, unit1, unit2, unit3, unit4, unit5, unit6
!     real(dp), allocatable :: phases(:,:)

!     open(newunit=unit1, file='new_pp_obs.dat', status='unknown')
!     open(newunit=unit2, file='new_pp_d_obs.dat', status='unknown')
!     open(newunit=unit3, file='new_np_obs.dat', status='unknown')
!     open(newunit=unit4, file='new_np_d_obs.dat', status='unknown')
!     open(newunit=unit5, file='new_pp_aps.dat', status='unknown')
!     open(newunit=unit6, file='new_np_aps.dat', status='unknown')

!     ! pp
!     do i = 50, 350, 50
!         momentum = real(i, kind=dp)
!         do j = 10, 180, 10
!             angle = real(j, dp)
!             do k = 1, size(obs_types)
!                 type = lower(obs_types(k))
!                 call observable(momentum, pre, angle, type, 'pp', obs, d_obs)
!                 write(unit1, *) i, j, obs_types(k), obs
!                 write(unit2, *) d_obs
!                 pre = momentum
!             end do
!         end do
!     end do

!     ! np
!     do i = 50, 350, 50
!         momentum = real(i, kind=dp)
!         do j = 10, 180, 10
!             angle = real(j, kind=dp)
!             do k = 1, size(obs_types)
!                 type = lower(obs_types(k))
!                 call observable(momentum, pre, angle, type, 'np', obs, d_obs)
!                 write(unit3, *) i, j, obs_types(k), obs
!                 write(unit4, *) d_obs
!                 pre = momentum
!             end do
!         end do
!     end do

!     do i = 50, 350, 50
!         momentum = real(i, kind=dp)
!         call just_phases(momentum, 'np', phases)
!         write(unit6, *) phases
!     end do
! end subroutine print_observables
end module read_write
