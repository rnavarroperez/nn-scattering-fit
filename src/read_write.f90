!!
!> @brief Output and input handling
!!
!! These subroutines will deal with reading in data and outputting results
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
module read_write

use precisions, only: dp
use amplitudes, only: em_amplitudes
use observables, only: observable, kinematics
use nn_phaseshifts, only: all_phaseshifts
use delta_shell, only: nn_model
use constants, only: pi
implicit none

private

public :: print_em_amplitudes, print_observables, write_phases

contains

!!
!> @brief      Test the output of em_amplitudes
!!
!! Makes a grid of t_lab from 1-350 in steps of 50 and angle from
!! 10pi/180 - 180pi/180 in steps of 10
!!
!! This subroutine is for benchmarking purposes only and should not be 
!! used in production runs!
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine print_em_amplitudes(channel)
    implicit none
    character(len=2), intent(in) :: channel !< reaction channel (pp or np)
    real(dp) :: angle, t_lab, theta
    complex(dp), dimension(1:5) :: em_amplitude
    integer :: i, j, unit
    character(len=128) :: fmt ! format specifier

    if (channel/='pp' .and. channel/='np') then
        print*, channel
        stop 'invalid reaction channel in print_em_amplitudes'
    endif
    ! open file for output
    open(newunit=unit, file='new_out_'//channel//'.dat', status='unknown')
    ! format
    fmt = '(I3, 4x, F21.15, 5(5x, E21.15,SP, E21.15, "j"))'
    ! file header
    write(unit, '(2(A,4x), 5(40x, A))') 't_lab', 'theta', 'a', 'b', 'c', 'd', 'e'

    do i = 50, 350, 50
        t_lab = real(i, kind=dp)
        do j = 10, 180, 10
            angle = real(j, kind=dp)
            theta = angle*pi/180.0_dp
            em_amplitude = em_amplitudes(t_lab, angle, channel)
            write(unit, fmt) i, theta, em_amplitude
        end do
    end do
end subroutine print_em_amplitudes

!!
!> @brief      Test the output of observables
!!
!! Makes a grid of t_lab from 1-350 in steps of 50 and angle from
!! 10pi/180 - 180pi/180 in steps of 10
!!
!! This subroutine is for benchmarking purposes only and should not be 
!! used in production runs!
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine print_observables(parameters, model, channel)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    type(nn_model), intent(in) :: model !< nn scattering model
    character(len=2), intent(in) :: channel !< reaction channel (pp or np)

    integer, parameter :: n_observables = 26
    character(len=4), dimension(1:n_observables), parameter :: &
           obs_types = ['dsg ','dt  ','ayy ','d   ','p   ','azz ','r   '&
             ,'rt  ','rpt ','at  ','d0sk','nskn','nssn','nnkk','a   '&
             ,'axx ','ckp ','rp  ','mssn','mskn','azx ','ap  ','dtrt'&
             ,'sgt ','sgtt','sgtl'] 
    !---------------------------------------------------------------------------
    type(kinematics) :: kinematic
    real(dp) :: obs
    real(dp), allocatable :: d_obs(:)
    integer :: i, j, k, unit1, unit2

    open(newunit=unit1, file='new_'//channel//'_obs.dat', status='unknown')
    open(newunit=unit2, file='new_'//channel//'_d_obs.dat', status='unknown')
    kinematic%em_amplitude = 0
    kinematic%channel = channel
    do i = 50, 350, 50
        kinematic%t_lab = real(i, kind=dp)
        do j = 10, 180, 10
            kinematic%angle = real(j, kind=dp)
            do k = 1, size(obs_types)
                kinematic%type = obs_types(k)
                call observable(kinematic, parameters, model, obs, d_obs)
                write(unit1, *) i, j, obs_types(k), obs
                write(unit2, *) d_obs
            end do
        end do
    end do
    close(unit1)
    close(unit2)
end subroutine print_observables

!!
!> @brief      Test the output of observables
!!
!! Makes a grid of t_lab from 1-350 in steps of 50
!!
!! This subroutine is for benchmarking purposes only and should not be 
!! used in production runs!
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine write_phases(parameters, model, channel)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    type(nn_model), intent(in) :: model !< nn scattering model
    character(len=2), intent(in) :: channel !< reaction channel (pp or np)

    integer, parameter :: jmax = 20
    
    integer :: unit, i
    real(dp) :: t_lab
    real(dp), dimension(1:5, 1:jmax) :: phases
    real(dp), allocatable, dimension(:, :, :) :: d_phases
    
    open(newunit=unit, file='new_'//channel//'_aps.dat', status='unknown')
    do i = 50, 350, 50
        t_lab = real(i, kind=dp)
        call all_phaseshifts(model, parameters, t_lab, channel, phases, d_phases)
        write(unit, *) phases
    end do
    close(unit)  
end subroutine write_phases

end module read_write
