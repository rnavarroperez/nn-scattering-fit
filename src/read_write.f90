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
use utilities, only : double_2darray_allocation, trim_2d_array
implicit none

private

public :: print_em_amplitudes, print_observables, write_phases, read_montecarlo_parameters, write_montecarlo_phases

contains

!!
!> @brief      Writes Monte Carlo phases.
!!
!! Writes to a file the low angular momentum phaseshifts from a given
!! model, MC sample of potential parameter, reaction channel, and
!! laboratory energy in MeV
!!
!! Since the value of t_lab is used to determine the name of the file to which
!! the phases are written, t_lab has to be a positive integer of up to three digits
!! (the potential parameters are fitted to data up to 350 MeV anyway)
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine write_montecarlo_phases(model, mc_paramerters, channel, t_lab)
    implicit none
    type(nn_model), intent(in) :: model !< NN potential model
    real(dp), intent(in), dimension(:, :) :: mc_paramerters !< MC sample of potential parameters
    character(len=*), intent(in) :: channel !< reaction channel ('pp' or 'np')
    integer, intent(in) :: t_lab !< laboratory energy in (MeV)

    integer :: unit, imc
    character(len=3) :: t_lab_string
    character(len=50) :: file_name
    integer, parameter :: j_max = 4, n_waves = 5
    real(dp), dimension(n_waves, j_max) :: phases
    real(dp), allocatable, dimension(:, :, :) :: d_phases

    write(t_lab_string, '(i3.3)') t_lab
    file_name = 'mc_phases_tlab_'//t_lab_string//'.dat'
    open(newunit = unit, file = file_name, status = 'unknown')
    write(unit, '(a5,17a15)') 'imc', '1S0', '3P0', '1P1', '3P1', '3S1', 'EP1', '3D1', '1D2', '3D2', &
        '3P2', 'EP2', '3F2', '1F3', '3F3', '3D3', 'EP3', '3G3'
    do imc = 1, size(mc_paramerters, 2)
        call all_phaseshifts(model, mc_paramerters(:, imc), real(t_lab, kind=dp), channel, phases, d_phases)
        write(unit, '(i5,17f15.8)') imc, phases(1, 1), phases(5, 1), phases(:, 2:4)
    enddo
    close(unit)
end subroutine write_montecarlo_phases

!!
!> @brief      reads a set of samples of fitted parameters
!!
!! Given the name of a file containing a sample of fitted parameters
!! (generated with a legace code), reads and stores the parameters
!! in an rank 2 array.
!!
!! Paramters that were set (and fixed) to zero in the orginal fitting
!! in the legacy code were not sabed in the files to be read. In order
!! to knew where to put the parameters that are fixed to zero a mask
!! (array of logicals) needs to be given. The positions where the
!! paramter is meant to be zero are indicated by a .false. value in the
!! mask.
!!
!! Having the original set of parameters (with the zeros) the mask can
!! be created pretty quickly with mask = abs(parameters) > 1.e-12_dp.
!! See program in mc_phases.f90 for an example
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine read_montecarlo_parameters(param_mask, file_name, mc_params)
    use iso_fortran_env, only : iostat_end
    implicit none
    logical, intent(in), dimension(:) :: param_mask
    character(len=*), intent(in) :: file_name
    real(dp), intent(out), allocatable, dimension(:, :) :: mc_params

    integer :: counter, limit, unit, io_status, imc, n_np_data, n_pp_data, n_active_parameters
    integer :: p_counter, ip
    real(dp) :: pp_chi2, np_chi2
    real(dp), allocatable, dimension(:) :: mc_lambdas

    n_active_parameters = count(param_mask )
    allocate(mc_lambdas(1:n_active_parameters))
    counter = 0
    p_counter = 0
    limit = 1
    allocate(mc_params(size(param_mask), 1))
    open(newunit=unit, file=trim(file_name), status='old')
    do
        read(unit, *, iostat=io_status) imc, n_np_data, n_pp_data, pp_chi2, np_chi2, mc_lambdas
        if(imc == 0) cycle
        if (io_status == iostat_end) exit
        counter = counter + 1
        if (counter > limit) then
            call double_2darray_allocation(mc_params)
            limit = 2*limit
        endif
        do ip = 1, size(mc_params, 1)
            if (param_mask(ip)) then
                p_counter = p_counter + 1
                mc_params(ip, counter) = mc_lambdas(p_counter)
            else
                mc_params(ip, counter) = 0._dp
            endif
        enddo
        p_counter = 0
    enddo
    close(unit)
    call trim_2d_array(counter, mc_params)
end subroutine read_montecarlo_parameters

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

! subroutine write_chi(chi2, n_points, parameters)
!     implicit none
!     real(dp),
! end subroutine write_chi

end module read_write
