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
use av18, only : set_av18_potential
use delta_shell, only : set_ds_potential
use av18_compatibility, only : write_marias_format
use chiral_potential, only : vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9
implicit none

private

public :: print_em_amplitudes, print_observables, write_phases, read_montecarlo_parameters, &
    write_montecarlo_phases, print_phases, write_potential_setup, setup_from_namelist, &
    write_optimization_results, plot_potential_components, plot_potential_partial_waves, &
    write_chiral_kernals

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
!> @brief      Saves list of phases to a file
!!
!! Writes into a file 3 tables of phase-shifts with error bars.
!!
!! The phase-shifts are calculated at the 11 'canonical' energies
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine write_phases(potential, parameters, covariance, file_name)
    implicit none
    type(nn_model), intent(in) :: potential !< nn scattering model
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    real(dp), intent(in), dimension(:, :) :: covariance !< Covariance matrix
    character(len=*), intent(in) :: file_name

    integer, parameter :: n_energies = 11
    real(dp), parameter, dimension(1:n_energies) :: energies = [1._dp, 5._dp, 10._dp, 25._dp, &
        50._dp, 100._dp, 150._dp, 200._dp, 250._dp, 300._dp, 350._dp]
    integer, parameter :: jmax = 4
    real(dp), dimension(1:5, 1:jmax) :: phases, phases_errors
    real(dp), allocatable, dimension(:, :, :) :: d_phases

    integer :: unit, i, j, k
    real(dp) :: t_lab

    open(newunit=unit, file=file_name, status='unknown')
    write(unit, *) 'pp phase shifts'
    write(unit, '(9a13)') 'T_lab', '1S0', '1D2', '3P0', '3P1', '3P2', 'Ep2', '3F2', '3F3'
    do i = 1, n_energies
        t_lab = energies(i)
        call all_phaseshifts(potential, parameters, t_lab, 'pp', phases, d_phases)
        phases = phases*180/pi
        d_phases = d_phases*180/pi
        do j=1, size(phases, 2)
            do k=1, size(phases, 1)
                phases_errors(k, j) = propagated_error_bar(d_phases(:, k, j), covariance)
            enddo
        enddo
        write(unit, '(9f13.6)') energies(i), phases(1, 1), phases(1, 3), phases(5, 1), phases(2, 2), &
            phases(3, 3), phases(4, 3), phases(5, 3), phases(2, 4)
        write(unit, '(13x,8f13.6)') phases_errors(1, 1), phases_errors(1, 3), phases_errors(5, 1), &
        phases_errors(2, 2), phases_errors(3, 3), phases_errors(4, 3), phases_errors(5, 3), &
        phases_errors(2, 4)
    end do
    write(unit,*)
    write(unit, *) 'np T=1 phase shifts'
    write(unit, '(9a13)') 'T_lab', '1S0', '1D2', '3P0', '3P1', '3P2', 'Ep2', '3F2', '3F3'
    do i = 1, n_energies
        t_lab = energies(i)
        call all_phaseshifts(potential, parameters, t_lab, 'np', phases, d_phases)
        phases = phases*180/pi
        d_phases = d_phases*180/pi
        do j=1, size(phases, 2)
            do k=1, size(phases, 1)
                phases_errors(k, j) = propagated_error_bar(d_phases(:, k, j), covariance)
            enddo
        enddo
        write(unit, '(9f13.6)') energies(i), phases(1, 1), phases(1, 3), phases(5, 1), phases(2, 2), &
            phases(3, 3), phases(4, 3), phases(5, 3), phases(2, 4)
        write(unit, '(13x,8f13.6)') phases_errors(1, 1), phases_errors(1, 3), phases_errors(5, 1), &
        phases_errors(2, 2), phases_errors(3, 3), phases_errors(4, 3), phases_errors(5, 3), &
        phases_errors(2, 4)
    end do
    write(unit,*)
    write(unit, *) 'np T=0 phase shifts'
    write(unit, '(10a13)') 'T_lab', '1P1', '1F3', '3S1', 'Ep1', '3D1', '3D2', '3D3', 'Ep3', '3G3'
    do i = 1, n_energies
        t_lab = energies(i)
        call all_phaseshifts(potential, parameters, t_lab, 'np', phases, d_phases)
        phases = phases*180/pi
        d_phases = d_phases*180/pi
        do j=1, size(phases, 2)
            do k=1, size(phases, 1)
                phases_errors(k, j) = propagated_error_bar(d_phases(:, k, j), covariance)
            enddo
        enddo
        write(unit, '(10f13.6)') energies(i), phases(1, 2), phases(1, 4), phases(3, 2), phases(4, 2), &
            phases(5, 2), phases(2, 3), phases(3, 4), phases(4, 4), phases(5, 4)
        write(unit, '(13x,9f13.6)') phases_errors(1, 2), phases_errors(1, 4), phases_errors(3, 2), &
            phases_errors(4, 2), phases_errors(5, 2), phases_errors(2, 3), phases_errors(3, 4), &
            phases_errors(4, 4), phases_errors(5, 4)
    end do
    close(unit)
end subroutine write_phases

!!
!> @brief      Saves list of phases to a file
!!
!! Writes into a file 3 tables of phase-shifts with error bars.
!!
!! The phase-shifts are calculated at the 11 'canonical' energies
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine plot_phases(potential, parameters, covariance, output_name)
    implicit none
    type(nn_model), intent(in) :: potential !< nn scattering model
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    real(dp), intent(in), dimension(:, :) :: covariance !< Covariance matrix
    character(len=*), intent(in) :: output_name

    integer, parameter :: jmax = 4
    real(dp), parameter :: tlab_step = 1._dp
    real(dp), parameter :: tlab_max = 1000._dp
    real(dp), dimension(1:5, 1:jmax) :: phases, phases_errors
    real(dp), allocatable, dimension(:, :, :) :: d_phases

    integer :: unit, j, k
    real(dp) :: t_lab

    open(newunit=unit, file=trim(output_name)//'_pp_phases.dat', status='unknown')
    write(unit, '(17a13)') 'T_lab', '1S0', 'sig_1S0', '1D2', 'sig_1D2', '3P0', 'sig_3P0', '3P1', 'sig_3P1', '3P2', &
        'sig_3P2', 'Ep2', 'sig_Ep2', '3F2', 'sig_3F2', '3F3', 'sig_3F3'
    t_lab = 1._dp
    do
        if(t_lab > tlab_max) exit
        call all_phaseshifts(potential, parameters, t_lab, 'pp', phases, d_phases)
        phases = phases*180/pi
        d_phases = d_phases*180/pi
        do j=1, size(phases, 2)
            do k=1, size(phases, 1)
                phases_errors(k, j) = propagated_error_bar(d_phases(:, k, j), covariance)
            enddo
        enddo
        write(unit, '(17f13.6)') t_lab, phases(1, 1), phases_errors(1, 1), phases(1, 3), phases_errors(1, 3), &
            phases(5, 1), phases_errors(5, 1), phases(2, 2), phases_errors(2, 2), phases(3, 3), phases_errors(3, 3), &
            phases(4, 3), phases_errors(4, 3), phases(5, 3), phases_errors(5, 3), phases(2, 4), phases_errors(2, 4)
        t_lab = t_lab + tlab_step
    end do
    close(unit)

    open(newunit=unit, file=trim(output_name)//'_np_phases.dat', status='unknown')
    write(unit, '(35a13)') 'T_lab', '1S0', 'sig_1S0', '3P0', 'sig_3P0', '1P1', 'sig_1P1', '3P1', 'sig_3P1', '3S1', &
        'sig_3S1', 'Ep1', 'sig_Ep1', '3D1', 'sig_3D1', '1D2', 'sig_1D2', '3D2', 'sig_3D2', '3P2', 'sig_3P2', 'Ep2', &
        'sig_Ep2', '3F2', 'sig_3F2', '1F3', 'sig_1F3', '3F3', 'sig_3F3', '3D3', 'sig_3D3', 'Ep3', 'sig_Ep3', '3G3', &
        'sig_3G3'
    t_lab = 1._dp
    do
        if(t_lab > tlab_max) exit
        call all_phaseshifts(potential, parameters, t_lab, 'np', phases, d_phases)
        phases = phases*180/pi
        d_phases = d_phases*180/pi
        do j=1, size(phases, 2)
            do k=1, size(phases, 1)
                phases_errors(k, j) = propagated_error_bar(d_phases(:, k, j), covariance)
            enddo
        enddo
        write(unit, '(35f13.6)') t_lab, phases(1, 1), phases_errors(1, 1), phases(5, 1), phases_errors(5, 1), &
            phases(1, 2), phases_errors(1, 2), phases(2, 2), phases_errors(2, 2), phases(3, 2), phases_errors(3, 2), &
            phases(4, 2), phases_errors(4, 2), phases(5, 2), phases_errors(5, 2), phases(1, 3), phases_errors(1, 3), &
            phases(2, 3), phases_errors(2, 3), phases(3, 3), phases_errors(3, 3), phases(4, 3), phases_errors(4, 3), &
            phases(5, 3), phases_errors(5, 3), phases(1, 4), phases_errors(1, 4), phases(2, 4), phases_errors(2, 4), &
            phases(3, 4), phases_errors(3, 4), phases(4, 4), phases_errors(4, 4), phases(5, 4), phases_errors(5, 4)
        t_lab = t_lab + tlab_step
    end do
    close(unit)
end subroutine plot_phases

!!
!> @brief      Prints phases in the format of the AV18 paper.
!!
!! Given a nn scattering model and a set of parameters appropriate for that model,
!! calculates pp, nn and np phase shifts and prints them in the format of tables 
!! IV through VII in Phys. Rev. C 51, 38
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine print_phases(parameters, model)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    type(nn_model), intent(in) :: model !< nn scattering model
    integer, parameter :: n_energies = 11
    real(dp), parameter, dimension(1:n_energies) :: energies = [1._dp, 5._dp, 10._dp, 25._dp, &
        50._dp, 100._dp, 150._dp, 200._dp, 250._dp, 300._dp, 350._dp]
    integer :: i
    integer, parameter :: jmax = 4
    real(dp), dimension(1:5, 1:jmax, 1:n_energies) :: phases
    real(dp), allocatable, dimension(:, :, :) :: d_phases

    do i = 1, size(energies)
        call all_phaseshifts(model, parameters, energies(i), 'pp', phases(:, :, i), d_phases)
    enddo
    phases = phases*180/pi
    print*, 'Reproducing pp phases (Table IV in WSS paper)'
    print '(9a9)', 'T_lab', '1S0', '1D2', '3P0', '3P1', '3P2', 'Ep2', '3F2', '3F3'
    do i = 1, size(energies)
        print '(9f9.2)', energies(i), phases(1, 1, i), phases(1, 3, i), phases(5, 1, i), phases(2, 2, i), &
            phases(3, 3, i), phases(4, 3, i), phases(5, 3, i), phases(2, 4, i)
    enddo
    print*, ''

    do i = 1, size(energies)
        call all_phaseshifts(model, parameters, energies(i), 'nn', phases(:, :, i), d_phases)
    enddo
    phases = phases*180/pi
    print*, 'Reproducing nn phases (Table V in WSS paper)'
    print '(9a9)', 'T_lab', '1S0', '1D2', '3P0', '3P1', '3P2', 'Ep2', '3F2', '3F3'
    do i = 1, size(energies)
        print '(9f9.2)', energies(i), phases(1, 1, i), phases(1, 3, i), phases(5, 1, i), phases(2, 2, i), &
            phases(3, 3, i), phases(4, 3, i), phases(5, 3, i), phases(2, 4, i)
    enddo
    print*, ''

    do i = 1, size(energies)
        call all_phaseshifts(model, parameters, energies(i), 'np', phases(:, :, i), d_phases)
    enddo
    phases = phases*180/pi
    print*, 'Reproducing np, T=1 phases (Table VI in WSS paper)'
    print '(9a9)', 'T_lab', '1S0', '1D2', '3P0', '3P1', '3P2', 'Ep2', '3F2', '3F3'
    do i = 1, size(energies)
        print '(9f9.2)', energies(i), phases(1, 1, i), phases(1, 3, i), phases(5, 1, i), phases(2, 2, i), &
            phases(3, 3, i), phases(4, 3, i), phases(5, 3, i), phases(2, 4, i)
    enddo
    print*, ''
    print*, 'Reproducing np, T=0 phases (Table VII in WSS paper)'
    print '(10a9)', 'T_lab', '1P1', '1F3', '3S1', 'Ep1', '3D1', '3D2', '3D3', 'Ep3', '3G3'
    do i = 1, size(energies)
        print '(10f9.2)', energies(i), phases(1, 2, i), phases(1, 4, i), phases(3, 2, i), phases(4, 2, i), &
            phases(5, 2, i), phases(2, 3, i), phases(3, 4, i), phases(4, 4, i), phases(5, 4, i)
    enddo
end subroutine print_phases

subroutine write_potential_setup(potential, parameters, mask, unit)
    implicit none
    type(nn_model), intent(in) :: potential
    real(dp), intent(in), dimension(:) :: parameters
    logical, intent(in), dimension(:) :: mask
    integer, intent(in) :: unit !< Unit where the output is sent to. Either and already opened file or output_unit from iso_fortran_env

    write(unit, *) 'Characteristics of the potential'
    write(unit, *) 'Name:', potential%name
    write(unit, *) 'Maximum integration radius:', potential%r_max
    select case(trim(potential%potential_type))
    case ('local')
        write(unit, *) 'Integration step:', potential%dr
    case ('delta_shell')
        write(unit, *) 'Number of internal delta shells:', potential%n_lambdas
        write(unit, *) 'Distance between internal delta shells:', potential%dr_core
        write(unit, *) 'Distance between external delta shells:', potential%dr_tail
    case default
        stop 'Unrecognized potential type in write_potential_setup'
    end select
    if (potential%relativistic_deuteron) then
        write(unit, *) 'The deuteron is calculated with RELATIVISITC kinematics'
    else
        write(unit, *) 'The deuteron is calculated with NON-RELATIVISITC kinematics'
    endif
    if (potential%full_em_wave) then
        write(unit, *) '1S0 pp phases are calculated against the full EM wave function'
    else
        write(unit, *) '1S0 pp phases are calculated against the Coulomb wave function'
    endif
    write(unit, *) 'Displaying potential parameters'
    write(unit, *) '* indicates that a parameter is kept fixed during the optimization'
    call potential%display_subroutine(parameters, mask, unit)
    write(unit, *) ''
end subroutine write_potential_setup

subroutine setup_from_namelist(namelist_file, potential, parameters, mask, database_file, &
    save_results, output_name)
    implicit none
    character(len=*), intent(in) :: namelist_file
    type(nn_model), intent(out) :: potential
    real(dp), intent(out), allocatable, dimension(:) :: parameters
    logical, intent(out), allocatable, dimension(:) :: mask
    character(len=*), intent(out) :: database_file
    logical, intent(out) :: save_results
    character(len=*), intent(out) :: output_name

    character(len=1024) :: name
    real(dp) :: r_max, delta_r, dr_core, dr_tail
    integer :: n_lambdas
    logical :: relativistic

    logical :: file_exists
    integer :: unit, ierror

    namelist /data_base/ database_file
    namelist /nn_potential/ name
    namelist /local_integration/ r_max, delta_r
    namelist /delta_shell_integration/ r_max, n_lambdas, dr_core, dr_tail
    namelist /deuteron/ relativistic
    namelist /potential_parameters/ parameters
    namelist /adjust_parameter/ mask
    namelist /output/ save_results, output_name


    database_file = 'database/granada_database.dat'
    name = 'AV18'
    save_results = .true.
    output_name = 'results'

    inquire(file=trim(namelist_file), exist = file_exists)
    if (file_exists) then
        print*, 'Setting up potential as specified in: ', trim(namelist_file)
        print*, ''
        open(newunit = unit, file = namelist_file )

        read(unit, nml = data_base, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading data_base namelist'
        endif
        rewind(unit)

        read(unit, nml = output, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading output namelist'
        endif
        rewind(unit)

        read(unit, nml = nn_potential, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading nn_potential namelist'
        endif
        rewind(unit)

        select case(trim(name))
            ! These subroutines setup the default versions and parameters,
            ! including allocating the correct size for the parameters array,
            ! those are later replaced by whatever is read in the namelist file
        case ('AV18')
            call set_av18_potential(potential, parameters)
            r_max = 12.5_dp
            delta_r = 1/128._dp
            relativistic = .false.
            n_lambdas = 0
            dr_core = 0.0_dp
            dr_tail = 0.0_dp
        case ('ds_ope30')
            call set_ds_potential(name, potential, parameters)
            name = 'ds_ope30'
            r_max = 13.0_dp
            n_lambdas = 5
            dr_core = 0.6_dp
            dr_tail = 0.5_dp
            relativistic = .true.
            delta_r = 0.0_dp
        case ('ds_ope30_fff')
            call set_ds_potential(name, potential, parameters)
            name = 'ds_ope30_fff'
            r_max = 13.0_dp
            n_lambdas = 5
            dr_core = 0.6_dp
            dr_tail = 0.5_dp
            relativistic = .true.
            delta_r = 0.0_dp
        case default
            print*, 'Unrecognized name ', trim(name), ' in nn_potential namelist'
            print*, 'stopping the program'
            stop 
        end select

        select case(trim(potential%potential_type))
        case ('local')
            read(unit, nml = local_integration, iostat = ierror)
            if(ierror /= 0) then
                stop 'Error reading local_integration namelist'
            endif
        case ('delta_shell')
            read(unit, nml = delta_shell_integration, iostat = ierror)
            if(ierror /= 0) then
                stop 'Error reading delta_shell_integration namelist'
            endif
        case default
            print*, 'Unrecognized type ', trim(potential%potential_type),' in potential namelist'
            print*, 'stopping the program'
            stop 
        end select
        rewind(unit)

        read(unit, nml = deuteron, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading deuteron namelist'
        endif
        rewind(unit)
        
        read(unit, nml = potential_parameters, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading potential_parameters namelist'
        endif
        rewind(unit)
        
        allocate(mask(1:size(parameters)))
        where (parameters == 0)
            mask = .false.
        else where
            mask = .true.
        end where
        
        read(unit, nml = adjust_parameter, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading adjust_parameter namelist'
        endif
        close(unit)
    else
        print*, 'Namelist file, ', trim(namelist_file)
        print*, 'does not exist. Ending program'
        stop
    endif

    ! Updating values if those where present in the namelist file
    potential%r_max = r_max
    potential%dr = delta_r

    potential%n_lambdas = n_lambdas
    potential%dr_core = dr_core
    potential%dr_tail = dr_tail

    potential%relativistic_deuteron = relativistic
  
end subroutine setup_from_namelist


subroutine write_optimization_results(model, initial_parameters, parameters, mask, chi2, n_points, &
        covariance, output_name)
    implicit none
    type(nn_model), intent(in) :: model
    real(dp), intent(in), dimension(:) :: initial_parameters
    real(dp), intent(in), dimension(:) :: parameters
    logical, intent(in), dimension(:) :: mask
    real(dp), intent(in) :: chi2
    integer, intent(in) :: n_points
    real(dp), intent(in), dimension(:, :) :: covariance
    character(len=*), intent(in) :: output_name

    character(len=31), parameter :: format = '(1x,a,f15.8,2x,a,i5,2x,a,f13.8)'
    real(dp), parameter :: r_min = 0.0_dp, r_max = 2.0_dp, r_step = 0.0078125_dp

    integer :: unit

    open(newunit=unit, file=trim(output_name)//'_parameters.txt')
        write(unit, *) 'Potential setup and initial parameters:'
        call write_potential_setup(model, initial_parameters, mask, unit)
        write(unit, *)
        write(unit, *) 'Final paramters:'
        call model%display_subroutine(parameters, mask, unit, covariance)
        write(unit, format) 'chi^2:', chi2, 'N_data:', n_points, 'chi^2/N_data:', chi2/n_points

    if(model%name == 'AV18') then
        call plot_potential_components(model, parameters, covariance, r_min, r_max, r_step, trim(output_name)//'_plots.dat')
        write(unit, *) 'AV18 components plots saved in: ', trim(output_name)//'_plots.dat'
        call write_marias_format(trim(output_name)//'.in', parameters)
        write(unit, *) 'Parameters in Marias input format saved in: ', trim(output_name)//'.in'
    endif
    if(model%potential_type == 'local') then
        call plot_potential_partial_waves(model, parameters, covariance, r_min, r_max, r_step, output_name)
        write(unit, *) 'Potential in partial waves plots saved in: ', trim(output_name)//'_pp_v_partial_wave.dat and ', &
        trim(output_name)//'_np_v_partial_wave.dat'
    endif
    call write_phases(model, parameters, covariance, trim(output_name)//'_phases.txt')
    write(unit, *) 'Phases listed in: ', trim(output_name)//'_phases.txt'
    call plot_phases(model, parameters, covariance, output_name)
    write(unit, *) 'Phaseshift plots in: ', trim(output_name)//'_pp_phases.dat and ', &
        trim(output_name)//'_np_phases.dat'
    close(unit)

end subroutine write_optimization_results

subroutine write_chiral_kernals(r, file_name)
    implicit none
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: file_name

    real(dp) :: u, u_max
    integer :: unit

    character(len=31), parameter :: format = '(f15.8,38e19.8e3)'

    u = 0.1_dp
    u_max = 20.0_dp

    open(newunit=unit, file=file_name)
    write(unit, *) 'mu', 'vf_1', 'vf_2', 'vf_3', 'vf_4', 'vf_5', 'vf_6', 'vf_7', 'vf_8', 'vf_9'

    do
        if (u > u_max) exit
        u = u + 0.1_dp
    end do
    close(unit)

end subroutine write_chiral_kernals

subroutine plot_potential_components(potential, parameters, covariance, r_min, r_max, r_step, file_name)
    implicit none
    type(nn_model), intent(in) :: potential
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in), dimension(:, :) :: covariance
    real(dp), intent(in) :: r_min
    real(dp), intent(in) :: r_max
    real(dp), intent(in) :: r_step
    character(len=*), intent(in) :: file_name

    real(dp) :: r
    real(dp), allocatable, dimension(:) ::  v
    real(dp), allocatable, dimension(:, :) :: dv
    real(dp), allocatable, dimension(:) :: v_error
    integer :: unit, i

    character(len=31), parameter :: format = '(f15.8,38e19.8e3)'

    r = r_min
    allocate(v(1:potential%n_components))
    allocate(v_error, mold=v)
    open(newunit=unit, file=trim(file_name))
    write(unit, '(a15, 38a19)') 'radius', 'v_c', 'sig_v_c', 'v_tau', 'sig_v_tau', 'v_sigma', 'sig_v_sigma', &
        'v_sigma_tau', 'sig_v_sigma_tau', 'v_t', 'sig_v_t', 'v_t_tau', 'sig_v_t_tau', 'v_ls', 'sig_v_ls', 'v_ls_tau', &
        'sig_v_ls_tau', 'v_l2', 'sig_v_l2', 'v_l2_tau', 'sig_v_l2_tau', 'v_l2_sigma', 'sig_v_l2_sigma', &
        'v_l2_sigma_tau', 'sig_v_l2_sigma_tau', 'v_ls2', 'sig_v_ls2', 'v_ls2_tau', 'sig_v_ls2_tau', 'v_T', 'sig_v_T', &
        'v_sigma_T', 'sig_v_sigma_T', 'v_t_T', 'sig_v_t_T', 'v_ls_T', 'sig_v_ls_T', 'v_tau_z', 'sig_v_tau_z'
    do
        if(r > r_max) exit
        call potential%potential_components(parameters, r, v, dv)
        do i = 1, size(v_error)
            v_error(i) = propagated_error_bar(dv(i, :), covariance)
        enddo
        write(unit, format) r, (v(i), v_error(i), i = 1, size(v))
        r = r + r_step
    enddo
    close(unit)
    
end subroutine plot_potential_components

subroutine plot_potential_partial_waves(potential, parameters, covariance, r_min, r_max, r_step, output_name)
    implicit none
    type(nn_model), intent(in) :: potential
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in), dimension(:, :) :: covariance
    real(dp), intent(in) :: r_min
    real(dp), intent(in) :: r_max
    real(dp), intent(in) :: r_step
    character(len=*), intent(in) :: output_name

    integer, parameter :: n_waves = 5
    integer, parameter :: j_max = 4

    real(dp) :: r
    real(dp), dimension(1:n_waves, 1:j_max) :: v_pw
    real(dp), dimension(1:n_waves, 1:j_max) :: v_pw_errors
    real(dp), allocatable, dimension(:, :, :) :: dv_pw
    integer :: unit, i, j

    character(len=2) :: channel

    channel = 'pp'
    r = r_min
    open(newunit=unit, file=trim(output_name)//'_pp_v_partial_wave.dat', status='unknown')
    write(unit, '(17a13)') 'radius', '1S0', 'sig_1S0', '1D2', 'sig_1D2', '3P0', 'sig_3P0', '3P1', 'sig_3P1', '3P2', &
        'sig_3P2', 'Ep2', 'sig_Ep2', '3F2', 'sig_3F2', '3F3', 'sig_3F3'
    do
        if(r > r_max) exit
        call potential%potential(parameters, r, channel, v_pw, dv_pw)
        do i = 1, size(v_pw_errors, 2)
            do j=1, size(v_pw_errors, 1)
                v_pw_errors(j, i) = propagated_error_bar(dv_pw(j, i, :), covariance)
            enddo
        enddo
        write(unit, '(17f13.6)') r, v_pw(1, 1), v_pw_errors(1, 1), v_pw(1, 3), v_pw_errors(1, 3), &
            v_pw(5, 1), v_pw_errors(5, 1), v_pw(2, 2), v_pw_errors(2, 2), v_pw(3, 3), v_pw_errors(3, 3), &
            v_pw(4, 3), v_pw_errors(4, 3), v_pw(5, 3), v_pw_errors(5, 3), v_pw(2, 4), v_pw_errors(2, 4)
        r = r + r_step
    enddo
    close(unit)

    channel = 'np'
    r = r_min
    open(newunit=unit, file=trim(output_name)//'_np_v_partial_wave.dat', status='unknown')
    write(unit, '(35a13)') 'radius', '1S0', 'sig_1S0', '3P0', 'sig_3P0', '1P1', 'sig_1P1', '3P1', 'sig_3P1', '3S1', &
        'sig_3S1', 'Ep1', 'sig_Ep1', '3D1', 'sig_3D1', '1D2', 'sig_1D2', '3D2', 'sig_3D2', '3P2', 'sig_3P2', 'Ep2', &
        'sig_Ep2', '3F2', 'sig_3F2', '1F3', 'sig_1F3', '3F3', 'sig_3F3', '3D3', 'sig_3D3', 'Ep3', 'sig_Ep3', '3G3', &
        'sig_3G3'
    do
        if(r > r_max) exit
        call potential%potential(parameters, r, channel, v_pw, dv_pw)
        do i = 1, size(v_pw_errors, 2)
            do j=1, size(v_pw_errors, 1)
                v_pw_errors(j, i) = propagated_error_bar(dv_pw(j, i, :), covariance)
            enddo
        enddo
        write(unit, '(35f13.6)') r, v_pw(1, 1), v_pw_errors(1, 1), v_pw(5, 1), v_pw_errors(5, 1), &
            v_pw(1, 2), v_pw_errors(1, 2), v_pw(2, 2), v_pw_errors(2, 2), v_pw(3, 2), v_pw_errors(3, 2), &
            v_pw(4, 2), v_pw_errors(4, 2), v_pw(5, 2), v_pw_errors(5, 2), v_pw(1, 3), v_pw_errors(1, 3), &
            v_pw(2, 3), v_pw_errors(2, 3), v_pw(3, 3), v_pw_errors(3, 3), v_pw(4, 3), v_pw_errors(4, 3), &
            v_pw(5, 3), v_pw_errors(5, 3), v_pw(1, 4), v_pw_errors(1, 4), v_pw(2, 4), v_pw_errors(2, 4), &
            v_pw(3, 4), v_pw_errors(3, 4), v_pw(4, 4), v_pw_errors(4, 4), v_pw(5, 4), v_pw_errors(5, 4)
        r = r + r_step
    enddo
    close(unit)

end subroutine plot_potential_partial_waves

real(dp) function propagated_error_bar(derivatives, covariance) result(error)
    implicit none
    real(dp), intent(in), dimension(:) :: derivatives
    real(dp), intent(in), dimension(:, :) :: covariance

    integer :: i, j

    error = 0._dp

    do i = 1, size(derivatives)
        do j=1, size(derivatives)
            error = error + derivatives(i)*derivatives(j)*covariance(j,i)
        enddo
    enddo

    error = sqrt(abs(error))
    
end function propagated_error_bar


end module read_write
