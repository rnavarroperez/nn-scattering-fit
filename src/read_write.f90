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
implicit none

private

public :: print_em_amplitudes, print_observables, write_phases, read_montecarlo_parameters, &
    write_montecarlo_phases, print_phases, print_potential_setup, setup_from_namelist

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

subroutine print_potential_setup(potential, parameters, mask)
    implicit none
    type(nn_model), intent(in) :: potential
    real(dp), intent(in), dimension(:) :: parameters
    logical, intent(in), dimension(:) :: mask

    print*, 'Characteristics of the potential being used'
    print*, 'Name:', potential%name
    print*, 'Maximum integration radius:', potential%r_max
    select case(trim(potential%potential_type))
    case ('local')
        print*, 'Integration step:', potential%dr
    case ('delta_shell')
        print*, 'Number of internal delta shells:', potential%n_lambdas
        print*, 'Distance between internal delta shells:', potential%dr_core
        print*, 'Distance between external delta shells:', potential%dr_tail
    case default
        stop 'Unrecognized potential type in print_potential_setup'
    end select
    if (potential%relativistic_deuteron) then
        print*, 'The deuteron will be calculated with RELATIVISITC kinematics'
    else
        print*, 'The deuteron will be calculated with NON-RELATIVISITC kinematics'
    endif
    print*, 'Displaying potential parameters'
    print*, '* indicates that a parameter is kept fixed during the optimization'
    call potential%display_subroutine(parameters, mask)
    print*, ''
end subroutine print_potential_setup

subroutine setup_from_namelist(namelist_file, potential, parameters, mask, database_file)
    implicit none
    character(len=*), intent(in) :: namelist_file
    type(nn_model), intent(out) :: potential
    real(dp), intent(out), allocatable, dimension(:) :: parameters
    logical, intent(out), allocatable, dimension(:) :: mask
    character(len=*), intent(out) :: database_file

    character(len=1024) :: name, type
    real(dp) :: r_max, delta_r, dr_core, dr_tail
    integer :: n_lambdas
    logical :: relativistic

    logical :: file_exists
    integer :: unit, ierror

    namelist /data_base/ database_file
    namelist /nn_potential/ name, type
    namelist /local_integration/ r_max, delta_r
    namelist /delta_shell_integration/ r_max, n_lambdas, dr_core, dr_tail
    namelist /deuteron/ relativistic
    namelist /potential_parameters/ parameters
    namelist /adjust_parameter/ mask


    !setting up default values in namelists
    name = 'AV18'
    type = 'local'
    r_max = 12.5_dp
    delta_r = 1/128._dp
    n_lambdas = 5
    dr_core = 0.6_dp
    dr_tail = 0.5_dp
    relativistic = .false.


    inquire(file=trim(namelist_file), exist = file_exists)
    if (file_exists) then
        print*, 'Setting up potential as specified in: ', trim(namelist_file)
        print*, ''
        open(newunit = unit, file = namelist_file )

        read(unit, nml = data_base, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading data_base namelist'
        endif

        read(unit, nml = nn_potential, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading nn_potential namelist'
        endif

        select case(trim(type))
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
            stop 'Unrecognized type in potential namelist'
        end select

        select case(trim(name))
            ! These subroutines setup the default versions and parameters,
            ! including allocating the correct size for the parameters array,
            ! those are later replaced by whatever is read in the namelist file
        case ('AV18')
            call set_av18_potential(potential, parameters)
        case ('ds_ope30')
            call set_ds_potential(name, potential, parameters)
        case ('ds_ope30_fff')
            call set_ds_potential(name, potential, parameters)
        case default
            stop 'Unrecognized name in potential namelist'
        end select
        
        read(unit, nml = deuteron, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading deuteron namelist'
        endif
        
        read(unit, nml = potential_parameters, iostat = ierror)
        if(ierror /= 0) then
            stop 'Error reading potential_parameters namelist'
        endif
        
        allocate(mask(1:size(parameters)))
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

    potential%name = trim(name)
    potential%potential_type = trim(type)

    potential%r_max = r_max
    potential%dr = delta_r

    potential%n_lambdas = n_lambdas
    potential%dr_core = dr_core
    potential%dr_tail = dr_tail

    potential%relativistic_deuteron = relativistic
  
end subroutine setup_from_namelist

end module read_write
