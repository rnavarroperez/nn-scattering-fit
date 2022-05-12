!!
!> @brief      Wrapper functions to test derivatives
!!
!! Set of wrapper functions that are used in conjuction with the
!! test_derivatives module to benchmark the analytic calculations
!! of derivatives of dependent variablies with respect of the
!! potential parameters (independent variables)
!!
!! What these function do is to take the arguments of the
!! generic type context and use the contents of that argument
!! to call on the corresponding subroutine to return either
!! the dependent variable or the dependent variables derivatives
!!
!! @author     Rodrigo Navarro Perez
!!
module derivatives

use precisions, only : dp
use num_recipes, only : context
use av18, only : av18_all_partial_waves, set_av18_potential
use delta_shell, only : nn_model
use observables, only : kinematics


implicit none

private

public :: f_scattering_length, df_scattering_length, f_av18, df_av18, f_av18_pw, &
    df_av18_pw, f_all_phaseshifts, df_all_phaseshifts, f_amplitudes, df_amplitudes, f_observable, df_observable, &
    f_deuteron_binding_energy, df_deuteron_binding_energy, f_all_phaseshifts_ds, df_all_phaseshifts_ds, &
    f_observable_ds, df_observable_ds, f_short_chiral_pot, df_short_chiral_pot

contains

real(dp) function f_short_chiral_pot(x, data) result(r)
    use short_range_chiral_potentials, only : short_range_potentials
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for short_range_potentials
    
    real(dp) :: radius
    real(dp), allocatable, dimension(:) :: parameters
    real(dp), allocatable, dimension(:) :: v_short
    real(dp), allocatable, dimension(:, :) :: d_v_short

    radius = data%a
    allocate(parameters, source=data%x)

    parameters(data%i) = x
    call short_range_potentials(radius, parameters, v_short, d_v_short)
    r = v_short(data%j)

end function f_short_chiral_pot

function df_short_chiral_pot(data) result(r)
    use short_range_chiral_potentials, only : short_range_potentials
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for av18_operator
    real(dp), allocatable :: r(:)

    real(dp) :: radius
    real(dp), allocatable, dimension(:) :: parameters
    real(dp), allocatable, dimension(:) :: v_short
    real(dp), allocatable, dimension(:, :) :: d_v_short

    allocate(parameters, source=data%x)
    radius = data%a
    call short_range_potentials(radius, parameters, v_short, d_v_short)

    r = d_v_short(data%j, :)

end function df_short_chiral_pot

!!
!> @brief      wrapper function for dueteron binding energy
!!
!! This wrapper function is used to test the derivatives of the binding_energy subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! binding_energy. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine.
!!
!! @returns    deuteron binding energy
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_deuteron_binding_energy(x, data) result(r)
    use deuteron, only : binding_energy
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for binding_energy

    real(dp), allocatable, dimension(:) :: parameters
    type(nn_model) :: model
    integer :: i_parameter
    real(dp) :: be
    real(dp), allocatable, dimension(:) :: dbe

    allocate(parameters, source = data%x)

    model%potential => av18_all_partial_waves
    model%dr = data%a
    model%r_max = data%b
    model%potential_type = 'local'
    model%relativistic_deuteron = .False.
    i_parameter = data%i

    parameters(i_parameter) = x

    call binding_energy(model, parameters, be, dbe)
    r = be

end function f_deuteron_binding_energy

!!
!> @brief      wrapper function for the derivatives of binding_energy
!!
!! This wrapper function is used to test the derivatives of the binding_energy subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! binding_energy.
!!
!! @returns    derivatives of the dueteron binding energy
!!
!! @author     Rodrigo Navarro Perez
!!
function df_deuteron_binding_energy(data) result(r)
    use deuteron, only : binding_energy
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for binding_energy
    real(dp), allocatable, dimension(:) :: r

    real(dp), allocatable, dimension(:) :: parameters
    type(nn_model) :: model
    real(dp) :: be
    real(dp), allocatable, dimension(:) :: dbe

    allocate(parameters, source = data%x)
    allocate(r, mold=data%x)
    model%potential => av18_all_partial_waves
    model%dr = data%a
    model%r_max = data%b
    model%potential_type = 'local'
    model%relativistic_deuteron = .False.

    call binding_energy(model, parameters, be, dbe)
    r = dbe

end function df_deuteron_binding_energy

!!
!> @brief      wrapper function for scattering_length
!!
!! This wrapper function is used to test the derivatives of the scattering_length subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! scattering_length. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which reaction channel will be used.
!!
!! @returns    \f$ ^1S_0 \f$ scattering lenght
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_scattering_length(x, data) result(r)
    use observables, only : scattering_length
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for scattering_length

    real(dp), allocatable, dimension(:) :: parameters
    type(nn_model) :: model
    integer :: i_target, i_parameter
    character(len=2), parameter, dimension(1:2) :: channels = ['np', 'nn']
    real(dp) :: a_length
    real(dp), allocatable, dimension(:) :: da_length

    allocate(parameters, source = data%x)

    model%potential => av18_all_partial_waves
    model%dr = data%a
    model%r_max = data%b
    model%potential_type = 'local'
    i_parameter = data%i
    i_target = data%j

    parameters(i_parameter) = x

    call scattering_length(model, parameters, channels(i_target), a_length, da_length)
    r = a_length

end function f_scattering_length

!!
!> @brief      wrapper function for the derivatives of scattering_length
!!
!! This wrapper function is used to test the derivatives of the scattering_length subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! scattering_length. The same data of type context is used to receive which reaction channel will be used.
!!
!! @returns    derivatives of the \f$^1S_0\f$ scattering length
!!
!! @author     Rodrigo Navarro Perez
!!
function df_scattering_length(data) result(r)
    use observables, only : scattering_length
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for scattering_length
    real(dp), allocatable, dimension(:) :: r

    real(dp), allocatable, dimension(:) :: parameters
    type(nn_model) :: model
    integer :: i_target
    character(len=2), parameter, dimension(1:2) :: channels = ['np', 'nn']
    real(dp) :: a_length
    real(dp), allocatable, dimension(:) :: da_length

    allocate(parameters, source = data%x)
    allocate(r, mold=data%x)
    model%potential => av18_all_partial_waves
    model%dr = data%a
    model%r_max = data%b
    model%potential_type = 'local'
    i_target = data%j
    call scattering_length(model, parameters, channels(i_target), a_length, da_length)
    r = da_length

end function df_scattering_length

!!
!> @brief      wrapper function for observable
!!
!! This wrapper function is used to test the derivatives of the observable subroutine using the AV18 potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! observable. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which type of observable will be calculated.
!!
!! @returns    A NN scattering observable at an specific lab frame energy and scattering angle
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_observable(x, data) result(r)
    use observables, only : observable
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for observable

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    type(kinematics) :: kinematic
    integer :: i_target, i_parameter
    real(dp) :: obs
    real(dp), allocatable :: d_obs(:)
    integer, parameter :: n_observables = 28
    character(len=4), dimension(1:n_observables), parameter :: &
    obs_types = ['dsg ', 'dt  ', 'ayy ', 'd   ', 'p   ', 'azz ', 'r   ', &
                 'rt  ', 'rpt ', 'at  ', 'd0sk', 'nskn', 'nssn', 'nnkk', 'a   ', &
                 'axx ', 'ckp ', 'rp  ', 'mssn', 'mskn', 'azx ', 'ap  ', 'dtrt', &
                 'sgt ', 'sgtt', 'sgtl', 'asl ', 'dbe ']

    call set_av18_potential(model, ap)
    ap = data%x
    kinematic%t_lab = data%a
    kinematic%angle = data%b
    kinematic%channel = trim(data%string)
    kinematic%em_amplitude = 0._dp
    i_parameter = data%i
    i_target = data%j

    kinematic%type = obs_types(i_target)
    ap(i_parameter) = x
    call observable(kinematic, ap, model, obs, d_obs)
    r = obs
end function f_observable

!!
!> @brief      wrapper function for the derivatives of observable
!!
!! This wrapper function is used to test the derivatives of the observable subroutine using the AV18 potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! observable. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which type of observable will be calculated.
!!
!! @returns    derivatives of an observable at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function df_observable(data) result(r)
    use observables, only : observable
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for observable
    real(dp), allocatable :: r(:)

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    type(kinematics) :: kinematic
    integer :: i_target, i_parameter
    real(dp) :: obs
    real(dp), allocatable :: d_obs(:)
    integer, parameter :: n_observables = 28
    character(len=4), dimension(1:n_observables), parameter :: &
    obs_types = ['dsg ', 'dt  ', 'ayy ', 'd   ', 'p   ', 'azz ', 'r   ', &
                 'rt  ', 'rpt ', 'at  ', 'd0sk', 'nskn', 'nssn', 'nnkk', 'a   ', &
                 'axx ', 'ckp ', 'rp  ', 'mssn', 'mskn', 'azx ', 'ap  ', 'dtrt', &
                 'sgt ', 'sgtt', 'sgtl', 'asl ', 'dbe ']

    call set_av18_potential(model, ap)
    ap = data%x
    kinematic%t_lab = data%a
    kinematic%angle = data%b
    kinematic%channel = trim(data%string)
    kinematic%em_amplitude = 0._dp
    i_parameter = data%i
    i_target = data%j

    kinematic%type = obs_types(i_target)
    call observable(kinematic, ap, model, obs, d_obs)
    r = d_obs
end function df_observable

!!
!> @brief      wrapper function for observable
!!
!! This wrapper function is used to test the derivatives of the observable subroutine using a DS potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! observable. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which type of observable will be calculated.
!!
!! @returns    A NN scattering observable at an specific lab frame energy and scattering angle
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_observable_ds(x, data) result(r)
    use observables, only : observable
    use delta_shell, only : set_ds_potential
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for observable

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    type(kinematics) :: kinematic
    integer :: i_target, i_parameter
    real(dp) :: obs
    real(dp), allocatable :: d_obs(:)
    integer, parameter :: n_observables = 28
    character(len=4), dimension(1:n_observables), parameter :: &
    obs_types = ['dsg ', 'dt  ', 'ayy ', 'd   ', 'p   ', 'azz ', 'r   ', &
                 'rt  ', 'rpt ', 'at  ', 'd0sk', 'nskn', 'nssn', 'nnkk', 'a   ', &
                 'axx ', 'ckp ', 'rp  ', 'mssn', 'mskn', 'azx ', 'ap  ', 'dtrt', &
                 'sgt ', 'sgtt', 'sgtl', 'asl ', 'dbe ']

    call set_ds_potential(data%string, model, ap)
    ap = data%x
    kinematic%t_lab = data%a
    kinematic%angle = data%b
    kinematic%channel = trim(data%string_2)
    kinematic%em_amplitude = 0._dp
    i_parameter = data%i
    i_target = data%j
    kinematic%type = obs_types(i_target)

    ap(i_parameter) = x
    call observable(kinematic, ap, model, obs, d_obs)
    r = obs
end function f_observable_ds

!!
!> @brief      wrapper function for the derivatives of observable
!!
!! This wrapper function is used to test the derivatives of the observable subroutine using a DS potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! observable. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which type of observable will be calculated.
!!
!! @returns    derivatives of an observable at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function df_observable_ds(data) result(r)
    use observables, only : observable
    use delta_shell, only : set_ds_potential
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for observable
    real(dp), allocatable :: r(:)

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    type(kinematics) :: kinematic
    integer :: i_target, i_parameter
    real(dp) :: obs
    real(dp), allocatable :: d_obs(:)
    integer, parameter :: n_observables = 28
    character(len=4), dimension(1:n_observables), parameter :: &
    obs_types = ['dsg ', 'dt  ', 'ayy ', 'd   ', 'p   ', 'azz ', 'r   ', &
                 'rt  ', 'rpt ', 'at  ', 'd0sk', 'nskn', 'nssn', 'nnkk', 'a   ', &
                 'axx ', 'ckp ', 'rp  ', 'mssn', 'mskn', 'azx ', 'ap  ', 'dtrt', &
                 'sgt ', 'sgtt', 'sgtl', 'asl ', 'dbe ']

    call set_ds_potential(data%string, model, ap)
    ap = data%x
    kinematic%t_lab = data%a
    kinematic%angle = data%b
    kinematic%channel = trim(data%string_2)
    kinematic%em_amplitude = 0._dp
    i_parameter = data%i
    i_target = data%j
    kinematic%type = obs_types(i_target)

    call observable(kinematic, ap, model, obs, d_obs)
    r = d_obs
end function df_observable_ds

!!
!> @brief      wrapper function for saclay_amplitudes
!!
!! This wrapper function is used to test the derivatives of the saclay_amplitudes subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! saclay_amplitudes. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    NN phase-shift at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_amplitudes(x, data) result(r)
    use nn_phaseshifts, only : all_phaseshifts, momentum_cm
    use amplitudes, only : saclay_amplitudes
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for saclay_amplitudes

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    real(dp) :: t_lab, theta, k_cm, phases(1:5,1:20)
    real(dp), allocatable :: d_phases(:, :, :)
    integer :: i_target, i_parameter
    character(len=2) :: reaction
    complex(dp) :: a, b, c, d, e
    complex(dp), allocatable, dimension(:) :: d_a, d_b, d_c, d_d, d_e

    allocate(ap, source = data%x)
    model%potential => av18_all_partial_waves
    t_lab = data%a
    model%r_max = data%b
    model%dr = data%c
    model%potential_type = 'local'
    theta = data%d
    reaction = trim(data%string)
    i_parameter = data%i
    i_target = data%j

    ap(i_parameter) = x
    call all_phaseshifts(model, ap, t_lab, reaction, phases, d_phases)
    k_cm = momentum_cm(t_lab, reaction)
    call saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, d_c, d_d, d_e)

    select case (i_target)
    case (1)
        r = real(a)
    case (2)
        r = aimag(a)
    case (3)
        r = real(b)
    case (4)
        r = aimag(b)
    case (5)
        r = real(c)
    case (6)
        r = aimag(c)
    case (7)
        r = real(d)
    case (8)
        r = aimag(d)
    case (9)
        r = real(e)
    case (10)
        r = aimag(e)
    case default
        r = 0
    end select
end function f_amplitudes

!!
!> @brief      wrapper function for the derivatives of saclay_amplitudes
!!
!! This wrapper function is used to test the derivatives of the saclay_amplitudes subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! saclay_amplitudes. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    derivatives of a NN phase-shift at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function df_amplitudes(data) result(r)
    use nn_phaseshifts, only : all_phaseshifts, momentum_cm
    use amplitudes, only : saclay_amplitudes
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for saclay_amplitudes
    real(dp), allocatable :: r(:)

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    real(dp) :: t_lab, theta, k_cm, phases(1:5, 1:20)
    real(dp), allocatable :: d_phases(:, :, :)
    integer :: i_target, i_parameter
    character(len=2) :: reaction
    complex(dp) :: a, b, c, d, e
    complex(dp), allocatable, dimension(:) :: d_a, d_b, d_c, d_d, d_e

    allocate(ap, source = data%x)
    model%potential => av18_all_partial_waves
    t_lab = data%a
    model%r_max = data%b
    model%dr = data%c
    model%potential_type = 'local'
    theta = data%d
    reaction = trim(data%string)
    i_parameter = data%i
    i_target = data%j

    call all_phaseshifts(model, ap, t_lab, reaction, phases, d_phases)
    k_cm = momentum_cm(t_lab, reaction)
    call saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, d_c, d_d, d_e)

    allocate(r, mold = ap)

    select case (i_target)
    case (1)
        r = real(d_a)
    case (2)
        r = aimag(d_a)
    case (3)
        r = real(d_b)
    case (4)
        r = aimag(d_b)
    case (5)
        r = real(d_c)
    case (6)
        r = aimag(d_c)
    case (7)
        r = real(d_d)
    case (8)
        r = aimag(d_d)
    case (9)
        r = real(d_e)
    case (10)
        r = aimag(d_e)
    case default
        r = 0
    end select
end function df_amplitudes


!!
!> @brief      wrapper function for all_phaseshifts
!!
!! This wrapper function is used to test the derivatives of the all_phaseshifts subroutine using the AV18 potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! all_phaseshifts. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    NN phase-shift at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_all_phaseshifts(x, data) result(r)
    use delta_shell, only: nn_model
    use nn_phaseshifts, only : all_phaseshifts
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for all_phaseshifts

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    real(dp) :: t_lab
    real(dp), allocatable :: phases(:, :), d_phases(:, :, :)
    integer :: i_target, i_parameter, ic, ij
    character(len=2) :: reaction

    allocate(ap, source = data%x)
    model%potential => av18_all_partial_waves
    t_lab = data%a
    model%r_max = data%b
    model%dr = data%c
    model%potential_type = 'local'
    reaction = trim(data%string)
    i_parameter = data%i
    i_target = data%j

    ap(i_parameter) = x

    ic = mod(i_target - 1, 5) + 1
    ij = 1 + (i_target - 1)/5

    allocate(phases(1:5, ij))
    call all_phaseshifts(model, ap, t_lab, reaction, phases, d_phases)
    r = phases(ic, ij)
end function f_all_phaseshifts

!!
!> @brief      wrapper function for the derivatives of all_phaseshifts
!!
!! This wrapper function is used to test the derivatives of the all_phaseshifts subroutine using the AV18 potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! all_phaseshifts. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    derivatives of a NN phase-shift at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function df_all_phaseshifts(data) result(r)
    use delta_shell, only : nn_model
    use nn_phaseshifts, only : all_phaseshifts
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for all_phaseshifts
    real(dp), allocatable :: r(:)

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    real(dp) :: t_lab
    real(dp), allocatable :: phases(:, :), d_phases(:, :, :)
    integer :: i_target, i_parameter, ic, ij
    character(len=2) :: reaction

    allocate(ap, source = data%x)
    model%potential => av18_all_partial_waves
    t_lab = data%a
    model%r_max = data%b
    model%dr = data%c
    model%potential_type = 'local'
    reaction = trim(data%string)
    i_parameter = data%i
    i_target = data%j

    ic = mod(i_target - 1, 5) + 1
    ij = 1 + (i_target - 1)/5

    allocate(phases(1:5, ij))
    call all_phaseshifts(model, ap, t_lab, reaction, phases, d_phases)
    r = d_phases(:, ic, ij)
end function df_all_phaseshifts

!!
!> @brief      wrapper function for all_phaseshifts
!!
!! This wrapper function is used to test the derivatives of the all_phaseshifts subroutine using a DS potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! all_phaseshifts. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    NN phase-shift at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_all_phaseshifts_ds(x, data) result(r)
    use delta_shell, only: nn_model
    use nn_phaseshifts, only : all_phaseshifts
    use pion_exchange, only : ope_all_partial_waves
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for all_phaseshifts

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    real(dp) :: t_lab
    real(dp), allocatable :: phases(:, :), d_phases(:, :, :)
    integer :: i_target, i_parameter, ic, ij
    character(len=2) :: reaction

    allocate(ap, source = data%x)
    model%potential => ope_all_partial_waves
    t_lab = data%a
    model%r_max = data%b
    model%dr_core = data%c
    model%dr_tail = data%d
    model%n_lambdas = data%k
    model%potential_type = 'delta_shell'
    reaction = trim(data%string)
    i_parameter = data%i
    i_target = data%j

    ap(i_parameter) = x

    ic = mod(i_target - 1, 5) + 1
    ij = 1 + (i_target - 1)/5

    allocate(phases(1:5, ij))
    call all_phaseshifts(model, ap, t_lab, reaction, phases, d_phases)
    r = phases(ic, ij)
end function f_all_phaseshifts_ds

!!
!> @brief      wrapper function for the derivatives of all_phaseshifts
!!
!! This wrapper function is used to test the derivatives of the all_phaseshifts subroutine using a DS potential.
!! The generic data of type context is used to receive all the arguments necessary to call
!! all_phaseshifts. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    derivatives of a NN phase-shift at an specific lab energy and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function df_all_phaseshifts_ds(data) result(r)
    use delta_shell, only : nn_model
    use nn_phaseshifts, only : all_phaseshifts
    use pion_exchange, only : ope_all_partial_waves
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for all_phaseshifts
    real(dp), allocatable :: r(:)

    real(dp), allocatable :: ap(:)
    type(nn_model) :: model
    real(dp) :: t_lab
    real(dp), allocatable :: phases(:, :), d_phases(:, :, :)
    integer :: i_target, i_parameter, ic, ij
    character(len=2) :: reaction

    allocate(ap, source = data%x)
    model%potential => ope_all_partial_waves
    t_lab = data%a
    model%r_max = data%b
    model%dr_core = data%c
    model%dr_tail = data%d
    model%n_lambdas = data%k
    model%potential_type = 'delta_shell'
    reaction = trim(data%string)
    i_parameter = data%i
    i_target = data%j

    ic = mod(i_target - 1, 5) + 1
    ij = 1 + (i_target - 1)/5

    allocate(phases(1:5, ij))
    call all_phaseshifts(model, ap, t_lab, reaction, phases, d_phases)
    r = d_phases(:, ic, ij)
end function df_all_phaseshifts_ds

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
    use av18, only : av18_operator, n_operators
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for av18_operator

    real(dp), allocatable, dimension(:) :: ap
    real(dp) :: radius
    real(dp) :: v_nn(1:n_operators)
    real(dp), allocatable, dimension(:, :) :: dv_nn
    integer :: i_target, i_parameter
    integer:: n_parameters

    allocate(ap, source=data%x)
    radius = data%a
    i_parameter = data%i
    i_target = data%j

    ap(i_parameter) = x
    n_parameters = size(ap)
    allocate(dv_nn(1:n_parameters, 1:n_operators))
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
    use av18, only : av18_operator, n_operators
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for av18_operator
    real(dp), allocatable :: r(:)

    real(dp), allocatable, dimension(:) :: ap
    real(dp) :: radius
    real(dp) :: v_nn(1:n_operators)
    real(dp), allocatable, dimension(:, :) :: dv_nn
    integer :: i_target, n_parameters

    allocate(ap, source=data%x)
    radius = data%a
    i_target = data%j
    n_parameters = size(ap)
    allocate(dv_nn(1:n_operators, 1:n_parameters))
    allocate(r, mold=data%x)
    call av18_operator(ap, radius, v_nn, dv_nn)

    r = dv_nn(i_target, :)

end function df_av18

!!
!> @brief      wrapper function for av18_all_partial_waves
!!
!! This wrapper function is used to test the derivatives of the av18_all_partial_waves subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! av18_all_partial_waves. The same data of type context is used to receive which parameter will
!! be varied by the dfridr subroutine and which partial wave will be returned.
!!
!! @returns    av18 potential at an specific radius and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function f_av18_pw(x, data) result(r)
    implicit none
    real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
    type(context), intent(in) :: data !< data structure with all the arguments for av18_all_partial_waves

    real(dp), allocatable, dimension(:) :: ap
    real(dp) :: radius
    real(dp), allocatable :: v_pw(:, :)
    real(dp), allocatable :: dv_pw(:, :, :)
    integer :: i_target, i_parameter, n_waves, ic, ij
    character(len=2) :: reaction

    allocate(ap, source=data%x)
    radius = data%a
    i_parameter = data%i
    i_target = data%j
    n_waves = data%k
    reaction = trim(data%string)

    allocate(v_pw(5,n_waves))

    ap(i_parameter) = x

    ic = mod(i_target - 1, 5) + 1
    ij = 1 + (i_target - 1)/5

    call av18_all_partial_waves(ap, radius, reaction, v_pw, dv_pw)

    r = v_pw(ic, ij)

end function f_av18_pw


!!
!> @brief      wrapper function for the derivatives of av18_all_partial_waves
!!
!! This wrapper function is used to test the derivatives of the av18_all_partial_waves subroutine.
!! The generic data of type context is used to receive all the arguments necessary to call
!! av18_all_partial_waves. The same data of type context is used to receive which operator
!! will be returned
!!
!! @returns    the derivatives of the av18 potential at an specific radius and partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function df_av18_pw(data) result(r)
    implicit none
    type(context), intent(in) :: data !< data structure with all the arguments for av18_all_partial_waves
    real(dp), allocatable :: r(:)

    real(dp), allocatable, dimension(:) :: ap
    real(dp) :: radius
    real(dp), allocatable :: v_pw(:, :)
    real(dp), allocatable :: dv_pw(:, :, :)
    integer :: i_target, n_waves, ic, ij
    character(len=2) :: reaction

    allocate(ap, source=data%x)
    radius = data%a
    i_target = data%j
    n_waves = data%k
    reaction = trim(data%string)

    allocate(v_pw(5,n_waves))

    ic = mod(i_target - 1, 5) + 1
    ij = 1 + (i_target - 1)/5

    call av18_all_partial_waves(ap, radius, reaction, v_pw, dv_pw)

    r = dv_pw(:, ic, ij)

end function df_av18_pw

end module derivatives
