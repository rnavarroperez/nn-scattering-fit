module observables
use nn_phaseshifts, only: all_phaseshifts, momentum_cm, nn_potential
use amplitudes, only: saclay_amplitudes
use precisions
use constants
implicit none
public observable, f_observable, df_observable!, just_phases
private


! Valid observables
! character(len=4), dimension(1:n_obs), parameter :: &
!        obs_types = ['DSG ','DT  ','AYY ','D   ','P   ','AZZ ','R   '&
!          ,'RT  ','RPT ','AT  ','D0SK','NSKN','NSSN','NNKK','A   '&
!          ,'AXX ','CKP ','RP  ','MSSN','MSKN','AZX ','AP  ','DTRT'&
!          ,'SGT ','SGTT','SGTL']

contains

!!
!> @brief Calculates a NN scattering observable for a given laboratory energy
!! in MeV and scattering angle in d_egrees.
!!
!! The observable to be calculated is d_etermined by the type integer
!! and corresponds to the list of observable labels in the strObs
!! array.
!!
!! To avoid recalculating phase-shifts (the more time consuming part
!! of the calculation) the phases are only calculated if t and tprev
!! are different. At the end of the subroutine tprev is upd_ated to
!! the value of t
!!
!! The d_erivative of the observable with respect of the parameters is stored
!! on the ap array
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine observable(model, params, type, t_lab, angle, reaction, r_max, dr, obs, d_obs)
    implicit none
    procedure(nn_potential) :: model
    real(dp), intent(in) :: params(:)
    character(len=*), intent(in) :: type !< ind_ex to indicate the type of observable
    real(dp), intent(in) :: t_lab !< laboratory energy
    real(dp), intent(in) :: angle !< scattering angle in d_egrees
    character(len=*), intent(in) :: reaction !< reaction channel
    real(dp), intent(in) :: r_max !< 
    real(dp), intent(in) :: dr !< integration radius and step in fm
    real(dp), intent(out) :: obs !< NN scattering observable
    real(dp), allocatable, intent(out) :: d_obs(:) !< derivative of the NN scattering observble
    real(dp), save :: pre_t_lab = -1._dp
    real(dp), save, allocatable :: pre_parameters(:)
    real(dp), save :: k_cm
    real(dp), save, allocatable :: phases(:,:)
    real(dp), save, allocatable :: d_phases(:,:,:)
    complex(dp) :: a, b, c, d, e
    complex(dp), allocatable :: d_a(:), d_b(:), d_c(:), d_d(:), d_e(:)
    integer :: n_parameters 
    real(dp) :: theta, sg, num, denom 
    real(dp), allocatable :: d_sg(:), d_num(:), d_denom(:) 
    integer, parameter :: j_max = 20

    ! Set number of parameters
    n_parameters = size(params)

    ! allocate all arrays
    allocate(d_obs(1:n_parameters))
    if(.not. allocated(phases)) allocate(phases(1:5, 1:j_max))
    if (.not. allocated(pre_parameters)) then
        allocate(pre_parameters, mold=params)
        pre_parameters = 0._dp
    end if
    allocate(d_sg(1:n_parameters))
    allocate(d_num(1:n_parameters))
    allocate(d_denom(1:n_parameters))

    if(t_lab /= pre_t_lab .or. .not. all(params == pre_parameters)) then
        call all_phaseshifts(model, params, t_lab, reaction, r_max, dr, phases, d_phases)
        k_cm = momentum_cm(t_lab, reaction)
        pre_t_lab = t_lab
        pre_parameters = params
    end if
    theta = angle*pi/180.0_dp ! angle in d_egrees to radians
    call saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, d_c, d_d, d_e)

    ! Initialize values for calculation observable
    obs = 0.0_dp
    d_obs = 0.0_dp
    sg = (abs(a)**2 + abs(b)**2 + abs(c)**2 + abs(d)**2 + abs(e)**2)*0.5_dp
    d_sg = real(a)*real(d_a) + aimag(a)*aimag(d_a) &
         + real(b)*real(d_b) + aimag(b)*aimag(d_b) &
         + real(c)*real(d_c) + aimag(c)*aimag(d_c) &
         + real(d)*real(d_d) + aimag(d)*aimag(d_d) &
         + real(e)*real(d_e) + aimag(e)*aimag(d_e)

    ! switch case for all observable
    select case (trim(type))
    case ('dsg')
        obs = sg*10.0_dp
        d_obs = d_sg*10.0_dp
    case ('dt')
        obs = 0.5_dp*(abs(a)**2 - abs(b)**2 + abs(c)**2 - abs(d)**2 + abs(e)**2)/sg
        d_num = real(a)*real(d_a) + aimag(a)*aimag(d_a) &
                  -real(b)*real(d_b) - aimag(b)*aimag(d_b) &
                  +real(c)*real(d_c) + aimag(c)*aimag(d_c) &
                  -real(d)*real(d_d) - aimag(d)*aimag(d_d) &
                  +real(e)*real(d_e) + aimag(e)*aimag(d_e)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('ayy')
       obs = 0.5_dp*(abs(a)**2-abs(b)**2-abs(c)**2+abs(d)**2+abs(e)**2)/sg
        d_num = real(a)*real(d_a) + aimag(a)*aimag(d_a) &
                     - real(b)*real(d_b) - aimag(b)*aimag(d_b) &
                     - real(c)*real(d_c) - aimag(c)*aimag(d_c) &
                     + real(d)*real(d_d) + aimag(d)*aimag(d_d) &
                     + real(e)*real(d_e) + aimag(e)*aimag(d_e)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('d')
        obs = 0.5_dp*(abs(a)**2+abs(b)**2-abs(c)**2-abs(d)**2+abs(e)**2)/sg
        d_num = real(a)*real(d_a) + aimag(a)*aimag(d_a) &
               +real(b)*real(d_b) + aimag(b)*aimag(d_b) &
               -real(c)*real(d_c) - aimag(c)*aimag(d_c) &
               -real(d)*real(d_d) - aimag(d)*aimag(d_d) &
               +real(e)*real(d_e) + aimag(e)*aimag(d_e)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('p')
        obs = real(conjg(a)*e)/sg
        d_num = real(conjg(d_a)*e + conjg(a)*d_e)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('azz')
        obs = (-real(conjg(a)*d)*cos(theta) + real(conjg(b)*c) + aimag(conjg(d)*e)*sin(theta))/sg
        d_num = - real(conjg(d_a)*d + conjg(a)*d_d)*cos(theta) &
                + real(conjg(d_b)*c + conjg(b)*d_c) &
                +aimag(conjg(d_d)*e + conjg(d)*d_e)*sin(theta)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('r')
        obs = (cos(0.5_dp*theta)*( real(conjg(a)*b + conjg(c)*d)) &
              -sin(0.5_dp*theta)*(aimag(conjg(b)*e)))/sg
        d_num = cos(0.5_dp*theta)*( real(conjg(d_a)*b + conjg(a)*d_b &
                                        +conjg(d_c)*d + conjg(c)*d_d)) &
               -sin(0.5_dp*theta)*(aimag(conjg(d_b)*e + conjg(b)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('rt')
        obs = (-cos(0.5_dp*(pi + theta))*( real(conjg(a)*c)) &
               -cos(0.5_dp*(pi - theta))*( real(conjg(b)*d)) &
               +sin(0.5_dp*(pi + theta))*(aimag(conjg(c)*e)))/sg
        d_num = -cos(0.5_dp*(pi + theta))*( real(conjg(d_a)*c + conjg(a)*d_c)) &
                -cos(0.5_dp*(pi - theta))*( real(conjg(d_b)*d + conjg(b)*d_d)) &
                +sin(0.5_dp*(pi + theta))*(aimag(conjg(d_c)*e + conjg(c)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('rpt')
        obs = ( sin(0.5_dp*(pi + theta))*( real(conjg(a)*c)) &
               +sin(0.5_dp*(pi - theta))*( real(conjg(b)*d)) &
               +cos(0.5_dp*(pi + theta))*(aimag(conjg(c)*e)))/sg
        d_num = sin(0.5_dp*(pi + theta))*( real(conjg(d_a)*c + conjg(a)*d_c)) &
               +sin(0.5_dp*(pi - theta))*( real(conjg(d_b)*d + conjg(b)*d_d)) &
               +cos(0.5_dp*(pi + theta))*(aimag(conjg(d_c)*e + conjg(c)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('at')
        obs = -( sin(0.5_dp*(pi + theta))*( real(conjg(a)*c)) &
                -sin(0.5_dp*(pi - theta))*( real(conjg(b)*d))&
                +cos(0.5_dp*(pi + theta))*(aimag(conjg(c)*e)))/sg
        d_num = -sin(0.5_dp*(pi + theta))*( real(conjg(d_a)*c + conjg(a)*d_c))&
                +sin(0.5_dp*(pi - theta))*( real(conjg(d_b)*d + conjg(b)*d_d))&
                -cos(0.5_dp*(pi + theta))*(aimag(conjg(d_c)*e + conjg(c)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('d0sk')
        obs = ( sin(0.5_dp*(pi + theta))*( real(conjg(a)*b)) &
               -sin(0.5_dp*(pi - theta))*( real(conjg(c)*d)) &
               +cos(0.5_dp*(pi + theta))*(aimag(conjg(b)*e)))/sg
        d_num = sin(0.5_dp*(pi + theta))*( real(conjg(d_a)*b + conjg(a)*d_b)) &
               -sin(0.5_dp*(pi - theta))*( real(conjg(d_c)*d + conjg(c)*d_d)) &
               +cos(0.5_dp*(pi + theta))*(aimag(conjg(d_b)*e + conjg(b)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('nskn')
        obs = (sin(0.5_dp*(pi + theta))*( real(conjg(c)*e)) &
              -cos(0.5_dp*(pi + theta))*(aimag(conjg(a)*c)) &
              +cos(0.5_dp*(pi - theta))*(aimag(conjg(b)*d)))/sg
        d_num = sin(0.5_dp*(pi + theta))*( real(conjg(d_c)*e + conjg(c)*d_e)) &
               -cos(0.5_dp*(pi + theta))*(aimag(conjg(d_a)*c + conjg(a)*d_c)) &
               +cos(0.5_dp*(pi - theta))*(aimag(conjg(d_b)*d + conjg(b)*d_d))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('nssn')
        obs = (-cos(0.5_dp*(pi + theta))*( real(conjg(c)*e)) &
               -sin(0.5_dp*(pi + theta))*(aimag(conjg(a)*c)) &
               -sin(0.5_dp*(pi - theta))*(aimag(conjg(b)*d)))/sg
        d_num = -cos(0.5_dp*(pi + theta))*( real(conjg(d_c)*e + conjg(c)*d_e)) &
                -sin(0.5_dp*(pi + theta))*(aimag(conjg(d_a)*c + conjg(a)*d_c)) &
                -sin(0.5_dp*(pi - theta))*(aimag(conjg(d_b)*d + conjg(b)*d_d))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('nnkk')
        obs = (-cos(theta)*(real(conjg(d)*e)) - sin(theta)*(aimag(conjg(a)*d)))/sg
        d_num = -cos(theta)*( real(conjg(d_d)*e + conjg(d)*d_e)) &
                -sin(theta)*(aimag(conjg(d_a)*d + conjg(a)*d_d))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('a')
        obs = (-sin(0.5_dp*theta)*( real(conjg(a)*b + conjg(c)*d)) &
               -cos(0.5_dp*theta)*(aimag(conjg(b)*e)))/sg
        d_num = -sin(0.5_dp*theta)*( real(conjg(d_a)*b + conjg(a)*d_b &
                                         +conjg(d_c)*d + conjg(c)*d_d)) &
                -cos(0.5_dp*theta)*(aimag(conjg(d_b)*e + conjg(b)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('axx')
        obs = (cos(theta)* real(conjg(a)*d) + real(conjg(b)*c) &
              -sin(theta)*aimag(conjg(d)*e))/sg
        d_num = cos(theta)*real(conjg(d_a)*d + conjg(a)*d_d) &
                          +real(conjg(d_b)*c + conjg(b)*d_c) &
               -sin(theta)*aimag(conjg(d_d)*e + conjg(d)*d_e)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('ckp')
        obs = aimag(conjg(d)*e)/sg
        d_num = aimag(conjg(d_d)*e + conjg(d)*d_e)
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('rp')
        obs = (sin(0.5_dp*theta)*( real(conjg(a)*b-conjg(c)*d)) &
              +cos(0.5_dp*theta)*(aimag(conjg(b)*e)))/sg
        d_num = sin(0.5_dp*theta)*( real(conjg(d_a)*b + conjg(a)*d_b &
                                        -conjg(d_c)*d - conjg(c)*d_d)) &
               +cos(0.5_dp*theta)*(aimag(conjg(d_b)*e+conjg(b)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('mssn')
        obs = (cos(0.5_dp*theta)*( real(conjg(b)*e)) &
              +sin(0.5_dp*theta)*(aimag(conjg(a)*b-conjg(c)*d)))/sg
        d_num = cos(0.5_dp*theta)*( real(conjg(d_b)*e + conjg(b)*d_e))&
               +sin(0.5_dp*theta)*(aimag(conjg(d_a)*b + conjg(a)*d_b &
                                        -conjg(d_c)*d - conjg(c)*d_d))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('mskn')
        obs = (-sin(0.5_dp*theta)*(real(conjg(b)*e)) &
               +cos(0.5_dp*theta)*(aimag(conjg(a)*b-conjg(c)*d)))/sg
        d_num= -sin(0.5_dp*theta)*( real(conjg(d_b)*e + conjg(b)*d_e)) &
               +cos(0.5_dp*theta)*(aimag(conjg(d_a)*b + conjg(a)*d_b &
                                        -conjg(d_c)*d - conjg(c)*d_d))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('azx')
        obs = (sin(theta)*( real(conjg(a)*d)) &
              +cos(theta)*(aimag(conjg(d)*e)))/sg
        d_num = sin(theta)*( real(conjg(d_a)*d + conjg(a)*d_d)) &
               +cos(theta)*(aimag(conjg(d_d)*e + conjg(d)*d_e))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('ap')
        obs = (-sin(0.5_dp*theta)*(aimag(conjg(b)*e)) &
               +cos(0.5_dp*theta)*( real(conjg(a)*b - conjg(c)*d)))/sg
        d_num = -sin(0.5_dp*theta)*(aimag(conjg(d_b)*e + conjg(b)*d_e)) &
                +cos(0.5_dp*theta)*( real(conjg(d_a)*b + conjg(a)*d_b &
                                         -conjg(d_c)*d - conjg(c)*d_d))
        num = obs*sg
        d_obs = (d_num*sg - num*d_sg)/sg**2
    case ('dtrt')
        num = 0.5_dp*(abs(a)**2 - abs(b)**2 + abs(c)**2 - abs(d)**2 + abs(e)**2)
        denom = -cos(0.5_dp*(pi + theta))*( real(conjg(a)*c)) &
                -cos(0.5_dp*(pi - theta))*( real(conjg(b)*d)) &
                +sin(0.5_dp*(pi + theta))*(aimag(conjg(c)*e))
        d_num = real(a)*real(d_a) + aimag(a)*aimag(d_a) &
               -real(b)*real(d_b) - aimag(b)*aimag(d_b) &
               +real(c)*real(d_c) + aimag(c)*aimag(d_c) &
               -real(d)*real(d_d) - aimag(d)*aimag(d_d) &
               +real(e)*real(d_e) + aimag(e)*aimag(d_e)
        d_denom = -cos(0.5_dp*(pi + theta))*( real(conjg(d_a)*c + conjg(a)*d_c)) &
                  -cos(0.5_dp*(pi - theta))*( real(conjg(d_b)*d + conjg(b)*d_d)) &
                  +sin(0.5_dp*(pi + theta))*(aimag(conjg(d_c)*e + conjg(c)*d_e))
        obs = num/denom
        d_obs = (d_num*denom - num*d_denom)/denom**2
    case ('sgt')
        obs = 20*pi*(aimag(a + b))/k_cm
        d_obs = 20*pi*(aimag(d_a + d_b))/k_cm
    case ('sgtt')
        obs = -40*pi*(aimag(a - b))/k_cm
        d_obs = -40*pi*(aimag(d_a - d_b))/k_cm
    case ('sgtl')
        obs = -40*pi*(aimag(c - d))/(k_cm)
        d_obs = -40*pi*(aimag(d_c - d_d))/k_cm
    case default
        stop 'INVALID OBSERVABLE in observable subroutine:'
    end select
end subroutine observable

! The functions below were written to test the derivatives of observable against numerical calculations
! and require the av18 module. All analytic calculations of the derivatives match the numerical ones.
! since the observables module should not depend on a specific module for the NN interaction the functions
! are commented and left here for reference.

! !!
! !! @brief      wrapper function for observable
! !!
! !! This wrapper function is used to test the derivatives of the observable subroutine.
! !! The generic data of type context is used to receive all the arguments necessary to call
! !! observable. The same data of type context is used to receive which parameter will
! !! be varied by the dfridr subroutine and which type of observable will be calculated.
! !!
! !! @returns    A NN scattering observable at an specific lab frame energy and scattering angle
! !!
! !! @author     Rodrigo Navarro Perez
! !!
! real(dp) function f_observable(x, data) result(r)
!     use num_recipes, only : context
!     use av18, only : av18_all_partial_waves
!     implicit none
!     real(dp), intent(in) :: x !< parameter that will be varied by the dfridr subroutine
!     type(context), intent(in) :: data !< data structure with all the arguments for av18_operator

!     real(dp), allocatable :: ap(:)
!     real(dp) :: t_lab, theta, r_max, dr
!     integer :: i_target, i_parameter
!     character(len=2) :: reaction
!     character(len=4) :: type
!     real(dp) :: obs
!     real(dp), allocatable :: d_obs(:)
!     integer, parameter :: n_observables = 26
!     character(len=4), dimension(1:n_observables), parameter :: &
!     obs_types = ['dsg ','dt  ','ayy ','d   ','p   ','azz ','r   ', &
!                  'rt  ','rpt ','at  ','d0sk','nskn','nssn','nnkk','a   ', &
!                  'axx ','ckp ','rp  ','mssn','mskn','azx ','ap  ','dtrt', &
!                  'sgt ','sgtt','sgtl']

!     allocate(ap, source = data%x)
!     t_lab = data%a
!     r_max = data%b
!     dr = data%c
!     theta = data%d
!     reaction = trim(data%string)
!     i_parameter = data%i
!     i_target = data%j

!     type = obs_types(i_target)
!     ap(i_parameter) = x
!     call observable(av18_all_partial_waves, ap, type, t_lab, theta, reaction, r_max, dr, obs, d_obs)
!     r = obs
! end function f_observable

! ! !!
! ! !> @brief      wrapper function for the derivatives of observable
! ! !!
! ! !! This wrapper function is used to test the derivatives of the observable subroutine.
! ! !! The generic data of type context is used to receive all the arguments necessary to call
! ! !! observable. The same data of type context is used to receive which parameter will
! ! !! be varied by the dfridr subroutine and which type of observable will be calculated.
! ! !!
! ! !! @returns    derivatives of an observable at an specific lab energy and partial wave
! ! !!
! ! !! @author     Rodrigo Navarro Perez
! ! !!
! function df_observable(data) result(r)
!     use num_recipes, only : context
!     use av18, only : av18_all_partial_waves
!     implicit none
!     type(context), intent(in) :: data !< data structure with all the arguments for av18_operator
!     real(dp), allocatable :: r(:)
    
!     real(dp), allocatable :: ap(:)
!     real(dp) :: t_lab, theta, r_max, dr
!     integer :: i_target, i_parameter
!     character(len=2) :: reaction
!     character(len=4) :: type
!     real(dp) :: obs
!     real(dp), allocatable :: d_obs(:)
!     integer, parameter :: n_observables = 26
!     character(len=4), dimension(1:n_observables), parameter :: &
!     obs_types = ['dsg ','dt  ','ayy ','d   ','p   ','azz ','r   ', &
!                  'rt  ','rpt ','at  ','d0sk','nskn','nssn','nnkk','a   ', &
!                  'axx ','ckp ','rp  ','mssn','mskn','azx ','ap  ','dtrt', &
!                  'sgt ','sgtt','sgtl']

!     allocate(ap, source = data%x)
!     t_lab = data%a
!     r_max = data%b
!     dr = data%c
!     theta = data%d
!     reaction = trim(data%string)
!     i_parameter = data%i
!     i_target = data%j

!     type = obs_types(i_target)
!     call observable(av18_all_partial_waves, ap, type, t_lab, theta, reaction, r_max, dr, obs, d_obs)
!     r = d_obs
! end function df_observable

end module observables
