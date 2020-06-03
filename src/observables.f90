!!
!> @brief      NN scattering observables
!!
!! Module to calculate nn scattering observables with a given model (nn local potential),
!! model parameters and kinematics (lab energy, scattering angle, reaction channel, etc)
!!
!! @author     Rodrigo Navarro Perez
!!
module observables

use nn_phaseshifts, only: all_phaseshifts, momentum_cm, nn_local_model, add_coulomb
use amplitudes, only: saclay_amplitudes
use precisions, only: dp
use constants, only: pi, hbar_c, m_p=>proton_mass, m_n=>neutron_mass

implicit none

public observable, kinematics, scattering_length
private

!!
!> @brief      kinematic variables for observables
!!
!! This can be considered the independent variables, which 
!! change from one observable to the next in a \f$ \chi^2 \f$ calculation/optimization
!!
!! Although the electromagnetic amplitudes are technically not "independent varibles"
!! we include them here since they do not dependent on the fitting parameters and do
!! chage from one observable to the next when calculating a \f$ \chi^2 \f$.
!!
!! @author     Rodrigo Navarro Perez
!!
type :: kinematics
    real(dp) :: t_lab !< Laboratory energy in MeV
    real(dp) :: angle !< scattering angles in degrees
    character(len=2) :: channel !< reaction channel, either 'pp' or 'np'
    character(len=4) :: type !< the type of observable ('dsg', 'sgt', ...)
    complex(dp), dimension(1:5) :: em_amplitude !< electromagnetic amplitude
end type kinematics

contains

!!
!> @brief Calculates a NN observable
!!
!! The type of observable is determined by the kinematic%type argument 
!!
!! This is a wrapper subroutine that base on the type of observable
!! calls the appropriate subroutine. Either the \f$ ^1S_0 \f$ np 
!! scattering, the dueteron binding energy, or one of the 26 observables
!! types that require to calculate all phase-shifts and scattering 
!! amplitudes parameters
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine observable(kinematic, params, model, obs, d_obs)
    implicit none
    type(kinematics), intent(in) :: kinematic !< kinematic variables
    real(dp), intent(in) :: params(:) !< adjustable parameters
    type(nn_local_model), intent(in) :: model !< nn scattering model
    real(dp), intent(out) :: obs !< NN scattering observable
    real(dp), allocatable, intent(out) :: d_obs(:) !< derivative of the NN scattering observble

    select case (trim(kinematic%type))
    case('asl')
        call scattering_length(model, params, kinematic%channel, obs, d_obs)
    case default
        call scattering_obs(kinematic, params, model, obs, d_obs)
    end select
end subroutine observable

!!
!> @brief Calculates a NN scattering observable
!!
!! The type of observable is determined by the kinematic%type argument 
!!
!! To avoid recalculating phase-shifts (the more time consuming part
!! of the calculation) the phases are only calculated if: kinematic%t_lab
!! is different from the previous call, any of the fitting parameters is 
!! different from the previews call, or the kinematic%channel is different 
!! from the previous call.
!!
!! The derivative of the observable with respect of the parameters is stored
!! on the d_obs array
!!
!! @author     Raul L Bernal-Gonzalez
!! @author     Rodrigo Navarro Perez
!!
subroutine scattering_obs(kinematic, params, model, obs, d_obs)
    implicit none
    type(kinematics), intent(in) :: kinematic !< kinematic variables
    real(dp), intent(in) :: params(:) !< adjustable parameters
    type(nn_local_model), intent(in) :: model !< nn scattering model
    real(dp), intent(out) :: obs !< NN scattering observable
    real(dp), allocatable, intent(out) :: d_obs(:) !< derivative of the NN scattering observble
    
    real(dp), save :: pre_t_lab = -1._dp
    real(dp), save, allocatable :: pre_parameters(:)
    real(dp), save :: k_cm
    real(dp), save, allocatable :: phases(:,:)
    real(dp), save, allocatable :: d_phases(:,:,:)
    character(len=2), save :: pre_channel = '  '
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

    if(kinematic%t_lab/=pre_t_lab .or. (.not. all(params==pre_parameters)) .or. &
       kinematic%channel/=pre_channel) then
        call all_phaseshifts(model, params, kinematic%t_lab, kinematic%channel, phases, d_phases)
        k_cm = momentum_cm(kinematic%t_lab, kinematic%channel)
        pre_t_lab = kinematic%t_lab
        pre_parameters = params
        pre_channel = kinematic%channel
    end if
    theta = kinematic%angle*pi/180.0_dp ! angle in d_egrees to radians
    call saclay_amplitudes(k_cm, theta, kinematic%channel, phases, d_phases, a, b, c, d, e, d_a, d_b, d_c, d_d, d_e)
    a = a + kinematic%em_amplitude(1)
    b = b + kinematic%em_amplitude(2)
    c = c + kinematic%em_amplitude(3)
    d = d + kinematic%em_amplitude(4)
    e = e + kinematic%em_amplitude(5)
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
    select case (trim(kinematic%type))
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
end subroutine scattering_obs

!!
!> @brief     Calulate the 1S0 scattering length
!!
!! Uses a delta shell representation and the variable phase equation
!! of the effective range expansion to calculate the \f$ ^1S_0 \f$ np 
!! scattering lenght
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine scattering_length(model, parameters, channel, a_length, da_length)
    implicit none
    type(nn_local_model), intent(in) :: model !< nn scattering model
    real(dp), intent(in), dimension(:) :: parameters !< potetial parameters
    character(len=2), intent(in) :: channel !< reaction channel ('pp' or 'np')
    real(dp), intent(out) :: a_length !< \f$ ^1S_0 \f$ scattering_length
    real(dp), intent(out), allocatable, dimension(:) :: da_length !< derivatives of the \f$ ^1S_0 \f$ scattering_length

    integer, parameter :: n_waves = 5, j_max = 1
    real(dp), parameter :: mu = 2*m_p*m_n/(m_p + m_n)
    real(dp) :: r, lambda
    real(dp), dimension(1:n_waves, 1:j_max) :: v_pw
    real(dp), allocatable, dimension(:, :, :) :: dv_pw
    real(dp), allocatable, dimension(:) :: d_lambda
    real(dp) :: diff, denominator, numerator
    real(dp), allocatable, dimension(:) :: d_denominator, d_numerator

    allocate(da_length, d_lambda, d_denominator, d_numerator, mold=parameters)
    
    if (channel == 'pp') then
        stop 'scattering length for pp channel is not implemented yet'
    endif

    r = model%dr/2._dp
    a_length = 0
    da_length = 0
    do
        if(r > model%r_max) exit
        call model%potential(parameters, r, channel, v_pw, dv_pw)
        lambda = v_pw(1, 1)*mu*model%dr/(hbar_c**2)
        d_lambda = dv_pw(:, 1, 1)*mu*model%dr/(hbar_c**2)
        diff = r - a_length
        numerator = a_length + r*lambda*diff
        denominator = 1 + lambda*diff
        d_numerator = da_length + r*(d_lambda*diff - lambda*da_length)
        d_denominator = d_lambda*diff - lambda*da_length
        a_length = numerator/denominator
        da_length = (d_numerator*denominator - numerator*d_denominator)/denominator**2
        r = r + model%dr
    enddo
    
end subroutine scattering_length


end module observables
