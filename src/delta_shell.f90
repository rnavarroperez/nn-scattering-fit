!!
!> @brief      delta shell interactions
!!
!! Subroutines and functions to create delta shell representations
!! of fully local potentials (AV18, SOG, etc) or purely delta shell
!! potentials (DS-OPE, DS-TPE, etc) 
!!
!! @author     Rodrigo Navarro Perez
!!
module delta_shell

use precisions, only : dp
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi, alpha
implicit none

private

public :: nn_local_model, all_delta_shells

!!
!> @brief      interface of nn local potentials
!!
!! @author     Rodrigo Navarro Perez
!!
interface
    subroutine local_potential(ap, r, reaction, v_pw, dv_pw)
        use precisions, only : dp
        implicit none
        real(dp), intent(in) :: ap(:) !< potential parameters
        real(dp), intent(in) :: r !< radius (in fm) at which the potential is evaluated
        character(len=2), intent(in) :: reaction !< reaction channel. 'pp' or 'np'
        real(dp), intent(out) :: v_pw(:, :) !< local potential in all partial waves
        real(dp), allocatable, intent(out) :: dv_pw(:, :, :) !< derivatives of the potential with respect of the parameters
    end subroutine local_potential
end interface

!!
!> @brief      nn model for local interactions
!!
!! potential and fixed parameters to calculate all nuclear phase shifts
!!
!! @author     Rodrigo Navarro Perez
!!
type :: nn_local_model
    procedure(local_potential), pointer, nopass :: potential !< local NN potential
    real(dp) :: r_max !< maximum intetgration radius
    real(dp) :: dr !< radial integration step
end type nn_local_model

contains

!!
!> @brief      delta shells in all partial waves
!!
!! Given a nn model (local potential, maximum integration radius and integration step),
!! returns a delta shell representation of said model.
!!
!! The delta shell representation includes an array of the concentration radii where 
!! the delta shells are located, the strength coefficients that multiply the delta shells,
!! an the derivatives of the strength coefficients with respect to the potential parameters.
!!
!! In the case of the proton-proton channel, a energy dependent Coulomb interaction is
!! added to nn potential. The energy dependence is determined by the center of mass momemtum.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine all_delta_shells(model, parameters, channel, k_cm, j_max, radii, v_pw, dv_pw)
    implicit none
    type(nn_local_model), intent(in) :: model !< nn model
    real(dp), intent(in), dimension(:) :: parameters !< fitting parameters
    character(len=*), intent(in) :: channel !< reaction channel (pp, np, or nn)
    real(dp), intent(in) :: k_cm !< center of mass momentum. In fm\f$^{-1}\f$
    integer, intent(in) :: j_max !< maximum j quantum number
    real(dp), intent(out), allocatable, dimension(:) :: radii !< concentration radii of the delta shells. In fm
    real(dp), intent(out), allocatable, dimension(:, :, :) :: v_pw !< delta shell strength parameters. In fm\f$^{-1}\f$
    real(dp), intent(out), allocatable, dimension(:, :, :, :) :: dv_pw !< derivatives of strength parameters with respect of potential parameters
    
    integer, parameter :: n_waves = 5
    integer :: n_radii, n_parameters, i
    real(dp) :: mu
    real(dp), allocatable, dimension(:, :, :) :: my_dv_pw

    select case (channel)
    case ('pp')
        mu = m_p
    case ('np')
        mu = 2*m_p*m_n/(m_p + m_n)
    case ('nn')
        mu = m_n
    case default
        stop 'incorrect reaction channel in all_delta_shells'
    end select

    n_radii = int(model%r_max/model%dr)
    n_parameters = size(parameters)
    allocate(radii(1:n_radii))
    allocate(v_pw(1:n_waves, 1:j_max, 1:n_radii))
    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max, 1:n_radii))

    do i = 1, n_radii
        radii(i) = (i - 0.5_dp)*model%dr
        call model%potential(parameters, radii(i), channel, v_pw(:, :, i), my_dv_pw)
        dv_pw(:, :, :, i) = my_dv_pw
        if (channel == 'pp') then
            call add_coulomb(radii(i), k_cm, v_pw(:, :, i))
        endif
    enddo
    v_pw = v_pw*mu*model%dr/(hbar_c**2)
    dv_pw = dv_pw*mu*model%dr/(hbar_c**2)
end subroutine all_delta_shells

!!
!> @brief      Adds energy dependent Coulomb term to pp potential
!!
!! Adds an energy dependent Coulomb term to the pp potential in all partial waves
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine add_coulomb(r, k, v_pw)
    implicit none
    real(dp), intent(in) :: r !< potential radius in fm
    real(dp), intent(in) :: k !< center of mass momentum in fm\f$^{-1}\f$
    real(dp), intent(out) :: v_pw(:, :) !< pp potential for all partial waves in MeV
    integer :: i
    real(dp) :: v_coul, fcoulr, br, kmev, alphap
    real(dp), parameter :: b = 4.27_dp, small = 0.e-5_dp
    kmev = k*hbar_c
    alphap = alpha*(1 + 2*kmev**2/m_p**2)/sqrt(1 + kmev**2/m_p**2)
    br = b*r
    if (r < small) then
       fcoulr = 5*b/16
    else
        fcoulr = (1 - (1 +   11*br/16   + 3*br**2/16 + br**3/48)*exp(-br))/r
    end if

    v_coul = alphap*hbar_c*fcoulr

    do i = 1, size(v_pw, 2)
        v_pw(1:3, i) = v_pw(1:3, i) + v_coul
        v_pw(5, i) = v_pw(5, i) + v_coul
    enddo
end subroutine add_coulomb

end module delta_shell