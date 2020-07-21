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
use num_recipes, only : kronecker_delta
implicit none

integer, parameter :: n_waves = 5
integer, parameter :: n_active_waves = 15
integer, parameter :: n_ds_parameters = 81
real(dp), parameter, dimension(1:n_ds_parameters) :: ds_ope30_params = &
    [ 1.31304950_dp, -0.71836850_dp, -0.19077786_dp,  0.00000000_dp, -0.02111031_dp, & !1S0 pp
     -0.15213746_dp, -0.05274411_dp,  0.03818241_dp,  0.00000000_dp, -0.00298297_dp, & !1S0(np - pp)
      0.00000000_dp,  0.95038079_dp, -0.32392986_dp, -0.05939468_dp, -0.02387808_dp, & !3P0
      0.00000000_dp,  1.19613616_dp,  0.00000000_dp,  0.07545449_dp,  0.00000000_dp, & !1P1
      0.00000000_dp,  1.35473545_dp,  0.00000000_dp,  0.05784667_dp,  0.00000000_dp, & !3P1
      1.79352488_dp, -0.47504116_dp,  0.00000000_dp, -0.07391203_dp,  0.00000000_dp, & !3S1
      0.00000000_dp, -1.64436516_dp, -0.33262344_dp, -0.22945708_dp, -0.02140123_dp, & !EP1
      0.00000000_dp,  0.00000000_dp,  0.40821749_dp,  0.06725153_dp,  0.02237340_dp, & !3D1
      0.00000000_dp, -0.20284762_dp, -0.20418276_dp,  0.00000000_dp, -0.01918031_dp, & !1D2
      0.00000000_dp, -1.01484880_dp, -0.17088142_dp, -0.23607873_dp, -0.01882201_dp, & !3D2
      0.00000000_dp, -0.48384594_dp, -0.00000000_dp, -0.02803913_dp, -0.00414759_dp, & !3P2
      0.00000000_dp,  0.29466313_dp,  0.19428354_dp,  0.04806034_dp,  0.01333910_dp, & !EP2
      0.00000000_dp,  3.45507794_dp, -0.22580425_dp,  0.00000000_dp, -0.01421141_dp, & !3F2
      0.00000000_dp,  0.00000000_dp,  0.11350103_dp,  0.09267627_dp,  0.00000000_dp, & !1F3
      0.00000000_dp,  0.53601375_dp,  0.00000000_dp,  0.00000000_dp,  0.00000000_dp, & !3D3
      0.00000000_dp,  0.00000000_dp,  0.00000000_dp, & !c1, c3, c4
      sqrt(0.075_dp), sqrt(0.075_dp), -sqrt(0.075_dp)  & !fc fp fn
    ]

    real(dp), parameter, dimension(1:6) :: l_coulomb6 = [0.02566424_dp, 0.01374824_dp, 0.01351364_dp, &
        0.00685448_dp, 0.00860307_dp, 0.00214887_dp]
    real(dp), parameter, dimension(1:5) :: l_coulomb5 = [0.02091441_dp, 0.01816750_dp, 0.00952244_dp, &
        0.01052224_dp, 0.00263887_dp]
    real(dp), parameter, dimension(1:4) :: l_coulomb4 = [0.02530507_dp, 0.01398493_dp, 0.01358732_dp, &
        0.00338383_dp]
    real(dp), parameter, dimension(1:3) :: l_coulomb3 = [0.02069940_dp, 0.01871309_dp, 0.00460163_dp]
    real(dp), parameter, dimension(1:2) :: l_coulomb2 = [0.02662183_dp, 0.00663334_dp]
private

public :: nn_model, all_delta_shells, ds_ope30_params

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
type :: nn_model
    procedure(local_potential), pointer, nopass :: potential !< local NN potential
    real(dp) :: r_max !< maximum intetgration radius
    real(dp) :: dr !< radial integration step
    character(len=24) :: potential_type !< Typr of nn potential ('local' or 'delta_shell')
    integer :: n_lambdas
    real(dp) :: dr_core !< Distance between the internal lambdas 
    real(dp) :: dr_tail 
end type nn_model

type :: active_lambdas
    real(dp) :: l1s0
    real(dp) :: l3p0
    real(dp) :: l1p1
    real(dp) :: l3p1
    real(dp) :: l3s1
    real(dp) :: lep1
    real(dp) :: l3d1
    real(dp) :: l1d2
    real(dp) :: l3d2
    real(dp) :: l3p2
    real(dp) :: lep2
    real(dp) :: l3f2
    real(dp) :: l1f3
    real(dp) :: l3d3
    real(dp), allocatable, dimension(:) :: d1s0
    real(dp), allocatable, dimension(:) :: d3p0
    real(dp), allocatable, dimension(:) :: d1p1
    real(dp), allocatable, dimension(:) :: d3p1
    real(dp), allocatable, dimension(:) :: d3s1
    real(dp), allocatable, dimension(:) :: dep1
    real(dp), allocatable, dimension(:) :: d3d1
    real(dp), allocatable, dimension(:) :: d1d2
    real(dp), allocatable, dimension(:) :: d3d2
    real(dp), allocatable, dimension(:) :: d3p2
    real(dp), allocatable, dimension(:) :: dep2
    real(dp), allocatable, dimension(:) :: d3f2
    real(dp), allocatable, dimension(:) :: d1f3
    real(dp), allocatable, dimension(:) :: d3d3
end type active_lambdas

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
    type(nn_model), intent(in) :: model !< nn model
    real(dp), intent(in), dimension(:) :: parameters !< fitting parameters
    character(len=*), intent(in) :: channel !< reaction channel (pp, np, or nn)
    real(dp), intent(in) :: k_cm !< center of mass momentum. In fm\f$^{-1}\f$
    integer, intent(in) :: j_max !< maximum j quantum number
    real(dp), intent(out), allocatable, dimension(:) :: radii !< concentration radii of the delta shells. In fm
    real(dp), intent(out), allocatable, dimension(:, :, :) :: v_pw !< delta shell strength parameters. In fm\f$^{-1}\f$
    real(dp), intent(out), allocatable, dimension(:, :, :, :) :: dv_pw !< derivatives of strength parameters with respect of potential parameters

    ! integer :: i

    call sample_radii(model, radii)
    call sample_local_potential(model, parameters, channel, k_cm, j_max, radii, v_pw, dv_pw)

    if (trim(model%potential_type) == 'delta_shell') then
        call inner_lambdas(model%n_lambdas, parameters, channel, v_pw, dv_pw)
    endif

    ! do i = 1, size(radii)
    !     print'(13e13.6)', radii(i), v_pw(1, 1, i), v_pw(5,1,i), v_pw(:, 2,i), v_pw(:,3,i)  
    ! enddo

    ! select case ()
    ! case ('local')
        
    ! case ('delta_shell')
        
    !     call set_strength_coefficients(n_lambdas, parameters, channel, j_max, radii, v_pw, dv_pw)
    ! case default
    !     stop 'Incorrect potential_type in all_delta_shells'
    ! end select

end subroutine all_delta_shells

! subroutine set_strength_coefficients(n_lambdas, parameters, channel, j_max, radii, v_pw, dv_pw)
!     implicit none
!     integer, intent(in) :: n_lambdas
!     real(dp), intent(in), dimension(:) :: parameters
!     character(len=*), intent(in) :: channel
!     integer, intent(in) :: j_max
!     real(dp), intent(in), dimension(:) :: radii
!     real(dp), intent(out), allocatable, dimension(:, :, :) :: v_pw !< delta shell strength parameters. In fm\f$^{-1}\f$
!     real(dp), intent(out), allocatable, dimension(:, :, :, :) :: dv_pw !< derivatives of strength parameters with respect of potential parameters

!     integer :: n_radii, n_parameters
!     n_radii = size(radii)
!     n_parameters = size(parameters)
!     allocate(v_pw(1:n_waves, 1:j_max, 1:n_radii))
!     allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max, 1:n_radii))
!     v_pw = 0._dp
!     dv_pw = 0._dp
    
!     ! call tail_lambdas(n_lambdas, parameters, channel, radii, v_pw, dv_pw)
    
! end subroutine set_strength_coefficients



type(active_lambdas) function extract_active_lambdas(index, n_lambdas, channel, parameters) result(r)
    implicit none
    integer, intent(in) :: index
    integer, intent(in) :: n_lambdas
    character(len=*), intent(in) :: channel
    real(dp), intent(in), dimension(:) :: parameters

    integer :: n_parameters

    n_parameters = size(parameters)

    if (allocated(r%d1s0)) deallocate(r%d1s0, r%d3p0, r%d1p1, r%d3p1, r%d3s1, r%dep1, r%d3d1, r%d1d2, &
        r%d3d2, r%d3p2, r%dep2, r%d3f2, r%d1f3, r%d3d3)

    allocate(r%d1s0(1:n_parameters))
    r%d1s0 = 0
    allocate(r%d3p0, r%d1p1, r%d3p1, r%d3s1, r%dep1, r%d3d1, r%d1d2, r%d3d2, r%d3p2, r%dep2, &
        r%d3f2, r%d1f3, r%d3d3, source = r%d1s0)

    select case(trim(channel))
    case('pp')
        r%l1s0 = parameters(index)
        r%d1s0(index) = 1._dp
    case('np')
        r%l1s0 = parameters(index) + parameters(n_lambdas + index)
        r%d1s0(index) = 1._dp
        r%d1s0(n_lambdas + index) = 1._dp
    case default
        stop 'Incorrect channel in extract_active_lambdas'
    end select

    r%l3p0 = parameters( 2*n_lambdas + index)
    r%l1p1 = parameters( 3*n_lambdas + index)
    r%l3p1 = parameters( 4*n_lambdas + index)
    r%l3s1 = parameters( 5*n_lambdas + index)
    r%lep1 = parameters( 6*n_lambdas + index)
    r%l3d1 = parameters( 7*n_lambdas + index)
    r%l1d2 = parameters( 8*n_lambdas + index)
    r%l3d2 = parameters( 9*n_lambdas + index)
    r%l3p2 = parameters(10*n_lambdas + index)
    r%lep2 = parameters(11*n_lambdas + index)
    r%l3f2 = parameters(12*n_lambdas + index)
    r%l1f3 = parameters(13*n_lambdas + index)
    r%l3d3 = parameters(14*n_lambdas + index)

    r%d3p0( 2*n_lambdas + index) = 1._dp
    r%d1p1( 3*n_lambdas + index) = 1._dp
    r%d3p1( 4*n_lambdas + index) = 1._dp
    r%d3s1( 5*n_lambdas + index) = 1._dp
    r%dep1( 6*n_lambdas + index) = 1._dp
    r%d3d1( 7*n_lambdas + index) = 1._dp
    r%d1d2( 8*n_lambdas + index) = 1._dp
    r%d3d2( 9*n_lambdas + index) = 1._dp
    r%d3p2(10*n_lambdas + index) = 1._dp
    r%dep2(11*n_lambdas + index) = 1._dp
    r%d3f2(12*n_lambdas + index) = 1._dp
    r%d1f3(13*n_lambdas + index) = 1._dp
    r%d3d3(14*n_lambdas + index) = 1._dp

end function extract_active_lambdas

subroutine partial_wave_lamba(s, t, j, l1, l2, channel, lwaves, lambda, dlambda)
    implicit none
    integer, intent(in) :: s
    integer, intent(in) :: t
    integer, intent(in) :: j
    integer, intent(in) :: l1
    integer, intent(in) :: l2
    character(len=*), intent(in) :: channel
    type(active_lambdas), intent(in) :: lwaves
    real(dp), intent(out) :: lambda
    real(dp), intent(out), dimension(:) :: dlambda

    real(dp), dimension(0:1, 0:1) :: vc, vt, vl2, vls, vls2
    real(dp), allocatable, dimension(: ,:, :) :: dvc, dvt, dvl2, dvls, dvls2
    real(dp) :: mu
    integer :: n_parameters

    if ((j == 0 .and. s == 1 .and. l1 /= 1) .or. (t == 0 .and. channel == 'pp')) then
        lambda = 0._dp
        dlambda = 0._dp
    else
        n_parameters = size(dlambda)
        allocate(dvc(1:n_parameters, 0:1, 0:1))
        dvc = 0._dp
        allocate(dvt, dvl2, dvls, dvls2, mold = dvc)
        vc(0, 0) = (6*lwaves%l1p1 - lwaves%l1f3)/5
        vc(1, 0) = lwaves%l3s1
        vc(0, 1) = lwaves%l1s0
        vc(1, 1) = 2*lwaves%l3p0/3 - lwaves%l3f2/5 + 8*lwaves%l3p2/15 + sqrt(2/3._dp)*16*lwaves%lep2/15
        dvc(:, 0, 0) = (6*lwaves%d1p1 - lwaves%d1f3)/5
        dvc(:, 1, 0) = lwaves%d3s1
        dvc(:, 0, 1) = lwaves%d1s0
        dvc(:, 1, 1) = 2*lwaves%d3p0/3 - lwaves%d3f2/5 + 8*lwaves%d3p2/15 + sqrt(2/3._dp)*16*lwaves%dep2/15

        vt(0, 0) = 0._dp
        vt(1, 0) = 0.5_dp*lwaves%lep1/sqrt(2._dp)
        vt(0, 1) = 0._dp
        vt(1, 1) = 5*lwaves%lep2/(6*sqrt(6._dp))
        dvt(:, 1, 0) = 0.5_dp*lwaves%dep1/sqrt(2._dp)
        dvt(:, 1, 1) = 5*lwaves%dep2/(6*sqrt(6._dp))

        vl2(0, 0) = (lwaves%l1f3 - lwaves%l1p1)/10
        vl2(1, 0) = ((lwaves%l3d3 - lwaves%l3d1)/10 + &
                       (lwaves%l3d2 - lwaves%l3s1)/2 - 2*sqrt(2._dp)*lwaves%lep1/7)/3
        vl2(0, 1) = (lwaves%l1d2 - lwaves%l1s0)/6
        vl2(1, 1) = (lwaves%l3f2 - lwaves%l3p2)/10 + (lwaves%l3p1 - lwaves%l3p0)/2 - 2*sqrt(6._dp)*lwaves%lep2/5
        dvl2(:, 0, 0) = (lwaves%d1f3 - lwaves%d1p1)/10
        dvl2(:, 1, 0) = ((lwaves%d3d3 - lwaves%d3d1)/10 + &
                       (lwaves%d3d2 - lwaves%d3s1)/2 - 2*sqrt(2._dp)*lwaves%dep1/7)/3
        dvl2(:, 0, 1) = (lwaves%d1d2 - lwaves%d1s0)/6
        dvl2(:, 1, 1) = (lwaves%d3f2 - lwaves%d3p2)/10 + (lwaves%d3p1 - lwaves%d3p0)/2 - 2*sqrt(6._dp)*lwaves%dep2/5

        vls(0, 0) = 0._dp
        vls(1, 0) = -lwaves%l3d1/10 - lwaves%l3d2/6 + 4*lwaves%l3d3/15 + lwaves%lep1/(7*sqrt(2._dp))
        vls(0, 1) = 0._dp
        vls(1, 1) = 0.5_dp*(lwaves%l3p2 - lwaves%l3p1) + lwaves%lep2/sqrt(6._dp)
        dvls(:, 1, 0) = -lwaves%d3d1/10 - lwaves%d3d2/6 + 4*lwaves%d3d3/15 + lwaves%dep1/(7*sqrt(2._dp))
        dvls(:, 1, 1) = 0.5_dp*(lwaves%d3p2 - lwaves%d3p1) + lwaves%dep2/sqrt(6._dp)
        
        vls2(0, 0) = 0._dp
        vls2(1, 0) = lwaves%l3d1/10 - lwaves%l3d2/6 + lwaves%l3d3/15 + sqrt(2._dp)*lwaves%lep1/7
        vls2(0, 1) = 0._dp
        vls2(1, 1) = lwaves%l3p0/3 - 0.5_dp*lwaves%l3p1 + lwaves%l3p2/6 + sqrt(2/3._dp)*lwaves%lep2
        dvls2(:, 1, 0) = lwaves%d3d1/10 - lwaves%d3d2/6 + lwaves%d3d3/15 + sqrt(2._dp)*lwaves%dep1/7
        dvls2(:, 1, 1) = lwaves%d3p0/3 - 0.5_dp*lwaves%d3p1 + lwaves%d3p2/6 + sqrt(2/3._dp)*lwaves%dep2

        lambda = ((1 - (-1)**(l1 + s + t))*(kronecker_delta(l1, l2)*(vc(s, t) + l1*(l1 + 1)*vl2(s, t) + &
                 ((j*(j + 1) - l1*(l1 + 1) - s*(s + 1))**2*kronecker_delta(l1, l2)**2*vls2(s, t))/4._dp + &
                 ((j*(j + 1) - l1*(l1 + 1) - s*(s + 1))*kronecker_delta(l1, l2)*vls(s, t))/2._dp) + &
                 kronecker_delta(1, s)*(2*kronecker_delta(j, l1)*kronecker_delta(j, l2) - &
                 (2*(j - 1)*kronecker_delta(j, 1 + l1)*kronecker_delta(j, 1 + l2))/(1._dp + 2*J) - &
                 (2*(j + 2)*kronecker_delta(1 + j, l1)*kronecker_delta(1 + j, l2))/(1._dp + 2*J) + &
                 (6*sqrt(j*(1._dp + j))* (kronecker_delta(j - 1, l2)*kronecker_delta(1 + j, l1) + &
                 kronecker_delta(j - 1, l1)*kronecker_delta(1 + j, l2)))/(1._dp + 2*j))*vt(s, t)))/2._dp
        dlambda = ((1 - (-1)**(l1 + s + t))*(kronecker_delta(l1, l2)*(dvc(: ,s, t) + l1*(l1 + 1)*dvl2(: ,s, t) + &
                  ((j*(j + 1) - l1*(l1 + 1) - s*(s + 1))**2*kronecker_delta(l1, l2)**2*dvls2(: ,s, t))/4._dp + &
                  ((j*(j + 1) - l1*(l1 + 1) - s*(s + 1))*kronecker_delta(l1, l2)*dvls(: ,s, t))/2._dp) + &
                  kronecker_delta(1, s)*(2*kronecker_delta(j, l1)*kronecker_delta(j, l2) - &
                  (2*(j - 1)*kronecker_delta(j, 1 + l1)*kronecker_delta(j, 1 + l2))/(1._dp + 2*J) - &
                  (2*(j + 2)*kronecker_delta(1 + j, l1)*kronecker_delta(1 + j, l2))/(1._dp + 2*J) + &
                  (6*sqrt(j*(1._dp + j))* (kronecker_delta(j - 1, l2)*kronecker_delta(1 + j, l1) + &
                  kronecker_delta(j - 1, l1)*kronecker_delta(1 + j, l2)))/(1._dp + 2*j))*dvt(: ,s, t)))/2._dp

    endif
    if (channel == 'np') then
        mu = reduced_mass(channel)
        lambda = mu/m_p*lambda
        dlambda = mu/m_p*dlambda
    endif
end subroutine partial_wave_lamba

subroutine inner_lambdas(n_lambdas, parameters, channel, v_pw, dv_pw)
    implicit none
    integer, intent(in) :: n_lambdas
    real(dp), intent(in), dimension(:) :: parameters
    character(len=*), intent(in) :: channel
    real(dp), intent(inout), dimension(:, :, :) :: v_pw
    real(dp), intent(inout), dimension(:, :, :, :) :: dv_pw
    

    type(active_lambdas) :: low_waves
    integer :: i, ij, j, l1, l2, s, t

    do i = 1, n_lambdas
        low_waves = extract_active_lambdas(i, n_lambdas, channel, parameters)
        do ij = 1, size(v_pw, 2)
            j = ij - 1
            l1 = j
            l2 = l1
            s = 0
            t = 1 - MOD(l1+s,2)
            call partial_wave_lamba(s, t, j, l1, l2, channel, low_waves, v_pw(1, ij, i), dv_pw(:, 1, ij, i))
            s = 1
            t = 1 - MOD(l1+s,2)
            call partial_wave_lamba(s, t, j, l1, l2, channel, low_waves, v_pw(2, ij, i), dv_pw(:, 2, ij, i))
            l1 = j-1
            l2 = l1
            t = 1 - MOD(j,2)
            call partial_wave_lamba(s, t, j, l1, l2, channel, low_waves, v_pw(3, ij, i), dv_pw(:, 3, ij, i))
            l1 = j-1
            l2 = j+1
            call partial_wave_lamba(s, t, j, l1, l2, channel, low_waves, v_pw(4, ij, i), dv_pw(:, 4, ij, i))
            l1 = j+1
            l2 = l1
            call partial_wave_lamba(s, t, j, l1, l2, channel, low_waves, v_pw(5, ij, i), dv_pw(:, 5, ij, i))
        enddo
    enddo
    if (trim(channel) == 'pp') then
        call add_delta_shell_coulomb(n_lambdas, v_pw)
    endif
end subroutine inner_lambdas

subroutine add_delta_shell_coulomb(n_lambdas, v_pw)
    implicit none
    integer, intent(in) :: n_lambdas
    real(dp), intent(inout), dimension(:, :, :) :: v_pw

    real(dp), allocatable, dimension(:) :: l_coulomb
    integer :: i, j

    allocate(l_coulomb(1:n_lambdas))

    select case(n_lambdas)
    case(2)
        l_coulomb = l_coulomb2
    case(3)
        l_coulomb = l_coulomb3
    case(4)
        l_coulomb = l_coulomb4
    case(5)
        l_coulomb = l_coulomb5
    case(6)
        l_coulomb = l_coulomb6
    case default
        stop 'delta shell models only have between 2 and 6 inner delta shells'
    end select

    do i = 1, n_lambdas
        do j = 1, size(v_pw, 2)
            v_pw(1:3, j, i) = v_pw(1:3, j, i) + l_coulomb(i)
            v_pw(5, j, i) = v_pw(5, j, i) + l_coulomb(i)
        enddo
    enddo
end subroutine add_delta_shell_coulomb

subroutine sample_radii(model, radii)
    implicit none
    type(nn_model), intent(in) :: model
    real(dp), intent(out), allocatable, dimension(:) :: radii

    integer :: i, n_radii, n_tail
    real(dp) :: r_cut

    select case (trim(model%potential_type))
    case ('local')
        n_radii = int(model%r_max/model%dr)
        allocate(radii(1:n_radii))
        do i = 1, n_radii
            radii(i) = (i - 0.5_dp)*model%dr
        enddo
    case ('delta_shell')
        n_tail = int((model%r_max - model%n_lambdas*model%dr_core)/model%dr_tail)
        allocate(radii(1:model%n_lambdas + n_tail))
        do i = 1, model%n_lambdas
            radii(i) = i*model%dr_core
        enddo
        r_cut = radii(model%n_lambdas) - 0.5_dp*model%dr_tail
        do i = 1, n_tail
            radii(i + model%n_lambdas) = r_cut + i*model%dr_tail
        enddo
    case default
        stop 'Incorrect potential_type in sample_radii'
    end select
end subroutine sample_radii

subroutine sample_local_potential(model, parameters, channel, k_cm, j_max, radii, v_pw, dv_pw)
    implicit none
    type(nn_model), intent(in) :: model !< nn model
    real(dp), intent(in), dimension(:) :: parameters !< fitting parameters
    character(len=*), intent(in) :: channel !< reaction channel (pp, np, or nn)
    real(dp), intent(in) :: k_cm !< center of mass momentum. In fm\f$^{-1}\f$
    integer, intent(in) :: j_max !< maximum j quantum number
    real(dp), intent(in), dimension(:) :: radii !< concentration radii of the delta shells. In fm
    real(dp), intent(out), allocatable, dimension(:, :, :) :: v_pw !< delta shell strength parameters. In fm\f$^{-1}\f$
    real(dp), intent(out), allocatable, dimension(:, :, :, :) :: dv_pw !< derivatives of strength parameters with respect of potential parameters
    
    
    integer :: n_radii, n_parameters, i, n_lambdas
    real(dp) :: mu, r, dr
    real(dp), allocatable, dimension(:, :, :) :: my_dv_pw

    mu = reduced_mass(channel)

    n_radii = size(radii)
    n_parameters = size(parameters)
    allocate(v_pw(1:n_waves, 1:j_max, 1:n_radii))
    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max, 1:n_radii))
    v_pw = 0._dp
    dv_pw = 0._dp

    select case (trim(model%potential_type))
    case ('local')
        n_lambdas = 0
        dr = model%dr
    case ('delta_shell')
        n_lambdas = model%n_lambdas
        dr = model%dr_tail
    case default
        stop 'Incorrect potential_type in sample_local_potential'
    end select

    do i = n_lambdas + 1, n_radii
        r = radii(i)
        call model%potential(parameters, r, channel, v_pw(:, :, i), my_dv_pw)
        dv_pw(:, :, :, i) = my_dv_pw
        if (trim(channel) == 'pp' .and. trim(model%potential_type) == 'local') then
            call add_coulomb(r, k_cm, v_pw(:, :, i))
        endif
    enddo
    v_pw = v_pw*mu*dr/(hbar_c**2)
    dv_pw = dv_pw*mu*dr/(hbar_c**2)    

end subroutine sample_local_potential

real(dp) function reduced_mass(channel) result(mu)
    implicit none
    character(len=*) :: channel
    select case (trim(channel))
    case ('pp')
        mu = m_p
    case ('np')
        mu = 2*m_p*m_n/(m_p + m_n)
    case ('nn')
        mu = m_n
    case default
        stop 'incorrect reaction channel in reduced_mass'
    end select
end function reduced_mass

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