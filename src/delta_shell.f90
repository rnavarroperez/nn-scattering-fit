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
use utilities, only : kronecker_delta
use pion_exchange, only : ope_all_partial_waves
use em_nn_potential, only : n_em_terms, em_potential
use string_functions, only : mask_to_string
implicit none

integer, parameter :: n_waves = 5 !< number of waves per angular momentum quantum number 
integer, parameter :: n_ds_parameters = 81 !< number of parameters in a DS potential
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
      0.00000000_dp, -0.48384594_dp,  0.00000000_dp, -0.02803913_dp, -0.00414759_dp, & !3P2
      0.00000000_dp,  0.29466313_dp,  0.19428354_dp,  0.04806034_dp,  0.01333910_dp, & !EP2
      0.00000000_dp,  3.45507794_dp, -0.22580425_dp,  0.00000000_dp, -0.01421141_dp, & !3F2
      0.00000000_dp,  0.00000000_dp,  0.11350103_dp,  0.09267627_dp,  0.00000000_dp, & !1F3
      0.00000000_dp,  0.53601375_dp,  0.00000000_dp,  0.00000000_dp,  0.00000000_dp, & !3D3
      0.00000000_dp,  0.00000000_dp,  0.00000000_dp, & !c1, c3, c4
      0.00000000_dp,  0.00000000_dp,  0.00000000_dp  & !fc fp fn
    ] !< Parameters  in the DS-OPE-30 potential

real(dp), parameter, dimension(1:n_ds_parameters) :: ds_ope30fff_params = &
    [ 1.30805321_dp, -0.71599779_dp, -0.19227977_dp,  0.00000000_dp, -0.02051403_dp, & !1S0 pp
     -0.14317826_dp, -0.05430227_dp,  0.03783094_dp,  0.00000000_dp, -0.00316311_dp, & !1S0 (np - pp)
      0.00000000_dp,  0.94270616_dp, -0.31863305_dp, -0.06229178_dp, -0.02266658_dp, & !3P0
      0.00000000_dp,  1.20156539_dp,  0.00000000_dp,  0.07454345_dp,  0.00000000_dp, & !1P1
      0.00000000_dp,  1.35425731_dp,  0.00000000_dp,  0.05697965_dp,  0.00000000_dp, & !3P1
      1.79212928_dp, -0.47487625_dp,  0.00000000_dp, -0.07201932_dp,  0.00000000_dp, & !3S1
      0.00000000_dp, -1.64879496_dp, -0.32628229_dp, -0.23328971_dp, -0.01824648_dp, & !EP1
      0.00000000_dp,  0.00000000_dp,  0.40446224_dp,  0.07026704_dp,  0.02084076_dp, & !3D1
      0.00000000_dp, -0.19647879_dp, -0.20564374_dp,  0.00000000_dp, -0.01870414_dp, & !1D2
      0.00000000_dp, -1.01391101_dp, -0.17034218_dp, -0.23731203_dp, -0.01603702_dp, & !3D2
      0.00000000_dp, -0.48237361_dp,  0.00000000_dp, -0.02886842_dp, -0.00370122_dp, & !3P2
      0.00000000_dp,  0.31926483_dp,  0.18971047_dp,  0.04959567_dp,  0.01266249_dp, & !EP2
      0.00000000_dp,  3.50354961_dp, -0.22935164_dp,  0.00000000_dp, -0.01402715_dp, & !3F3
      0.00000000_dp,  0.00000000_dp,  0.12395971_dp,  0.08901852_dp,  0.00000000_dp, & !1F3
      0.00000000_dp,  0.53571542_dp,  0.00000000_dp,  0.00000000_dp,  0.00000000_dp, & !3D3
      0.00000000_dp,  0.00000000_dp,  0.00000000_dp, & !c1, c3, c4
      0.27535448_dp,  0.27641819_dp, -0.28177242_dp  & !fc fp fn
    ] !< Parameters  in the DS-OPE-30-FFF potential

    real(dp), parameter, dimension(1:6) :: l_coulomb6 = [0.02566424_dp, 0.01374824_dp, 0.01351364_dp, &
        0.00685448_dp, 0.00860307_dp, 0.00214887_dp] !< Coulomb strength coefficients for a DS potential with 6 inner lambdas
    real(dp), parameter, dimension(1:5) :: l_coulomb5 = [0.02091441_dp, 0.01816750_dp, 0.00952244_dp, &
        0.01052224_dp, 0.00263887_dp] !< Coulomb strength coefficients for a DS potential with 5 inner lambdas
    real(dp), parameter, dimension(1:4) :: l_coulomb4 = [0.02530507_dp, 0.01398493_dp, 0.01358732_dp, &
        0.00338383_dp] !< Coulomb strength coefficients for a DS potential with 4 inner lambdas
    real(dp), parameter, dimension(1:3) :: l_coulomb3 = [0.02069940_dp, 0.01871309_dp, 0.00460163_dp] !< Coulomb strength coefficients for a DS potential with 3 inner lambdas
    real(dp), parameter, dimension(1:2) :: l_coulomb2 = [0.02662183_dp, 0.00663334_dp] !< Coulomb strength coefficients for a DS potential with 2 inner lambdas

private

public :: nn_model, all_delta_shells, ds_ope30_params, ds_ope30fff_params, set_ds_potential

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
!> @brief      interface for subroutines that display potential parameters
!!
!! @author     Rodrigo Navarro Perez
!!
interface
    subroutine display_parameters(ap, mask, cv)
        use precisions, only : dp
        implicit none
        real(dp), intent(in), dimension(:) :: ap !< potential parameters
        logical, intent(in), dimension(:) :: mask !< which parameters are kept fixed during the optimization
        real(dp), intent(in), optional, dimension(:, :) :: cv !< covariance matrix
    end subroutine display_parameters
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
    procedure(display_parameters), pointer, nopass :: display_subroutine !< local NN potential
    real(dp) :: r_max !< maximum intetgration radius
    real(dp) :: dr !< radial integration step
    character(len=24) :: potential_type !< Typr of nn potential ('local' or 'delta_shell')
    character(len=24) :: name !< potential name
    integer :: n_lambdas !< number of internal delta shells in a DS potential
    real(dp) :: dr_core !< Distance between the internal lambdas 
    real(dp) :: dr_tail !< Distance between the external lambdas (usually pion exchange)
    logical :: relativistic_deuteron !< Should relativistic kinematics be used when calculating the deuteron binding energy?
    logical :: full_em_wave !< Should the full EM wave (with 2 photon exchange and vacuum polarization) calculated for the 1S0 partial wave
end type nn_model

!!
!> @brief      strength coefficients set directly by potential parameters
!!
!! When decomposing a DS potential into partial wave we make a distinction between
!! 'active' partial waves (those set directly by the potential parameters) and
!! 'derived' partial waves (those obtained by converting the active ones to operator
!! basis and decomposing back to a desired partial wave).
!!
!! This type collects all the strength coefficients at a particular interaction radius
!! as well as it's derivatives (either zero or one) with respect of the potential
!! parameters
!!
!! @author     Rodrigo Navarro Perez
!!
type :: active_lambdas
    real(dp) :: l1s0 !< \f$ ^1S_0 \f$ strength coefficient
    real(dp) :: l3p0 !< \f$ ^3P_0 \f$ strength coefficient
    real(dp) :: l1p1 !< \f$ ^1P_1 \f$ strength coefficient
    real(dp) :: l3p1 !< \f$ ^3P_1 \f$ strength coefficient
    real(dp) :: l3s1 !< \f$ ^3S_1 \f$ strength coefficient
    real(dp) :: lep1 !< \f$ \epsilon_1 \f$ strength coefficient
    real(dp) :: l3d1 !< \f$ ^3D_1 \f$ strength coefficient
    real(dp) :: l1d2 !< \f$ ^1D_2 \f$ strength coefficient
    real(dp) :: l3d2 !< \f$ ^3D_2 \f$ strength coefficient
    real(dp) :: l3p2 !< \f$ ^3P_2 \f$ strength coefficient
    real(dp) :: lep2 !< \f$ \epsilon_2 \f$ strength coefficient
    real(dp) :: l3f2 !< \f$ ^3F_2 \f$ strength coefficient
    real(dp) :: l1f3 !< \f$ ^1F_3 \f$ strength coefficient
    real(dp) :: l3d3 !< \f$ ^3D_3 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d1s0 !< derivatives of \f$ ^1S_0 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3p0 !< derivatives of \f$ ^3P_0 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d1p1 !< derivatives of \f$ ^1P_1 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3p1 !< derivatives of \f$ ^3P_1 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3s1 !< derivatives of \f$ ^3S_1 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: dep1 !< derivatives of \f$ \epsilon_1 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3d1 !< derivatives of \f$ ^3D_1 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d1d2 !< derivatives of \f$ ^1D_2 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3d2 !< derivatives of \f$ ^3D_2 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3p2 !< derivatives of \f$ ^3P_2 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: dep2 !< derivatives of \f$ \epsilon_2 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3f2 !< derivatives of \f$ ^3F_2 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d1f3 !< derivatives of \f$ ^1F_3 \f$ strength coefficient
    real(dp), allocatable, dimension(:) :: d3d3 !< derivatives of \f$ ^3D_3 \f$ strength coefficient
end type active_lambdas

contains

!!
!> @brief      Sets the delta shell potential.
!!
!! Giventhe name of a delta shell potential, returns an object 
!! of type nn_model with the corresponding values set.
!!
!! The set of 'default' parameters to the corresponding potential
!! are also returned. These default parameters where determined
!! by fitting to the granada data base with a legacy code. 
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine set_ds_potential(name, ds_potential, parameters)
    implicit none
    character(len=*), intent(in) :: name
    type(nn_model), intent(out) :: ds_potential
    real(dp), intent(out), allocatable, dimension(:) :: parameters
    allocate(parameters(1:n_ds_parameters))
    select case(trim(name))
    case('ds_ope30')
        parameters = ds_ope30_params
        ds_potential%potential => ope_all_partial_waves
        ds_potential%display_subroutine => display_ds_parameters
        ds_potential%r_max = 13.0_dp
        ds_potential%potential_type = 'delta_shell'
        ds_potential%name = name
        ds_potential%n_lambdas = 5
        ds_potential%dr_core = 0.6_dp
        ds_potential%dr_tail = 0.5_dp
        ds_potential%relativistic_deuteron = .True.
        ds_potential%full_em_wave = .false.
    case('ds_ope30_fff')
        parameters = ds_ope30fff_params
        ds_potential%potential => ope_all_partial_waves
        ds_potential%display_subroutine => display_ds_parameters
        ds_potential%r_max = 13.0_dp
        ds_potential%potential_type = 'delta_shell'
        ds_potential%name = name
        ds_potential%n_lambdas = 5
        ds_potential%dr_core = 0.6_dp
        ds_potential%dr_tail = 0.5_dp
        ds_potential%relativistic_deuteron = .True.
        ds_potential%full_em_wave = .false.
    case default
        stop 'unrecognized ds potential name in set_ds_potential'
    end select

    ! Property of a local potential. We set it to zero
    ds_potential%dr = 0._dp

end subroutine set_ds_potential

!!
!> @brief      Display the DS parameters
!!
!! Given a set of parameters for a delta shell potential with
!! the corresponding mask to indicate fixed parameters, displays
!! the parameters to screen.
!!
!! When the optional argument for the covariance matrix is given, 
!! the uncertainty for each parameter is display right below the 
!! parameter.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine display_ds_parameters(ap, mask, cv)
    implicit none
    real(dp), intent(in), dimension(:) :: ap !< parameters for the Delta-Shell potential
    logical, intent(in), dimension(:) :: mask !< Indicates which parameters are optimized
    real(dp), intent(in), optional, dimension(:, :) :: cv !< Covariance matrix of the parameters

    character(len=size(ap)) :: s1
    integer :: i

    s1 = mask_to_string(mask, ' ', '*')

    print*, ' '
    print'(a15,4a16)', 'lambda_1', 'lambda_2', 'lambda_3', 'lambda_4', 'lambda_5'
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i = 1, 5), ' 1S0_pp'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 1, 5)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i = 6,10), ' 1S0_np - 1S0_pp'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 6, 10)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =11,15), ' 3P0'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 11, 15)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =16,20), ' 1P1'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 16, 20)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =21,25), ' 3P1'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 21, 25)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =26,30), ' 3S1'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 26, 30)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =31,35), ' Ep1'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 31, 35)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =36,40), ' 3D1'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 36, 40)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =41,45), ' 1D2'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 41, 45)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =46,50), ' 3D2'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 46, 50)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =51,55), ' 3P2'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 51, 55)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =56,60), ' Ep2'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 56, 60)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =61,65), ' 3F2'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 61, 65)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =66,70), ' 1F3'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 66, 70)
    print'(5(f15.8,a1),a)', (ap(i), s1(i:i), i =71,75), ' 3D3'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 71, 75)
    print'(3(f15.8,a1),a)', (ap(i), s1(i:i), i =76,78), ' c_1, c_3, c_4'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 76, 78)
    print'(3(f15.8,a1),a)', (ap(i), s1(i:i), i =79,81), ' f_c, f_p, f_n'
    if(present(cv)) print'(5(f15.8,1x))', (cv(i,i), i= 79, 81)

end subroutine display_ds_parameters

!!
!> @brief      delta shells in all partial waves
!!
!! Given a nn_model (local potential, maximum integration radius and integration step, etc),
!! returns a delta shell representation of said model.
!!
!! The delta shell representation includes an array of the concentration radii where 
!! the delta shells are located, the strength coefficients that multiply the delta shells,
!! an the derivatives of the strength coefficients with respect to the potential parameters.
!!
!! In the case of the proton-proton channel and a fully local model, a energy dependent Coulomb 
!! interaction is added to strength coefficients. The energy dependence is determined by the center
!! of mass momentum.
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

    call sample_radii(model, radii)
    call sample_local_potential(model, parameters, channel, k_cm, j_max, radii, v_pw, dv_pw)

    if (trim(model%potential_type) == 'delta_shell') then
        call inner_lambdas(model%n_lambdas, parameters, channel, v_pw, dv_pw)
    endif
end subroutine all_delta_shells

!!
!> @brief      get active lambdas from parameters array
!!
!! Given the array of DS potential parameters and the index of the interaction radius,
!! returns all the strength coefficient of all the 'active' partial waves in the
!! indicated interaction radius. It also returns the derivatives of the strength coefficients
!! with respect of the potential parameters (either 0 or 1)
!!
!! @return     object of type active_lambdas with all the inner strength coefficients
!!
!! @author     Rodrigo Navarro Perez
!!
type(active_lambdas) function extract_active_lambdas(index, n_lambdas, channel, parameters) result(r)
    implicit none
    integer, intent(in) :: index !< position of the concentration radius (between 1 and n_lambdas)
    integer, intent(in) :: n_lambdas !< number of inner delta shells
    character(len=*), intent(in) :: channel !< reaction channel ('pp' or 'np')
    real(dp), intent(in), dimension(:) :: parameters !< DS potential parameters

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

!!
!> @brief      strength coefficient of a given partial wave
!!
!! Given the quantum numbers of a specific partial wave, reaction channel and 
!! set of 'active' strength coefficients, returns the strength coefficient
!! for the corresponding partial wave. The derivatives of the strength coefficient
!! with respect of the potential parameters are also returned
!!
!! The strength coefficient is calculate by transforming the 'active' strength
!! coefficients to the spin-isospin basis and decomposing back to the partial
!! wave indicated by the quantum numbers
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine partial_wave_lamba(s, t, j, l1, l2, channel, lwaves, lambda, dlambda)
    implicit none
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: t !< isospin quantum number
    integer, intent(in) :: j !< total angular momentum quantum number
    integer, intent(in) :: l1 !< first orbital angular momentum quantum number
    integer, intent(in) :: l2 !< second orbital angular momentum quantum number
    character(len=*), intent(in) :: channel !< reaction channel ('pp' or 'np')
    type(active_lambdas), intent(in) :: lwaves !< set of 'active' strength coefficients (and their derivatives)
    real(dp), intent(out) :: lambda !< strength coefficient in the given partial wave
    real(dp), intent(out), dimension(:) :: dlambda !< derivatives of the strength coefficient in the given partial wave

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
        allocate(dvt, dvl2, dvls, dvls2, source = dvc)
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

!!
!> @brief      All the inner strength coefficients of a DS potential
!!
!! Given the DS potential parameters, returns the 'inner' strength coefficients
!! in all the partial waves. The strength coefficients are stored in the 
!! v_pw array. The dv_pw array stores the derivatives of the strength
!! coefficients with respect of the potential parameters
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine inner_lambdas(n_lambdas, parameters, channel, v_pw, dv_pw)
    implicit none
    integer, intent(in) :: n_lambdas !< number of inner lambdas in the DS potential
    real(dp), intent(in), dimension(:) :: parameters !< DS potential parameters
    character(len=*), intent(in) :: channel !< reaction channel ('pp' or 'np')
    real(dp), intent(inout), dimension(:, :, :) :: v_pw !< delta shell strength parameters. In fm\f$^{-1}\f$
    real(dp), intent(inout), dimension(:, :, :, :) :: dv_pw !< derivatives of the strength coefficient in the given partial wave
    

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

!!
!> @brief      Adds a coulomb strengths to the delta shell strength coefficients.
!!
!! For 'pp' potentials. Adds a set of strengths to the 'inner' strength
!! coefficients that correspond to a coarse grained representation of the coulomb
!! potential.
!!
!! The strengths added where calculated (with a legacy code) by adjusting a 
!! specific number of delta-shells (between 2 and 6) every 0.6 fm to reproduce
!! the nuclear phase shifts that result from the coulomb potential in the 
!! same range.
!!
!! A delta shell potential with a concentration radii with a distance different
!! from 0.6 fm (or more than 6 delta shells) can NOT be used with the strengths
!! used in this subroutine
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine add_delta_shell_coulomb(n_lambdas, v_pw)
    implicit none
    integer, intent(in) :: n_lambdas !< number of inner lambdas in the DS potential
    real(dp), intent(inout), dimension(:, :, :) :: v_pw !< delta shell strength parameters. In fm\f$^{-1}\f$

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

!!
!> @brief      Samples the concentration radii for a DS representation
!!
!! Given a nn_model type object, samples the concentration radii where
!! the DS strength coefficients will be located
!!
!! If the potential type is local, the radii are set to be equidistant
!! with the distance set in dr and starting with the first radii at 0.5*dr
!!
!! If the potential type is delta_shell, the a number of inner radii are
!! set to be equidistant starting with dr_core and ending at n_lambdas*dr_core.
!! The rest of the radii are also equidistant with the distance set by dr_tail
!! and starting at  n_lambdas*dr_core + 0.5*dr_tail.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine sample_radii(model, radii)
    implicit none
    type(nn_model), intent(in) :: model !< nn model
    real(dp), intent(out), allocatable, dimension(:) :: radii !< concentration radii in fm

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

!!
!> @brief      samples a local potential in a delta shell representation
!!
!! Given a nn model with a local potential, calculates a delta shell representation
!! of said potential in the given interaction radii
!!
!! If the potential type is local, the sampling is done for all given radii.
!! 
!! If the potential type is delta_shell, the sampling is done only for the 'outer' radii
!!
!! For a potential type local and a 'pp' reaction channel, an energy dependent coulomb 
!! contribution is added to the DS strength coefficients.
!!
!! @author     Rodrigo Navarro Perez
!!
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

!!
!> @brief      Reduced mass of a two nucleon system
!!
!! Given a reaction channel, calculates the corresponding reduced mass
!!
!! @returns    reduced mass in MeV
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function reduced_mass(channel) result(mu)
    implicit none
    character(len=*) :: channel !< reaction channel ('pp', 'nn', or 'np')
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
    real(dp), intent(inout) :: v_pw(:, :) !< pp potential for all partial waves in MeV

    real(dp), dimension(1:n_em_terms) :: v_em
    real(dp) :: kmev, relativistic_correction
    integer :: i
    kmev = k*hbar_c
    relativistic_correction = (1 + 2*kmev**2/m_p**2)/sqrt(1 + kmev**2/m_p**2)
    v_em = em_potential(r)
    v_em(1) = relativistic_correction*v_em(1)
    v_em(3:4) = relativistic_correction*v_em(3:4)
    v_pw(1, 1) = v_pw(1, 1) + v_em(3) + v_em(4)
    do i = 1, size(v_pw, 2)
        v_pw(1:3, i) = v_pw(1:3, i) + v_em(1)
        v_pw(5, i) = v_pw(5, i) + v_em(1)
    enddo
end subroutine add_coulomb

end module delta_shell