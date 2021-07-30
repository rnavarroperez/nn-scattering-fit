!!
!> @brief      av18 potential
!!
!! Module to calculate the AV18 potential in operator basis. See derivations and details in:
!!
!! Phys.Rev. C51 (1995) 38-51
!!
!! @author Rodrigo Navarro Perez
!!
module av18


use precisions, only : dp
use constants, only : hbar_c, mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, f2=>f_pi_n_2, &
    pi
use em_nn_potential, only : n_em_terms, em_potential, em_potential_in_st_basis
use st_basis_2_partial_waves, only : n_st_terms, uncoupled_pot, coupled_pot, st_2_pw_basis
use delta_shell, only : nn_model
use string_functions, only : mask_to_string

implicit none

private 

public :: n_parameters, default_params, av18_all_partial_waves, av18_operator, n_operators, display_parameters, set_av18_potential

integer, parameter :: n_parameters = 44 !< Number of phenomenological parameters
integer, parameter :: n_operators = 18  !< Number of operators in the AV18 basis
real(dp), parameter :: r_max = 12.5_dp !< Maximum integration radius for phase-shifts. In units of fm
real(dp), parameter :: delta_r = 1/128._dp !< Integration step for phases and the dueteron. In units of fm

real(dp), parameter, dimension(1:n_parameters) :: default_params = &
    [  -7.627010_dp,  1815.492000_dp, 1847.805900_dp, & ! S=1, T=1 c
        1.079850_dp,  -190.094900_dp, -811.204000_dp, & ! S=1, T=1 t
       -0.626970_dp,  -570.557100_dp,  819.122200_dp, & ! S=1, T=1 ls
        0.067090_dp,   342.066900_dp, -615.233900_dp, & ! S=1, T=1 l2
        0.741290_dp,     9.341800_dp, -376.438400_dp, & ! S=1, T=1 ls2
       -8.627700_dp,  2605.268200_dp,  441.973300_dp, & ! S=1, T=0 c
        1.485601_dp, -1126.835900_dp,  370.132400_dp, & ! S=1, T=0 t
        0.101800_dp,    86.065800_dp, -356.517500_dp, & ! S=1, T=0 ls
       -0.132010_dp,   253.435000_dp,   -1.007600_dp, & ! S=1, T=0 l2
        0.073570_dp,  -217.579100_dp,   18.393500_dp, & ! S=1, T=0 ls2
      -11.270280_dp,  3346.687400_dp, & ! S=0, T=1 c pp
      -10.667880_dp,  3126.554200_dp, & ! S=0, T=1 c np
        0.124720_dp,    16.778000_dp, & ! S=0, T=1 l2
       -2.099710_dp,  1204.430100_dp, & ! S=0, T=0 c
       -0.314520_dp,   217.455900_dp, & ! S=0, T=1 l2
       -3.921000_dp, & ! P CD c term 
        2.100000_dp,     0.500000_dp,    0.200000_dp & ! c_pi, r_ws, a_ws
    ] !< default parameters in the AV18 potential
contains

subroutine set_av18_potential(potential, parameters)
    implicit none
    type(nn_model), intent(out) :: potential
    real(dp), intent(out), allocatable, dimension(:) :: parameters

    allocate(parameters, source = default_params)
    potential%potential => av18_all_partial_waves
    potential%potential_components => av18_operator
    potential%display_subroutine => display_parameters
    potential%r_max = r_max
    potential%dr = delta_r
    potential%potential_type = 'local'
    potential%name = 'AV18'
    potential%relativistic_deuteron = .False.
    potential%full_em_wave = .true.
    potential%n_components = n_operators

    ! Properties of a delta-shell potential. We set them to zero
    potential%n_lambdas = 0
    potential%dr_core = 0._dp
    potential%dr_tail = 0._dp
end subroutine set_av18_potential

!!
!> @brief      av18 potential in all partial waves
!!
!! Given a set of parameters and a radius (in units of fm) returns the nuclear part of the 
!! AV18 potential (in units of MeV) in all partial waves and the derivative of each potential with
!! respect of the given parameters.
!! 
!! The role of every parameter and the derivation of the potential can be found
!!
!! Phys.Rev. C51 (1995) 38-51
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine av18_all_partial_waves(ap, r, reaction, v_pw, dv_pw)
    implicit none
    real(dp), intent(in) :: ap(:) !< Phenomenological parameters
    real(dp), intent(in) :: r !< radius in units of fm
    character(len=2), intent(in) :: reaction !< reaction channel: np, pp, or nn
    real(dp), intent(out) :: v_pw(:, :) !< AV18 potential in all partial waves, units of MeV
    real(dp), allocatable, intent(out) :: dv_pw(:, :, :) !< derivatives of v_nn with respect to the parameters in ap

    real(dp) :: v_nn(1:n_operators) 
    real(dp), allocatable :: dv_nn(:, :)
    real(dp) :: v_00(1:n_st_terms), v_01(1:n_st_terms), v_10(1:n_st_terms), v_11(1:n_st_terms), &
                v_em_st(1:n_st_terms)
    real(dp) :: dv_00(1:n_st_terms, 1:n_parameters), dv_01(1:n_st_terms, 1:n_parameters), &
                dv_10(1:n_st_terms, 1:n_parameters), dv_11(1:n_st_terms, 1:n_parameters)
    real(dp) :: v_em(1:n_em_terms)
    integer :: n_jwaves, n_waves
    integer :: tz1, tz2, l, s, j

    if (size(ap) /= n_parameters) stop 'incorrect number of parameters for v_pw in av18_all_partial_waves'

    call av18_operator(ap, r, v_nn, dv_nn)
    v_em = em_potential(r)

    n_waves  = size(v_pw, 1)
    n_jwaves = size(v_pw, 2)

    if (n_waves /= 5) stop 'incorrect number of waves for v_pw in av18_all_partial_waves'

    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:n_jwaves))
    v_pw = 0
    dv_pw = 0

    select case (trim(reaction))
    case ('pp')
        tz1 = 1
        tz2 = 1
    case ('np')
        tz1 = -1
        tz2 = 1
    case ('nn')
        tz1 = -1
        tz2 = -1
    case default
        stop 'incorrect reaction channel in av18_all_partial_waves'
    end select

    v_01 = operator_2_st_basis(tz1, tz2, 0, 1, v_nn)
    v_11 = operator_2_st_basis(tz1, tz2, 1, 1, v_nn)
    dv_01 = d_operator_2_st_basis(tz1, tz2, 0, 1, dv_nn)
    dv_11 = d_operator_2_st_basis(tz1, tz2, 1, 1, dv_nn)

    v_00 = 0._dp
    v_10 = 0._dp
    dv_00 = 0._dp
    dv_10 = 0._dp
    if (tz1*tz2 == -1) then
        v_00 = operator_2_st_basis(tz1, tz2, 0, 0, v_nn)
        v_10 = operator_2_st_basis(tz1, tz2, 1, 0, v_nn)
        dv_00 = d_operator_2_st_basis(tz1, tz2, 0, 0, dv_nn)
        dv_10 = d_operator_2_st_basis(tz1, tz2, 1, 0, dv_nn)
    endif

    call st_2_pw_basis(reaction, v_00, v_01, v_10, v_11, dv_00, dv_01, dv_10, dv_11, v_pw, dv_pw)

    ! Adding EM terms to the 1S0 partial wave
    v_em = em_potential(r)
    l = 0
    s = 0
    j = 0
    call em_potential_in_st_basis(reaction, s, v_em, v_em_st)
    v_pw(1, 1) = v_pw(1, 1) + uncoupled_pot(l, s, j, v_em_st)

    if (n_jwaves > 1) then
        ! Adding EM terms to the 3S1-3D1 coupled channel
        s = 1
        j = 1
        v_10 = 0._dp
        call em_potential_in_st_basis(reaction, s, v_em, v_em_st)
        v_pw(3:5, 2) = v_pw(3:5, 2) + coupled_pot(j, v_em_st)
    endif

end subroutine av18_all_partial_waves

!!
!> @brief      transform potential from operator to st basis
!!
!! Given a potential in the AV18 operator basis and the corresponding spin and isospin quantum
!! numbers, returns the potential in the spin-isospin basis.
!!
!! The five terms in the basis are central, tensor, spin-orbit, l squared, and spin-orbit squared
!!
!! @return     potential in spin-isospin basis
!!
!! @author     Rodrigo Navarro Perez
!!
function operator_2_st_basis(tz1, tz2, s, t, v_op) result(v_st)
    implicit none
    integer, intent(in) :: tz1 !< isospin projected in the z direction for the first particle
    integer, intent(in) :: tz2 !< isospin projected in the z direction for the second particle
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: t !< isospin quantum number 
    real(dp), intent(in) :: v_op(1:n_operators) !< potential in the AV18 operator basis
    real(dp) :: v_st(1:n_st_terms) !< potential in the spin-isospin basis
    integer :: s1ds2, t1dt2,t12
    s1ds2 = 4*s - 3
    t1dt2 = 4*t - 3
    t12 = 3*tz1*tz2 - t1dt2
    ! central term, c
    v_st(1) = v_op(1) + t1dt2*v_op(2) + s1ds2*v_op(3) + s1ds2*t1dt2*v_op(4) + t12*v_op(15) & 
        + s1ds2*t12*v_op(16) + (tz1+tz2)*v_op(18)
    ! tensor term, t
    v_st(2) = v_op(5) + t1dt2*v_op(6) + t12*v_op(17)
    ! spin-orbit term, ls
    v_st(3) = v_op(7) + t1dt2*v_op(8)
    ! l squared term, l2
    v_st(4) = v_op(9) + t1dt2*v_op(10) + s1ds2*v_op(11) + s1ds2*t1dt2*v_op(12)
    ! spin-orbit squared term, ls2
    v_st(5) = v_op(13) + t1dt2*v_op(14)
end function operator_2_st_basis

!!
!> @brief      transform derivatives potential from operator to st basis
!!
!! Given the derivatives of potential in the AV18 operator basis with respect to phenomenological
!! parameters and the corresponding spin and isospin quantum numbers, returns the derivatives in the
!! spin-isospin basis.
!!
!! The five terms in the basis are central, tensor, spin-orbit, l squared, and spin-orbit squared
!!
!! @return     derivatives of the potential in spin-isospin basis
!!
!! @author     Rodrigo Navarro Perez
!!
function d_operator_2_st_basis(tz1, tz2, s, t, dv_op) result(dv_st)
    implicit none
    integer, intent(in) :: tz1 !< isospin projected in the z direction for the first particle
    integer, intent(in) :: tz2 !< isospin projected in the z direction for the second particle
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: t !< isospin quantum number 
    real(dp), intent(in) :: dv_op(:, :) !< derivatives of the potential in the AV18 operator basis
    real(dp), allocatable :: dv_st(:, :) !< derivatives of the potential in the spin-isospin basis
    integer :: n_p, n_o, i
    n_p = size(dv_op,2)
    n_o = size(dv_op,1)
    if (n_o /= n_operators) stop 'incorrect number of operators in d_operator_2_st_basis'
    allocate(dv_st(1:n_st_terms, 1:n_p))
    do i = 1, n_p
        dv_st(:, i) = operator_2_st_basis(tz1, tz2, s, t, dv_op(:, i))
    enddo

end function d_operator_2_st_basis

!!
!> @brief      av18 potential in operator basis
!!
!! Given a set of parameters and a radius (in units of fm) returns the nuclear part of the 
!! AV18 potential (in units of MeV) in operator basis and the derivative of each potential with
!! respect of the given parameters.
!! 
!! The role of every parameter and the derivation of the potential can be found
!!
!! Phys.Rev. C51 (1995) 38-51
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine av18_operator(ap, r, v_nn, dv_nn)
    implicit none
    real(dp), intent(in)  :: ap(:) !< Phenomenological parameters
    real(dp), intent(in)  :: r !< radius in units of fm
    real(dp), intent(out) :: v_nn(:) !< AV18 potential in operator basis, units of MeV
    real(dp), allocatable, intent(out) :: dv_nn(:, :) !< derivatives of v_nn with respect to the parameters in ap

    real(dp), parameter :: &
        mu0 = mpi0/hbar_c, &
        muc = mpic/hbar_c, &
        mu  = mpi/hbar_c,  &
        small = 1.0e-4_dp

    real(dp) :: cpi, rws, aiws

    real(dp) :: x, x0, xc, ypi, tpi, ypi0, tpi0, ypic, tpic, rcut, ws, ws0, wsp, &
        wsx, wsx2, dypi00, dypic0, ypi0p, ypicp, tpi2, p11pp, p11np, p11nn, &
        pt1pp, pt1np, pt1nn, pls1, pl211, pls21, p10, pt0, pls0, pl210, pls20, &
        p01pp, p01np, p01nn, pl201, p00, pl200, p11, p11cd, p11cs, pt1, pt1cd, &
        pt1cs, p01, p01cd, p01cs

    real(dp), dimension(1:n_parameters) :: d_p11pp, d_p11np, d_p11nn, d_pt1pp, d_pt1np, d_pt1nn, &
        d_pls1, d_pl211, d_pls21, d_p10, d_pt0, d_pls0, d_pl210, d_pls20, d_p01pp, d_p01np, &
        d_p01nn, d_pl201, d_p00, d_pl200, d_p11, d_p11cd, d_p11cs, d_pt1, d_pt1cd, d_pt1cs, &
        d_p01, d_p01cd, d_p01cs

    real(dp), dimension(1:n_parameters) :: d_cpi, d_rws, d_aiws, d_rcut, d_ypi, d_tpi, d_ypi0, &
        d_tpi0, d_ypic, d_tpic, d_tpi2, d_ws , d_ws0, d_wsp, d_wsx, d_wsx2, d_dypi00, d_dypic0, &
        d_ypi0p, d_ypicp


    integer :: ip

    allocate(dv_nn(1:n_operators, 1:n_parameters))

    v_nn = 0
    dv_nn = 0

    d_p11pp =  0
    d_p11np =  0
    d_p11nn =  0
    d_pt1pp =  0
    d_pt1np =  0
    d_pt1nn =  0
    d_pls1  =  0
    d_pl211 =  0
    d_pls21 =  0
    d_p10   =  0
    d_pt0   =  0
    d_pls0  =  0
    d_pl210 =  0
    d_pls20 =  0
    d_p01pp =  0
    d_p01np =  0
    d_p01nn =  0
    d_pl201 =  0
    d_p00   =  0
    d_pl200 =  0
    d_cpi   =  0
    d_rws   =  0
    d_aiws  =  0
    d_rcut  =  0
    d_ypi   =  0
    d_tpi   =  0
    d_ypi0  =  0
    d_tpi0  =  0
    d_ypic  =  0
    d_tpic  =  0
    d_tpi2  =  0
    d_ws    =  0
    d_ws0   =  0
    d_wsp   =  0
    d_wsx   =  0
    d_wsx2  =  0
    d_dypi00 = 0
    d_dypic0 = 0
    d_ypi0p =  0
    d_ypicp =  0
    d_ypi   =  0
    d_tpi   =  0

    cpi = ap(42)
    rws = ap(43)
    aiws = 1._dp/ap(44)


    d_cpi(42) = 1._dp
    d_rws(43) = 1._dp
    d_aiws(44) = -1._dp/ap(44)**2

    x = mu*r
    x0 = mu0*r
    xc = muc*r

    if (r.le.small) then
        tpi = 3*cpi**2*r/mu**3
        ypi0 = (mpi0/mpic)**2*(mpi0/3)*cpi*r/mu0
        tpi0 = 3*cpi*ypi0/mu0**2
        ypic = (mpic/3)*cpi*r/muc
        tpic = 3*cpi*ypic/muc**2

        d_tpi = 6*cpi*d_cpi*r/mu**3
        d_ypi0 = (mpi0/mpic)**2*(mpi0/3)*d_cpi*r/mu0
        d_tpi0 = 3*d_cpi*ypi0/mu0**2 + 3*cpi*d_ypi0/mu0**2
        d_ypic = (mpic/3)*d_cpi*r/muc
        d_tpic = 3*d_cpi*ypic/muc**2 + 3*cpi*d_ypic/muc**2
    else
        rcut = 1-exp(-cpi*r*r) 
        ypi = exp(-x)*rcut/x   
        tpi = (1 + (3 + 3/x)/x)*ypi*rcut 
        ypi0 = (mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        tpi0 = (1 + (3 + 3/x0)/x0)*ypi0*rcut
        ypic = (mpic/3)*exp(-xc)*rcut/xc
        tpic = (1 + (3 + 3/xc)/xc)*ypic*rcut

        d_rcut = r*r*exp(-cpi*r*r)*d_cpi
        d_ypi = exp(-x)*d_rcut/x   
        d_tpi = (1 + (3 + 3/x)/x)*(d_ypi*rcut + ypi*d_rcut) 
        d_ypi0 = (mpi0/mpic)**2*(mpi0/3)*exp(-x0)*d_rcut/x0
        d_tpi0 = (1 + (3 + 3/x0)/x0)*(d_ypi0*rcut + ypi0*d_rcut)
        d_ypic = (mpic/3)*exp(-xc)*d_rcut/xc
        d_tpic = (1 + (3 + 3/xc)/xc)*(d_ypic*rcut + ypic*d_rcut)
    end if

    ypi0 = f2*ypi0
    ypic = f2*ypic
    tpi0 = f2*tpi0
    tpic = f2*tpic
    tpi2 = tpi*tpi 

    d_ypi0 = f2*d_ypi0
    d_ypic = f2*d_ypic
    d_tpi0 = f2*d_tpi0
    d_tpic = f2*d_tpic
    d_tpi2 = 2*tpi*d_tpi 

    ws = 1/(1 + exp((r - rws)*aiws)) 
    ws0 = 1/(1 + exp(-rws*aiws))   
    wsp = ws*(1 + aiws*exp(-rws*aiws)*ws0*r) 
    wsx = ws*x
    wsx2 = wsx*x

    d_ws = -1/(1 + exp((r - rws)*aiws))**2*(-aiws*d_rws + (r - rws)*d_aiws)*exp((r - rws)*aiws)
    d_ws0 = 1/(1 + exp(-rws*aiws))**2*(d_rws*aiws + rws*d_aiws)*exp(-rws*aiws)
    d_wsp = d_ws*(1 + aiws*exp(-rws*aiws)*ws0*r) + &
        ws*r*exp(-rws*aiws)*(d_ws0*aiws + ws0*d_aiws - ws0*aiws*(d_rws*aiws + rws*d_aiws)) 
    d_wsx = d_ws*x
    d_wsx2 = d_wsx*x

    dypi00 = (mpi0/mpic)**2*(mpi0/3)*cpi/mu0 
    dypic0 = (mpic/3)*cpi/muc
    ypi0p = ypi0 - f2*dypi00*ws*r/ws0
    ypicp = ypic - f2*dypic0*ws*r/ws0
    ypi = (ypi0 + 2*ypic)/3
    tpi = (tpi0 + 2*tpic)/3

    d_dypi00 = (mpi0/mpic)**2*(mpi0/3)*d_cpi/mu0 
    d_dypic0 = (mpic/3)*d_cpi/muc
    d_ypi0p = d_ypi0 - f2*r*((d_dypi00*ws + dypi00*d_ws)/ws0 - dypi00*ws/ws0**2*d_ws0)
    d_ypicp = d_ypic - f2*r*((d_dypic0*ws + dypic0*d_ws)/ws0 - dypic0*ws/ws0**2*d_ws0)
    d_ypi = (d_ypi0 + 2*d_ypic)/3
    d_tpi = (d_tpi0 + 2*d_tpic)/3

    p11pp =  ap( 1)*tpi2 + ap( 2)*wsp + ap( 3)*wsx2 + ypi0p
    p11np =  ap( 1)*tpi2 + (ap( 2) + 0.5_dp*ap(41))*wsp + ap( 3)*wsx2 - ypi0p + 2*ypicp
    p11nn =  ap( 1)*tpi2 + (ap( 2) + ap(41))*wsp + ap( 3)*wsx2 + ypi0p
    pt1pp =  ap( 4)*tpi2 + ap( 5)*wsx + ap( 6)*wsx2 + tpi0
    pt1np =  ap( 4)*tpi2 + ap( 5)*wsx + ap( 6)*wsx2 - tpi0  + 2*tpic
    pt1nn =  ap( 4)*tpi2 + ap( 5)*wsx + ap( 6)*wsx2 + tpi0
    pls1  =  ap( 7)*tpi2 + ap( 8)*wsp + ap( 9)*wsx2
    pl211 =  ap(10)*tpi2 + ap(11)*wsp + ap(12)*wsx2
    pls21 =  ap(13)*tpi2 + ap(14)*wsp + ap(15)*wsx2
    p10   =  ap(16)*tpi2 + ap(17)*wsp + ap(18)*wsx2 - ypi0p - 2*ypicp
    pt0   =  ap(19)*tpi2 + ap(20)*wsx + ap(21)*wsx2 - tpi0  - 2*tpic
    pls0  =  ap(22)*tpi2 + ap(23)*wsp + ap(24)*wsx2
    pl210 =  ap(25)*tpi2 + ap(26)*wsp + ap(27)*wsx2
    pls20 =  ap(28)*tpi2 + ap(29)*wsp + ap(30)*wsx2
    p01pp =  ap(31)*tpi2 + ap(32)*wsp - 3*ypi0p
    p01np =  ap(33)*tpi2 + ap(34)*wsp - 3*(-ypi0p + 2*ypicp)
    p01nn =  ap(31)*tpi2 + (ap(32) + ap(41))*wsp - 3*ypi0p
    pl201 =  ap(35)*tpi2 + ap(36)*wsp
    p00   =  ap(37)*tpi2 + ap(38)*wsp - 3*(-ypi0p - 2*ypicp)
    pl200 =  ap(39)*tpi2 + ap(40)*wsp

    d_p11pp =  ap( 1)*d_tpi2 + ap( 2)*d_wsp + ap( 3)*d_wsx2 + d_ypi0p
    d_p11np =  ap( 1)*d_tpi2 + (ap( 2) + 0.5_dp*ap(41))*d_wsp + ap( 3)*d_wsx2 - d_ypi0p + 2*d_ypicp
    d_p11nn =  ap( 1)*d_tpi2 + (ap( 2) + ap(41))*d_wsp + ap( 3)*d_wsx2 + d_ypi0p
    d_pt1pp =  ap( 4)*d_tpi2 + ap( 5)*d_wsx + ap( 6)*d_wsx2 + d_tpi0
    d_pt1np =  ap( 4)*d_tpi2 + ap( 5)*d_wsx + ap( 6)*d_wsx2 - d_tpi0  + 2*d_tpic
    d_pt1nn =  ap( 4)*d_tpi2 + ap( 5)*d_wsx + ap( 6)*d_wsx2 + d_tpi0
    d_pls1  =  ap( 7)*d_tpi2 + ap( 8)*d_wsp + ap( 9)*d_wsx2
    d_pl211 =  ap(10)*d_tpi2 + ap(11)*d_wsp + ap(12)*d_wsx2
    d_pls21 =  ap(13)*d_tpi2 + ap(14)*d_wsp + ap(15)*d_wsx2
    d_p10   =  ap(16)*d_tpi2 + ap(17)*d_wsp + ap(18)*d_wsx2 - d_ypi0p - 2*d_ypicp
    d_pt0   =  ap(19)*d_tpi2 + ap(20)*d_wsx + ap(21)*d_wsx2 - d_tpi0  - 2*d_tpic
    d_pls0  =  ap(22)*d_tpi2 + ap(23)*d_wsp + ap(24)*d_wsx2
    d_pl210 =  ap(25)*d_tpi2 + ap(26)*d_wsp + ap(27)*d_wsx2
    d_pls20 =  ap(28)*d_tpi2 + ap(29)*d_wsp + ap(30)*d_wsx2
    d_p01pp =  ap(31)*d_tpi2 + ap(32)*d_wsp - 3*d_ypi0p
    d_p01np =  ap(33)*d_tpi2 + ap(34)*d_wsp - 3*(-d_ypi0p + 2*d_ypicp)
    d_p01nn =  ap(31)*d_tpi2 + (ap(32) + ap(41))*d_wsp - 3*d_ypi0p
    d_pl201 =  ap(35)*d_tpi2 + ap(36)*d_wsp
    d_p00   =  ap(37)*d_tpi2 + ap(38)*d_wsp - 3*(-d_ypi0p - 2*d_ypicp)
    d_pl200 =  ap(39)*d_tpi2 + ap(40)*d_wsp

    d_p11pp(1) =  tpi2
    d_p11pp(2) =  wsp
    d_p11pp(3) =  wsx2

    d_p11np( 1) =  tpi2
    d_p11np( 2) =  wsp
    d_p11np(41) =  0.5_dp*wsp
    d_p11np( 3) =  wsx2

    d_p11nn( 1) =  tpi2
    d_p11nn( 2) =  wsp
    d_p11nn(41) =  wsp
    d_p11nn( 3) =  wsx2

    d_pt1pp(4) =  tpi2
    d_pt1pp(5) =  wsx
    d_pt1pp(6) =  wsx2

    d_pt1np(4) =  tpi2
    d_pt1np(5) =  wsx
    d_pt1np(6) =  wsx2

    d_pt1nn(4) =  tpi2
    d_pt1nn(5) =  wsx
    d_pt1nn(6) =  wsx2

    d_pls1(7) =  tpi2
    d_pls1(8) =  wsp
    d_pls1(9) =  wsx2

    d_pl211(10) =  tpi2
    d_pl211(11) =  wsp
    d_pl211(12) =  wsx2

    d_pls21(13) =  tpi2
    d_pls21(14) =  wsp
    d_pls21(15) =  wsx2

    d_p10(16) =  tpi2
    d_p10(17) =  wsp
    d_p10(18) =  wsx2

    d_pt0(19) =  tpi2
    d_pt0(20) =  wsx
    d_pt0(21) =  wsx2

    d_pls0(22) =  tpi2
    d_pls0(23) =  wsp
    d_pls0(24) =  wsx2

    d_pl210(25) =  tpi2
    d_pl210(26) =  wsp
    d_pl210(27) =  wsx2

    d_pls20(28) =  tpi2
    d_pls20(29) =  wsp
    d_pls20(30) =  wsx2

    d_p01pp(31) =  tpi2
    d_p01pp(32) =  wsp

    d_p01np(33) =  tpi2
    d_p01np(34) =  wsp

    d_p01nn(31) =  tpi2
    d_p01nn(32) =  wsp
    d_p01nn(41) =  wsp

    d_pl201(35) =  tpi2
    d_pl201(36) =  wsp

    d_p00(37) =  tpi2
    d_p00(38) =  wsp

    d_pl200(39) =  tpi2
    d_pl200(40) =  wsp

    ! p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p
    ! p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp
    ! p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p
    ! pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
    ! pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic
    ! pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
    ! pls1=    -.62697*tpi2 -570.5571*wsp +819.1222*wsx2
    ! pl211=    .06709*tpi2 +342.0669*wsp -615.2339*wsx2
    ! pls21=    .74129*tpi2   +9.3418*wsp -376.4384*wsx2
    ! p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2*ypicp
    ! pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2*tpic
    ! pls0=     .10180*tpi2  +86.0658*wsp -356.5175*wsx2
    ! pl210=   -.13201*tpi2 +253.4350*wsp   -1.0076*wsx2
    ! pls20=    .07357*tpi2 -217.5791*wsp  +18.3935*wsx2
    ! p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p
    ! p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp)
    ! p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p
    ! pl201=    .12472*tpi2  +16.7780*wsp
    ! p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2*ypicp)
    ! pl200=   -.31452*tpi2 +217.4559*wsp

    p11   = (p11pp + p11nn + p11np)/3
    pt1   = (pt1pp + pt1nn + pt1np)/3
    p01   = (p01pp + p01nn + p01np)/3

    p11cd = ((p11pp + p11nn)/2 - p11np)/6
    p01cd = ((p01pp + p01nn)/2 - p01np)/6
    pt1cd = ((pt1pp + pt1nn)/2 - pt1np)/6
    
    p11cs = (p11pp - p11nn)/4
    pt1cs = (pt1pp - pt1nn)/4
    p01cs = (p01pp - p01nn)/4

    d_p11   = (d_p11pp + d_p11nn + d_p11np)/3
    d_pt1   = (d_pt1pp + d_pt1nn + d_pt1np)/3
    d_p01   = (d_p01pp + d_p01nn + d_p01np)/3

    d_p11cd = ((d_p11pp + d_p11nn)/2 - d_p11np)/6
    d_p01cd = ((d_p01pp + d_p01nn)/2 - d_p01np)/6
    d_pt1cd = ((d_pt1pp + d_pt1nn)/2 - d_pt1np)/6
    

    d_p11cs = (d_p11pp - d_p11nn)/4
    d_pt1cs = (d_pt1pp - d_pt1nn)/4
    d_p01cs = (d_p01pp - d_p01nn)/4

    v_nn( 1) = (9*p11 + 3*p10 + 3*p01 + p00)/16 ! v_c
    v_nn( 2) = (3*p11 - 3*p10 +   p01 - p00)/16 ! v_tau
    v_nn( 3) = (3*p11 +   p10 - 3*p01 - p00)/16 ! v_sigma
    v_nn( 4) = (  p11 -   p10 -   p01 + p00)/16 ! v_sigma_tau
    v_nn( 5) = (3*pt1  + pt0 )/4 ! v_t
    v_nn( 6) = (  pt1  - pt0 )/4 ! v_t_tau
    v_nn( 7) = (3*pls1 + pls0)/4 ! v_ls
    v_nn( 8) = (  pls1 - pls0)/4 ! v_ls_tau
    v_nn( 9) = (9*pl211 + 3*pl210 + 3*pl201 + pl200)/16 ! v_l2
    v_nn(10) = (3*pl211 - 3*pl210 +   pl201 - pl200)/16 ! v_l2_tau
    v_nn(11) = (3*pl211 +   pl210 - 3*pl201 - pl200)/16 ! v_l2_sigma
    v_nn(12) = (  pl211 -   pl210 -   pl201 + pl200)/16 ! v_l2_sigma_tau
    v_nn(13) = (3*pls21 + pls20)/4 ! v_ls2
    v_nn(14) = (  pls21 - pls20)/4 ! v_ls2_tau
    v_nn(15) = (3*p11cd + p01cd)/4 ! v_T
    v_nn(16) = (  p11cd - p01cd)/4 ! v_sigma_T
    v_nn(17) = pt1cd ! v_t_T
    v_nn(18) = p01cs ! v_tau_z

    do ip = 1, size(dv_nn,2)
        dv_nn( 1, ip) = (9*d_p11(ip) + 3*d_p10(ip) + 3*d_p01(ip) + d_p00(ip))/16
        dv_nn( 2, ip) = (3*d_p11(ip) - 3*d_p10(ip) +   d_p01(ip) - d_p00(ip))/16
        dv_nn( 3, ip) = (3*d_p11(ip) +   d_p10(ip) - 3*d_p01(ip) - d_p00(ip))/16
        dv_nn( 4, ip) = (  d_p11(ip) -   d_p10(ip) -   d_p01(ip) + d_p00(ip))/16
        dv_nn( 5, ip) = (3*d_pt1(ip)  + d_pt0(ip) )/4
        dv_nn( 6, ip) = (  d_pt1(ip)  - d_pt0(ip) )/4
        dv_nn( 7, ip) = (3*d_pls1(ip) + d_pls0(ip))/4
        dv_nn( 8, ip) = (  d_pls1(ip) - d_pls0(ip))/4
        dv_nn( 9, ip) = (9*d_pl211(ip) + 3*d_pl210(ip) + 3*d_pl201(ip) + d_pl200(ip))/16
        dv_nn(10, ip) = (3*d_pl211(ip) - 3*d_pl210(ip) +   d_pl201(ip) - d_pl200(ip))/16
        dv_nn(11, ip) = (3*d_pl211(ip) +   d_pl210(ip) - 3*d_pl201(ip) - d_pl200(ip))/16
        dv_nn(12, ip) = (  d_pl211(ip) -   d_pl210(ip) -   d_pl201(ip) + d_pl200(ip))/16
        dv_nn(13, ip) = (3*d_pls21(ip) + d_pls20(ip))/4
        dv_nn(14, ip) = (  d_pls21(ip) - d_pls20(ip))/4
        dv_nn(15, ip) = (3*d_p11cd(ip) + d_p01cd(ip))/4
        dv_nn(16, ip) = (  d_p11cd(ip) - d_p01cd(ip))/4
        dv_nn(17, ip) = d_pt1cd(ip)
        dv_nn(18, ip) = d_p01cs(ip)
    enddo

end subroutine av18_operator

!!
!> @brief      Display the AV18 parameters
!!
!! Given a set of parameters for the AV18 potential with
!! the corresponding mask to indicate fixed parameters, displays
!! the parameters to screen.
!!
!! When the optional argument for the covariance matrix is given, 
!! the uncertainty for each parameter is display right below the 
!! parameter.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine display_parameters(ap, mask, unit, cv)
    implicit none
    real(dp), intent(in), dimension(:) :: ap !< parameters for the AV18 potential
    logical, intent(in), dimension(:) :: mask !< Indicates which parameters are optimized
    integer, intent(in) :: unit !< Unit where the output is sent to. Either and already opened file or output_unit from iso_fortran_env
    real(dp), intent(in), optional, dimension(:, :) :: cv !< Covariance matrix of the parameters

    character(len=18), parameter :: format_1 = '(1(2x,a,f15.8,a1))'
    character(len=18), parameter :: format_2 = '(2(2x,a,f15.8,a1))'
    character(len=18), parameter :: format_3 = '(3(2x,a,f15.8,a1))'
    character(len=17), parameter :: format_4 = '(3(12x,f15.8,1x))'

    character(len=size(ap)) :: s1
    integer :: i

    s1 = mask_to_string(mask, ' ', '*')
    write(unit, format_3) 'I_11    c:', ap( 1), s1(1:1), 'P_11    c:', ap( 2), s1(2:2), 'R_11    c:', ap( 3), s1(3:3)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=1, 3)
    write(unit, format_3) 'I_11    t:', ap( 4), s1(4:4), 'Q_11    t:', ap( 5), s1(5:5), 'R_11    t:', ap( 6), s1(6:6)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=4, 6)
    write(unit, format_3) 'I_11   ls:', ap( 7), s1(7:7), 'P_11   ls:', ap( 8), s1(8:8), 'R_11   ls:', ap( 9), s1(9:9)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=7, 9)
    write(unit, format_3) 'I_11   l2:', ap(10), s1(10:10), 'P_11   l2:', ap(11), s1(11:11), 'R_11   l2:', ap(12), s1(12:12)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=10, 12)
    write(unit, format_3) 'I_11  ls2:', ap(13), s1(13:13), 'P_11  ls2:', ap(14), s1(14:14), 'R_11  ls2:', ap(15), s1(15:15)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=13, 15)
    write(unit, format_3) 'I_10    c:', ap(16), s1(16:16), 'P_10    c:', ap(17), s1(17:17), 'R_10    c:', ap(18), s1(18:18)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=16, 18)
    write(unit, format_3) 'I_10    t:', ap(19), s1(19:19), 'Q_10    t:', ap(20), s1(20:20), 'R_10    t:', ap(21), s1(21:21)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=19, 21)
    write(unit, format_3) 'I_10   ls:', ap(22), s1(22:22), 'P_10   ls:', ap(23), s1(23:23), 'R_10   ls:', ap(24), s1(24:24)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=22, 24)
    write(unit, format_3) 'I_10   l2:', ap(25), s1(25:25), 'P_10   l2:', ap(26), s1(26:26), 'R_10   l2:', ap(27), s1(27:27)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=25, 27)
    write(unit, format_3) 'I_10  ls2:', ap(28), s1(28:28), 'P_10  ls2:', ap(29), s1(29:29), 'R_10  ls2:', ap(30), s1(30:30)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=28, 30)
    write(unit, format_2) 'I_01 c pp:', ap(31), s1(31:31), 'P_01 c pp:', ap(32), s1(32:32)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=31, 32)
    write(unit, format_2) 'I_01 c np:', ap(33), s1(33:33), 'P_01 c np:', ap(34), s1(34:34)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=33, 34)
    write(unit, format_2) 'I_01   l2:', ap(35), s1(35:35), 'P_01   l2:', ap(36), s1(36:36)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=35, 36)
    write(unit, format_2) 'I_00    c:', ap(37), s1(37:37), 'P_00    c:', ap(38), s1(38:38)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=37, 38)
    write(unit, format_2) 'I_00   l2:', ap(39), s1(39:39), 'P_00   l2:', ap(40), s1(40:40)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=39, 40)
    write(unit, format_1) 'P CD    c:', ap(41), s1(41:41) 
    if(present(cv)) write(unit, format_4) sqrt(cv(41,41))
    write(unit, format_3) 'c_pi     :', ap(42), s1(42:42), 'r_ws     :', ap(43), s1(43:43), 'a_ws     :', ap(44), s1(44:44)
    if(present(cv)) write(unit, format_4) (sqrt(cv(i,i)), i=42, 44)

end subroutine display_parameters

end module av18