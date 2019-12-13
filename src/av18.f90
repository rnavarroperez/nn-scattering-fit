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
    pi, mu_p=>mu_proton, mu_n=>mu_neutron, alpha, m_e=>electron_mass, m_p=>proton_mass, &
    m_n=>neutron_mass, gamma=>euler_mascheroni


implicit none

private 

public :: av18_operator, n_parameters, n_operators, default_params, av18_all_partial_waves

integer, parameter :: n_parameters = 44 !< Number of phenomenological parameters
integer, parameter :: n_operators = 18  !< Number of operators in the AV18 basis
integer, parameter :: n_st_terms = 5 !< Number of terms in the spin-isospin basis

real(dp), parameter, dimension(1:n_parameters) :: default_params = &
    [  -7.62701_dp, 1815.49200_dp, 1847.80590_dp,  1813.53150_dp, 1811.57100_dp,    1.07985_dp, &
     -190.09490_dp, -811.20400_dp,   -0.62697_dp,  -570.55710_dp,  819.12220_dp,    0.06709_dp, &
      342.06690_dp, -615.23390_dp,    0.74129_dp,     9.34180_dp, -376.43840_dp,   -8.62770_dp, &
     2605.26820_dp,  441.97330_dp,    1.48560_dp, -1126.83590_dp,  370.13240_dp,    0.10180_dp, &
       86.06580_dp, -356.51750_dp,   -0.13201_dp,   253.43500_dp,   -1.00760_dp,    0.07357_dp, &
     -217.57910_dp,   18.39350_dp,  -11.27028_dp,  3346.68740_dp,  -10.66788_dp, 3126.55420_dp, &
      -11.27028_dp, 3342.76640_dp,    0.12472_dp,    16.77800_dp,   -2.09971_dp, 1204.43010_dp, &
       -0.31452_dp,  217.45590_dp] !< default parameters in the AV18 potential
contains

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
    real(dp), intent(in) :: ap(1:n_parameters) !< Phenomenological parameters
    real(dp), intent(in) :: r !< radius in units of fm
    character(len=2), intent(in) :: reaction !< reaction channel: np, pp, or nn
    real(dp), intent(out) :: v_pw(:, :) !< AV18 potential in all partial waves, units of MeV
    real(dp), allocatable, intent(out) :: dv_pw(:, :, :) !< derivatives of v_nn with respect to the parameters in ap

    real(dp) :: v_nn(1:n_operators), dv_nn(1:n_parameters, 1:n_operators)
    real(dp) :: v_00(1:n_st_terms), v_01(1:n_st_terms), v_10(1:n_st_terms), v_11(1:n_st_terms)
    real(dp) :: dv_00(1:n_parameters, 1:n_st_terms), dv_01(1:n_parameters, 1:n_st_terms), &
                dv_10(1:n_parameters, 1:n_st_terms), dv_11(1:n_parameters, 1:n_st_terms)
    integer :: n_jwaves, n_waves
    integer :: tz1, tz2, ij, l, s, j, t, ip

    call av18_operator(ap, r, v_nn, dv_nn)

    n_waves  = size(v_pw, 1)
    n_jwaves = size(v_pw, 2)

    if (n_waves /= 5) stop 'incorrect number of waves for v_pw in av18_all_partial_waves'

    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:n_jwaves))
    dv_pw = 0

    select case (reaction)
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

    if (tz1*tz2 == -1) then
        v_00 = operator_2_st_basis(tz1, tz2, 0, 0, v_nn)
        v_10 = operator_2_st_basis(tz1, tz2, 1, 0, v_nn)
        dv_00 = d_operator_2_st_basis(tz1, tz2, 0, 0, dv_nn)
        dv_10 = d_operator_2_st_basis(tz1, tz2, 1, 0, dv_nn)
    endif

    ! 1s0
    l = 0
    s = 0
    j = 0
    v_pw(1, 1) = uncoupled_pot(l, s, j, v_01)

    ! 3p0
    l = 1
    s = 1
    j = 0
    v_pw(5, 1) = uncoupled_pot(l, s, j, v_11)

    ! everything with j >= 1
    do ij = 2, n_jwaves
        j = ij - 1
        l = j

        ! singlets
        s = 0
        t = 1 - mod(l+s, 2)
        if (t == 1) then
            v_pw(1, ij) = uncoupled_pot(l, s, j, v_01)
            do ip = 1, n_parameters
                dv_pw(ip, 1, ij) = uncoupled_pot(l, s, j, dv_01(ip,:))
            enddo
        elseif (tz1*tz2 == -1) then ! only present in np
            v_pw(1, ij) = uncoupled_pot(l, s, j, v_00)
            do ip = 1, n_parameters
                dv_pw(ip, 1, ij) = uncoupled_pot(l, s, j, dv_00(ip,:))
            enddo
        endif

        ! triplets
        s = 1
        t = 1 - mod(l+s, 2)
        if (t == 1) then
            v_pw(2, ij) = uncoupled_pot(l, s, j, v_11)
            do ip = 1, n_parameters
                dv_pw(ip, 2, ij) = uncoupled_pot(l, s, j, dv_11(ip,:))
            enddo
        elseif (tz1*tz2 == -1) then !only present in np
            v_pw(2, ij) = uncoupled_pot(l, s, j, v_10)
            do ip = 1, n_parameters
                dv_pw(ip, 2, ij) = uncoupled_pot(l, s, j, dv_10(ip,:))
            enddo
        endif

        ! coupled channels
        t = 1 - mod(j,2)
        if (t == 1) then
            v_pw(3:5, ij) = coupled_pot(j, v_11)
            do ip = 1, n_parameters
                dv_pw(ip, 3:5, ij) = coupled_pot(j, dv_11(ip,:))
            enddo
        elseif (tz1*tz2 == -1) then !only present in np
            v_pw(3:5, ij) = coupled_pot(j, v_10)
            do ip = 1, n_parameters
                dv_pw(ip, 3:5, ij) = coupled_pot(j, dv_10(ip,:))
            enddo
            if (j == 1 ) then
            ! add em part 
            endif
        endif
    enddo
   
end subroutine av18_all_partial_waves

!!
!> @brief      electromagnetic potential in NN scattering
!!
!! EM potential in NN scattering as calculated in the AV18 potential
!!
!!
!! For details see Phys.Rev. C51 (1995) 38-51
!!
!! @return     electromagnetic potential 
!!
!! @author     Rodrigo Navarro Perez
!!
function em_potential(r) result(v_em)
    implicit none
    real(dp), intent(in) :: r !< radius at which the potential is evaluated
    integer, parameter :: n_em_terms = 14
    real(dp) :: v_em(1:n_em_terms) !< electromagnetic potential

    real(dp), parameter :: b=4.27_dp, beta=0.0189_dp, small=0.e-5_dp
    real(dp) :: br, mr, fcoulr, ftr3, flsr3, kr, fivp, fdelta, fnpr

    v_em = 0
    br = b*r
    mr = m_p*m_n/(m_p + m_n)
    if (r < small) then
       fcoulr = 5*b/16
       ftr3 = b**3*br**2/720
       flsr3 = b**3/48
       kr = m_e*small/hbar_c
    else
        fcoulr = (1 - (1 +   11*br/16   + 3*br**2/16 + br**3/48)*exp(-br))/r
        ftr3 =   (1 - (1 + br + br**2/2 +   br**3/6  + br**4/24 + br**5/144)*exp(-br))/r**3
        flsr3 =  (1-  (1 + br + br**2/2 + 7*br**3/48 + br**4/48)*exp(-br))/r**3
        kr = m_e*r/hbar_c
    end if
    fivp = -gamma - 5/6._dp + abs(log(kr)) + 6*pi*kr/8
    fdelta = b**3*(1 + br + br**2/3)*exp(-br)/16
    fnpr = b**3*(15 + 15*br + 6*br**2 + br**3)*exp(-br)/384
    v_em(1) =  alpha*hbar_c*fcoulr
    v_em(2) = -alpha*hbar_c**3*fdelta/(4*m_p**2)
    v_em(3) = -v_em(1)**2/m_p
    v_em(4) = 2*alpha*v_em(1)*fivp/(3*pi)
    v_em(5) = alpha*hbar_c*beta*fnpr
    v_em(6) = -alpha*hbar_c**3*mu_p**2*fdelta/(6*m_p**2)
    v_em(7) = -alpha*hbar_c**3*mu_n**2*fdelta/(6*m_n**2)
    v_em(8) = -alpha*hbar_c**3*mu_p*mu_n*fdelta/(6*m_n*m_p)
    v_em(9) = -alpha*hbar_c**3*mu_p**2*ftr3/(4*m_p**2)
    v_em(10) = -alpha*hbar_c**3*mu_n**2*ftr3/(4*m_n**2)
    v_em(11) = -alpha*hbar_c**3*mu_p*mu_n*ftr3/(4*m_p*m_n)
    v_em(12) = -alpha*hbar_c**3*(4*mu_p - 1)*flsr3/(2*m_p**2)
    v_em(13) = 0
    v_em(14) = -alpha*hbar_c**3*mu_n*flsr3/(2*m_n*mr)
end function em_potential

!!
!> @brief      potential in uncoupled partial wave
!!
!! Given a potential in spin-isospin basis and a set of
!! angular momentum quantum numbers (l, s, j), calculates
!! the uncoupled potential in the corresponding partial wave
!!
!! @return     uncoupled potential in the given partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function uncoupled_pot(l, s, j, v_st) result(v_pw)
    implicit none
    integer, intent(in) :: l !< orbital angular momentum quantum number
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: j !< total angular momentum quantum number
    real(dp), intent(in) :: v_st(1:n_st_terms) !< potential in spin-isospin basis
    real(dp) :: v_pw !< uncoupled potential in the given partial wave
    integer :: ls, l2, s_12


    if (s == 0 .and. l == j) then ! singlets
        s_12 = 0
    elseif (s == 1 .and. l == j) then !triplet
        s_12 = 2
    elseif (s == 1 .and. l == 1 .and. j == 0) then !3p0
        s_12 = -4
    else
        stop 'incorrect quantum numbers in uncoupled_pot'
    endif

    ls = (j*(j + 1) - l*(l + 1) - s*(s + 1))/2 !numerator is always even, no need to worry about integer division
    l2 = l*(l + 1)
    v_pw =v_st(1) + s_12*v_st(2) + ls*v_st(3) + l2*v_st(4) + ls**2*v_st(5)
end function uncoupled_pot

!!
!> @brief      potential in uncoupled partial wave
!!
!! Given a potential in spin-isospin basis and the total angular momentum j, calculates the
!! coupled potential in the corresponding channel
!!
!! @return     coupled potential in the given j channel
!!
!! @author     Rodrigo Navarro Perez
!!
function coupled_pot(j, v_st) result(v_pw)
    implicit none
    integer, intent(in) :: j !< total angular momentum quantum number
    real(dp), intent(in) :: v_st(1:n_st_terms) !< potential in spin-isospin basis
    real(dp) :: v_pw(1:3) !< coupled potential in the j channel
    real(dp) :: s_12m, s_12, s_12p
    integer :: lm, lp, lm2, lp2, lsm, lsp

    lm = j - 1
    lp = j + 1
    lm2 = lm*(lm + 1)
    lp2 = lp*(lp + 1)
    lsm = j - 1    ! (j*(j + 1) - l*(l + 1) - s*(s + 1))/2 with l = j - 1 and s = 1
    lsp = -(j + 2) ! (j*(j + 1) - l*(l + 1) - s*(s + 1))/2 with l = j + 1 and s = 1
    s_12m = -2*(j - 1)/(2*j + 1._dp)
    s_12  = sqrt(36._dp*j*(j+1))/(2*j + 1)
    s_12p = -2*(j + 2)/(2*j + 1._dp)
    v_pw(1) = v_st(1) + s_12m*v_st(2) + lsm*v_st(3) + lm2*v_st(4) + lsm**2*v_st(5)
    v_pw(2) = s_12*v_st(2)
    v_pw(3) = v_st(1) + s_12p*v_st(2) + lsp*v_st(3) + lp2*v_st(4) + lsp**2*v_st(5)

end function coupled_pot

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
    ! central term
    v_st(1) = v_op(1) + t1dt2*v_op(2) + s1ds2*v_op(3) + s1ds2*t1dt2*v_op(4) + t12*v_op(15) & 
        + s1ds2*t12*v_op(16) + (tz1+tz2)*v_op(18)
    ! tensor term
    v_st(2) = v_op(5) + t1dt2*v_op(6) + t12*v_op(17)
    ! spin-orbit term
    v_st(3) = v_op(7) + t1dt2*v_op(8)
    ! l squared term
    v_st(4) = v_op(9) + t1dt2*v_op(10) + s1ds2*v_op(11) + s1ds2*t1dt2*v_op(12)
    ! spin-orbit squared term
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
    n_p = size(dv_op,1)
    n_o = size(dv_op,2)
    if (n_o /= n_operators) stop 'incorrect number of operators in doperator_2_st_basis'
    allocate(dv_st(1:n_p, 1:n_st_terms))
    do i = 1, n_p
        dv_st(i, :) = operator_2_st_basis(tz1, tz2, s, t, dv_op(i, :))
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
    real(dp), intent(in)  :: ap(1:n_parameters) !< Phenomenological parameters
    real(dp), intent(in)  :: r !< radius in units of fm
    real(dp), intent(out) :: v_nn(1:n_operators) !< AV18 potential in operator basis, units of MeV
    real(dp), intent(out) :: dv_nn(1:n_parameters, 1:n_operators) !< derivatives of v_nn with respect to the parameters in ap

    real(dp), parameter :: &
        mu0 = mpi0/hbar_c, &
        muc = mpic/hbar_c, &
        mu  = mpi/hbar_c,  &
        cpi = 2.1_dp, &
        rws = 0.5_dp, &
        aiws = 5._dp, &
        small = 1.0e-4_dp

    real(dp) :: x, x0, xc, ypi, tpi, ypi0, tpi0, ypic, tpic, rcut, ws, ws0, wsp, &
        wsx, wsx2, dypi00, dypic0, ypi0p, ypicp, tpi2, p11pp, p11np, p11nn, &
        pt1pp, pt1np, pt1nn, pls1, pl211, pls21, p10, pt0, pls0, pl210, pls20, &
        p01pp, p01np, p01nn, pl201, p00, pl200, p11, p11cd, p11cs, pt1, pt1cd, &
        pt1cs, p01, p01cd, p01cs

    real(dp), dimension(1:n_parameters) :: d_p11pp, d_p11np, d_p11nn, d_pt1pp, d_pt1np, d_pt1nn, &
        d_pls1, d_pl211, d_pls21, d_p10, d_pt0, d_pls0, d_pl210, d_pls20, d_p01pp, d_p01np, &
        d_p01nn, d_pl201, d_p00, d_pl200, d_p11, d_p11cd, d_p11cs, d_pt1, d_pt1cd, d_pt1cs, &
        d_p01, d_p01cd, d_p01cs

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

    x = mu*r
    x0 = mu0*r
    xc = muc*r

    if (r.le.small) then
        tpi = 3*cpi**2*r/mu**3
        ypi0 = (mpi0/mpic)**2*(mpi0/3)*cpi*r/mu0
        tpi0 = 3*cpi*ypi0/mu0**2
        ypic = (mpic/3)*cpi*r/muc
        tpic = 3*cpi*ypic/muc**2
    else
        rcut = 1-exp(-cpi*r*r) 
        ypi = exp(-x)*rcut/x   
        tpi = (1 + (3 + 3/x)/x)*ypi*rcut 
        ypi0 = (mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        tpi0 = (1 + (3 + 3/x0)/x0)*ypi0*rcut
        ypic = (mpic/3)*exp(-xc)*rcut/xc
        tpic = (1 + (3 + 3/xc)/xc)*ypic*rcut
    end if

    ypi0 = f2*ypi0
    ypic = f2*ypic
    tpi0 = f2*tpi0
    tpic = f2*tpic
    tpi2 = tpi*tpi 

    ws = 1/(1 + exp((r - rws)*aiws)) 
    ws0 = 1/(1 + exp(-rws*aiws))   
    wsp = ws*(1 + aiws*exp(-rws*aiws)*ws0*r) 
    wsx = ws*x
    wsx2 = wsx*x

    dypi00 = (mpi0/mpic)**2*(mpi0/3)*cpi/mu0 
    dypic0 = (mpic/3)*cpi/muc
    ypi0p = ypi0 - f2*dypi00*ws*r/ws0
    ypicp = ypic - f2*dypic0*ws*r/ws0
    ypi = (ypi0 + 2*ypic)/3
    tpi = (tpi0 + 2*tpic)/3

    p11pp =  ap( 1)*tpi2 + ap( 2)*wsp + ap( 3)*wsx2 + ypi0p
    p11np =  ap( 1)*tpi2 + ap( 4)*wsp + ap( 3)*wsx2 - ypi0p + 2*ypicp
    p11nn =  ap( 1)*tpi2 + ap( 5)*wsp + ap( 3)*wsx2 + ypi0p
    pt1pp =  ap( 6)*tpi2 + ap( 7)*wsx + ap( 8)*wsx2 + tpi0
    pt1np =  ap( 6)*tpi2 + ap( 7)*wsx + ap( 8)*wsx2 - tpi0 + 2*tpic
    pt1nn =  ap( 6)*tpi2 + ap( 7)*wsx + ap( 8)*wsx2 + tpi0
    pls1  =  ap( 9)*tpi2 + ap(10)*wsp + ap(11)*wsx2
    pl211 =  ap(12)*tpi2 + ap(13)*wsp + ap(14)*wsx2
    pls21 =  ap(15)*tpi2 + ap(16)*wsp + ap(17)*wsx2
    p10   =  ap(18)*tpi2 + ap(19)*wsp + ap(20)*wsx2 - ypi0p - 2*ypicp
    pt0   =  ap(21)*tpi2 + ap(22)*wsx + ap(23)*wsx2 - tpi0 - 2*tpic
    pls0  =  ap(24)*tpi2 + ap(25)*wsp + ap(26)*wsx2
    pl210 =  ap(27)*tpi2 + ap(28)*wsp + ap(29)*wsx2
    pls20 =  ap(30)*tpi2 + ap(31)*wsp + ap(32)*wsx2
    p01pp =  ap(33)*tpi2 + ap(34)*wsp - 3*ypi0p
    p01np =  ap(35)*tpi2 + ap(36)*wsp - 3*(-ypi0p + 2*ypicp)
    p01nn =  ap(37)*tpi2 + ap(38)*wsp - 3*ypi0p
    pl201 =  ap(39)*tpi2 + ap(40)*wsp
    p00   =  ap(41)*tpi2 + ap(42)*wsp - 3*(-ypi0p - 2*ypicp)
    pl200 =  ap(43)*tpi2 + ap(44)*wsp

    d_p11pp(1) =  tpi2
    d_p11pp(2) =  wsp
    d_p11pp(3) =  wsx2

    d_p11np(1) =  tpi2
    d_p11np(4) =  wsp
    d_p11np(3) =  wsx2

    d_p11nn(1) =  tpi2
    d_p11nn(5) =  wsp
    d_p11nn(3) =  wsx2

    d_pt1pp(6) =  tpi2
    d_pt1pp(7) =  wsx
    d_pt1pp(8) =  wsx2

    d_pt1np(6) =  tpi2
    d_pt1np(7) =  wsx
    d_pt1np(8) =  wsx2

    d_pt1nn(6) =  tpi2
    d_pt1nn(7) =  wsx
    d_pt1nn(8) =  wsx2

    d_pls1( 9) =  tpi2
    d_pls1(10) =  wsp
    d_pls1(11) =  wsx2

    d_pl211(12) =  tpi2
    d_pl211(13) =  wsp
    d_pl211(14) =  wsx2

    d_pls21(15) =  tpi2
    d_pls21(16) =  wsp
    d_pls21(17) =  wsx2

    d_p10(18) =  tpi2
    d_p10(19) =  wsp
    d_p10(20) =  wsx2

    d_pt0(21) =  tpi2
    d_pt0(22) =  wsx
    d_pt0(23) =  wsx2

    d_pls0(24) =  tpi2
    d_pls0(25) =  wsp
    d_pls0(26) =  wsx2

    d_pl210(27) =  tpi2
    d_pl210(28) =  wsp
    d_pl210(29) =  wsx2

    d_pls20(30) =  tpi2
    d_pls20(31) =  wsp
    d_pls20(32) =  wsx2

    d_p01pp(33) =  tpi2
    d_p01pp(34) =  wsp

    d_p01np(35) =  tpi2
    d_p01np(36) =  wsp

    d_p01nn(37) =  tpi2
    d_p01nn(38) =  wsp

    d_pl201(39) =  tpi2
    d_pl201(40) =  wsp

    d_p00(41) =  tpi2
    d_p00(42) =  wsp

    d_pl200(43) =  tpi2
    d_pl200(44) =  wsp

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

    v_nn( 1) = (9*p11 + 3*p10 + 3*p01 + p00)/16
    v_nn( 2) = (3*p11 - 3*p10 +   p01 - p00)/16
    v_nn( 3) = (3*p11 +   p10 - 3*p01 - p00)/16
    v_nn( 4) = (  p11 -   p10 -   p01 + p00)/16
    v_nn( 5) = (3*pt1  + pt0 )/4
    v_nn( 6) = (  pt1  - pt0 )/4
    v_nn( 7) = (3*pls1 + pls0)/4
    v_nn( 8) = (  pls1 - pls0)/4
    v_nn( 9) = (9*pl211 + 3*pl210 + 3*pl201 + pl200)/16
    v_nn(10) = (3*pl211 - 3*pl210 +   pl201 - pl200)/16
    v_nn(11) = (3*pl211 +   pl210 - 3*pl201 - pl200)/16
    v_nn(12) = (  pl211 -   pl210 -   pl201 + pl200)/16
    v_nn(13) = (3*pls21 + pls20)/4
    v_nn(14) = (  pls21 - pls20)/4
    v_nn(15) = (3*p11cd + p01cd)/4
    v_nn(16) = (  p11cd - p01cd)/4
    v_nn(17) = pt1cd
    v_nn(18) = p01cs

    dv_nn(:,  1) = (9*d_p11 + 3*d_p10 + 3*d_p01 + d_p00)/16
    dv_nn(:,  2) = (3*d_p11 - 3*d_p10 +   d_p01 - d_p00)/16
    dv_nn(:,  3) = (3*d_p11 +   d_p10 - 3*d_p01 - d_p00)/16
    dv_nn(:,  4) = (  d_p11 -   d_p10 -   d_p01 + d_p00)/16
    dv_nn(:,  5) = (3*d_pt1  + d_pt0 )/4
    dv_nn(:,  6) = (  d_pt1  - d_pt0 )/4
    dv_nn(:,  7) = (3*d_pls1 + d_pls0)/4
    dv_nn(:,  8) = (  d_pls1 - d_pls0)/4
    dv_nn(:,  9) = (9*d_pl211 + 3*d_pl210 + 3*d_pl201 + d_pl200)/16
    dv_nn(:, 10) = (3*d_pl211 - 3*d_pl210 +   d_pl201 - d_pl200)/16
    dv_nn(:, 11) = (3*d_pl211 +   d_pl210 - 3*d_pl201 - d_pl200)/16
    dv_nn(:, 12) = (  d_pl211 -   d_pl210 -   d_pl201 + d_pl200)/16
    dv_nn(:, 13) = (3*d_pls21 + d_pls20)/4
    dv_nn(:, 14) = (  d_pls21 - d_pls20)/4
    dv_nn(:, 15) = (3*d_p11cd + d_p01cd)/4
    dv_nn(:, 16) = (  d_p11cd - d_p01cd)/4
    dv_nn(:, 17) = d_pt1cd
    dv_nn(:, 18) = d_p01cs

end subroutine av18_operator


end module av18