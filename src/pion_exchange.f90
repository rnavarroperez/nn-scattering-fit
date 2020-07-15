module pion_exchange
use precisions, only : dp
use constants, only : hbar_c, mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, &
    pi, mu_p=>mu_proton, mu_n=>mu_neutron, alpha, m_e=>electron_mass, m_p=>proton_mass, &
    m_n=>neutron_mass, gamma=>euler_mascheroni

implicit none

private

public :: ope_all_partial_waves

integer, parameter :: n_st_terms = 5 !< Number of terms in the spin-isospin basis
integer, parameter :: n_em_terms = 14 !< Number of terms in the EM potential
integer, parameter :: n_fpis = 3

contains

subroutine ope_all_partial_waves(parameters, r, reaction, v_pw, dv_pw)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: reaction
    real(dp), intent(out), dimension(:, :) :: v_pw
    real(dp), intent(out), allocatable, dimension(:, :, :) :: dv_pw

    real(dp), dimension(1:n_st_terms) :: v_00, v_10, v_01, v_11, v_01_p_em
    real(dp), allocatable, dimension(:, :) :: dv_00, dv_10, dv_01, dv_11

    integer :: n_waves, j_max, n_parameters, s, t, l, j, ij, ip
    real(dp) :: v_em(1:n_em_terms)


    n_waves = size(v_pw, 1)
    j_max = size(v_pw, 2)
    n_parameters = size(parameters)

    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max))

    s = 0
    t = 0
    call ope_st_basis(parameters, r, reaction, s, t, v_00, dv_00)
    s = 0
    t = 1
    call ope_st_basis(parameters, r, reaction, s, t, v_01, dv_01)
    s = 1
    t = 0
    call ope_st_basis(parameters, r, reaction, s, t, v_10, dv_10)
    s = 1
    t = 1
    call ope_st_basis(parameters, r, reaction, s, t, v_11, dv_11)

    v_em = em_potential(r)
    v_01_p_em = v_01
    call add_em_potential(reaction, 0, v_em, v_01_p_em)
    ! 1s0
    l = 0
    s = 0
    j = 0
    v_pw(1, 1) = uncoupled_pot(l, s, j, v_01_p_em)
    do ip = 1, n_parameters
        dv_pw(ip, 1, 1) = uncoupled_pot(l, s, j, dv_01(:, ip))
    enddo
    ! 3p0
    l = 1
    s = 1
    j = 0
    v_pw(5, 1) = uncoupled_pot(l, s, j, v_11)
    do ip = 1, n_parameters
        dv_pw(ip, 5, 1) = uncoupled_pot(l, s, j, dv_11(:, ip))
    enddo
    ! everything with j >= 1
    do ij = 2, j_max
        j = ij - 1
        l = j
        ! singlets
        s = 0
        t = 1 - mod(l+s, 2)
        if (t == 1) then
            v_pw(1, ij) = uncoupled_pot(l, s, j, v_01)
            do ip = 1, n_parameters
                dv_pw(ip, 1, ij) = uncoupled_pot(l, s, j, dv_01(:, ip))
            enddo
        elseif (trim(reaction) == 'np') then ! only present in np
            v_pw(1, ij) = uncoupled_pot(l, s, j, v_00)
            do ip = 1, n_parameters
                dv_pw(ip, 1, ij) = uncoupled_pot(l, s, j, dv_00(:, ip))
            enddo
        endif
        ! triplets
        s = 1
        t = 1 - mod(l+s, 2)
        if (t == 1) then
            v_pw(2, ij) = uncoupled_pot(l, s, j, v_11)
            do ip = 1, n_parameters
                dv_pw(ip, 2, ij) = uncoupled_pot(l, s, j, dv_11(:, ip))
            enddo
        elseif (trim(reaction) == 'np') then !only present in np
            v_pw(2, ij) = uncoupled_pot(l, s, j, v_10)
            do ip = 1, n_parameters
                dv_pw(ip, 2, ij) = uncoupled_pot(l, s, j, dv_10(:, ip))
            enddo
        endif
        ! coupled channels
        t = 1 - mod(j,2)
        if (t == 1) then
            v_pw(3:5, ij) = coupled_pot(j, v_11)
            do ip = 1, n_parameters
                dv_pw(ip, 3:5, ij) = coupled_pot(j, dv_11(:, ip))
            enddo
        elseif (trim(reaction) == 'np') then !only present in np
            v_pw(3:5, ij) = coupled_pot(j, v_10)
            do ip = 1, n_parameters
                dv_pw(ip, 3:5, ij) = coupled_pot(j, dv_10(:, ip))
            enddo
        endif
    enddo
end subroutine ope_all_partial_waves

subroutine ope_st_basis(parameters, r, reaction, s,t, v_st, dv_st)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: reaction
    integer, intent(in) :: s
    integer, intent(in) :: t
    real(dp), intent(out), dimension(1:n_st_terms) :: v_st
    real(dp), intent(out), allocatable, dimension(:, :) ::  dv_st
    real(dp) :: y_0, t_0, y_c, t_c, y_total, t_total
    real(dp) :: f2c, f2pp, f2nn, f2np
    real(dp), dimension(1:n_fpis) :: fpi, df2c, df2pp, df2nn, df2np, dy_total, dt_total
    
    integer :: s1ds2, n_parameters, i
    s1ds2 = 4*s - 3
    n_parameters = size(parameters)
    allocate(dv_st(1:n_st_terms, 1:n_parameters))
    v_st = 0._dp
    dv_st = 0._dp
    fpi = parameters(n_parameters-2:n_parameters)
    call set_pion_couplings(fpi, f2c, f2pp, f2nn, f2np, df2c, df2pp, df2nn, df2np)
    y_0 = v_ope_s(mpi0, r)
    t_0 = v_ope_t(mpi0, r)
    select case(trim(reaction))
    case('pp')
        y_total = f2pp*y_0
        t_total = f2pp*t_0
        dy_total = df2pp*y_0
        dt_total = df2pp*t_0
    case('nn')
        y_total = f2nn*y_0 
        t_total = f2nn*t_0
        dy_total = df2nn*y_0 
        dt_total = df2nn*t_0
    case('np')
        y_c = v_ope_s(mpic, r)
        t_c = v_ope_t(mpic, r)
        y_total = f2np*y_0 + (-1)**(t+1)*2*f2c*y_c
        t_total = f2np*t_0 + (-1)**(t+1)*2*f2c*t_c
        dy_total = df2np*y_0 + (-1)**(t+1)*2*df2c*y_c
        dt_total = df2np*t_0 + (-1)**(t+1)*2*df2c*t_c
    case default
        stop 'incorrect reaction channel in ope_st_basis'
    end select
    v_st(1)  = s1ds2*y_total
    v_st(2)  = t_total
    do i = 1, n_fpis
        dv_st(1, n_parameters - 3 + i) = s1ds2*dy_total(i)
        dv_st(2, n_parameters - 3 + i) = dt_total(i)
    enddo
end subroutine ope_st_basis

subroutine set_pion_couplings(fpi, f2c, f2pp, f2nn, f2np, df2c, df2pp, df2nn, df2np)
    implicit none
    real(dp), intent(in), dimension(1:n_fpis) :: fpi
    real(dp), intent(out) :: f2c
    real(dp), intent(out) :: f2pp
    real(dp), intent(out) :: f2nn
    real(dp), intent(out) :: f2np
    real(dp), intent(out), dimension(1:n_fpis) :: df2c
    real(dp), intent(out), dimension(1:n_fpis) :: df2pp
    real(dp), intent(out), dimension(1:n_fpis) :: df2nn
    real(dp), intent(out), dimension(1:n_fpis) :: df2np

    real(dp), parameter :: small = 1.e-12_dp
    integer :: i, n_active

    df2c = 0._dp
    df2pp = 0._dp
    df2nn = 0._dp
    df2np = 0._dp

    n_active = 0
    do i = 1, n_fpis
        if (abs(fpi(i)) > epsilon(fpi(i))) n_active = n_active + 1
    enddo

    select case(n_active)
    case(0)
        f2c  =  0.075_dp
        f2pp =  0.075_dp
        f2nn =  0.075_dp
        f2np = -0.75_dp
    case(1)
        f2c  =  fpi(1)**2
        f2pp =  f2c
        f2nn =  f2c
        f2np = -f2c
        df2c(1) = 2*fpi(1)
        df2pp =  df2c
        df2nn =  df2c
        df2np = -df2c
    case(2)
        f2c  =  fpi(1)**2
        f2pp =  fpi(2)**2
        f2nn =  f2pp
        f2np = -f2c
        df2c(1)  = 2*fpi(1)
        df2pp(2) = 2*fpi(2)
        df2nn =  df2pp
        df2np = -df2c
    case(3)
        f2c  = fpi(1)**2
        f2pp = fpi(2)**2
        f2nn = fpi(3)**2
        f2np = fpi(2)*fpi(3)
        df2c(1)  = 2*fpi(1)
        df2pp(2) = 2*fpi(2)
        df2nn(3) = 2*fpi(3)
        df2np(2) = fpi(3)
        df2np(3) = fpi(2)
    end select
end subroutine set_pion_couplings

real(dp) function v_ope_s(m, r) result(v)
    implicit none
    real(dp) :: m, r
    real(dp) :: xi
    xi = m/mpic
    v = xi**2*m*yukawa_y(m/hbar_c, r)/3
end function v_ope_s

real(dp) function v_ope_t(m, r) result(v)
    implicit none
    real(dp) :: m, r
    real(dp) :: xi
    xi = m/mpic
    v = xi**2*m*yukawa_t(m/hbar_c, r)/3
end function v_ope_t


real(dp) function yukawa_y(m, r) result(y)
    implicit none
    real(dp), intent(in) :: m
    real(dp), intent(in) :: r
    y = exp(-m*r)/(m*r)
end function yukawa_y

real(dp) function yukawa_t(m, r) result(t)
    implicit none
    real(dp), intent(in) :: m
    real(dp), intent(in) :: r
    t = yukawa_y(m, r)*(1 + 3/(m*r) + 3/(m*r)**2)
end function yukawa_t


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
!> @brief      Adds EM potential to strong potential in basis
!!
!! Given an already calculated EM potential, reaction channel and spin quantum number, adds the 
!! EM potential to the strong potential in the spin-isospin channel
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine add_em_potential(reaction, s, v_em, v_st)
    implicit none
    character(len=2), intent(in) :: reaction !< reaction channel, 'pp', 'np' or 'nn'
    integer, intent(in) :: s !< spin quantum number
    real(dp), intent(in) :: v_em(1:n_em_terms) !< EM potential as defined in AV18 paper
    real(dp), intent(inout) :: v_st(1:n_st_terms) !< Strong potential in spin-isospin basis

    integer :: s1ds2

    s1ds2 = 4*s - 3
    select case (trim(reaction))
    case ('pp')
        v_st(1) = v_st(1) + v_em(1)*0 + v_em(2) + v_em(3) + v_em(4) + s1ds2*v_em(6)
        v_st(2) = v_st(2) + v_em(9)
        v_st(3) = v_st(3) + v_em(12)
    case ('np')
        v_st(1) = v_st(1) + v_em(5) + s1ds2*v_em(8)
        v_st(2) = v_st(2) + v_em(11)
        v_st(3) = v_st(3) + v_em(14)
    case ('nn')
        v_st(1) = v_st(1) + s1ds2*v_em(7)
        v_st(2) = v_st(2) + v_em(10)
    case default
        stop 'incorrect reaction channel in add_em_potential'
    end select

end subroutine add_em_potential

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

end module pion_exchange
