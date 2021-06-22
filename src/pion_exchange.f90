!!
!> @brief      Pion exchange potentials
!!
!! Subroutines and functions to calculate the one (and later multiple)
!! pion exchange potential and its derivatives with respect of the pion
!! coupling constants.
!!
!! @author     Rodrigo Navarro Perez
!!
module pion_exchange
use precisions, only : dp
use constants, only : hbar_c, mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, &
    pi
use em_nn_potential, only : n_em_terms, em_potential, add_em_potential
use st_basis_2_partial_waves, only : n_st_terms, uncoupled_pot, coupled_pot

implicit none

private

public :: ope_all_partial_waves

integer, parameter :: n_fpis = 3 !< number of nucleon-pion coupling constants

contains

!!
!> @brief      One pion exchange in all partial waves
!!
!! Given a set of parameters (whose last three elements are the 
!! nucleon-pion coupling constants, see set_pion_couplings for details), 
!! calculates the one pion exchange potential for all partial waves in the 
!! v_pw array. The derivatives of the potential with respect of the pion
!! couplings are also calculated.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine ope_all_partial_waves(parameters, r, reaction, v_pw, dv_pw)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    real(dp), intent(in) :: r !< Radius at which the potential is calculated
    character(len=*), intent(in) :: reaction !< reaction channel ('pp', 'nn', or 'np')
    real(dp), intent(out), dimension(:, :) :: v_pw !< OPE potential in all partial waves
    real(dp), intent(out), allocatable, dimension(:, :, :) :: dv_pw !< derivatives of the OPE potential with respect of the pion couplings

    real(dp), dimension(1:n_st_terms) :: v_00, v_10, v_01, v_11, v_01_p_em
    real(dp), allocatable, dimension(:, :) :: dv_00, dv_10, dv_01, dv_11

    integer :: n_waves, j_max, n_parameters, s, t, l, j, ij, ip
    real(dp) :: v_em(1:n_em_terms)

    logical, parameter :: full_pp = .true.


    n_waves = size(v_pw, 1)
    j_max = size(v_pw, 2)
    n_parameters = size(parameters)

    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max))
    dv_pw = 0._dp

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
    call add_em_potential(reaction, 0, v_em, full_pp, v_01_p_em)
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

!!
!> @brief      One pion exchange potential in spin isospin basis
!!
!! Given a set of parameters (whose last three elements are the 
!! nucleon-pion coupling constants, see set_pion_couplings for details), 
!! calculates the one pion exchange potential in the spin-isospin basis.
!! The derivatives of the potential with respect of the pion
!! couplings are also calculated.
!!
!! The structure and order of the spin-isospin basis is the same
!! that was adopted in the av18 module. (see st_basis_2_partial_waves module
!! for details)
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine ope_st_basis(parameters, r, reaction, s, t, v_st, dv_st)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    real(dp), intent(in) :: r !< Radius at which the potential is calculated
    character(len=*), intent(in) :: reaction !< reaction channel ('pp', 'nn', or 'np')
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: t !< isospin quantum number
    real(dp), intent(out), dimension(1:n_st_terms) :: v_st !< OPE potential in spin-isospin basis
    real(dp), intent(out), allocatable, dimension(:, :) ::  dv_st !< derivatives of OPE potential in spin-isospin basis
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

!!
!> @brief      Calculates the nucleon-pion coupling constants
!!
!! Given an array with the three coupling constants \f$ f_c \f$,
!! \f$ f_p \f$, and\f$ f_n \f$, calculates the nucleon-pion coupling 
!! constants \f$ f_c^2 \f$, \f$ f_{pp}^2 \f$, \f$ f_{nn}^2 \f$, 
!! and \f$ f_{np}^2 \f$. The derivatives with respect to the given constants are also
!! calculated
!!
!! If the three given constants are numerically zero (their absulte value smaller
!! than \f$ 10^{-12} \f$), the recommended default value of 
!! \f$ f_c^2 = f_{pp}^2 = f_{nn}^2 = -f_{np}^2 = 0.075 \f$ is used.
!!
!! If only the first element of the fpi array is numerically nonzero, then  
!! \f$ f_c^2 = f_{pp}^2 = f_{nn}^2 = -f_{np}^2 = \f$ fpi(1) is returned.
!!
!! If only the first and second elements of the fpi array are numerically nonzero, then
!! \f$ f_c^2 = -f_{np}^2 = \f$ fpi(1) and \f$ f_{pp}^2 = f_{nn}^2 = \f$ fpi(2) is returned.
!!
!! If all three elements of the fpi array are numerically nonzero, then
!! \f$ f_c^2 = \f$ fpi(1), \f$ f_{pp}^2 = \f$ fpi(2), \f$ f_{nn}^2 = \f$ fpi(3)
!! and \f$ f_{np}^2 = \f$ fpi(2)*fpi(3). Notice that for appropriate coupling constants
!! fpi(3) should be negative
!!
!! See Phys. Rev. C 95(2017) 064001 for detail
!!
! @return     { description_of_the_return_value }
!!
subroutine set_pion_couplings(fpi, f2c, f2pp, f2nn, f2np, df2c, df2pp, df2nn, df2np)
    implicit none
    real(dp), intent(in), dimension(1:n_fpis) :: fpi !< array containing the \f$f_c\f$, \f$f_{p}\f$, and \f$f_{n}\f$ coupling
    real(dp), intent(out) :: f2c !< \f$f_{c}^2\f$ coupling constant, dimensionless
    real(dp), intent(out) :: f2pp !< \f$f_{pp}^2\f$ coupling constant, dimensionless
    real(dp), intent(out) :: f2nn !< \f$f_{nn}^2\f$ coupling constant, dimensionless
    real(dp), intent(out) :: f2np !< \f$f_{np}^2\f$ coupling constant, dimensionless
    real(dp), intent(out), dimension(1:n_fpis) :: df2c !< derivatives of \f$f_{c}^2\f$ coupling constant with respect of \f$f_c\f$, \f$f_{p}\f$, and \f$f_{n}\f$
    real(dp), intent(out), dimension(1:n_fpis) :: df2pp !< derivatives of \f$f_{pp}^2\f$ coupling constant with respect of \f$f_c\f$, \f$f_{p}\f$, and \f$f_{n}\f$
    real(dp), intent(out), dimension(1:n_fpis) :: df2nn !< derivatives of \f$f_{nn}^2\f$ coupling constant with respect of \f$f_c\f$, \f$f_{p}\f$, and \f$f_{n}\f$
    real(dp), intent(out), dimension(1:n_fpis) :: df2np !< derivatives of \f$f_{np}^2\f$ coupling constant with respect of \f$f_c\f$, \f$f_{p}\f$, and \f$f_{n}\f$

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
        f2np = -0.075_dp
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

!!
!> @brief      spin component of one pion exchange potential
!!
!! Calculates the charge dependent (via the pion mass) spin component of the
!! one pion exchange potential
!!
!! See Phys. Rev. C 95(2017) 064001 for detail
!!
!! @returns    \f$ \left(\frac{m}{m_{\pi^\pm}}\right)^2 \frac{m}{3} Y_m(r)\f$
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function v_ope_s(m, r) result(v)
    implicit none
    real(dp), intent(in) :: m !< pion mass, in MeV
    real(dp), intent(in) :: r !< radius, in fm
    real(dp) :: xi
    xi = m/mpic
    v = xi**2*m*yukawa_y(m/hbar_c, r)/3
end function v_ope_s

!!
!> @brief      tensor component of one pion exchange potential
!!
!! Calculates the charge dependent (via the pion mass) tensor component of the
!! one pion exchange potential
!!
!! See Phys. Rev. C 95(2017) 064001 for detail
!!
!! @returns    \f$ \left(\frac{m}{m_{\pi^\pm}}\right)^2 \frac{m}{3} T_m(r)\f$ in units of MeV
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function v_ope_t(m, r) result(v)
    implicit none
    real(dp), intent(in) :: m !< pion mass, in MeV
    real(dp), intent(in) :: r !< radius, in fm
    real(dp) :: xi
    xi = m/mpic
    v = xi**2*m*yukawa_t(m/hbar_c, r)/3
end function v_ope_t

!!
!> @brief      Y Yukawa function
!!
!! The usual Y Yukawa function
!!
!! @return     \f$ \frac{e^{-mr}}{mr} \f$, dimensionless
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function yukawa_y(m, r) result(y)
    implicit none
    real(dp), intent(in) :: m !< pion mass, in fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< radius, in fm
    y = exp(-m*r)/(m*r)
end function yukawa_y

!!
!> @brief      T Yukawa function
!!
!! The usual T Yukawa function
!!
!! @return     \f$ \frac{e^{-mr}}{mr}\left[1 + \frac{3}{mr} +  \frac{3}{(mr)^2} \right] \f$, dimensionless
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function yukawa_t(m, r) result(t)
    implicit none
    real(dp), intent(in) :: m !< pion mass, in fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< radius, in fm
    t = yukawa_y(m, r)*(1 + 3/(m*r) + 3/(m*r)**2)
end function yukawa_t

end module pion_exchange
