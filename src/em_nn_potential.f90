!!
!> @brief      electromagnetic NN potential
!!
!! Module for calculating the electromagnetic nucleon nucleon potential (WITHOUT Coulomb)
!! that is added to the local NN potential used to calculate phasehifts 
!!
!! @author     Rodrigo Navarro Peerez
!!
module em_nn_potential
use precisions, only : dp
use constants, only : hbar_c, mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, f2=>f_pi_n_2, &
    pi, mu_p=>mu_proton, mu_n=>mu_neutron, alpha, m_e=>electron_mass, m_p=>proton_mass, &
    m_n=>neutron_mass, gamma=>euler_mascheroni
use st_basis_2_partial_waves, only : n_st_terms
use quadrature, only : booles_quadrature
implicit none

private

public :: n_em_terms, add_em_potential, em_potential, vacuum_polarization_integral, compare_Ivps

integer, parameter :: n_em_terms = 14 !< Number of terms in the EM potential
contains

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
        v_st(1) = v_st(1) + v_em(1)*0 + v_em(2) + v_em(3)*0 + v_em(4)*0 + s1ds2*v_em(6)
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

    real(dp), parameter :: b=4.27_dp, beta=0.0189_dp, small=1.e-5_dp, vp_cutoff=0.025_dp
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
    if (2*kr < vp_cutoff) then
        fivp = -gamma - 5/6._dp + abs(log(kr)) + 6*pi*kr/8
    else
        fivp = vacuum_polarization_integral(r)
    endif
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

subroutine compare_Ivps(r_min, r_max, dr)
    implicit none
    real(dp), intent(in) :: r_min, r_max, dr

    real(dp) :: r, Ivp, Ivp_small, Ivp_large, kr

    r = r_min
    do
        if(r > r_max) exit
        kr = m_e*r/hbar_c
        Ivp = vacuum_polarization_integral(r)
        Ivp_small = -gamma - 5/6._dp + abs(log(kr)) + 6*pi*kr/8
        Ivp_large = 2*sqrt(2*pi)/4._dp*exp(-2*kr)/((2*kr)**1.5_dp)
        print*, r, 2*kr, Ivp, Ivp_small, Ivp_large
        r = r + dr
    enddo
    
end subroutine compare_Ivps

real(dp) function vacuum_polarization_integral(r) result(Ivp)
    implicit none
    real(dp), intent(in) :: r

    integer, parameter :: n_points = 401
    real(dp), parameter :: t_min = 0._dp
    real(dp), parameter :: t_max = 1._dp

    real(dp), dimension(1:n_points) :: ft
    real(dp) :: t, delta_t
    integer :: i

    ft = 0._dp
    delta_t = (t_max - t_min)/(n_points - 1._dp)

    do i = 2, n_points - 1
        t = (i-1)*delta_t
        ft(i) = inverse_vp_kernel(r, t)
    enddo
    Ivp = booles_quadrature(ft, delta_t)


end function vacuum_polarization_integral

real(dp) function vacuum_polarization_kernel(r, x) result(vpk)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(in) :: x

    vpk = exp(-2*m_e*r*x/hbar_c)*(1 + 1._dp/(2*x**2))*sqrt(x**2 - 1)/x**2
    
end function vacuum_polarization_kernel

real(dp) function inverse_vp_kernel(r, t) result(ivpk)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(in) :: t

    ivpk = vacuum_polarization_kernel(r, 1/t)/(t**2)
    
end function inverse_vp_kernel
    
end module em_nn_potential