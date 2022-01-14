!!
!! Chiral Potential
!!
!! Functions from "Minimally nonlocal nucleon-nucleon potentials with chiral two-pion exchange including delta resonances"
!!
!! Phys.Rev. C91 (2015)
!!
!! @author      Ky Putnam
!!
module chiral_potential
use precisions, only : dp
use constants, only : pi, gA, hA, Fpi=>pion_decay_amplitude, c1, c2, c3, c4, b38=>b3_b8, &
    mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, hbar_c, delta_nucleon_mass_difference
    ! Fpi=2fpi
use special_functions, only : bessel_k0, bessel_k1  
use quadrature, only : booles_quadrature
implicit none

private

public :: vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9, vf_integral


interface

    real(dp) function chiral_kernel(mu, r) result(ck)
        use precisions, only : dp
        implicit none
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: r
    end function chiral_kernel
end interface

contains

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (1) at leading order
!!
!! Corresponds to (A1) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_sigma_tau(r) result(vlst)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: Y_n, Y_p

    Y_n = Y_pion(mpi0,r) !< Y for neutral pion
    Y_p = Y_pion(mpic,r) !< Y for (positively) charged pion

    vlst = (Y_n + 2*Y_p)/3

end function v_lo_sigma_tau

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (2) at leading order
!!
!! Corresponds to (A2) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_t_tau(r) result(vltt)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: T_n, T_p

    T_n = T_pion(mpi0,r) !< T for neutral pion
    T_p = T_pion(mpic,r) !< T for (positively) charged pion

    vltt = (T_n + 2*T_p)/3

end function v_lo_t_tau

!!
!< @brief       For the vst and vtt functions above
!!
!! Y, a function of pion mass and radius for use in the preceeding 
!! OPE LO functions (v_lo_sigma_tau and v_lo_t_tau)
!!
!! Corresponds to (A3) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function Y_pion(mpi, r) result(Y)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp), intent(in) :: mpi !< pion mass
    real(dp) :: x

    x = mpi * r / hbar_c
    
    Y = gA**2 * mpi**3 * exp(-x) / (12*pi * Fpi**2 * x)

end function Y_pion

!!
!! @        For the vtt function above
!!
!! T, a function of radius and pion mass, for use in the preceeding
!! OPE LO functions (v_lo_t_tau, specifically)
!!
!! Corresponds to (A4) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function T_pion(mpi,r) result(T)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp), intent(in) :: mpi !< pion mass
    real(dp) :: x

    x = mpi * r / hbar_c
    
    T = Y_pion(mpi,r) * (1 + 3/x + 3/x**2)
    
end function T_pion

!!
!> @brief       TPE at NLO
!!
!! Two pion exchange potential contribution (1) at next leading order
!!
!! Corresponds to (A5) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_tau(r) result(vntau)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vntau = mpi*(x * (1 + 10*gA**2 - gA**4 * (23 + 4*x**2))*bessel_k0(2*x) &
        + (1 + 2*gA**2 * (5 + 2*x**2) - gA**4 * (23 + 12*x**2))*bessel_k1(2*x)) &
        / (2*pi**3 * r**4 * Fpi**4)
    
end function v_nlo_tau

!!
!> @brief       TPE at NLO
!!
!! Two pion exchange potential contribution (2) at next leading order
!!
!! Corresponds to (A6) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigma(r) result(vns)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vns = gA**4 * mpi * (3*x * bessel_k0(2*x) + (3 + 2*x**2) * bessel_k1(2*x)) &
        / (2*pi**3 * r**4 * Fpi**4)

end function v_nlo_sigma

!!
!> @brief       TPE at NLO
!!
!! Two pion exchange potential contribution (3) at next leading order
!!
!! Corresponds to (A7) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_t(r) result(vnt)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vnt = -gA**4 * mpi * (12*x * bessel_k0(2*x) + (15 + 4*x**2) * bessel_k1(2*x))

end function v_nlo_t

!!
!> @brief       TPE at NLO
!!
!! Two pion exchange potential contribution (4) at next leading order
!!
!! Corresponds to (A8) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_c_d(r) result(vncd)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vncd = - gA**2 * hA**2 * exp(-2*x) * (6 + 12*x + 10*x**2 + 4*x**3 + x**4) &
        / ( 6 * pi **2 * r**5 * y * Fpi**4)

end function v_nlo_c_d

!!
!> @brief       TPE at NLO
!!
!! Two pion exchange potential contribution (7) at next leading order
!!
!! Corresponds to (A11) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_st_d(r) result(vnstd)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vnstd = gA**2 * hA**2 * exp(-2*x) * (1 + x) * (3 + 3*x + x**2) &
        / ( 54 * pi **2 * r**5 * y * Fpi**4)

end function v_nlo_st_d

!!
!> @brief       TPE at NLO
!!
!! Two pion exchange potential contribution (9) at next leading order
!!
!! Corresponds to (A13) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_tt_d(r) result(vnttd)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vnttd = - gA**2 * hA**2 * exp(-2*x) * (1 + x) * (3 + 3*x + 2*x**2) &
        / (54*pi **2 * r**5 * y * Fpi**4)

end function v_nlo_tt_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (1) at N2LO
!!
!! Corresponds to (A20) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_c(r) result(vn2c)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vn2c = 3*gA**2 * exp(-2*x) * (2*c1*x**2 *(1 + x)**2 + c3*(6 + 12*x + 10*x**2 + 4*x**3 + x**4)) &
        / (2*pi**2 * r**6 * Fpi**4)
end function v_n2lo_c

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (2) at N2LO
!!
!! Corresponds to (A21) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigma_tau(r) result(vn2st)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vn2st = 3*gA**2 * c4*exp(-2*x) * (1 + x)*(3 + 3*x + 2*x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_sigma_tau

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (3) at N2LO
!!
!! Corresponds to (A22) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_t_tau(r) result(vn2tt)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vn2tt = - 3*gA**2 * c4*exp(-2*x) * (1 + x)*(3 + 3*x + 2*x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_t_tau

function vf_integral(vf, r) result(i_vf1)
    implicit none
    procedure(chiral_kernel) :: vf
    real(dp), intent(in) :: r
    real(dp) :: i_vf1
    integer, parameter :: n_points = 725
    real(dp), parameter :: mu_max = 30
    real(dp), parameter :: mu_min = 0
    real(dp), dimension(1:n_points) :: f_mu
    real(dp) :: delta_mu, mu
    integer :: i

    delta_mu = (mu_max - mu_min)/(n_points - 1)   
    
    do i=1, n_points
        mu = mu_min + (i-1) * delta_mu
        f_mu(i) = vf(mu, r)
    enddo

    i_vf1 = booles_quadrature(f_mu, delta_mu)
end function vf_integral

!!
!!
!! THE FOLLOWING FUNCTIONS ARE INTEGRANDS OF POTENTIALS, AND ARE IMPLEMENTED SPECIFICALLY FOR PLOTTING
!!
!! Integrands are labeled based on their locations in the "Organizing Integrals" spreadsheet:
!! https://docs.google.com/spreadsheets/d/1dLe5_UPNaTz_hzaVfSETQZh9XZ9dlar5mD3eGOGdOVw/edit?usp=sharing
!!


!! potential integrand function 1: vf1
!!
!! labeled "1" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!!
!! @author      Ky Putnam
!!
real(dp) function vf_1(u , r) result(vf1)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf1 = u**4 * exp(-sqrt(u**2 + 4*x**2)) &
        /sqrt(u**2 + 4*x**2)
end function vf_1

!! potential integrand function 2: vf2
!!
!! labeled "2" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_2(u , r) result(vf2)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf2 = u**2 * exp(-sqrt(u**2 + 4*x**2)) &
        /(sqrt(u**2 + 4*x**2))
end function vf_2

!! potential integrand function 3: vf3
!!
!! labeled "3" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_3(u, r) result(vf3)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf3 = u**2 * exp(-sqrt(u**2 + 4*x**2))
end function vf_3

!! potential integrand function 4: vf4
!!
!! labeled "4" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_4(u , r) result(vf4)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf4 = u**2 * (2*x**2 + u**2 + 2*y**2)**2 * exp(-sqrt(u**2 + 4*x**2)) &
        /(sqrt(u**2 + 4*x**2) * (u**2 +4*y**2))
end function vf_4

!! potential integrand function 5: vf5
!!
!! labeled "5" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_5(u , r) result(vf5)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf5 = u * atan(u/(2*y)) * exp(-sqrt(u**2 + 4*x**2)) &
        /(sqrt(u**2 + 4*x**2))
end function vf_5

!! potential integrand function 6: vf6
!!
!! labeled "6" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_6(u , r) result(vf6)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf6 = u**3 * atan(u/(2*y)) * exp(-sqrt(u**2 + 4*x**2)) &
        /(sqrt(u**2 + 4*x**2))
end function vf_6

!! potential integrand function 7: vf7
!!
!! labeled "7" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_7(u , r) result(vf7)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf7 = u**5 * atan(u/(2*y)) * exp(-sqrt(u**2 + 4*x**2)) &
        /(sqrt(u**2 + 4*x**2))
end function vf_7

!! potential integrand function 8: vf8
!!
!! labeled "8" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_8(u , r) result(vf8)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf8 = u**3 * atan(u/(2*y)) * exp(-sqrt(u**2 + 4*x**2))
end function vf_8

!! potential integrand function 9: vf9
!!
!! labeled "9" in "organizing integrals" spreadsheet
!! currently unspecified value of r
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_9(u , r) result(vf9)
    implicit none
    real(dp), intent(in) :: u , r !< u=mu, parametric parameter, r: radius at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf9 = u * atan(u/(2*y)) * exp(-sqrt(u**2 + 4*x**2))
end function vf_9

end module chiral_potential