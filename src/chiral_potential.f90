!!
!! Chiral Potential
!!
!! Functions from "Minimally nonlocal nucleon-nucleon potentials with chiral two-pion exchange including delta resonances"
!!
!! Phys.Rev. C91 (2015)
!!
!! @author Ky Putnam
!!
module chiral_potential
use precisions, only : dp
use constants, only : pi, gA, hA, Fpi=>pion_decay_amplitude, c1, c2, c3, c4, b38=>b3_b8, &
    mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, hbar_c, delta_nucleon_mass_difference
    ! Fpi=2fpi
use special_functions, only : bessel_k0, bessel_k1   
implicit none

private

public :: v_lo_sigma_tau, v_lo_t_tau

contains

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (1) at leading order
!!
!! Corresponds to (A1) in the appendix
!!
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
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
!! @author Ky Putnam
!!
real(dp) function v_n2lo_t_tau(r) result(vn2tt)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: x

    x = mpi * r / hbar_c

    vn2tt = - 3*gA**2 * c4*exp(-2*x) * (1 + x)*(3 + 3*x + 2*x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_t_tau

end module chiral_potential