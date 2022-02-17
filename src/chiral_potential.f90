!!
!! Chiral Potential
!!
!! Functions from "Minimally nonlocal nucleon-nucleon potentials with chiral two-pion exchange including delta resonances"
!!
!! Phys.Rev. C91 (2015)
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
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

public :: vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9, vf_integral, chiral_integrals, write_all_potential_functions


!!
!> @brief       Interface of chiral_kernels
!!
!! [I'm not sure what to write about this one as I do not really understand what it does. - Ky Putnam]
!!
!! @author      Rodrigo Navarro Pérez
!!

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
!> @brief x as a function of r
!!
!! x = (avg. pion mass)(r in fentometers)/(hbar_c to make quantity unitless) : for use in the chiral potentials
!!
!! @author      Ky Putnam
!!
real(dp) function pion_mass_r(r) result(x)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm

    x = mpi * r / hbar_c

end function pion_mass_r

!!
!> @brief y as a function of r
!!
!! y = (delta-nucleon mass difference))(r in fentometers)/(hbar_c to make quantity unitless) : for use in the chiral potentials
!!
!! @author      Ky Putnam
!!
real(dp) function delta_nucleon_mass_difference_r(r) result(y)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm

    y = delta_nucleon_mass_difference * r / hbar_c

end function delta_nucleon_mass_difference_r

!!
!< @brief       For the vst and vtt functions
!!
!! Y, a function of pion mass and distance for use in the preceeding 
!! OPE LO functions (v_lo_sigmatau and v_lo_ttau)
!!
!! Corresponds to (A3) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function Y_pion(mpi, r) result(Y)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp), intent(in) :: mpi !< pion mass
    real(dp) :: x

    x = pion_mass_r(r)

    Y = gA**2 * mpi**3 * exp(-x) / (12*pi * Fpi**2 * x)

end function Y_pion

!!
!! @        For the vtt function
!!
!! T, a function of distance and pion mass, for use in the preceeding
!! OPE LO functions (v_lo_t_tau, specifically)
!!
!! Corresponds to (A4) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function T_pion(mpi,r) result(T)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp), intent(in) :: mpi !< pion mass
    real(dp) :: x

    x = pion_mass_r(r)
    
    T = Y_pion(mpi,r) * (1 + 3/x + 3/x**2)
    
end function T_pion

!! LEADING ORDER POTENTIALS

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (1) at leading order
!!
!! Corresponds to (A1) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_sigmatau(r) result(vlstau)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: Y_n, Y_p

    Y_n = Y_pion(mpi0,r) !< Y for neutral pion
    Y_p = Y_pion(mpic,r) !< Y for (positively) charged pion

    vlstau = (Y_n + 2*Y_p)/3

end function v_lo_sigmatau

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (2) at leading order
!!
!! Corresponds to (A2) in the appendix
!!
!! @author      Ky Putnam
!!

real(dp) function v_lo_ttau(r) result(vlttau)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: T_n, T_p

    T_n = T_pion(mpi0,r) !< T for neutral pion
    T_p = T_pion(mpic,r) !< T for (positively) charged pion

    vlttau = (T_n + 2*T_p)/3

end function v_lo_ttau

!! NEXT-TO-LEADING ORDER POTENTIALS

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
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x

    x = pion_mass_r(r)

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
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x

    x = pion_mass_r(r)

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
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x

    x = pion_mass_r(r)

    vnt = -gA**4 * mpi * (12*x * bessel_k0(2*x) + (15 + 4*x**2) * bessel_k1(2*x))

end function v_nlo_t

!!
!> @brief       TPE at NLO with a single Delta intermediate state
!!
!! Two pion exchange potential contribution (4) at next leading order
!!
!! Corresponds to (A8) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_c_d(r) result(vncd)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vncd = - gA**2 * hA**2 * exp(-2*x) * (6 + 12*x + 10*x**2 + 4*x**3 + x**4) &
        / ( 6 * pi **2 * r**5 * y * Fpi**4)

end function v_nlo_c_d

!!
!> @brief       TPE at NLO with a single Delta intermediate state
!!
!! Two pion exchange potential contribution (5) at next leading order
!!
!! Corresponds to (A9) in the appendix
!!
!! @author      Ky Putnam
!!

real(dp) function v_nlo_tau_d(r, vf1, vf2, vf5, vf6, vf7) result(vntaud)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vntaud = -hA**2 * ( (5 - 11*gA**2)*vf1 + (x**2 * (12 - 24*gA**2) + y**2 * (12 - 12*gA**2))*vf2 &
        + (-12*y*(2*x**2 + 2*y**2) + 6*gA**2*(4*x**2 + 4*y**2 + 8*x**2 * y**2)/y)*vf5 &
        + (-12*y + 6*gA**2*(4*x**2 + 4*y**2)/y)*vf6 + 6*gA**2*vf7/y ) &
        / (216 * pi**3 * r**5 * Fpi**4)
end function v_nlo_tau_d

!!
!> @brief       TPE at NLO with a single Delta intermediate state
!!
!! Two pion exchange potential contribution (6) at next leading order
!!
!! Corresponds to (A10) in the appendix
!!
!! @author      Ky Putnam
!!

real(dp) function v_nlo_sigma_d(r, vf1, vf2, vf5, vf6, vf7) result(vnsd)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnsd = -gA**2 * hA**2 * (2*vf1 + 8*x**2*vf2 - 16*y*x**2*vf5 - (4*y**2 + 4*x**2)*vf6/y - vf7/y) & 
        / (72 * pi**3 * r**5 * Fpi**4)
end function v_nlo_sigma_d

!!
!> @brief       TPE at NLO with a single Delta intermediate state
!!
!! Two pion exchange potential contribution (7) at next leading order
!!
!! Corresponds to (A11) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigmatau_d(r) result(vnstaud)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnstaud = gA**2 * hA**2 * exp(-2*x) * (1 + x) * (3 + 3*x + x**2) &
        / ( 54 * pi **2 * r**5 * y * Fpi**4)

end function v_nlo_sigmatau_d

!!
!> @brief       TPE at NLO with a single Delta intermediate state
!!
!! Two pion exchange potential contribution (8) at next leading order
!!
!! Corresponds to (A12) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vntd)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vntd = gA**2 * hA**2 * (2*vf1 + 2*(3 + 4*x**2)*vf2 + 6*vf3 -(12*y**2 + 16*y**2*x**2)*vf5/y &
            - (4*y**2 + 4*x**2 + 3)*vf6/y - vf7/y - 3*vf8/y -12*y*vf9) &
            / (144 * pi**3 * r**5 * Fpi**4)

end function v_nlo_t_d

!!
!> @brief       TPE at NLO with a single Delta intermediate state
!!
!! Two pion exchange potential contribution (9) at next leading order
!!
!! Corresponds to (A13) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_ttau_d(r) result(vnttaud)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnttaud = - gA**2 * hA**2 * exp(-2*x) * (1 + x) * (3 + 3*x + 2*x**2) &
        / (54*pi **2 * r**5 * y * Fpi**4)

end function v_nlo_ttau_d

!!
!> @brief       TPE at NLO with 2 Delta intermediate states
!!
!! Two pion exchange potential contribution (10) at next leading order
!!
!! Corresponds to (A14) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_c_2d(r, vf2, vf4, vf5, vf6, vf7) result(vnc2d)
    implicit none
    real(dp), intent(in) :: r, vf2, vf4, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnc2d = -hA**4 * (4*y**2*vf2 + 2*vf4 + (4*x**4 - 12*y**4 - 8*x**2*y**2)*vf5 &
        + (4*x**2 - 4*y**2)*vf6/y + vf7/y) &
        / (108 * pi**3 * r**5 * Fpi**4)

end function v_nlo_c_2d

!!
!> @brief       TPE at NLO with 2 Delta intermediate states
!!
!! Two pion exchange potential contribution (11) at next leading order
!!
!! Corresponds to (A15) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) result(vntau2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf4, vf5, vf6, vf7
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vntau2d = -hA**4 * (11*vf1 + (24*x**2 + 24*y**2)*vf2 + 6*vf4 - 3*(24*x**2*y**2 + 4*x**4 +20*y**4)*vf5/y &
        - 3*(4*x**2 + 12*y**2)*vf6/y - 3*vf7/y) &
        / (1944 * pi**3 * r**5 * Fpi**4)

end function v_nlo_tau_2d

!!
!> @brief       TPE at NLO with 2 Delta intermediate states
!!
!! Two pion exchange potential contribution (12) at next leading order
!!
!! Corresponds to (A16) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7) result(vnsigma2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnsigma2d = -hA**4 * (-6*vf1 - 24*x**2*vf2 + 48*y*x**2*vf5 + (4*x**2 + 12*y**2)*vf6/y + vf7/y) &
        / (1296 * pi**3 * r**5 * Fpi**4)  

end function v_nlo_sigma_2d

!!
!> @brief       TPE at NLO with 2 Delta intermediate states
!!
!! Two pion exchange potential contribution (13) at next leading order
!!
!! Corresponds to (A17) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7) result(vnstau2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7 !< distance at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnstau2d = -hA**4 * (-2*vf1 - 8*x**2*vf2 + 16*y*x**2*vf5 + (4*y**2 - 4*x**2)*vf6/y - vf7/y) &
        / (7776 * pi**3 * r**5 * Fpi**4)

end function v_nlo_sigmatau_2d

!!
!> @brief       TPE at NLO with 2 Delta intermediate states
!!
!! Two pion exchange potential contribution (14) at next leading order
!!
!! Corresponds to (A18) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vnt2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnt2d = hA**4 * (-6*vf1 - 6*(3 + 4*x**2)*vf2 - 18*vf3 + (36*y**2 + 48*x**2*y**2)*vf5/y &
        + (3 + 4*x**2 + 12*y**2)*vf6/y + vf7/y + 3*vf8/y + 36*y*vf9) &
        / (2592 * pi**3 * r**5 * Fpi**4)

end function v_nlo_t_2d

!!
!> @brief       TPE at NLO with 2 Delta intermediate states
!!
!! Two pion exchange potential contribution (15) at next leading order
!!
!! Corresponds to (A19) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vnttau2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x , y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnttau2d = hA**4 * (-2*vf1 - 2*(3 + 4*x**2)*vf2 - 6*vf3 + (12*y**2 + 16*x**2*y**2)*vf5/y &
        + (4*y**2 - 4*x**2 - 3)*vf6/y - vf7/y -3*vf8/y + 12*y*vf9) &
        / (15552 * pi**3 * r**5 * Fpi**4)  

end function v_nlo_ttau_2d

!! NEXT-TO-NEXT-TO LEADING ORDER POTENTIALS

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
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x

    x = pion_mass_r(r)

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
real(dp) function v_n2lo_sigmatau(r) result(vn2stau)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x

    x = pion_mass_r(r)

    vn2stau = 3*gA**2 * c4*exp(-2*x) * (1 + x)*(3 + 3*x + 2*x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_sigmatau

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (3) at N2LO
!!
!! Corresponds to (A22) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_ttau(r) result(vn2ttau)
    implicit none
    real(dp), intent(in) :: r !< distance at which the function will be evaluated, in fm
    real(dp) :: x

    x = pion_mass_r(r)

    vn2ttau = - 3*gA**2 * c4*exp(-2*x) * (1 + x)*(3 + 3*x + 2*x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_ttau

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (4) at N2LO
!!
!! Corresponds to (A23) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_c_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2cd)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2cd = hA**2 * y * ((5*c2-6*c3)*vf1 + ((-24*c1 + 12*c2 - 12*c3)*x**2 + 12*c2*y**2)*vf2 &
            + 6*((8*c1 + 4*c3)*x**4 + (8*c1 - 4*c2 + 4*c3)*x**2*y**2 -4*c2*y**4)*vf5/y &
            + 6*((4*c1 + 4*c3)*x**2 + (-2*c2 + 2*c3)*y**2)*vf6 +6*c3*vf7/y) &
            / (18*pi**3 * r**6 * Fpi**4)

end function v_n2lo_c_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (5) at N2LO
!!
!! Corresponds to (A24) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_tau_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2taud)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2taud = -b38 * hA * y * ((5 - 11*gA**2)*vf1 + (12*x**2 + 12*y**2 - gA**2 * (24*x**2 + 12*y**2))*vf2 &
            + (-24*y*(x**2 + y**2) +24*gA**2 * (x**4 + 2*x**2*y**2 + y**4))*vf5/y &
            + 6*(-12*y  + 24*gA**2 * (x**2 + y**2))*vf6/y + 6*gA**2*vf7/y) &
            / (54*pi**3 * r**6 * Fpi**4)

end function v_n2lo_tau_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (6) at N2LO
!!
!! Corresponds to (A25) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigma_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2sigmad)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2sigmad = -b38 * hA * gA**2 * y * (2*vf1 + 8*x**2*vf2 - 16*x**2*y*vf5 - 4*vf6/y - vf7/y) &
            / (18*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigma_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (7) at N2LO
!!
!! Corresponds to (A26) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2staud)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2staud = -c4* hA**2 * y * (2*vf1 + 8*x**2*vf2 - 16*x**2*y*vf5 - 4*vf6/y - vf7/y) &
            / (108*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigmatau_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (8) at N2LO
!!
!! Corresponds to (A27) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2td)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2td = -c4* hA**2 * y * (2*vf1 + (6+8*x**2)*vf2 + 6*vf3 - (12*y**2 + 16*x**2 * y**2)*vf5/y &
            - (3 + 4*x**2 + 4*y**2)*vf6/y - vf7/y - 3*vf8/y - 12*y*vf9) &
            / (36*pi**3 * r**6 * Fpi**4)

end function v_n2lo_t_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (9) at N2LO
!!
!! Corresponds to (A28) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_ttau_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2ttaud)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2ttaud = c4* hA**2 * y * (2*vf1 + (6+8*x**2)*vf2 + 6*vf3 - (12*y**2 + 16*x**2 * y**2)*vf5/y &
            - (3 + 4*x**2 + 4*y**2)*vf6/y - vf7/y - 3*vf8/y - 12*y*vf9) &
            / (216*pi**3 * r**6 * Fpi**4)

end function v_n2lo_ttau_d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (10) at N2LO
!!
!! Corresponds to (A29) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_c_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) result(vn2c2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf4, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2c2d = -2*b38 * hA**3 * y * (11*vf1 + (24*x**2 + 12*y**2)*vf2 + 6*vf4 - 3*(24*x**2*y**2 + 4*x**2 + 20*y**2)*vf5/y &
            - 3*(4*x**2 + 12*y**2)*vf6/y - 3*vf7/y) &
            / (81*pi**3 * r**6 * Fpi**4)

end function v_n2lo_c_2d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (11) at N2LO
!!
!! Corresponds to (A30) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) result(vn2tau2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf4, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2tau2d = -b38 * hA**3 * y * (11*vf1 + (24*x**2 + 12*y**2)*vf2 + 6*vf4 - 3*(24*x**2*y**2 + 4*x**2 + 20*y**2)*vf5/y &
            - 3*(4*x**2 + 12*y**2)*vf6/y - 3*vf7/y) &
            / (243*pi**3 * r**6 * Fpi**4)

end function v_n2lo_tau_2d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (12) at N2LO
!!
!! Corresponds to (A31) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7) result(vn2sigma2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2sigma2d = -b38 * hA**3 * y * (-6*vf1 - 24*x**2*vf2 + 48*x**2*y*vf5 + (4*x**2 + 12*y**2)*vf6/y + vf7/y) &
            / (162*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigma_2d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (13) at N2LO
!!
!! Corresponds to (A32) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7) result(vn2stau2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf5, vf6, vf7
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2stau2d = -b38 * hA**3 * y * (-6*vf1 - 24*x**2*vf2 + 48*x**2*y*vf5 + (4*x**2 + 12*y**2)*vf6/y + vf7/y) &
            / (972*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigmatau_2d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (13) at N2LO
!!
!! Corresponds to (A33) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2t2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2t2d = b38 * hA**3 * y * (-6*vf1 - 6*(3 + 4*x**2)*vf2 - 24*vf3 + (36*y**2 + 48*x**2*y**2)*vf5/y &
            + (3 + 4*x**2 + 12*y**2)*vf6/y + vf7/y + 3*vf8/y + 36*y*vf9) &
            / (324*pi**3 * r**6 * Fpi**4)

end function v_n2lo_t_2d

!!
!> @brief       N2LO loop correction
!!
!! Loop correction term (13) at N2LO
!!
!! Corresponds to (A33) in the appendix
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2ttau2d)
    implicit none
    real(dp), intent(in) :: r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9
    real(dp) :: x, y

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2ttau2d = b38 * hA**3 * y * (-6*vf1 - 6*(3 + 4*x**2)*vf2 - 24*vf3 + (36*y**2 + 48*x**2*y**2)*vf5/y &
            + (3 + 4*x**2 + 12*y**2)*vf6/y + vf7/y + 3*vf8/y + 36*y*vf9) &
            / (1944*pi**3 * r**6 * Fpi**4)

end function v_n2lo_ttau_2d

!!
!> @brief       Integration function for potentials
!!
!! Function that can be called to calculate potential integrals.
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!

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
!> @brief       Calculates potential integrals
!!
!! Calls the above function (vf_integral) and integrates the 9 potential integrands for a value of
!! r, then sends the results to an array.
!!
!! @author      Ky Putnam
!!

subroutine chiral_integrals(r, mu_integrals)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(out), dimension(1:9) :: mu_integrals

    mu_integrals(1) = vf_integral(vf_1, r)
    mu_integrals(2) = vf_integral(vf_2, r)
    mu_integrals(3) = vf_integral(vf_3, r)
    mu_integrals(4) = vf_integral(vf_4, r)
    mu_integrals(5) = vf_integral(vf_5, r)
    mu_integrals(6) = vf_integral(vf_6, r)
    mu_integrals(7) = vf_integral(vf_7, r)
    mu_integrals(8) = vf_integral(vf_8, r)
    mu_integrals(9) = vf_integral(vf_9, r)

end subroutine

!> LO potential functions, Delta-less
subroutine leading_order_potentials(r, v_lo)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(out), dimension(1:2) :: v_lo
    
    v_lo(1) = v_lo_sigmatau(r)
    v_lo(2) = v_lo_ttau(r)

end subroutine

!> NLO potential functions, Delta-less
subroutine next_to_leading_order_potentials_deltaless(r)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: vntau, vns, vnt

    vntau = v_nlo_tau(r)
    vns = v_nlo_sigma(r)
    vnt = v_nlo_t(r)

end subroutine

!NLO potential functions with 1 Delta intermediate state
subroutine next_to_leading_order_potentials_1delta(r, mu_integrals)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: vncd, vntaud, vnsd, vnstaud, vntd, vnttaud
    real(dp), intent(in), dimension(:) :: mu_integrals
    real(dp) :: vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8, vf9

    vncd = v_nlo_c_d(r)
    vntaud = v_nlo_tau_d(r, vf1, vf2, vf5, vf6, vf7)
    vnsd = v_nlo_sigma_d(r, vf1, vf2, vf5, vf6, vf7)
    vnstaud = v_nlo_sigmatau_d(r)
    vntd = v_nlo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    vnttaud = v_nlo_ttau_d(r)

end subroutine

!NLO potential functions with 2 Delta intermediate states
subroutine next_to_leading_order_potentials_2delta(r, mu_integrals, v_nlo_2delta)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(in), dimension(:) :: mu_integrals
    real(dp), intent(out), dimension(:) :: v_nlo_2delta
    real(dp) :: vnc2d, vntau2d, vnsigma2d, vnstau2d, vnt2d, vnttau2d
    real(dp) :: vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8, vf9

    vnc2d = v_nlo_c_2d(r, vf2, vf4, vf5, vf6, vf7)
    vntau2d = v_nlo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7)
    vnsigma2d = v_nlo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7)
    vnstau2d = v_nlo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7)
    vnt2d = v_nlo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    vnttau2d = v_nlo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)

end subroutine

!N2LO potential functions, Delta-less
subroutine next_to_next_to_leading_order_potentials_deltaless(r)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: vn2c, vn2stau, vn2ttau

    vn2c = v_n2lo_c(r)
    vn2stau = v_n2lo_sigmatau(r)
    vn2ttau = v_n2lo_ttau(r)

end subroutine

!NLO potential functions with 1 Delta intermediate state
subroutine next_to_next_to_leading_order_potentials_1delta(r, mu_integrals)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: vn2cd, vn2taud, vn2sigmad, vn2staud, vn2td, vn2ttaud
    real(dp), intent(in), dimension(:) :: mu_integrals
    real(dp) :: vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8, vf9

    vn2cd = v_n2lo_c_d(r, vf1, vf2, vf5, vf6, vf7)
    vn2taud = v_n2lo_tau_d(r, vf1, vf2, vf5, vf6, vf7)
    vn2sigmad = v_n2lo_sigma_d(r, vf1, vf2, vf5, vf6, vf7)
    vn2staud = v_n2lo_sigmatau_d(r, vf1, vf2, vf5, vf6, vf7)
    vn2td = v_n2lo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    vn2ttaud = v_n2lo_ttau_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)

end subroutine

!NLO potential functions with 1 Delta intermediate state
subroutine next_to_next_to_leading_order_potentials_2delta(r, mu_integrals)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: vn2c2d, vn2tau2d, vn2sigma2d, vn2stau2d, vn2t2d, vn2ttau2d
    real(dp), intent(in), dimension(:) :: mu_integrals

    vn2c2d = v_n2lo_c_2d(r, vf1, vf2, vf4, vf5, vf6, vf7)
    vn2tau2d = v_n2lo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7)
    vn2sigma2d = v_n2lo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7)
    vn2stau2d = v_n2lo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7)
    vn2t2d = v_n2lo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    vn2ttau2d = v_n2lo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    
end subroutine

subroutine unpack_chiral_integrals(r, v_lo, v_nlo_deltaless, .....)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(out),  dimension(1:2) :: v_lo

    real(dp), dimension(1:9) :: mu_integrals

    call chiral_integrals(r, mu_integrals)

    call leading_order_potentials(r, v_lo)
    call next_to_leading_order_potentials_deltaless(r, v_nlo_deltaless)
    call next_to_leading_order_potentials_1delta(r, mu_integrals, v_nlo_1delta)
    call next_to_leading_order_potentials_2delta(r, mu_integrals, v_nlo_2delta)
    call next_to_next_to_leading_order_potentials_deltaless(r, v_n2lo_deltaless)
    call next_to_next_to_leading_order_potentials_1delta(r, mu_integrals, v_n2lo_1delta)
    call next_to_next_to_leading_order_potentials_2delta(r, mu_integrals, v_n2lo_2delta)


end subroutine

subroutine write_all_potential_functions()
    implicit none
    real(dp) :: r, r_max
    real(dp) :: vlstau, vlttau ! LO potentials
    !real(dp) :: vntau, vns, vnt ! NLO delta-less potentials
    integer :: unit1!, unit2, unit3, unit4, unit5, unit6, unit7
    real(dp), dimension(1:2) :: v_lo
    ! all the others

    !call potential subroutines

    !opens data files
    open(newunit = unit1, file = "lo.dat", status = 'unknown')
    ! open(newunit = unit2, file = "nlo.dat", status = 'unknown')
    ! open(newunit = unit3, file = "nlo_d.dat", status = 'unknown')
    ! open(newunit = unit4, file = "nlo_2d.dat", status = 'unknown')
    ! open(newunit = unit5, file = "n2lo.dat", status = 'unknown')
    ! open(newunit = unit6, file = "n2lo_d.dat", status = 'unknown')
    ! open(newunit = unit7, file = "n2lo_2d.dat", status = 'unknown')

    r = 0.1_dp
    r_max = 12._dp
    do
        if (r > r_max) exit
        call unpack_chiral_integrals(r, v_lo, vnlo_deltaless ..... )

        !write data files
        write(unit1,*) r, v_lo
        !write(unit2,*) r, vntau, vns, vnt

        r = r + 0.1_dp
    end do

    !closes data files
    close(unit1)
    ! close(unit2)
    ! close(unit3)
    ! close(unit4)
    ! close(unit5)
    ! close(unit6)
    ! close(unit7)

end subroutine

!!
!> @brief       Calculates Next-to-Leading Order potentials
!!
!! Calls the chiral_integrals subroutine and uses the integrals to calculate the 9 Next-to-Leading Order (NLO) chiral potentials 
!! for a value of r, then sends the results to an array.
!!
!! @author      Ky Putnam
!!
! subroutine next_to_leading_order(r, nlo_potentials)
!     implicit none
!     real(dp), intent(in) :: r
!     real(dp), dimension(1:9) :: mu_integrals
!     real(dp), intent(out), dimension(1:15) :: nlo_potentials

!     real(dp) :: vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8, vf9

!     call chiral_integrals(r, mu_integrals)

! ! assigning mu_integrals shorter names: vfi = mu_integrals(i)
!     vf1 = mu_integrals(1)
!     vf2 = mu_integrals(2)
!     vf3 = mu_integrals(3)
!     vf4 = mu_integrals(4)
!     vf5 = mu_integrals(5)
!     vf6 = mu_integrals(6)
!     vf7 = mu_integrals(7)
!     vf8 = mu_integrals(8)
!     vf9 = mu_integrals(9)

! ! sending next-to-leading order potentials to nlo_potentials array
!     nlo_potentials(1) = v_nlo_tau(r) ! corresponds to A5 in the appendix
!     nlo_potentials(2) = v_nlo_sigma(r) ! corresponds to A6 in the appendix
!     nlo_potentials(3) = v_nlo_t(r) ! corresponds to A7 in the appendix
!     nlo_potentials(4) = v_nlo_c_d(r) ! corresponds to A8 in the appendix
!     nlo_potentials(5) = v_nlo_tau_d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A9 in the appendix
!     nlo_potentials(6) = v_nlo_sigma_d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A10 in the appendix
!     nlo_potentials(7) = v_nlo_sigmatau_d(r) ! corresponds to A11 in the appendix
!     nlo_potentials(8) = v_nlo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A12 in the appendix
!     nlo_potentials(9) = v_nlo_ttau_d(r) ! corresponds to A13 in the appendix
!     nlo_potentials(10) = v_nlo_c_2d(r, vf2, vf4, vf5, vf6, vf7) ! corresponds to A14 in the appendix
!     nlo_potentials(11) = v_nlo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) ! corresponds to A15 in the appendix
!     nlo_potentials(12) = v_nlo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A16 in the appendix
!     nlo_potentials(13) = v_nlo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A17 in the appendix
!     nlo_potentials(14) = v_nlo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A18 in the appendix
!     nlo_potentials(15) = v_nlo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A19 in the appendix
      
! end subroutine

!!
!> @brief       Calculates Next-to-Next-to-Leading Order potentials
!!
!! Calls the chiral_integrals subroutine and uses the integrals to calculate the 11 Next-to-Next-to-Leading Order (N2LO) chiral potentials 
!! for a value of r, then sends the results to an array.
!!
!! @author      Ky Putnam
!!

! subroutine next_to_next_to_leading_order(r, mu_integrals, n2lo_potentials)
!     implicit none
!     real(dp), intent(in) :: r
!     real(dp), intent(out), dimension(1:15) :: n2lo_potentials
!     real(dp), intent(in), dimension(:) :: mu_integrals
!     real(dp) :: vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8, vf9

!     ! call chiral_integrals(r, mu_integrals)

! ! assigning mu_integrals shorter names: vfi = mu_integrals(i)
!     vf1 = mu_integrals(1)
!     vf2 = mu_integrals(2)
!     vf3 = mu_integrals(3)
!     vf4 = mu_integrals(4)
!     vf5 = mu_integrals(5)
!     vf6 = mu_integrals(6)
!     vf7 = mu_integrals(7)
!     vf8 = mu_integrals(8)
!     vf9 = mu_integrals(9)

! ! sending next-to-next-to-leading order potentials to n2lo_potentials array
!     n2lo_potentials(1) = v_n2lo_c(r) ! corresponds to A20 in the appendix
!     n2lo_potentials(2) = v_n2lo_sigmatau(r) ! corresponds to A21 in the appendix
!     n2lo_potentials(3) = v_n2lo_ttau(r) ! corresponds to A22 in the appendix
!     n2lo_potentials(4) = v_n2lo_c_d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A23 in the appendix
!     n2lo_potentials(5) = v_n2lo_tau_d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A24 in the appendix
!     n2lo_potentials(6) = v_n2lo_sigma_d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A25 in the appendix
!     n2lo_potentials(7) = v_n2lo_sigmatau_d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A26 in the appendix
!     n2lo_potentials(8) = v_n2lo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A27 in the appendix
!     n2lo_potentials(9) = v_n2lo_ttau_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A28 in the appendix
!     n2lo_potentials(10) = v_n2lo_c_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) ! corresponds to A29 in the appendix
!     n2lo_potentials(11) = v_n2lo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) ! corresponds to A30 in the appendix
!     n2lo_potentials(12) = v_n2lo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A31 in the appendix
!     n2lo_potentials(13) = v_n2lo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7) ! corresponds to A32 in the appendix
!     n2lo_potentials(14) = v_n2lo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A33 in the appendix
!     n2lo_potentials(15) = v_n2lo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) ! corresponds to A34 in the appendix

! end subroutine

!!
!!
!! POTENTIAL INTEGRANDS
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu:, parametric parameter, r: distance at which the function will be evaluated, in fm
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
    real(dp), intent(in) :: u , r !< u=mu, parametric parameter, r: distance at which the function will be evaluated, in fm
    real(dp) :: x , y

    x = mpi * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop
    y = delta_nucleon_mass_difference * r / hbar_c !< r needs to be replaced with value of r that we are using for the do loop

    vf9 = u * atan(u/(2*y)) * exp(-sqrt(u**2 + 4*x**2))
end function vf_9

end module chiral_potential