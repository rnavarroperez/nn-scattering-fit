!!
!> @brief       long range chiral potentials
!!
!! Module to calculate long range components of potential from Phys. Rev. C 91, 024003(2015).
!! 
!! The potential components are calculated separately then summed. Includes linear sums of
!! 9 integrands that the potential functions are combinations of.
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
module long_range_chiral_potentials
use precisions, only : dp
use constants, only : pi, gA, hA, Fpi=>pion_decay_amplitude, c1, c2, c3, c4, b38=>b3_b8, &
    mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, hbar_c, delta_nucleon_mass_difference
    ! Fpi=2fpi
use special_functions, only : bessel_k0, bessel_k1  
use quadrature, only : booles_quadrature
implicit none

real(dp), parameter :: b_small = 0.045_dp

private

public :: vf_1, vf_2, vf_3, vf_4, vf_5, vf_6, vf_7, vf_8, vf_9, vf_integral, chiral_integrals, &
        calculate_chiral_potentials, long_range_potentials


!!
!> @brief       interface of chiral_kernels
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
!> @brief       subroutine that fills an array with long range potential components
!!
!! Linear combinations of potential components; from pg. 16 of Phys. Rev. C 91, 024003(2015). Corresponds to A35-A42 in the appendix.
!! Output is an array with all long-range potential components.
!!
!! Charge-independent components:
!! \f[
!!      v^{c}_{\mathrm{L}}(r) = v^{\mathrm{NLO}}_{c}(r;\Delta) + v^{\mathrm{NLO}}_{c}(r;2\Delta) + v^{\mathrm{N2LO}}_{c}(r;\mathrm{no }\Delta) + v^{\mathrm{N2LO}}_{c}(r;\Delta) + v^{\mathrm{N2LO}}_{c}(r;2\Delta) \\
!!      v^{\tau}_{\mathrm{L}}(r) = v^{\mathrm{NLO}}_{\tau}(r;\mathrm{no }\Delta) + v^{\mathrm{NLO}}_{\tau}(r;\Delta) + v^{\mathrm{NLO}}_{\tau}(r;2\Delta) + v^{\mathrm{N2LO}}_{\tau}(r;\Delta) + v^{\mathrm{N2LO}}_{\tau}(r;2\Delta) \\
!!      v^{\sigma}_{\mathrm{L}}(r) = v^{\mathrm{NLO}}_{\sigma}(r;\mathrm{no }\Delta) + v^{\mathrm{NLO}}_{\sigma}(r;\Delta) + v^{\mathrm{NLO}}_{\sigma}(r;2\Delta) + v^{\mathrm{NLO}}_{\sigma}(r;\Delta) + v^{\mathrm{NLO}}_{\sigma}(r;2\Delta)
!!      v^{\sigma \tau}_{\mathrm{L}}(r) = v^{\mathrm{LO}}_{\sigma \tau}(r) + v^{\mathrm{NLO}}_{\sigma \tau}(r;\Delta) + v^{\mathrm{NLO}}_{\sigma \tau}(r;2\Delta) + v^{\mathrm{N2LO}}_{\sigma \tau}(r;\mathrm{no }\Delta) + v^{\mathrm{N2LO}}_{\sigma \tau}(r;\Delta) + v^{\mathrm{N2LO}}_{\sigma \tau}(r;2\Delta) \\
!!      v^{t}_{\mathrm{L}}(r) = v^{\mathrm{NLO}}_{t}(r;\mathrm{no }\Delta) + v^{\mathrm{NLO}}_{t}(r;\Delta) + v^{\mathrm{NLO}}_{t}(r;2\Delta) + v^{\mathrm{N2LO}}_{t}(r;\Delta) + v^{\mathrm{N2LO}}_{t}(r;2\Delta) \\
!!      v^{t \tau}_{\mathrm{L}}(r) = v^{\mathrm{LO}}_{t \tau}(r) + v^{\mathrm{NLO}}_{t \tau}(r;\Delta) + v^{\mathrm{NLO}}_{t \tau}(r;2\Delta) + v^{\mathrm{N2LO}}_{t \tau}(r;\mathrm{no }\Delta) + v^{\mathrm{N2LO}}_{t \tau}(r;\Delta) + v^{\mathrm{N2LO}}_{t \tau}(r;2\Delta) \\
!! \f]
!!
!! Charge-dependent components:
!! \f[
!!      v^{\sigma T}_{\mathrm{L}}(r) = \frac{Y_{0}(r) - Y_{+}(r)}{3} \\
!!      v^{t T}_{\mathrm{L}}(r) = \frac{T_{0}(r) - T_{+}(r)}{3}
!! \f]
!!
!! @author      Ky Putnam
!!
subroutine long_range_potentials(r, R_L, a_L, v_long)
    implicit none
    real(dp), dimension(1:19), intent(out) :: v_long
    real(dp), intent(in) :: r                       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L                     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L                     !< parameter of long range regulator (unitless constant)
    real(dp), dimension(1:4) :: v_lo                !< array of leading order potential terms
    real(dp), dimension(1:3) :: v_nlo_deltaless     !< array of next-to-leading order potential terms, deltaless
    real(dp), dimension(1:6) :: v_nlo_1delta        !< array of next-to-leading order potential terms, 1 delta
    real(dp), dimension(1:6) :: v_nlo_2delta        !< array of next-to-leading order potential terms, 2 deltas
    real(dp), dimension(1:3) :: v_n2lo_deltaless    !< array of next-to-next-to-leading order potential terms, deltaless
    real(dp), dimension(1:6) :: v_n2lo_1delta       !< array of next-to-next-to-leading order potential terms, 1 delta
    real(dp), dimension(1:6) :: v_n2lo_2delta       !< array of next-to-next-to-leading order potential terms, 2 deltas
    real(dp) :: c_RL                                !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    call calculate_chiral_potentials(r, R_L, a_L, v_lo, v_nlo_deltaless, v_nlo_1delta, v_nlo_2delta, v_n2lo_deltaless,&
    v_n2lo_1delta, v_n2lo_2delta)

    ! initialize all elements to 0
    v_long = 0 

    ! fill array with charge-independent components
    v_long( 1) = v_nlo_1delta(1) + v_nlo_2delta(1) + v_n2lo_deltaless(1) + v_n2lo_1delta(1) + v_n2lo_2delta(1)    ! v_c
    v_long( 2) = v_nlo_deltaless(1) + v_nlo_1delta(2) + v_nlo_2delta(2) + v_n2lo_1delta(2) + v_n2lo_2delta(2)     ! v_tau
    v_long( 3) = v_nlo_deltaless(2) + v_nlo_1delta(3) + v_nlo_2delta(3) + v_n2lo_1delta(3) + v_n2lo_2delta(3)     ! v_sigma
    v_long( 4) = v_lo(1) + v_nlo_1delta(4) + v_nlo_2delta(4) + v_n2lo_deltaless(2) + v_n2lo_1delta(4) + v_n2lo_2delta(4)    ! v_sigma_tau
    v_long( 5) = v_nlo_deltaless(3) + v_nlo_1delta(5) + v_nlo_2delta(5) + v_n2lo_1delta(5) + v_n2lo_2delta(5)     ! v_t
    v_long( 6) = v_lo(2) + v_nlo_1delta(6) + v_nlo_2delta(6) + v_n2lo_deltaless(3) + v_n2lo_1delta(6) + v_n2lo_2delta(6)    ! v_t_tau

    ! fill array with charge-depenent components
    v_long(16) = v_lo(3) ! v_sigma_T
    v_long(17) = v_lo(4) ! v_t_T

end subroutine long_range_potentials

!!
!> @brief       integration function
!!
!! Function that performs integration for a given input. Uses booles_quadrature for integration.
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
function vf_integral(vf, r) result(i_vf1)
    implicit none
    procedure(chiral_kernel) :: vf
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
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
!> @brief       Integration subroutine for the 9 common potential integrals
!!
!! Subroutine that integrates the 9 integrands that are common to the long range potential components (input) and sends the results to an array (output).
!! Integrands are defined in functions vf_i(u,r) where i=1,2,...9.
!!
!! @author      Ky Putnam
!!
subroutine chiral_integrals(r, mu_integrals)
    implicit none
    real(dp), intent(in) :: r                               !< point at which potential will be evaluated (in fm)
    real(dp), intent(out), dimension(1:9) :: mu_integrals   !< array of the 9 repeating integrated functions

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

!!
!> @brief       LO potential functions
!!
!! Subroutine that calculates delta-less leading order potential functions and stores them in an array (output).
!!
!! @author      Ky Putnam
!!
subroutine lo_potentials(r, R_L, a_L, v_lo)
    implicit none
    real(dp), intent(in) :: r                       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L                     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L                     !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:4) :: v_lo   !< array of LO potential components
    real(dp) :: c_RL                                !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)
    
    ! fill v_lo array
    if (r < b_small) then
        v_lo(1) = v_lo_sigmatau_small_r(r, R_L, a_L)
        v_lo(2) = v_lo_ttau_small_r(r, R_L, a_L)
        v_lo(3) = v_lo_sigmaT_small_r(r, R_L, a_L)
        v_lo(4) = v_lo_tT_small_r(r, R_L, a_L)
    else        
        v_lo(1) = v_lo_sigmatau(r)
        v_lo(2) = v_lo_ttau(r)
        v_lo(3) = v_lo_sigmaT(r)
        v_lo(4) = v_lo_tT(r)

        v_lo = c_RL*v_lo
    endif

end subroutine

!!
!> @brief       NLO potential functions, \f$ \Delta \f$-less
!!
!! Subroutine that calculates \f$ \Delta \f$-less next-to-leading order potential functions and stores them in an array (output).
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
subroutine nlo_potentials_deltaless(r, R_L, a_L, v_nlo_deltaless)
    use precisions, only : sp
    implicit none
    real(dp), intent(in) :: r                                   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L                                 !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L                                 !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:3) :: v_nlo_deltaless    !< array of LO potential components, no \f$ \Delta \f$
    real(dp) :: c_RL                                            !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    !fill v_nlo_deltaless array
    if (r < tiny(1._sp)) then
        v_nlo_deltaless = 0._dp
    elseif (r < b_small) then
        v_nlo_deltaless(1) = v_nlo_tau_small_r(r, R_L, a_L)
        v_nlo_deltaless(2) = v_nlo_sigma_small_r(r, R_L, a_L)
        v_nlo_deltaless(3) = v_nlo_t_small_r(r, R_L, a_L)
    else
        v_nlo_deltaless(1) = v_nlo_tau(r)
        v_nlo_deltaless(2) = v_nlo_sigma(r)
        v_nlo_deltaless(3) = v_nlo_t(r)

        v_nlo_deltaless = c_RL*v_nlo_deltaless
    endif

end subroutine

!!
!> @brief       NLO potential functions, 1 \f$ \Delta \f$
!!
!! Subroutine that calculates next-to-leading order potential functions with one delta intermediate state and stores them in an array.
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
subroutine nlo_potentials_1delta(r, R_L, a_L, mu_integrals, v_nlo_1delta)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:9) :: mu_integrals   !< array of the 9 repeating integrated functions
    real(dp), intent(out), dimension(1:6) :: v_nlo_1delta   !< array of NLO potential components, 1 \f$ \Delta \f$
    real(dp) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: c_RL                !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    !unpack mu integrals
    vf1 = mu_integrals(1)
    vf2 = mu_integrals(2)
    vf3 = mu_integrals(3)
    vf5 = mu_integrals(5)
    vf6 = mu_integrals(6)
    vf7 = mu_integrals(7)
    vf8 = mu_integrals(8)
    vf9 = mu_integrals(9)
    
    !fill v_nlo_1delta array
    if (r < b_small) then
        v_nlo_1delta(1) = v_nlo_c_d_small_r(r, R_L, a_L)
        v_nlo_1delta(2) = v_nlo_tau_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_nlo_1delta(3) = v_nlo_sigma_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_nlo_1delta(4) = v_nlo_sigmatau_d_small_r(r, R_L, a_L)
        v_nlo_1delta(5) = v_nlo_t_d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_nlo_1delta(6) = v_nlo_ttau_d_small_r(r, R_L, a_L)        
    else
        v_nlo_1delta(1) = v_nlo_c_d(r)
        v_nlo_1delta(2) = v_nlo_tau_d(r, vf1, vf2, vf5, vf6, vf7)
        v_nlo_1delta(3) = v_nlo_sigma_d(r, vf1, vf2, vf5, vf6, vf7)
        v_nlo_1delta(4) = v_nlo_sigmatau_d(r)
        v_nlo_1delta(5) = v_nlo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_nlo_1delta(6) = v_nlo_ttau_d(r)

        v_nlo_1delta = c_RL*v_nlo_1delta
    endif

end subroutine

!!
!> @brief       NLO potential functions, 2 \f$ \Delta \f$
!!
!! Subroutine that calculates next-to-leading order potential functions with 2 \f$ \Delta \f$ intermediate states, and stores them in an array.
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
subroutine nlo_potentials_2delta(r, R_L, a_L, mu_integrals, v_nlo_2delta)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:9) :: mu_integrals   !< array of the 9 repeating integrated functions
    real(dp), intent(out), dimension(1:6) :: v_nlo_2delta   !< array of NLO potential components, 2 \f$ \Delta \f$
    real(dp) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: c_RL    !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    !unpack mu_integrals
    vf1 = mu_integrals(1)
    vf2 = mu_integrals(2)
    vf3 = mu_integrals(3)
    vf4 = mu_integrals(4)
    vf5 = mu_integrals(5)
    vf6 = mu_integrals(6)
    vf7 = mu_integrals(7)
    vf8 = mu_integrals(8)
    vf9 = mu_integrals(9)

    !fill v_nlo_2delta array
    if (r < b_small) then
        v_nlo_2delta(1) = v_nlo_c_2d_small_r(r, R_L, a_L, vf2, vf4, vf5, vf6, vf7)
        v_nlo_2delta(2) = v_nlo_tau_2d_small_r(r, R_L, a_L, vf1, vf2, vf4, vf5, vf6, vf7)
        v_nlo_2delta(3) = v_nlo_sigma_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_nlo_2delta(4) = v_nlo_sigmatau_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_nlo_2delta(5) = v_nlo_t_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_nlo_2delta(6) = v_nlo_ttau_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)        
    else
        v_nlo_2delta(1) = v_nlo_c_2d(r, vf2, vf4, vf5, vf6, vf7)
        v_nlo_2delta(2) = v_nlo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7)
        v_nlo_2delta(3) = v_nlo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7)
        v_nlo_2delta(4) = v_nlo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7)
        v_nlo_2delta(5) = v_nlo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_nlo_2delta(6) = v_nlo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)

        v_nlo_2delta = c_RL*v_nlo_2delta
    endif

end subroutine

!!
!> @brief       N2LO potential functions, \f$ \Delta \f$-less
!!
!! Subroutine that calculates next-to-next-to-leading order \f$ \Delta \f$-less potential functions and stores them in an array.
!!
!! @author      Ky Putnam
!!
subroutine n2lo_potentials_deltaless(r, R_L, a_L, v_n2lo_deltaless)
    implicit none
    real(dp), intent(in) :: r                                   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L                                 !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L                                 !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:3) :: v_n2lo_deltaless   !< array of N2LO potential components, no \f$ \Delta \f$
    real(dp) :: c_RL                                            !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    !fill v_n2lo_deltaless array
    if (r < b_small) then
        v_n2lo_deltaless(1) = v_n2lo_c_small_r(r, R_L, a_L)
        v_n2lo_deltaless(2) = v_n2lo_sigmatau_small_r(r, R_L, a_L)
        v_n2lo_deltaless(3) = v_n2lo_ttau_small_r(r, R_L, a_L)
    else        
        v_n2lo_deltaless(1) = v_n2lo_c(r)
        v_n2lo_deltaless(2) = v_n2lo_sigmatau(r)
        v_n2lo_deltaless(3) = v_n2lo_ttau(r)

        v_n2lo_deltaless = c_RL*v_n2lo_deltaless
    endif

end subroutine

!!
!> @brief       N2LO potential functions, 1 \f$ \Delta \f$
!!
!! Subroutine that calculates next-to-next-to-leading order potential functions with 1 \f$ \Delta \f$ intermediate state, and stores them in an array.
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
subroutine n2lo_potentials_1delta(r, R_L, a_L, mu_integrals, v_n2lo_1delta)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:9) :: mu_integrals   !< array of the 9 repeating integrated functions
    real(dp), intent(out), dimension(1:6) :: v_n2lo_1delta  !< array of N2LO potential components, 1 \f$ \Delta \f$
    real(dp) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: c_RL    !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    !unpack mu_integrals
    vf1 = mu_integrals(1)
    vf2 = mu_integrals(2)
    vf3 = mu_integrals(3)
    vf5 = mu_integrals(5)
    vf6 = mu_integrals(6)
    vf7 = mu_integrals(7)
    vf8 = mu_integrals(8)
    vf9 = mu_integrals(9)

    !fill v_n2lo_1delta array
    if (r < b_small) then
        v_n2lo_1delta(1) = v_n2lo_c_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(2) = v_n2lo_tau_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(3) = v_n2lo_sigma_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(4) = v_n2lo_sigmatau_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(5) = v_n2lo_t_d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_n2lo_1delta(6) = v_n2lo_ttau_d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    else        
        v_n2lo_1delta(1) = v_n2lo_c_d(r, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(2) = v_n2lo_tau_d(r, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(3) = v_n2lo_sigma_d(r, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(4) = v_n2lo_sigmatau_d(r, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_1delta(5) = v_n2lo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_n2lo_1delta(6) = v_n2lo_ttau_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)

        v_n2lo_1delta = c_RL*v_n2lo_1delta
    endif

end subroutine

!!
!> @brief       N2LO potential functions, 2 \f$ \Delta \f$
!!
!! Subroutine that calculates next-to-next-to-leading order potential functions with 2 \f$ \Delta \f$ intermediate states, and stores them in an array.
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
subroutine n2lo_potentials_2delta(r, R_L, a_L, mu_integrals, v_n2lo_2delta)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:9) :: mu_integrals   !< array of the 9 repeating integrated functions
    real(dp), intent(out), dimension(1:6) :: v_n2lo_2delta  !< array of N2LO potential components, 2 \f$ \Delta \f$
    real(dp) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: c_RL    !< long range regulator (unitless constant)

    c_RL = long_range_regulator(r, R_L, a_L)

    !unpack mu_integrals
    vf1 = mu_integrals(1)
    vf2 = mu_integrals(2)
    vf3 = mu_integrals(3)
    vf4 = mu_integrals(4)
    vf5 = mu_integrals(5)
    vf6 = mu_integrals(6)
    vf7 = mu_integrals(7)
    vf8 = mu_integrals(8)
    vf9 = mu_integrals(9)

    !fill v_n2lo_2delta array
    if (r < b_small) then
        v_n2lo_2delta(1) = v_n2lo_c_2d_small_r(r, R_L, a_L, vf1, vf2, vf4, vf5, vf6, vf7)
        v_n2lo_2delta(2) = v_n2lo_tau_2d_small_r(r, R_L, a_L, vf1, vf2, vf4, vf5, vf6, vf7)
        v_n2lo_2delta(3) = v_n2lo_sigma_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_2delta(4) = v_n2lo_sigmatau_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_2delta(5) = v_n2lo_t_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_n2lo_2delta(6) = v_n2lo_ttau_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
    else        
        v_n2lo_2delta(1) = v_n2lo_c_2d(r, vf1, vf2, vf4, vf5, vf6, vf7)
        v_n2lo_2delta(2) = v_n2lo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7)
        v_n2lo_2delta(3) = v_n2lo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_2delta(4) = v_n2lo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7)
        v_n2lo_2delta(5) = v_n2lo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)
        v_n2lo_2delta(6) = v_n2lo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9)

        v_n2lo_2delta = c_RL*v_n2lo_2delta
    endif
    
end subroutine

!!
!> @brief       calls subroutines that calculate long-range potential components
!!
!! Calls subroutines that calculate long-range potential components.
!!
!! @author      Ky Putnam
!!
subroutine calculate_chiral_potentials(r, R_L, a_L, v_lo, v_nlo_deltaless, v_nlo_1delta, v_nlo_2delta, v_n2lo_deltaless,&
         v_n2lo_1delta, v_n2lo_2delta)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(out), dimension(1:4) :: v_lo               !< array of LO potential components
    real(dp), intent(out), dimension(1:3) :: v_nlo_deltaless    !< array of NLO potential components (no \f$ \Delta \f$)
    real(dp), intent(out), dimension(1:6) :: v_nlo_1delta       !< array of NLO potential components (1 \f$ \Delta \f$)
    real(dp), intent(out), dimension(1:6) :: v_nlo_2delta       !< array of NLO potential components (2 \f$ \Delta \f$)
    real(dp), intent(out), dimension(1:3) :: v_n2lo_deltaless   !< array of N2LO potential components (no \f$ \Delta \f$)
    real(dp), intent(out), dimension(1:6) :: v_n2lo_1delta      !< array of N2LO potential components (1 \f$ \Delta \f$)
    real(dp), intent(out), dimension(1:6) :: v_n2lo_2delta      !< array of N2LO potential components (2 \f$ \Delta \f$)
    real(dp), dimension(1:9) :: mu_integrals                    !< array of the 9 repeating integrated functions

    call chiral_integrals(r, mu_integrals)

    ! call subroutines containing LO, NLO, N2LO potential components
    call lo_potentials(r, R_L, a_L, v_lo)
    call nlo_potentials_deltaless(r, R_L, a_L, v_nlo_deltaless)
    call nlo_potentials_1delta(r, R_L, a_L, mu_integrals, v_nlo_1delta)
    call nlo_potentials_2delta(r, R_L, a_L, mu_integrals, v_nlo_2delta)
    call n2lo_potentials_deltaless(r, R_L, a_L, v_n2lo_deltaless)
    call n2lo_potentials_1delta(r, R_L, a_L, mu_integrals, v_n2lo_1delta)
    call n2lo_potentials_2delta(r, R_L, a_L, mu_integrals, v_n2lo_2delta)
    
end subroutine calculate_chiral_potentials

!!
!> @brief       long range regulator for potentials
!!
!! \f[ C_{R_{L}}(r)= 1-\frac{1}{(r/R_{L})^{6}e^{(r-R_{L})/a_{L}+1}} \f]
!!
!! @author      Ky Putnam
!!
real(dp) function long_range_regulator(r, R_L, a_L) result(c_RL)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)

    c_RL = 1 - 1/((r/R_L)**6 * exp((r-R_L)/a_L) + 1)

end function long_range_regulator

!!
!> @brief \f$ x = m_{\pi}r \f$ (unitless), average pion mass
!!
!! Calculates a value x(r), where r is the point at which a potential is being evaluated and \f$ m_{\pi} \f$ is the average pion mass.
!!
!! @author      Ky Putnam
!!
real(dp) function pion_mass_r(r) result(x)
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)

    x = mpi * r / hbar_c

end function pion_mass_r

!!
!> @brief \f$ x = m_{\pi}r \f$ (unitless), average pion mass
!!
!! Calculates a value x(r), where r is the point at which a potential is being evaluated and \f$ m_{\pi} \f$ is pion mass (type depends on input).
!!
!! @author      Ky Putnam
!!
real(dp) function variable_pion_mass_r(r, mpi) result(x)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: mpi     !< pion mass (in MeV)

    x = mpi * r / hbar_c

end function variable_pion_mass_r

!!
!> @brief   \f$ y = \Delta M r \f$ (unitless)
!!
!! \f$ y = \Delta \f$ is the \f$ \Delta \f$-nucleon mass difference.
!!
!! @author      Ky Putnam
!!
real(dp) function delta_nucleon_mass_difference_r(r) result(y)
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)

    y = delta_nucleon_mass_difference * r / hbar_c

end function delta_nucleon_mass_difference_r

!!
!< @brief  input for leading order (OPE) function (MeV), \f$ v^{LO}_{\sigma \tau} (r)\f$
!!
!! Corresponds to (A3) in the appendix
!!
!! \f[     Y_{\alpha}(r) = \frac{g_{A}^{2}}{12\pi} \frac{m^{3}_{\pi_{\alpha}}}{F^{2}_{\pi}} \frac{e^{-x_{\alpha}}}{x_{\alpha}} \f]
!!
!! @author      Ky Putnam
!!
real(dp) function Y_pion(mpi, r) result(Y)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: mpi     !< pion mass (in MeV)
    real(dp) :: x                   !< variable pion mass (unitless)

    x = variable_pion_mass_r(r,mpi) ! depends on which mass of the pion is received

    Y = gA**2 * mpi**3 * exp(-x) / (12*pi * Fpi**2 * x)

end function Y_pion

!!
!< @brief  inputs for leading order (OPE) function, \f$ v^{LO}_{\sigma \tau} (r)\f$, modified to deal with small r
!!
!! Corresponds to (A3) in the appendix
!!
!! \f[         Y_{\alpha}(r) = \frac{g^{2}_{A} m^{2}_{\pi_{\alpha}} e^{-x_{\alpha}}}
!!    {12\pi F^{2}_{\pi} } e^{-R_{L}/a^{6}_{L}}(1 + \frac{r}{a_{L}} + \frac{(r/a_{L})^{2}}{2} + \frac{(r/a_{L})^{3}}{6}) r^{5} \f]
!!
!! @author      Rodrigo Navarro Pérez
!!
real(dp) function Y_pion_small_r(mpi, r, R_L, a_L) result(Y)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: mpi     !< pion mass (in MeV)
    real(dp) :: x                   !< variable pion mass (unitless)

    x = variable_pion_mass_r(r,mpi) ! depends on which mass of the pion is received (unitless)

    Y = gA**2*mpi**2*exp(-x)*hbar_c/(12*pi*Fpi**2)*&
    (exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)*r**5 

end function Y_pion_small_r

!!
!< @brief       inputs for leading order (OPE) function (MeV), \f$ v^{LO}_{t \tau} (r)\f$
!!
!! Corresponds to (A4) in the appendix
!!
!! \f[ T_{\alpha}(r) = Y_{\alpha}\big(1+\frac{3}{x_{\alpha}} + \frac{3}{x^{3}_{\alpha}}\big) \f]
!!
!! @author      Ky Putnam
!!
real(dp) function T_pion(mpi,r) result(T)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: mpi     !< pion mass (in MeV)
    real(dp) :: x                   !< variable pion mass (unitless)

    x = variable_pion_mass_r(r, mpi)
    
    T = Y_pion(mpi,r) * (1 + 3/x + 3/x**2)
    
end function T_pion


!!
!< @brief       inputs for leading order (OPE) function, \f$ v^{LO}_{t \tau} (r)\f$, modified to deal with small r
!!
!! Corresponds to (A4) in the appendix
!!
!! \f[     T_{\alpha}(r) = \frac{g^{2}_{A} e^{-x_{\alpha}}}
!!    {12\pi F^{2}_{\pi} } (x^{2}_{\alpha} + 3x + 3) e^{-R_{L}/a^{6}_{L}}(1 + \frac{r}{a_{L}} + \frac{(r/a_{L})^{2}}{2} + \frac{(r/a_{L})^{3}}{6}) r^{3} \f]
!!
!! @author      Rodrigo Navarro Pérez
!!
real(dp) function T_pion_small_r(mpi, r, R_L, a_L) result(T)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: mpi     !< pion mass (in MeV)
    real(dp) :: x                   !< variable pion mass (unitless)

    x = variable_pion_mass_r(r, mpi)
    
    T = gA**2*exp(-x)*(hbar_c)**3/(12*pi*Fpi**2)*(x**2 + 3*x + 3)*&
    (exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)*r**3 
    
end function T_pion_small_r

!! LEADING ORDER POTENTIALS

!!
!> @brief       OPE at LO, \f$ v^{LO}_{\sigma \tau} (r) \f$
!!
!! One pion exchange potential contribution (1) at leading order
!! Corresponds to (A1) in the appendix
!!
!! \f[ v^{LO}_{\sigma \tau} (r) = \frac{Y_{0}(r) + 2Y_{+}(r)}{3}\f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_sigmatau(r) result(vlstau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: Y_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: Y_p                 !< function of positively charged pion and r (MeV)

    Y_n = Y_pion(mpi0,r)
    Y_p = Y_pion(mpic,r)

    vlstau = (Y_n + 2*Y_p)/3

end function v_lo_sigmatau

!!
!> @brief       OPE at LO, \f$ v^{LO}_{\sigma \tau} (r) \f$, modified to deal with small r
!!
!! One pion exchange potential contribution (1) at leading order
!! Corresponds to (A1) in the appendix
!!
!! \f[ v^{LO}_{\sigma \tau} (r) = \frac{Y_{0}(r, R_{L}, a_{L}) + 2Y_{+}(r, R_{L}, a_{L})}{3}\f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_lo_sigmatau_small_r(r, R_L, a_L) result(vlstau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: Y_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: Y_p                 !< function of positively charged pion and r (MeV)

    Y_n = Y_pion_small_r(mpi0, r, R_L, a_L)
    Y_p = Y_pion_small_r(mpic, r, R_L, a_L)

    vlstau = (Y_n + 2*Y_p)/3

end function v_lo_sigmatau_small_r

!!
!> @brief       OPE at LO, \f$ v^{LO}_{t \tau} (r) \f$
!!
!! One pion exchange potential contribution (2) at leading order
!! Corresponds to (A2) in the appendix
!!
!! \f[ v^{LO}_{t \tau} (r) = \frac{T_{0}(r, R_{L}, a_{L}) + 2T_{+}(r, R_{L}, a_{L})}{3}\f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_ttau(r) result(vlttau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: T_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: T_p                 !< function of positively charged pion and r (MeV)


    T_n = T_pion(mpi0,r)
    T_p = T_pion(mpic,r)

    vlttau = (T_n + 2*T_p)/3

end function v_lo_ttau

!!
!> @brief       OPE at LO, \f$ v^{LO}_{t \tau} (r) \f$, modified to deal with small r
!!
!! Corresponds to (A2) in the appendix
!!
!! One pion exchange potential contribution (2) at leading order
!! \f[ v^{LO}_{t \tau} (r) = \frac{T_{0}(r, R_{L}, a_{L}) + 2T_{+}(r, R_{L}, a_{L})}{3} \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_lo_ttau_small_r(r, R_L, a_L) result(vlttau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: T_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: T_p                 !< function of positively charged pion and r (MeV)

    T_n = T_pion_small_r(mpi0, r, R_L, a_L)
    T_p = T_pion_small_r(mpic, r, R_L, a_L)

    vlttau = (T_n + 2*T_p)/3

end function v_lo_ttau_small_r

!!
!> @brief       charge-dependent term of long-range potential, \f$ v^{\sigma T}_{L} (r) \f$
!!
!! Corresponds to (A41) in the appendix
!!
!! \f[ v^{\sigma T}_{L} (r) = \frac{Y_{0}(r) - 2Y_{+}(r)}{3}\f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_sigmaT(r) result(vlstau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: Y_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: Y_p                 !< function of positively charged pion and r (MeV)

    Y_n = Y_pion(mpi0,r)
    Y_p = Y_pion(mpic,r)

    vlstau = (Y_n - Y_p)/3

end function v_lo_sigmaT

!!
!> @brief       charge-dependent term of long-range potential, \f$ v^{\sigma T}_{L} (r) \f$, modified to deal with small r
!!
!! Corresponds to (A41) in the appendix
!!
!! \f[ v_{\sigma T}^{LO} (r) = \frac{Y_{0}(r, R_{L}, a_{L}) - Y_{+}(r, R_{L}, a_{L})}{3} \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_lo_sigmaT_small_r(r, R_L, a_L) result(vlstau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: Y_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: Y_p                 !< function of positively charged pion and r (MeV)

    Y_n = Y_pion_small_r(mpi0, r, R_L, a_L)
    Y_p = Y_pion_small_r(mpic, r, R_L, a_L)

    vlstau = (Y_n - Y_p)/3

end function v_lo_sigmaT_small_r

!!
!> @brief       charge-dependent term of long-range potential, \f$ v^{t T}_{L} (r) \f$, modified to deal with small r
!!
!! Corresponds to (A42) in the appendix
!!
!! \f[ v^{t T}_{LO} (r) = \frac{T_{0}(r, R_{L}, a_{L}) - 2T_{+}(r, R_{L}, a_{L})}{3}\f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_lo_tT(r) result(vlttau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: T_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: T_p                 !< function of positively charged pion and r (MeV)

    T_n = T_pion(mpi0,r)
    T_p = T_pion(mpic,r)

    vlttau = (T_n - T_p)/3

end function v_lo_tT

!!
!> @brief       charge-dependent term of long-range potential, \f$ v^{t T}_{L} (r) \f$, modified to deal with small r
!!
!! Corresponds to (A42) in the appendix
!!
!! \f[ v_{t T}^{LO} (r) = \frac{T_{0}(r, R_{L}, a_{L}) - T_{+}(r, R_{L}, a_{L})}{3} \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_lo_tT_small_r(r, R_L, a_L) result(vlttau)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: T_n                 !< function of neutral charged pion and r (MeV)
    real(dp) :: T_p                 !< function of positively charged pion and r (MeV)

    T_n = T_pion_small_r(mpi0, r, R_L, a_L)
    T_p = T_pion_small_r(mpic, r, R_L, a_L)

    vlttau = (T_n - T_p)/3

end function v_lo_tT_small_r

!! NEXT-TO-LEADING ORDER POTENTIALS

!> @brief       TPE at NLO: \f$ v_{\tau}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order, with no \$[ \Delta \$] intermediate state.
!! Corresponds to (A5) in the appendix
!!
!! \f[ v_{\tau}^{\mathrm{NLO}}(r) = \frac{m_\pi \hbar_c^4}{8 \pi^3 r^4 F_\pi^4} \left[ x (1 + 10g_A^2 - g_A^4 (23 + 4x^2))K_0(2x) + (1 + 2g_A^2 (5 + 2x^2) - g_A^4 (23 + 12x^2))K_1(2x) \right] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_tau(r) result(vntau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    vntau = mpi * hbar_c**4 * (x * (1 + 10*gA**2 - gA**4 * (23 + 4*x**2))*bessel_k0(2*x) &
        + (1 + 2*gA**2 * (5 + 2*x**2) - gA**4 * (23 + 12*x**2))*bessel_k1(2*x)) &
        / (8*pi**3 * r**4 * Fpi**4)
    
end function v_nlo_tau

!!
!> @brief       TPE at NLO: \f$ v_{\tau \mathrm{small r}}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution (1) at next leading order
!! Corresponds to (A5) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\tau} = \frac{m_{\pi} x }{8 \pi^3 F^{4}_{\pi} R_{L}^6} (1 + g^{2}_{A} - g^{4}_{A} (23 + 4x^{2})K_{0}(2x) + (1 + 2g^{2}_{A}(5 + 2x^{2}) - g^{4}_{A}(23 + 12x^{2}))K_{1}(2x))
!!    e^{-R_{L}/a_{L}}(1 + \frac{r}{a_{L}} + \frac{(r/a_{L})^{2}}{2} + \frac{(r/a_{L})^{3}}{6}) r^{2} \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_tau_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = mpi*hbar_c**4*(x*(1 + 10*gA**2 - gA**4 * (23 + 4*x**2))*bessel_k0(2*x) &
        + (1 + 2*gA**2*(5 + 2*x**2) - gA**4*(23 + 12*x**2))*bessel_k1(2*x)) &
        /(8*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)*r**2
    
end function v_nlo_tau_small_r

!> @brief       TPE at NLO: \f$ v_{\sigma}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order, with no \$[ \Delta \$] intermediate state.
!! Corresponds to (A6) in the appendix
!!
!! \f[ v_{\sigma}^{\mathrm{NLO}}(r) = \frac{g_A^4 m_\pi}{2\pi^3 r^4 F_\pi^4}\left(3x K_0(2x) + (3+2x^2)K_1(2x)\right) \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigma(r) result(vns)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    vns = gA**4 * mpi * hbar_c**4 * (3*x * bessel_k0(2*x) + (3 + 2*x**2) * bessel_k1(2*x)) &
        / (2*pi**3 * r**4 * Fpi**4)

end function v_nlo_sigma

!!
!> @brief       TPE at NLO: \f$ v_{\sigma \mathrm{small r}}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A6) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\sigma} = \frac{g_A^4 m_\pi}{2\pi^3 F_\pi^4} (3xK_0(2x)+(3+2x^2)K_1(2x))
!!  \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 +
!!  \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) r^2 \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_sigma_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = gA**4*mpi*hbar_c**4*(3*x*bessel_k0(2*x) + (3 + 2*x**2)*bessel_k1(2*x)) &
        / (2*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)*r**2

end function v_nlo_sigma_small_r

!> @brief       TPE at NLO: \f$ v_{t}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order, with no \$[ \Delta \$] intermediate state.
!! Corresponds to (A7) in the appendix
!!
!! \f[ v_{t}^{\mathrm{NLO}}(r) = -\frac{g_A^4 m_\pi}{8\pi^3 r^4 F_\pi^4} \left(12x K_0(2x) + (15+4x^2) K_1(2x)\right) \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_t(r) result(vnt)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    vnt = -gA**4 * mpi * hbar_c**4 * (12*x * bessel_k0(2*x) + (15 + 4*x**2) * bessel_k1(2*x))&
    / (8*pi**3 * r**4 * Fpi**4)

end function v_nlo_t


!!
!> @brief       TPE at NLO: \f$ v_{t \mathrm{small r}}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution (1) at next leading order
!! Corresponds to (A7) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{t} = -\frac{g_A^4 m_\pi}{8\pi^3 F_\pi^4} \left(12xK_0(2x) + (15+4x^2)K_1(2x)\right) 
!!  \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!  \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) r^2 \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_t_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = -gA**4*mpi*hbar_c**4*(12*x*bessel_k0(2*x) + (15 + 4*x**2)*bessel_k1(2*x)) &
    /(8*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)*r**2

end function v_nlo_t_small_r

!> @brief       TPE at NLO: \f$ v_{c}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order, with 1 \$[ \Delta \$] intermediate state.
!! Corresponds to (A8) in the appendix
!!
!! \f[ v_{c}^{\mathrm{NLO}}(r) = - \frac{g_A^2 h_A^2}{6\pi^2 r^5 y F_\pi^4} e^{-2x} \left(6 + 12x + 10x^2 + 4x^3 + x^4\right) \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_c_d(r) result(vncd)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vncd = - gA**2 * hA**2 * hbar_c**5 * exp(-2*x) * (6 + 12*x + 10*x**2 + 4*x**3 + x**4) &
        / ( 6 * pi **2 * r**5 * y * Fpi**4)

end function v_nlo_c_d

!!
!> @brief       TPE at NLO: \f$ v_{c \mathrm{small r}}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A8) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{c} = -\frac{e^{-\frac{R_L}{a_L}}}{R_L^6}\left(1+\frac{r}{a_L}+
!!  \frac{\left(\frac{r}{a_L}\right)^2}{2}+\frac{\left(\frac{r}{a_L}\right)^3}{6}\right)
!!  \frac{g_A^2 h_A^2 e^{-2x}(6+12x+10x^2+4x^3+x^4)}{6\pi^2\Delta M F_\pi^4} \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_c_d_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = -(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)*gA**2*hA**2*hbar_c**6*exp(-2*x)*&
        (6 + 12*x + 10*x**2 + 4*x**3 + x**4)/(6*pi**2*delta_nucleon_mass_difference*Fpi**4)
    
end function v_nlo_c_d_small_r

!> @brief       TPE at NLO: \f$ v_{\tau}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order, with 1 \$[ \Delta \$] intermediate state.
!! Corresponds to (A9) in the appendix
!!
!! \f[ v_{\tau}^{\mathrm{NLO}}(r;\Delta) =  -\frac{1}{216\pi^{3}r^{5}} \frac{h_{A}^{2}}{F_{\pi}^{4}}
!!    \big[ \left(5-11g_{A}^{2}\right)\mathrm{vf1} +  \left(x^{2}\left(12-24g_{A}^{2}\right) + y^{2}\left(12-12g_{A}^{2}\right)\right) \mathrm{vf2}   \\
!!    + \big(-12y\left(2x^{2}+2y^{2}\right) + \frac{6}{y}g_{A}^{2}\left(4x^{2} + 4y^{2} + 8x^{2}y^{2}\right) \big) \mathrm{vf5} + \big(-12y + \frac{6}{y}g_{A}^{2}\left(4x{2} + 4y^{2}\right)\big) \mathrm{vf6} \\
!!    + \big(\frac{6}{y}g_{A}^{2}\big)\mathrm{vf7} \big] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_tau_d(r, vf1, vf2, vf5, vf6, vf7) result(vntaud)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vntaud = -hA**2 * hbar_c**5 * ( (5 - 11*gA**2)*vf1 + (12*x**2 + 12*y**2 - gA**2*(24*x**2 + 12*y**2))*vf2 &
        + (-12*y*(2*x**2 + 2*y**2) + 6*gA**2*(4*x**4 + 4*y**4 + 8*x**2 * y**2)/y)*vf5 &
        + (-12*y + 6*gA**2*(4*x**2 + 4*y**2)/y)*vf6 + 6*gA**2*vf7/y ) &
        / (216 * pi**3 * r**5 * Fpi**4)
end function v_nlo_tau_d

!!
!> @brief       TPE at NLO: \f$ v_{\tau \small r}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A9) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\tau} = -\frac{h_A^2}{216\pi^3r^5F_\pi^4}\left((5-11g_A^2)y\mathrm{vf1}+(12x^2+12y^2-g_A^2(24x^2+12y^2))y\mathrm{vf2}+ \\
!!    (-12y(2x^2+2y^2)+6g_A^2(4x^4+4y^4+8x^2y^2)/y)\mathrm{vf5}+(-12y+6g_A^2(4x^2+4y^2))\mathrm{vf6}+6g_A^2\mathrm{vf7}\right)
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_tau_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = -hA**2*hbar_c**6*( (5-11*gA**2)*y*vf1 + (12*x**2 + 12*y**2 - gA**2*(24*x**2 + 12*y**2))*y*vf2 &
        + (-12*y**2*(2*x**2 + 2*y**2) + 6*gA**2*(4*x**4 + 4*y**4 + 8*x**2 * y**2) )*vf5 &
        + (-12*y**2 + 6*gA**2*(4*x**2 + 4*y**2))*vf6 + 6*gA**2*vf7 ) &
        / (216*pi**3*Fpi**4*delta_nucleon_mass_difference)*(exp(-R_L/a_L)/R_L**6)&
        * (1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)
end function v_nlo_tau_d_small_r

!> @brief       TPE at NLO: \f$ v_{\sigma}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order, with 1 \$[ \Delta \$] intermediate state.
!! Corresponds to (A10) in the appendix
!!
!! \f[ v_{\sigma}^{\mathrm{NLO}}(r;\Delta) = &
!!    -\frac{1}{72\pi^{3}r^{5}} \frac{g_{A}^{2} h_{A}^{2}}{F_{\pi}^{4}}
!!    \big[
!!    2 \mathrm{vf1} + 8x^{2} \cdot \mathrm{vf2} - (16x^{2}y) \mathrm{vf5}\\
!!    & - \frac{1}{y} \left(4y^{2}+4x^{2}\right)\mathrm{vf6} -\frac{1}{y}\mathrm{vf7}
!!    \big] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigma_d(r, vf1, vf2, vf5, vf6, vf7) result(vnsd)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnsd = -gA**2 * hA**2 * hbar_c**5 * (2*vf1 + 8*x**2*vf2 - (16*x**2*y**2*vf5 + (4*x**2 + 4*y**2)*vf6 + vf7)/y) & 
        / (72 * pi**3 * r**5 * Fpi**4)
end function v_nlo_sigma_d

!!
!> @brief       TPE at NLO: \f$ v_{\sigma \mathrm{small r}}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution (1) at next leading order
!! Corresponds to (A10) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\sigma} = -\frac{g_A^2 h_A^2}{72\pi^3 F_\pi^4 \Delta M}(2y\mathrm{vf1} + 8yx^2\mathrm{vf2}  -16x^2y^2\mathrm{vf5} - (4x^2 + 4y^2)\mathrm{vf6} - \mathrm{vf7})
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_sigma_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = -gA**2*hA**2*hbar_c**6*(2*y*vf1 + 8*y*x**2*vf2 - 16*x**2*y**2*vf5 - (4*x**2 + 4*y**2)*vf6 - vf7) & 
        / (72*pi**3*delta_nucleon_mass_difference*Fpi**4)*(exp(-R_L/a_L)/R_L**6)&
        * (1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)
end function v_nlo_sigma_d_small_r

!> @brief       TPE at NLO: \f$ v_{\sigma \tau}^{\mathrm{N2LO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order, with 1 \$[ \Delta \$] intermediate state.
!! Corresponds to (A11) in the appendix
!!
!! \f[ v_{\sigma \tau}^{\mathrm{N2LO}}(r) = \frac{g_A^2 h_A^2}{54 \pi^2 r^5 y F_\pi^4} e^{-2x} (1+x)(3+3x+x^2) \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigmatau_d(r) result(vnstaud)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnstaud = gA**2 * hA**2 * hbar_c**5 * exp(-2*x) * (1 + x) * (3 + 3*x + x**2) &
        / ( 54 * pi **2 * r**5 * y * Fpi**4)

end function v_nlo_sigmatau_d

!!
!> @brief       TPE at NLO: \f$ v_{\sigma \tau \mathrm{small r}}^{\mathrm{N2LO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A11) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\sigma \tau} = \frac{g_A^2 h_A^2}{54\pi^2 F_\pi^4 \Delta M} e^{-2x} (1 + x) (3 + 3x + x^2)
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 +
!!    \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_sigmatau_d_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = gA**2*hA**2*hbar_c**6*exp(-2*x)*(1 + x)*(3 + 3*x + x**2) &
        /(54*pi**2*delta_nucleon_mass_difference*Fpi**4)*(exp(-R_L/a_L)/R_L**6)&
        * (1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_sigmatau_d_small_r

!> @brief       TPE at NLO: \f$ v_{t}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order, with 1 \$[ \Delta \$] intermediate state.
!! Corresponds to (A12) in the appendix
!!
!! \f[ v_{t}^{\mathrm{NLO}}(r;\Delta) = &
!!    \frac{1}{144\pi^{3}r^{5}} \frac{g_{A}^{2} h_{A}^{2}}{F_{\pi}^{4}}
!!    \big[
!!    2 \mathrm{vf1} + (6+8x^{2})\mathrm{vf2} 
!!    + 6\mathrm{vf3} - (12y + 16yx^{2}) \mathrm{vf5}\\
!!    & - \frac{1}{y} \left(4y^{2}+4x^{2}+3\right)\mathrm{vf6} -\frac{1}{y}\mathrm{vf7} -\frac{3}{y}\mathrm{vf8} -12y \cdot \mathrm{vf9}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vntd)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vntd = gA**2 * hA**2 * hbar_c**5 * (2*vf1 + 2*(3 + 4*x**2)*vf2 + 6*vf3 -(12*y**2 + 16*y**2*x**2)*vf5/y &
            - (4*y**2 + 4*x**2 + 3)*vf6/y - vf7/y - 3*vf8/y -12*y*vf9) &
            / (144 * pi**3 * r**5 * Fpi**4)

end function v_nlo_t_d

!!
!> @brief       TPE at NLO: \f$ v_{t \mathrm{small r}}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A12) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{t} = \frac{g_A^2 h_A^2}{144\pi^3 F_\pi^4 \Delta M}(2y\mathrm{vf1} + 2(3+4x^2)y\mathrm{vf2} + 
!!    6y\mathrm{vf3} - (12y^2 + 16y^2x^2)\mathrm{vf5} - (4y^2 +4x^2 +3)\mathrm{vf6} -\mathrm{vf7} -3\mathrm{vf8} - 
!!    12y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!    \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_t_d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = gA**2*hA**2*hbar_c**6*(2*y*vf1 + 2*(3 + 4*x**2)*y*vf2 + 6*y*vf3 -(12*y**2 + 16*y**2*x**2)*vf5 &
        - (4*y**2 + 4*x**2 + 3)*vf6 - vf7 - 3*vf8 -12*y**2*vf9) &
        / (144*pi**3*delta_nucleon_mass_difference*Fpi**4)*(exp(-R_L/a_L)/R_L**6)&
        * (1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_t_d_small_r

!> @brief       TPE at NLO: \f$ v_{t \tau}^{\mathrm{NLO}}(r) \f$
!!
!! Two pion exchange potential contribution at next leading order, with 1 \$[ \Delta \$] intermediate state.
!! Corresponds to (A13) in the appendix
!!
!! \f[ v_{t \tau}^{\mathrm{N2LO}}(r) = -\frac{g_A^2 h_A^2}{54 \pi^2 r^5 y F_\pi^4} e^{-2x} (1+x)(3+3x+2x^2)  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_ttau_d(r) result(vnttaud)
    implicit none
    real(dp), intent(in) :: r !< point at which the function will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnttaud = - gA**2 * hA**2 * hbar_c**5 * exp(-2*x) * (1 + x) * (3 + 3*x + 2*x**2) &
        / (54*pi **2 * r**5 * y * Fpi**4)

end function v_nlo_ttau_d

!!
!> @brief       TPE at NLO: \f$ v_{t \tau \mathrm{small r}}^{\mathrm{NLO}}(r;\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A13) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{t \tau} = -\frac{g_A^2 h_A^2}{54\pi^2 F_\pi^4 \Delta M} e^{-2x} (1 + x) (3 + 3x + x^2)
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 +
!!    \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_ttau_d_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = -gA**2*hA**2*hbar_c**6*exp(-2*x)*(1 + x)*(3 + 3*x + 2*x**2) &
        /(54*pi**2*delta_nucleon_mass_difference*Fpi**4)*(exp(-R_L/a_L)/R_L**6)&
        *(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_ttau_d_small_r

!> @brief       NLO loop correction: \f$ v_{c}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A14) in the appendix
!!
!! \f[ v_{c}^{\mathrm{NLO}}(r;2\Delta) = &
!!    -\frac{1}{108\pi^{3}r^{5}} \frac{h_{A}^{4}}{F_{\pi}^{4}}
!!    \big[
!!    4y^{2} \cdot \mathrm{vf2} + 2\mathrm{vf4} + \frac{1}{y} (4x^{4}-12y^{4}-8x^{2}y^{2}) \mathrm{vf5}\\
!!    & \frac{1}{y} \left(4x^{2}-4y^{2}\right)\mathrm{vf6} +\frac{1}{y}\mathrm{vf7}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_c_2d(r, vf2, vf4, vf5, vf6, vf7) result(vnc2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)
    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnc2d = -hA**4 * hbar_c**5 * (4*y**2*vf2 + 2*vf4 + ((4*x**4 -8*x**2*y**2 - 12*y**4)*vf5 &
        + (4*x**2 - 4*y**2)*vf6 + vf7)/y) &
        / (108 * pi**3 * r**5 * Fpi**4)

end function v_nlo_c_2d

!!
!> @brief       TPE at NLO: \f$ v_{c \mathrm{small r}}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A14) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{c} = -\frac{h_A^2}{108\pi^3 F_\pi^4 \Delta M}(4y^3\mathrm{vf2} + 2y\mathrm{vf4} + 
!!    (4x^4 - 8x^2y^2 -12y^4)\mathrm{vf5} + (4x^2 - 4y^2)\mathrm{vf6} + \mathrm{vf7})
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!    \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_c_2d_small_r(r, R_L, a_L, vf2, vf4, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = -hA**4*hbar_c**6*(4*y**3*vf2 + 2*y*vf4 + (4*x**4 -8*x**2*y**2 - 12*y**4)*vf5 &
        + (4*x**2 - 4*y**2)*vf6 + vf7)/(108*pi**3*delta_nucleon_mass_difference*Fpi**4) &
        *(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_c_2d_small_r

!> @brief       NLO loop correction: \f$ v_{\tau}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order.
!! Corresponds to (A15) in the appendix
!!
!! \f[ v_{\tau}^{\mathrm{NLO}}(r;2\Delta) = &
!!    -\frac{1}{1944\pi^{3}r^{5}} \frac{h_{A}^{4}}{F_{\pi}^{4}}
!!    \big[
!!    11 \mathrm{vf1} + (24x^{2}+24y^{2})\mathrm{vf2} + 6\mathrm{vf4}\\
!!    & -\frac{3}{y} \left(24x^{2}y^{2}+4x^{4}+20y^{4}\right)\mathrm{vf5} -\frac{3}{y}(4x^{2}+12y^{2})\mathrm{vf6} -\frac{3}{y}\mathrm{vf7}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) result(vntau2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vntau2d = -hA**4 * hbar_c**5 * (11*vf1 + (24*x**2 + 24*y**2)*vf2 + 6*vf4 - 3*(24*x**2*y**2 + 4*x**4 +20*y**4)*vf5/y &
        - 3*(4*x**2 + 12*y**2)*vf6/y - 3*vf7/y) &
        / (1944 * pi**3 * r**5 * Fpi**4)

end function v_nlo_tau_2d

!!
!> @brief       TPE at NLO: \f$ v_{\tau \mathrm{small r}}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A15) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\tau} = -\frac{h_A^2}{1944\pi^3 F_\pi^4 \Delta M}(11y\mathrm{vf1} + 
!!    y(24x^2 + 24y^2)\mathrm{vf2} + 6y\mathrm{vf4} - 3(24x^4y^2 + 4x^4 + 20y^4)\mathrm{vf5} - 
!!    3(4x^2 + 12y^2)\mathrm{vf6} - 3\mathrm{vf7}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + 
!!    \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_tau_2d_small_r(r, R_L, a_L, vf1, vf2, vf4, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = -hA**4*hbar_c**6*(11*y*vf1 + y*(24*x**2 + 24*y**2)*vf2 + 6*y*vf4 - 3*(24*x**2*y**2 + 4*x**4 + 20*y**4)*vf5 &
        - 3*(4*x**2 + 12*y**2)*vf6 - 3*vf7)/(1944*pi**3*delta_nucleon_mass_difference*Fpi**4) &
        *(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_tau_2d_small_r

!> @brief       NLO loop correction: \f$ v_{\sigma}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order.
!! Corresponds to (A16) in the appendix
!!
!! \f[ v_{\sigma}^{\mathrm{NLO}}(r;2\Delta) = &
!!    -\frac{1}{1296\pi^{3}r^{5}} \frac{h_{A}^{4}}{F_{\pi}^{4}}
!!    \big[ -6\mathrm{vf1} -24x^{2} \cdot \mathrm{vf2} + 48x^{2}y\cdot\mathrm{vf5}\\
!!    & + \frac{1}{y}(4x^{2} + 12y^{2})\mathrm{vf6} + \frac{1}{y}\mathrm{vf7}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7) result(vnsigma2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
        real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnsigma2d = -hA**4 * hbar_c**5 * (-6*vf1 - 24*x**2*vf2 + 48*y*x**2*vf5 + (4*x**2 + 12*y**2)*vf6/y + vf7/y) &
        / (1296 * pi**3 * r**5 * Fpi**4)  

end function v_nlo_sigma_2d

!!
!> @brief       TPE at NLO: \f$ v_{\sigma \mathrm{small r}}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A16) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\sigma} = -\frac{h_A^2}{1296 \pi^3 F_\pi^4 \Delta M}(-6y\mathrm{vf1} - 24x^2y\mathrm{vf2} + 
!!  48y^2x^2\mathrm{vf5} + (4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7}) \frac{e^{-R_L/a_L}}{R_L^6}
!!  \left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_sigma_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = -hA**4*hbar_c**6*(-6*y*vf1 - 24*x**2*y*vf2 + 48*y**2*x**2*vf5 + (4*x**2 + 12*y**2)*vf6 + vf7) &
        /(1296*pi**3*delta_nucleon_mass_difference*Fpi**4) &
        *(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_sigma_2d_small_r

!!
!> @brief       NLO loop correction: \f$ v_{\sigma \tau}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order.
!! Corresponds to (A17) in the appendix
!!
!! \f[ v_{\sigma \tau}^{\mathrm{NLO}}(r;2\Delta) = &
!!    -\frac{1}{7776\pi^{3}r^{5}} \frac{h_{A}^{4}}{F_{\pi}^{4}}
!!    \big[
!!    -2\mathrm{vf1} -8x^{2} \cdot \mathrm{vf2} + \frac{1}{y}(16x^{2}y^{2})\mathrm{vf5}\\
!!    & + \frac{1}{y}(4y^{2} - 4x^{2})\mathrm{vf6} - \frac{1}{y}\mathrm{vf7}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7) result(vnstau2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnstau2d = -hA**4 * hbar_c**5 * (-2*vf1 - 8*x**2*vf2 + 16*y*x**2*vf5 + (4*y**2 - 4*x**2)*vf6/y - vf7/y) &
        / (7776 * pi**3 * r**5 * Fpi**4)

end function v_nlo_sigmatau_2d

!!
!> @brief       TPE at NLO: \f$ v_{\sigma \tau \mathrm{small r}}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A17) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{\sigma \tau} = -\frac{h_A^4}{7776 \pi^3 F_\pi^4 \Delta M}(-2y\mathrm{vf1} - 8x^2y\mathrm{vf2} + 
!!  16y^2x^2\mathrm{vf5} + (4y^2 - 4x^2)\mathrm{vf6} - \mathrm{vf7}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + 
!!  \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_sigmatau_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = -hA**4*hbar_c**6*(-2*y*vf1 - 8*x**2*y*vf2 + 16*y**2*x**2*vf5 + (4*y**2 - 4*x**2)*vf6 - vf7) &
        /(7776*pi**3*delta_nucleon_mass_difference*Fpi**4) &
        *(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)


end function v_nlo_sigmatau_2d_small_r


!!
!> @brief       NLO loop correction: \f$ v_{t}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order.
!! Corresponds to (A18) in the appendix
!!
!! \f[     v_{t}^{\mathrm{NLO}}(r;2\Delta) = &
!!    \frac{1}{2592\pi^{3}r^{5}} \frac{h_{A}^{4}}{F_{\pi}^{4}}
!!    \big[
!!    -6 \mathrm{vf1} - 6(3 + 4x^{2})\mathrm{vf2} - 18 \mathrm{vf3}\\
!!    & + \frac{1}{y}(36y^{2}+48x^{2}y^{2})\mathrm{vf5} + \frac{1}{y}(3+4x^{2}+12y^{2})\mathrm{vf6}\\
!!    & + \frac{1}{y}\mathrm{vf7} + \frac{3}{y}\mathrm{vf8} + 36y \cdot \mathrm{vf9}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vnt2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnt2d = hA**4 * hbar_c**5 * (-6*vf1 - 6*(3 + 4*x**2)*vf2 - 18*vf3 + (36*y**2 + 48*x**2*y**2)*vf5/y &
        + (3 + 4*x**2 + 12*y**2)*vf6/y + vf7/y + 3*vf8/y + 36*y*vf9) &
        / (2592 * pi**3 * r**5 * Fpi**4)

end function v_nlo_t_2d

!!
!> @brief       TPE at NLO: \f$ v_{t \mathrm{small r}}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A18) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{t} = \frac{h_A^4}{2592 \pi^3 F_\pi^4 \Delta M} (-6y\mathrm{vf1} - 6y(3 + 4x^2)\mathrm{vf2} - 
!!  18y\mathrm{vf3} + (36y^2 + 48x^2y^2)\mathrm{vf5} + (3 + 4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7} + 3\mathrm{vf8} + 
!!  36y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!  \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_t_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = hA**4*hbar_c**6*(-6*y*vf1 - 6*y*(3 + 4*x**2)*vf2 - 18*y*vf3 + (36*y**2 + 48*x**2*y**2)*vf5 &
        + (3 + 4*x**2 + 12*y**2)*vf6 + vf7 + 3*vf8 + 36*y**2*vf9) &
        / (2592*pi**3*delta_nucleon_mass_difference*Fpi**4) &
        *(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_t_2d_small_r

!!
!> @brief       NLO loop correction: \f$ v_{t\tau}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order.
!! Corresponds to (A19) in the appendix
!!
!! \f[ v_{t\tau}^{\mathrm{NLO}}(r;2\Delta) = &
!!    \frac{1}{15552\pi^{3}r^{5}} \frac{h_{A}^{4}}{F_{\pi}^{4}}
!!    \big[
!!    -2\mathrm{vf1} - 2(3+4x^{2})\mathrm{vf2} - 6 \mathrm{vf3}\\
!!    & + \frac{1}{y}(12y^{2}+16x^{2}y^{2})\mathrm{vf5} + \frac{1}{y}(4y^{2}-4x^{2}-3)\mathrm{vf6}\\
!!    & - \frac{1}{y}\mathrm{vf7} - \frac{3}{y}\mathrm{vf8} + 12y \cdot \mathrm{vf9}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_nlo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vnttau2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)


    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    vnttau2d = hA**4 * hbar_c**5 * (-2*vf1 - 2*(3 + 4*x**2)*vf2 - 6*vf3 + (12*y**2 + 16*x**2*y**2)*vf5/y &
        + (4*y**2 - 4*x**2 - 3)*vf6/y - vf7/y -3*vf8/y + 12*y*vf9) &
        / (15552 * pi**3 * r**5 * Fpi**4)  

end function v_nlo_ttau_2d

!!
!> @brief       TPE at NLO: \f$ v_{t \tau \mathrm{small r}}^{\mathrm{NLO}}(r;2\Delta) \f$
!!
!! Two pion exchange potential contribution at next leading order
!! Corresponds to (A19) in the appendix
!!
!! [\f v^{\mathrm{NLO}}_{t \tau} = \frac{h_A^4}{15552 \pi^3 F_\pi^4 \Delta M} (-2y\mathrm{vf1} - 2y(3 + 4x^2)\mathrm{vf2} - 
!!  6y\mathrm{vf3} + (12y^2 + 16x^2y^2)\mathrm{vf5} + (4y^2 - 4x^2 - 3)\mathrm{vf6} - \mathrm{vf7} - 3\mathrm{vf8} + 
!!  12y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!  \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_nlo_ttau_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)


    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = hA**4*hbar_c**6*(-2*y*vf1 - 2*y*(3 + 4*x**2)*vf2 - 6*y*vf3 + (12*y**2 + 16*x**2*y**2)*vf5 &
        + (4*y**2 - 4*x**2 - 3)*vf6 - vf7 - 3*vf8 + 12*y**2*vf9) &
        / (15552*pi**3*delta_nucleon_mass_difference*Fpi**4) &
        *(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_nlo_ttau_2d_small_r

!! NEXT-TO-NEXT-TO LEADING ORDER POTENTIALS

!!
!> @brief       N2LO loop correction: \f$ v_c^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A20) in the appendix
!!
!! \f[ v_c^{\mathrm{N2LO}}(r) = \frac{3g_A^2}{2\pi^2r^6F_\pi^4} e^{-2x} \left[2c_1x^2(1+x)^2 + c_3(6+12x+10x^2+4x^3+x^4)\right]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_c(r) result(vn2c)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    vn2c = 3*gA**2 * hbar_c**6 * exp(-2*x) * (2*c1*x**2 *(1 + x)**2 + c3*(6 + 12*x + 10*x**2 + 4*x**3 + x**4)) &
        / (2*pi**2 * r**6 * Fpi**4)
end function v_n2lo_c

!!
!> @brief       N2LO loop correction: \f$ v_{c \mathrm{small r}}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A20) in the appendix
!!
!! [\f v_c^{\mathrm{N2LO}}(r) = \frac{3g_A^2}{2\pi^2 F_\pi^4} e^{-2x} \left(2c_1x^2(1 + x)^2 + c_3(6 + 12x + 10x^2 + 4x^3 + x^4)\right)
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_c_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = 3*gA**2*hbar_c**6*exp(-2*x)*(2*c1*x**2*(1 + x)**2 + c3*(6 + 12*x + 10*x**2 + 4*x**3 + x**4)) &
        /(2*pi**2*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)
end function v_n2lo_c_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma \tau}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A21) in the appendix
!!
!! \f[ v_{\sigma \tau}^{\mathrm{N2LO}}(r) = \frac{g_A^2}{3 \pi^2 r^6 F_\pi^4} c_4 e^{-2x} (1+x)(3+3x+2x^2)\f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigmatau(r) result(vn2stau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    vn2stau = gA**2 * hbar_c**6 * c4*exp(-2*x) * (1 + x)*(3 + 3*x + 2*x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_sigmatau

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma \tau \mathrm{small r}}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A21) in the appendix
!!
!! [\f v_{\sigma \tau}^{\mathrm{N2LO}}(r) = \frac{g_A^2}{3 \pi^2 F_\pi^4} c_4 e^{-2x} (1 + x)(3 + 3x + 2x^2)
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!    \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = gA**2*hbar_c**6*c4*exp(-2*x)*(1 + x)*(3 + 3*x + 2*x**2) &
        /(3*pi**2*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)
end function v_n2lo_sigmatau_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{t \tau}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A22) in the appendix
!!
!! \f[ v_{t \tau}^{\mathrm{N2LO}}(r) = - \frac{g_A^2}{3 \pi^2 r^6 F_\pi^4} c_4 e^{-2x} (1+x)(3+3x+x^2) \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_ttau(r) result(vn2ttau)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    vn2ttau = - gA**2 * hbar_c**6 * c4 * exp(-2*x) * (1 + x)*(3 + 3*x + x**2) &
        / (3*pi**2 * r**6 * Fpi**4)
end function v_n2lo_ttau

!!
!> @brief       N2LO loop correction: \f$ v_{t \tau \mathrm{small r}}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A22) in the appendix
!!
!! [\f v_{\sigma \tau}^{\mathrm{N2LO}}(r) = \frac{- g_A^2}{3 \pi^2 F_\pi^4} c_4 e^{-2x} (1 + x)(3 + 3x + x^2)
!!    \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!    \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_ttau_small_r(r, R_L, a_L) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp) :: x                   !< pion mass (unitless)

    x = pion_mass_r(r)

    v = -gA**2*hbar_c**6*c4*exp(-2*x)*(1 + x)*(3 + 3*x + x**2) &
        /(3*pi**2*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)
end function v_n2lo_ttau_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{c}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A23) in the appendix
!!
!! \f[ v_{c}^{\mathrm{N2LO}}(r;\Delta) = &
!!    \frac{1}{18\pi^{3}r^{6}} \frac{h_{A}^{2}y}{F_{\pi}^{4}}
!!    \big[
!!    (5c_{2} - 6c_{3})\mathrm{vf1} \\
!!    & + \big( (-24c_{1} + 12c_{2} - 12c_{3})x^{2} + 12c_{2}y^{2}\big)\mathrm{vf2} \\
!!    & + \frac{6}{y}\big( (8c_{1}+4c_{3})x^{4} + (8c_{1}-4c_{2}+4c_{3})x^{2}y^{2} -4c_{2}y^{4}\big)\mathrm{vf5}\\
!!    & + \frac{6}{y}\big( (4c_{1}+4c_{3})x^{2} + (-2c_{2}+2c_{3})y^{2}\big)\mathrm{vf6} + \frac{6}{y}c_{3}\cdot\mathrm{vf7}
!!    \big]         \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_c_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2cd)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2cd = hA**2 * y * hbar_c**6 * ((5*c2-6*c3)*vf1 + ((-24*c1 + 12*c2 - 12*c3)*x**2 + 12*c2*y**2)*vf2 &
            + 6*(((8*c1 + 4*c3)*x**4 + (8*c1 - 4*c2 + 4*c3)*x**2*y**2 - 4*c2*y**4)*vf5 &
            + ((4*c1 + 4*c3)*x**2 + (-2*c2 + 2*c3)*y**2)*vf6 + c3*vf7)/y) &
            / (18*pi**3 * r**6 * Fpi**4)

end function v_n2lo_c_d

!!
!> @brief       N2LO loop correction: \f$ v_{c \mathrm{small r}}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A23) in the appendix
!!
!! [\f v^{\mathrm{N2LO}}_{c} = \frac{h_A^2}{18 \pi^3 F_\pi^4} (y(5c_2 - 6c_3)\mathrm{vf1} + y((-24c_1 + 12c_2 -12c_3)x^2y^2 + 
!!      (12c_2 y^2))\mathrm{vf2} +  6\left(((8c_1 + 4c_3)x^4 + (8c_1 - 4c_2 + 4c_3)x^2y^2 - 4c_2y^4)\mathrm{vf5} + 
!!      ((4c_1 + 4c_3)x^2 + (-2c_2 + 2c_3)y^2)\mathrm{vf6} + c_3\mathrm{vf7}\right)) \frac{e^{-R_L/a_L}}{R_L^6}
!!      \left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_c_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = hA**2*hbar_c**6*(y*(5*c2-6*c3)*vf1 + y*((-24*c1 + 12*c2 - 12*c3)*x**2 + 12*c2*y**2)*vf2 &
        + 6*(((8*c1 + 4*c3)*x**4 + (8*c1 - 4*c2 + 4*c3)*x**2*y**2 - 4*c2*y**4)*vf5 &
        + ((4*c1 + 4*c3)*x**2 + (-2*c2 + 2*c3)*y**2)*vf6 + c3*vf7)) &
        /(18*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_c_d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\tau}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A24) in the appendix
!!
!! \f[ v_{\tau}^{\mathrm{N2LO}}(r;\Delta) = &
!!    -\frac{1}{54\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}y}{F_{\pi}^{4}}
!!    \big[
!!    (5-11g_{A}^{2})\mathrm{vf1} \\
!!    & + \big( 12x^{2}+12y^{2}-g_{A}^{2}(24x^{2}+12y^{2} \big)\mathrm{vf2} \\
!!    & + \big( -24y(x^{2}+y^{2}) + \frac{24}{y}g_{A}^{2}(x^{4}+2x^{2}y^{2}+y^{4}) \big)\mathrm{vf5}\\
!!    & + \frac{6}{y}\big( -12y + \frac{24}{y}g_{A}^{2}(x^{2}+y^{2}) \big)\mathrm{vf6} + \frac{6}{y}g_{A}^{2}\cdot\mathrm{vf7}
!!    \big]       \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_tau_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2taud)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)


    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2taud = -b38 * hA * y * hbar_c**6 * ((5 - 11*gA**2)*vf1 + (12*x**2 + 12*y**2 - gA**2 * (24*x**2 + 12*y**2))*vf2 &
            + (-12*y*(2*x**2 + 2*y**2) + 6*gA**2*(4*x**4 + 8*x**2*y**2 + 4*y**4)/y)*vf5 &
            + (-12*y + 6*gA**2*(4*x**2 + 4*y**2)/y)*vf6 + 6*gA**2*vf7/y) &
            / (54*pi**3 * r**6 * Fpi**4)

end function v_n2lo_tau_d

!!
!> @brief       N2LO loop correction: \f$ v_{c \mathrm{small r}}^{\mathrm{N2LO}}(r) \f$
!!
!! Corresponds to (A24) in the appendix
!!
!! [\f v^{\mathrm{N2LO}}_{c} = \frac{-b_{38}h_A}{54 \pi^3 F_\pi^4} (y(5 - 11g_A^2)\mathrm{vf1} + 
!!      y(12x^2 + 12y^2 - g_A^2(24x^2 + 12y^2))\mathrm{vf2} + (-12y^2(2x^2 + 2y^2) + 6g_A^2(4x^4 + 8x^2y^2 + 4y^2))\mathrm{vf5} + 
!!      (-12y^2 + 6g_A^2(4x^2 + 12y^2))\mathrm{vf6} + 6g_A^2\mathrm{vf7}) \frac{e^{-R_L/a_L}}{R_L^6}
!!      +\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_tau_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)


    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -b38*hA*hbar_c**6*(y*(5 - 11*gA**2)*vf1 + y*(12*x**2 + 12*y**2 - gA**2 * (24*x**2 + 12*y**2))*vf2 &
        + (-12*y**2*(2*x**2 + 2*y**2) + 6*gA**2*(4*x**4 + 8*x**2*y**2 + 4*y**4))*vf5 &
        + (-12*y**2 + 6*gA**2*(4*x**2 + 4*y**2))*vf6 + 6*gA**2*vf7) &
        / (54*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_tau_d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A25) in the appendix
!!
!! \f[ v_{\sigma}^{\mathrm{N2LO}}(r;\Delta) = &
!!    -\frac{1}{18\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}g_{A}^{2}y}{F_{\pi}^{4}}
!!    \big[
!!    2\mathrm{vf1} + 8x^{2}\cdot\mathrm{vf2} -16x^{2}y\cdot\mathrm{vf5}\\
!!    & - \frac{4}{y}\big( x^{2} + y^{2} \big)\mathrm{vf6} - \frac{1}{y}\mathrm{vf7}
!!    \big]     \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigma_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2sigmad)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2sigmad = -b38 * hA * gA**2 * y * hbar_c**6 * (2*vf1 + 8*x**2*vf2 - (16*x**2*y**2*vf5 + 4*(x**2 + y**2)*vf6 + vf7)/y) &
            / (18*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigma_d

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma \mathrm{small r}}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A25) in the appendix
!!
!! [\f     v^{\mathrm{N2LO}}_{\sigma} = \frac{-b_{38} h_A g_A}{18 \pi^3 F_\pi^4} (2y\mathrm{vf1} 8x^2y\mathrm{vf2} - 
!!  16x^2y^2\mathrm{vf5} - 4(x^2 + y^2)\mathrm{vf6} - \mathrm{vf7}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + 
!!  \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_sigma_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -b38*hA*gA**2*hbar_c**6*(2*y*vf1 + 8*x**2*y*vf2 - 16*x**2*y**2*vf5 - 4*(x**2 + y**2)*vf6 - vf7) &
        /(18*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_sigma_d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma\tau}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A26) in the appendix
!!
!! \f[ v_{\sigma\tau}^{\mathrm{N2LO}}(r;\Delta) = &
!!    -\frac{1}{108\pi^{3}r^{6}} \frac{c_{4}h_{A}^{2}y}{F_{\pi}^{4}}
!!    2\mathrm{vf1} + 8x^{2}\cdot\mathrm{vf2} -16x^{2}y\cdot\mathrm{vf5}\\
!!    & - \frac{4}{y}\big( x^{2} + y^{2} \big)\mathrm{vf6} - \frac{1}{y}\mathrm{vf7}
!!    \big]   \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_d(r, vf1, vf2, vf5, vf6, vf7) result(vn2staud)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2staud = -c4* hA**2 * y * hbar_c**6 * (2*vf1 + 8*x**2*vf2 - (16*x**2*y**2*vf5 + 4*(x**2 + y**2)*vf6 + vf7)/y) &
            / (108*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigmatau_d

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma \tau \mathrm{small r}}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A26) in the appendix
!!
!! [\f     v^{\mathrm{N2LO}}_{\sigma \tau} = \frac{-c_4 h_A}{108 \pi^3 F_\pi^4} (2y\mathrm{vf1} + 8x^2y\mathrm{vf2} - 
!!      16x^2y^2\mathrm{vf5} - 4(x^2 + y^2)\mathrm{vf6} - \mathrm{vf7}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + 
!!      \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -c4*hA**2*hbar_c**6*(2*y*vf1 + 8*x**2*y*vf2 - 16*x**2*y**2*vf5 - 4*(x**2 + y**2)*vf6 - vf7) &
        /(108*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_sigmatau_d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{t}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A27) in the appendix
!!
!! \f[ v_{t}^{\mathrm{N2LO}}(r;\Delta) = &
!!   \frac{1}{36\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}g_{A}^{2}y}{F_{\pi}^{4}}
!!   \big[
!!   2\mathrm{vf1} + \big(6+8x^{2}\big)\mathrm{vf2} + 6\mathrm{vf3}\\
!!   & - \frac{1}{y}\big(12y^{2}+16x^{2}y^{2}\big)\mathrm{vf5} - \frac{1}{y}\big(3+4x^{2}+4y^{2}\big)\mathrm{vf6} - \frac{1}{y}\mathrm{vf7}\\
!!   & - \frac{3}{y}\mathrm{vf8} -12y\cdot\mathrm{vf9}
!!   \big]   \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_t_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2td)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2td = b38*hA*gA**2 * y*hbar_c**6 * (2*vf1 + (6+8*x**2)*vf2 + 6*vf3 - (12*y**2 + 16*x**2 * y**2)*vf5/y &
            - (3 + 4*x**2 + 4*y**2)*vf6/y - vf7/y - 3*vf8/y - 12*y*vf9) &
            / (36*pi**3 * r**6 * Fpi**4)

end function v_n2lo_t_d

!!
!> @brief       N2LO loop correction: \f$ v_{t \mathrm{small r}}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A27) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{t} = \frac{b_{38} h_A g_A^2}{36 \pi^3 F_\pi^4} (2y\mathrm{vf1} + y(6 + 8x^2)\mathrm{vf2} + 6y\mathrm{vf3} - 
!!      (12y^2 + 16x^2y^2)mathrm{vf5} -(3 + 4x^2 + 4y^2)\mathrm{vf6} - \mathrm{vf7} -3\mathrm{vf8} -12y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}
!!      \left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_t_d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = b38*hA*gA**2*hbar_c**6*(2*y*vf1 + y*(6 + 8*x**2)*vf2 + 6*y*vf3 - (12*y**2 + 16*x**2 * y**2)*vf5 &
        - (3 + 4*x**2 + 4*y**2)*vf6 - vf7 - 3*vf8 - 12*y**2*vf9) &
        /(36*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_t_d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{t\tau}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A28) in the appendix
!!
!! \f[ v_{t\tau}^{\mathrm{N2LO}}(r;\Delta) = &
!!    \frac{1}{216\pi^{3}r^{6}} \frac{c_{4}h_{A}^{2}y}{F_{\pi}^{4}}
!!    \big[
!!    2\mathrm{vf1} + \big(6+8x^{2}\big)\mathrm{vf2} + 6\mathrm{vf3}\\
!!    & - \frac{1}{y}\big(12y^{2}+16x^{2}y^{2}\big)\mathrm{vf5} - \frac{1}{y}\big(3+4x^{2}+4y^{2}\big)\mathrm{vf6} - \frac{1}{y}\mathrm{vf7}\\
!!    & - \frac{3}{y}\mathrm{vf8} -12y\cdot\mathrm{vf9}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_ttau_d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2ttaud)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2ttaud = c4* hA**2 * y * hbar_c**6 * (2*vf1 + (6+8*x**2)*vf2 + 6*vf3 - (12*y**2 + 16*x**2 * y**2)*vf5/y &
            - (3 + 4*x**2 + 4*y**2)*vf6/y - vf7/y - 3*vf8/y - 12*y*vf9) &
            / (216*pi**3 * r**6 * Fpi**4)

end function v_n2lo_ttau_d

!!
!> @brief       N2LO loop correction: \f$ v_{t \tau \mathrm{small r}}^{\mathrm{N2LO}}(r;\Delta) \f$
!!
!! Corresponds to (A28) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{t \tau} = \frac{c_4 h_A^2}{216 \pi^3 F_\pi^4} (2y\mathrm{vf1} + (6 + 8x^2)y\mathrm{vf2} + 
!!      6y\mathrm{vf3} - (12y^2 + 16x^2y^2)mathrm{vf5} - (3 + 4x^2 + 4y^2)\mathrm{vf6} - \mathrm{vf7} - 3\mathrm{vf8} - 
!!      12y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!      \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_ttau_d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)

    v = c4*hA**2*hbar_c**6*(2*y*vf1 + (6 + 8*x**2)*y*vf2 + 6*y*vf3 - (12*y**2 + 16*x**2 * y**2)*vf5 &
        - (3 + 4*x**2 + 4*y**2)*vf6 - vf7 - 3*vf8 - 12*y**2*vf9) &
        /(216*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_ttau_d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{c}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A29) in the appendix
!!
!! \f[ v_{c}^{\mathrm{N2LO}}(r;2\Delta) = &
!!    -\frac{2}{81\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}^{3}y}{F_{\pi}^{4}}
!!    \big[
!!    11\mathrm{vf1} + \big(24x^{2}+12y^{2}\big)\mathrm{vf2} + 6\mathrm{vf4}\\
!!    & -\frac{3}{y} \left(24x^{2}y^{2}+4x^{2}+20y^{2}\right)\mathrm{vf5} -\frac{3}{y}(4x^{2}+12y^{2})\mathrm{vf6} -\frac{3}{y}\mathrm{vf7}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_c_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) result(vn2c2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2c2d = -2*b38 * hA**3 * y * hbar_c**6 * (11*vf1 + (24*x**2 + 12*y**2)*vf2 + &
             6*vf4 - 3*((4*x**4 + 24*x**2*y**2 + 20*y**4)*vf5 + (4*x**2 + 12*y**2)*vf6 &
             + vf7)/y) &
            / (81*pi**3 * r**6 * Fpi**4)

end function v_n2lo_c_2d

!!
!> @brief       N2LO loop correction: \f$ v_{c \mathrm{small r}}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A29) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{t \tau} = \frac{-2 b_{38} h_A^3}{81 \pi^3 F_\pi^4} (11y\mathrm{vf1} + 
!!      (24x^2 + 12y^2)y\mathrm{vf2} + 6y\mathrm{vf4} - 3((4x^4 + 24x^2y^2 + 20y^4)mathrm{vf5} + 
!!      (4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7})) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + 
!!      \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_c_2d_small_r(r, R_L, a_L, vf1, vf2, vf4, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -2*b38*hA**3*hbar_c**6*(11*y*vf1 + (24*x**2 + 12*y**2)*y*vf2 + &
        6*y*vf4 - 3*((4*x**4 + 24*x**2*y**2 + 20*y**4)*vf5 + (4*x**2 + 12*y**2)*vf6 &
        + vf7))/(81*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_c_2d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\tau}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A30) in the appendix
!!
!! \f[ v_{\tau}^{\mathrm{N2LO}}(r;2\Delta) = &
!!    -\frac{1}{243\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}^{3}y}{F_{\pi}^{4}}
!!    \big[
!!    11\mathrm{vf1} + \big(24x^{2}+12y^{2}\big)\mathrm{vf2} + 6\mathrm{vf4}\\
!!    & -\frac{3}{y} \left(24x^{2}y^{2}+4x^{2}+20y^{2}\right)\mathrm{vf5} -\frac{3}{y}(4x^{2}+12y^{2})\mathrm{vf6} -\frac{3}{y}\mathrm{vf7}
!!    \big] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_tau_2d(r, vf1, vf2, vf4, vf5, vf6, vf7) result(vn2tau2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2tau2d = -b38 * hA**3 * y * hbar_c**6 * (11*vf1 + (24*x**2 + 12*y**2)*vf2 + &
        6*vf4 - 3*((4*x**4 + 24*x**2*y**2 + 20*y**4)*vf5 + (4*x**2 + 12*y**2)*vf6 &
        + vf7)/y) &
            / (243*pi**3 * r**6 * Fpi**4)

end function v_n2lo_tau_2d

!!
!> @brief       N2LO loop correction: \f$ v_{\tau \mathrm{small r}}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A30) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{\tau} = \frac{-b_{38} h_A^3}{243 \pi^3 F_\pi^4} (11y\mathrm{vf1} + (24x^2 + 12y^2)y\mathrm{vf2} + 
!!      6y\mathrm{vf4} - 3((4x^4 + 24x^2y^2 + 20y^4)mathrm{vf5} + (4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7})) \frac{e^{-R_L/a_L}}{R_L^6}
!!      \left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_tau_2d_small_r(r, R_L, a_L, vf1, vf2, vf4, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf4     !< evaluated integral 4; integrand defined in function vf_4(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -b38*hA**3*hbar_c**6*(11*y*vf1 + (24*x**2 + 12*y**2)*y*vf2 + &
        6*y*vf4 - 3*((4*x**4 + 24*x**2*y**2 + 20*y**4)*vf5 + (4*x**2 + 12*y**2)*vf6 &
        + vf7))/(243*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_tau_2d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A31) in the appendix
!!
!! \f[ v_{\sigma}^{\mathrm{N2LO}}(r;2\Delta) = &
!!    -\frac{1}{162\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}^{3}y}{F_{\pi}^{4}}
!!    \big[
!!    -6\mathrm{vf1} - 24x^{2}\cdot\mathrm{vf2} + 48x^{2}y\cdot\mathrm{vf5}\\
!!    & + \frac{1}{y}(4x^{2} + 12y^{2})\mathrm{vf6} + \frac{1}{y}\mathrm{vf7}
!!    \big]    \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigma_2d(r, vf1, vf2, vf5, vf6, vf7) result(vn2sigma2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2sigma2d = -b38 * hA**3 * y * hbar_c**6 * (-6*vf1 - 24*x**2*vf2 + 48*x**2*y*vf5 + (4*x**2 + 12*y**2)*vf6/y + vf7/y) &
            / (162*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigma_2d

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma \mathrm{small r}}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A31) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{\sigma} = \frac{-b_{38} h_A^3}{162 \pi^3 F_\pi^4} (-6y\mathrm{vf1} - 24x^2y\mathrm{vf2} + 48x^2y^2mathrm{vf5} + 
!!      (4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7})) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!      \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_sigma_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -b38*hA**3*hbar_c**6*(-6*y*vf1 - 24*x**2*y*vf2 + 48*x**2*y**2*vf5 + (4*x**2 + 12*y**2)*vf6 + vf7) &
        /(162*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_sigma_2d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma\tau}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A32) in the appendix
!!
!! \f[ v_{\sigma\tau}^{\mathrm{N2LO}}(r;2\Delta) = &
!!    -\frac{1}{972\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}^{3}y}{F_{\pi}^{4}}
!!    \big[
!!    -6\mathrm{vf1} - 24x^{2}\cdot\mathrm{vf2} + 48x^{2}y\cdot\mathrm{vf5}\\
!!    & + \frac{1}{y}(4x^{2} + 12y^{2})\mathrm{vf6} + \frac{1}{y}\mathrm{vf7}
!!    \big]  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_2d(r, vf1, vf2, vf5, vf6, vf7) result(vn2stau2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2stau2d = -b38 * hA**3 * y * hbar_c**6 * (-6*vf1 - 24*x**2*vf2 + 48*x**2*y*vf5 + (4*x**2 + 12*y**2)*vf6/y + vf7/y) &
            / (972*pi**3 * r**6 * Fpi**4)

end function v_n2lo_sigmatau_2d

!!
!> @brief       N2LO loop correction: \f$ v_{\sigma \tau \mathrm{small r}}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A32) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{\sigma \tau} = \frac{-b_{38} h_A^3}{972 \pi^3 F_\pi^4} (-6y\mathrm{vf1} - 24x^2y\mathrm{vf2} + 48x^2y^2mathrm{vf5} + 
!!      (4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7})) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!      \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_sigmatau_2d_small_r(r, R_L, a_L, vf1, vf2, vf5, vf6, vf7) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = -b38*hA**3*hbar_c**6*(-6*y*vf1 - 24*x**2*y*vf2 + 48*x**2*y**2*vf5 + (4*x**2 + 12*y**2)*vf6 + vf7) &
        /(972*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_sigmatau_2d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{t}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A33) in the appendix
!!
!! \f[ v_{t}^{\mathrm{N2LO}}(r;2\Delta) = &
!! \frac{1}{324\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}^{3}y}{F_{\pi}^{4}}
!! \big[
!! -6\mathrm{vf1} - 6\big(3+4x^{2}\big)\mathrm{vf2} - 24\mathrm{vf3}\\
!! & + \frac{1}{y}\big(36y^{2}+48x^{2}y^{2}\big)\mathrm{vf5} + \frac{1}{y}\big(3+4x^{2}+12y^{2}\big)\mathrm{vf6} + \frac{1}{y}\mathrm{vf7}\\
!! & + \frac{3}{y}\mathrm{vf8} + 36y\cdot\mathrm{vf9}
!! \big] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_t_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2t2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2t2d = b38 * hA**3 * y * hbar_c**6 * (-6*vf1 - 6*(3 + 4*x**2)*vf2 - 18*vf3 + ((36*y**2 + 48*x**2*y**2)*vf5 &
            + (3 + 4*x**2 + 12*y**2)*vf6 + vf7 + 3*vf8 + 36*y**2*vf9)/y) &
            / (324*pi**3 * r**6 * Fpi**4)

end function v_n2lo_t_2d

!!
!> @brief       N2LO loop correction: \f$ v_{t \mathrm{small r}}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A33) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{t} = \frac{b_{38} h_A^3}{324 \pi^3 F_\pi^4} (-6y\mathrm{vf1} - 6(3 + 4x^2)y\mathrm{vf2} - 
!!      18y\mathrm{vf3} + (36y^2 + 48x^2y^2)\mathrm{vf5} + (3 + 4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7} + 3\mathrm{vf8} + 
!!      36y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!      \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_t_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = b38*hA**3*hbar_c**6*(-6*y*vf1 - 6*(3 + 4*x**2)*y*vf2 - 18*y*vf3 + ((36*y**2 + 48*x**2*y**2)*vf5 &
        + (3 + 4*x**2 + 12*y**2)*vf6 + vf7 + 3*vf8 + 36*y**2*vf9)) &
        /(324*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_t_2d_small_r

!!
!> @brief       N2LO loop correction: \f$ v_{t\tau}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A34) in appendix
!!
!! \f[ v_{t\tau}^{\mathrm{N2LO}}(r;2\Delta) = &
!!    \frac{1}{1944\pi^{3}r^{6}} \frac{(b_{3}+b_{8})h_{A}^{3}y}{F_{\pi}^{4}}
!!    \big[
!!    -6\mathrm{vf1} - 6\big(3+4x^{2}\big)\mathrm{vf2} - 24\mathrm{vf3}\\
!!    & + \frac{1}{y}\big(36y^{2}+48x^{2}y^{2}\big)\mathrm{vf5} + \frac{1}{y}\big(3+4x^{2}+12y^{2}\big)\mathrm{vf6} + \frac{1}{y}\mathrm{vf7}\\
!!    & + \frac{3}{y}\mathrm{vf8} + 36y\cdot\mathrm{vf9}
!!    \big] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function v_n2lo_ttau_2d(r, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(vn2ttau2d)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    vn2ttau2d = b38 * hA**3 * y * hbar_c**6 * (-6*vf1 - 6*(3 + 4*x**2)*vf2 - 18*vf3 + ((36*y**2 + 48*x**2*y**2)*vf5 &
    + (3 + 4*x**2 + 12*y**2)*vf6 + vf7 + 3*vf8 + 36*y**2*vf9)/y) &
            / (1944*pi**3 * r**6 * Fpi**4)

end function v_n2lo_ttau_2d

!!
!> @brief       N2LO loop correction: \f$ v_{t \tau \mathrm{small r}}^{\mathrm{N2LO}}(r;2\Delta) \f$
!!
!! Corresponds to (A34) in the appendix
!!
!! [\f  v^{\mathrm{N2LO}}_{t \tau} = \frac{b_{38} h_A^3}{1944 \pi^3 F_\pi^4} (-6y\mathrm{vf1} - 6(3 + 4x^2)y\mathrm{vf2} - 
!!      18y\mathrm{vf3} + (36y^2 + 48x^2y^2)\mathrm{vf5} + (3 + 4x^2 + 12y^2)\mathrm{vf6} + \mathrm{vf7} + 3\mathrm{vf8} + 
!!      36y^2\mathrm{vf9}) \frac{e^{-R_L/a_L}}{R_L^6}\left(1 + \frac{r}{a_L} + \frac{1}{2}\left(\frac{r}{a_L}\right)^2 + 
!!      \frac{1}{6}\left(\frac{r}{a_L}\right)^3\right) \f]
!!
!! @author      Rodrigo Navarro Pérez, Ky Putnam
!!
real(dp) function v_n2lo_ttau_2d_small_r(r, R_L, a_L, vf1, vf2, vf3, vf5, vf6, vf7, vf8, vf9) result(v)
    implicit none
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: a_L     !< parameter of long range regulator (unitless constant)
    real(dp), intent(in) :: vf1     !< evaluated integral 1; integrand defined in function vf_1(u, r)
    real(dp), intent(in) :: vf2     !< evaluated integral 2; integrand defined in function vf_2(u, r)
    real(dp), intent(in) :: vf3     !< evaluated integral 3; integrand defined in function vf_3(u, r)
    real(dp), intent(in) :: vf5     !< evaluated integral 5; integrand defined in function vf_5(u, r)
    real(dp), intent(in) :: vf6     !< evaluated integral 6; integrand defined in function vf_6(u, r)
    real(dp), intent(in) :: vf7     !< evaluated integral 7; integrand defined in function vf_7(u, r)
    real(dp), intent(in) :: vf8     !< evaluated integral 8; integrand defined in function vf_8(u, r)
    real(dp), intent(in) :: vf9     !< evaluated integral 9; integrand defined in function vf_9(u, r)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = pion_mass_r(r)
    y = delta_nucleon_mass_difference_r(r)


    v = b38*hA**3*hbar_c**6*(-6*y*vf1 - 6*(3 + 4*x**2)*y*vf2 - 18*y*vf3 + ((36*y**2 + 48*x**2*y**2)*vf5 &
    + (3 + 4*x**2 + 12*y**2)*vf6 + vf7 + 3*vf8 + 36*y**2*vf9)) &
        /(1944*pi**3*Fpi**4)*(exp(-R_L/a_L)/R_L**6)*(1 + r/a_L + (r/a_L)**2/2 + (r/a_L)**3/6)

end function v_n2lo_ttau_2d_small_r

!!
!!
!! POTENTIAL INTEGRANDS
!!
!! Integrands are labeled based on their locations in the "Organizing Integrals" spreadsheet:
!! https://docs.google.com/spreadsheets/d/1dLe5_UPNaTz_hzaVfSETQZh9XZ9dlar5mD3eGOGdOVw/edit?usp=sharing
!!

!> @brief       potential integrand function 1: vf1
!!
!! output: vf1 (integrand "1" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf1} = \int_{0}^{\infty} d\mu \frac{\mu^{2}}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}}
!!   \left[\mu^{2}\right] \f]
!!
!! @author      Ky Putnam
!!
real(dp) function vf_1(u , r) result(vf1)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = mpi * r / hbar_c

    vf1 = u**4 * exp(-sqrt(u**2 + 4*x**2)) &
        /sqrt(u**2 + 4*x**2 + tiny(1._dp))
end function vf_1

!> @brief       potential integrand function 2: vf2
!!
!! output: vf2 (integrand "2" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf2} = \int_{0}^{\infty} d\mu \frac{\mu^{2}}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}} \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_2(u , r) result(vf2)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = mpi * r / hbar_c

    vf2 = u**2 * exp(-sqrt(u**2 + 4*x**2)) &
        /sqrt(u**2 + 4*x**2 + tiny(1._dp))
end function vf_2

!> @brief       potential integrand function 3: vf3
!!
!! output: vf3 (integrand "3" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf3} = \int_{0}^{\infty} d\mu \frac{\mu^{2}}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}}
!!    \left[\sqrt{\mu^{2}+4x^{2}}\right] \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_3(u, r) result(vf3)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)

    x = mpi * r / hbar_c

    vf3 = u**2 * exp(-sqrt(u**2 + 4*x**2))
end function vf_3

!> @brief       potential integrand function 4: vf4
!!
!! output: vf4 (integrand "4" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf4} = \int_{0}^{\infty} d\mu \frac{\mu^{2}}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}}
!!    \left[\frac{(2x^{2}+\mu^{2}+2y^{2})^{2}}{\mu^2+4y^{2}}\right] \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_4(u , r) result(vf4)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vf4 = u**2 * (2*x**2 + u**2 + 2*y**2)**2 * exp(-sqrt(u**2 + 4*x**2)) &
        /(sqrt(u**2 + 4*x**2) * (u**2 +4*y**2) + tiny(1._dp))
end function vf_4

!> @brief       potential integrand function 5: vf5
!!
!! output: vf5 (integrand "5" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf5} = \int_{0}^{\infty} d\mu \frac{\mu}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}} arctan{\left(\frac{\mu}{2y}\right)} \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_5(u , r) result(vf5)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vf5 = u * atan(u/(2*y + tiny(1._dp))) * exp(-sqrt(u**2 + 4*x**2)) &
        /sqrt(u**2 + 4*x**2 + tiny(1._dp))
end function vf_5

!> @brief       potential integrand function 6: vf6
!!
!! output: vf6 (integrand "6" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf6} = \int_{0}^{\infty} d\mu \frac{\mu}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}} arctan{\left(\frac{\mu}{2y}\right)}
!!    \left[\mu^{2}\right] \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_6(u , r) result(vf6)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vf6 = u**3 * atan(u/(2*y + tiny(1._dp))) * exp(-sqrt(u**2 + 4*x**2)) &
        /sqrt(u**2 + 4*x**2 + tiny(1._dp))
end function vf_6

!> @brief       potential integrand function 7: vf7
!!
!! output: vf7 (integrand "7" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf7} = \int_{0}^{\infty} d\mu \frac{\mu}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}} arctan{\left(\frac{\mu}{2y}\right)}
!!    \left[\mu^{4}\right] \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_7(u , r) result(vf7)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vf7 = u**5 * atan(u/(2*y + tiny(1._dp))) * exp(-sqrt(u**2 + 4*x**2)) &
        /sqrt(u**2 + 4*x**2 + tiny(1._dp))
end function vf_7

!> @brief       potential integrand function 8: vf8
!!
!! output: vf8 (integrand "8" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf8} = \int_{0}^{\infty} d\mu \frac{\mu}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}} arctan{\left(\frac{\mu}{2y}\right)}
!!    \left[\mu^{2} \sqrt{\mu^{2}+4x^{2}}\right] \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_8(u , r) result(vf8)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vf8 = u**3 * atan(u/(2*y + tiny(1._dp))) * exp(-sqrt(u**2 + 4*x**2))
end function vf_8

!> @brief       potential integrand function 9: vf9
!!
!! output: vf9 (integrand "9" from the "organizing integrals" sheet)
!! 
!! \f[ \mathrm{vf9} = \int_{0}^{\infty} d\mu \frac{\mu}{\sqrt{\mu^{2}+4x^{2}}} e^{-\sqrt{\mu^{2}+4x^{2}}} arctan{\left(\frac{\mu}{2y}\right)}
!!    \left[\sqrt{\mu^{2}+4x^{2}}\right] \f]
!! 
!! @author      Ky Putnam
!!
real(dp) function vf_9(u , r) result(vf9)
    implicit none
    real(dp), intent(in) :: u       !< integration variable, \f$ \mu \f$ (dimensionless)
    real(dp), intent(in) :: r       !< point at which potential will be evaluated (in fm)
    real(dp) :: x                   !< pion mass (unitless)
    real(dp) :: y                   !< delta nucleon mass difference (unitless)

    x = mpi * r / hbar_c
    y = delta_nucleon_mass_difference * r / hbar_c

    vf9 = u * atan(u/(2*y + tiny(1._dp))) * exp(-sqrt(u**2 + 4*x**2))
end function vf_9

end module long_range_chiral_potentials