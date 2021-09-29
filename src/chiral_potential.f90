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
use constants, only : gA, hA, Fpi=>pion_decay_amplitude, c1, c2, c3, c4, b38=>b3_b8, &
    mpi0=>pion_0_mass, mpic=>pion_c_mass
    ! Fpi=2fpi
                     
implicit none

private

contains

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (1) at leading order
!!
!! @author Ky Putnam
!!
real(dp) function v_lo_sigma_tau(r) result(vst)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: Y_n, Y_p

    Y_n = Y_pion(mpi0,r) !< Y for neutral pion
    Y_p = Y_pion(mpic,r) !< Y for (positively) charged pion

    vst = (Y_n + 2*Y_p)/3

end function v_lo_sigma_tau

!!
!> @brief       OPE at LO
!!
!! One pion exchange potential contribution (2) at leading order
!!
!! @author Ky Putnam
!!
real(dp) function v_lo_t_tau(r) result(vtt)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp) :: T_n, T_p

    T_n = T_pion(mpi0,r) !< T for neutral pion
    T_p = T_pion(mpic,r) !< T for (positively) charged pion

    vtt = (T_n + 2*T_p)/3

end function v_lo_t_tau

!!
!< @brief       For the vst and vtt functions above
!!
!! Y, a function of pion mass and radius for use in the preceeding 
!! OPE LO functions (v_lo_sigma_tau and v_lo_t_tau)
!!
!! @author Ky Putnam
!!
real(dp) function Y_pion(mpi, r) result(Y)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp), intent(in) :: mpi !< pion mass
    real(dp) :: x

    x = mpi * r

    Y = gA**2 * mpi**3 exp(-x) / (12*pi Fpi**2 * x)

end function Y_pion

!!
!! @        For the vtt function above
!!
!! T, a function of radius and pion mass, for use in the preceeding
!! OPE LO functions (v_lo_t_tau, specifically)
!!
!! @author Ky Putnam
!!
real(dp) function T_pion(mpi,r) result(T)
    implicit none
    real(dp), intent(in) :: r !< radius at which the function will be evaluated, in fm
    real(dp), intent(in) :: mpi !< pion mass
    real(dp) :: x, Y_n, Y_p

    x = mpi * r
    
    T = Y_pion(mpi,r) * (1 + 3/x + 3/x**2)
    
end function T_pion

end module chiral_potential