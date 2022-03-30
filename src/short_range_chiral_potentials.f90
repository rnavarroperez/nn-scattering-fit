!!
!! Short range chiral potentials
!!
!! Functions from "Minimally nonlocal nucleon-nucleon potentials with chiral two-pion exchange including delta resonances"
!!
!! Phys.Rev. C91 (2015)
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
module short_range_chiral_potentials
    use precisions, only : dp
    use constants, only : pi, gA, hA, Fpi=>pion_decay_amplitude, c1, c2, c3, c4, b38=>b3_b8, &
        mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, hbar_c, delta_nucleon_mass_difference
        ! Fpi=2fpi
    implicit none
    
    private

real(dp) function C_Rs(r) result(CRs)
    implicit none
    real(dp), intent(in) :: r ! radius

    CRs = hbar_c*exp(-(r/R_S)**2)/(pi)**(3/2)*R_S**3)

end function C_Rs
    
! real(dp) function d1_c_Rs(r, R_S) result(d1cRs)
!     implicit none
!     real(dp), intent(in) :: r ! r = radius
!     real(dp) :: d1cRs

!     d1cRs = -2*r*hbar_c*exp(-(r/R_S)**2)/(pi)**(3/2)*R_S**5)
    
! end function

! real(dp) function d2_c_Rs(r) result(d2cRs)
!     implicit none
!     real(dp), intent(in) :: r ! r = radius
!     real(dp) :: d2cRs

!     d2cRs = -2*hbar_c*exp(-(r/R_S)**2)*(R_S**2-2*r**2)/(pi)**3/2*R_S**7)
    
! end function

! real(dp) function d3_c_Rs(r) result(d3cRs)
!     implicit none
!     real(dp), intent(in) :: r ! r = radius
!     real(dp) :: d3cRs

!     d3cRs = 4*r*hbar_c*exp(-(r/R_S)**2)*(3*R_S**2-2*r**2)/(pi)**3/2*R_S**9)
    
! end function

! real(dp) function d4_c_Rs(r) result(d4cRs)
!     implicit none
!     real(dp), intent(in) :: r ! r = radius
!     real(dp) :: d4cRs

!     d4cRs = 4*hbar_c*exp(-(r/R_S)**2)*(4*r**4 - 12*r**2*R_S**2 + 3*R_S**4)/(pi)**3/2*R_S**11)
    
! end function

! ! !! following are the short-range potentials

! ! ! B11
! ! real(dp) function v_s_c(r) result(vsc)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius
! !     integer :: N 

! !     c_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vsc = C_s * C_rs + C_1 * (-C_Rs_2 - 2*C_Rs_1/r) &
! !         + D_1*(C_Rs_4 + 4*C_Rs_3/r)

! ! end function

! ! ! B12
! ! real(dp) function v_s_tau(r,R_s) result(vstau)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius
! !     integer :: N 

! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vstau = C_2 * (-C_Rs_2 - 2*C_Rs_1/r) + D_2*(C_Rs_4 + 4*C_Rs_3/r)

! ! end function

! ! ! B13
! ! real(dp) function v_s_sigma(r) result(vssigma)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vssigma = C_T*C_Rs + C_3 * (-C_Rs_2 - 2*C_Rs_1/r) + D_3*(C_Rs_4 + 4*C_Rs_3/r)


! ! end function

! ! ! B14
! ! real(dp) function v_s_sigma_tau(r) result(vssigmatau)
! ! implicit none
! ! real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vssigmatau = C_4 * (-C_Rs_2 - 2*C_Rs_1/r) + D_4*(C_Rs_4 + 4*C_Rs_3/r)

! ! end function

! ! ! B15
! ! real(dp) function v_s_t(r) result(vst)
! ! implicit none
! ! real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vst = -C_5 * (C_Rs_2 - C_Rs_1/r) + &
! !          D_5*(C_Rs_4 + C_Rs_3/r + 6*C_Rs_2/r**2 + 6*C_Rs_1/r**3)

! ! end function

! ! ! B16
! ! real(dp) function v_s_t_tau(r) result(vsttau)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vsttau = -C_6 * (C_Rs_2 - C_Rs_1/r) + &
! !          D_6*(C_Rs_4 + C_Rs_3/r + 6*C_Rs_2/r**2 + 6*C_Rs_1/r**3)


! ! end function

! ! ! B17
! ! real(dp) function v_s_b(r) result(vsb)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vsb = -C_7 * C_Rs_1/r + D_7*(C_Rs_3/r + 2*C_Rs_2/r**2 - 2*C_Rs_1/r**3)

! ! end function

! ! ! B18
! ! real(dp) function v_s_b_tau(r) result(vsbtau)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vsbtau = D_8*(C_Rs_3/r + 2*C_Rs_2/r**2 - 2*C_Rs_1/r**3)

! ! end function

! ! ! B19
! ! real(dp) function v_s_b_b(r) result(vsbb)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vsbb = -D_9*(C_Rs_2 - C_Rs_1/r)/r**2

! ! end function

! ! ! B20
! ! real(dp) function v_s_q(r) result(vsq)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)

! !     vsq = -D_10*(C_Rs_2 - C_Rs_1/r)/r**2

! ! end function

! ! ! B21
! ! real(dp) function v_s_q_sigma(r) result(vsqsigma)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)

! !     vsqsigma = -D_11*(C_Rs_2 - C_Rs_1/r)/r**2

! ! end function

! ! ! B25
! ! real(dp) function v_s_p_t_tau(r) result(vspttau)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vspttau = -D_15*(C_Rs_2 - C_Rs_1/r)

! ! end function

! ! ! B26
! ! real(dp) function v_s_upperT(r) result(vsupperT)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vsupperT = C_IV_0*C_Rs + C_IT_1*(-C_Rs_2 - 2*C_Rs_1/r)

! ! end function

! ! ! B27
! ! real(dp) function v_s_tau_z(r) result(vstauz)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vstauz = C_IV_0*C_Rs + C_IT_1*(-C_Rs_2 - 2*C_Rs_1/r)

! ! end function

! ! ! B28
! ! real(dp) function v_sigma_upperT_s(r) result(vssigmaupperT)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vssigmaupperT = C_IT_2*(-C_Rs_2 - 2*C_Rs_1/r)

! ! end function

! ! ! B30
! ! real(dp) function v_s_t_upperT(r)
! !     implicit none
! !     real(dp), intent(in) :: r ! r = radius

! !     C_Rs = C_Rs(r, 0)
! !     c_Rs_1 = C_Rs(r, 1)
! !     c_Rs_2 = C_Rs(r, 2)
! !     c_Rs_3 = C_Rs(r, 3)
! !     c_Rs_4 = C_Rs(r, 4)

! !     vssigmatauz = -C_IT_3*(C_Rs_2 - C_Rs_1/r)

! ! end function

! ! ! B32
! ! real(dp) function v_b_upperT_s(r) result(vsbupperT)
! !     implicit none
! !     real(dp) intent(in) :: r
! !     integer :: N

! !     C_Rs_1 = C_Rs(r, 1)

! !     vsbupperT = - C_IT_4 * C_Rs_1/r  

! ! end function

end module short_range_chiral_potentials