!!
!! Short range chiral potentials
!!
!! Functions from "Minimally nonlocal nucleon-nucleon potentials with chiral two-pion exchange including delta resonances"
!!
!! Phys.Rev. C91 (2015)
!!
!! @author      Ky Putnam, Rodrigo Navarro Pérez
!!
module short_range_chiral
    use precisions, only : dp
    use constants, only : pi, gA, hA, Fpi=>pion_decay_amplitude, c1, c2, c3, c4, b38=>b3_b8, &
        mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, hbar_c, delta_nucleon_mass_difference
        ! Fpi=2fpi
    implicit none
    
    private

    ! important! These functions depend on derivatives that will need to be implemented before they can be used

real(dp) function C_Rs(r, N)
    implicit none
    real(dp), intent(in) :: r ! r = radius
    integer :: N !I'm imagining have N, the # of derivatives, be an argument of the function
    ! that other functions can pass to it. Maybe using if/else statements? Unless there's a more
    ! efficient way?
    ! do I need to write these all in individually? What differentiation method should I use?


end function

!! following are the short-range potentials

! B11
real(dp) function v_s_c(r) result(vsc)
    implicit none
    real(dp), intent(in) :: r ! r = radius
    integer :: N 

    c_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsc = C_s * C_rs + C_1 * (-C_Rs_2 - 2*C_Rs_1/r) &
        + D_1*(C_Rs_4 + 4*C_Rs_3/r)

end function

! B12
real(dp) function v_s_tau(r) result(vstau)
    implicit none
    real(dp), intent(in) :: r ! r = radius
    integer :: N 

    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vstau = C_2 * (-C_Rs_2 - 2*C_Rs_1/r) + D_2*(C_Rs_4 + 4*C_Rs_3/r)

end function

! B13
real(dp) function v_s_sigma(r) result(vssigma)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vssigma = C_T*C_Rs + C_3 * (-C_Rs_2 - 2*C_Rs_1/r) + D_3*(C_Rs_4 + 4*C_Rs_3/r)


end function

! B14
real(dp) function v_s_sigma_tau(r) result(vssigmatau)
implicit none
real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vssigmatau = C_4 * (-C_Rs_2 - 2*C_Rs_1/r) + D_4*(C_Rs_4 + 4*C_Rs_3/r)

end function

! B15
real(dp) function v_s_t(r) result(vst)
implicit none
real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vst = -C_5 * (C_Rs_2 - C_Rs_1/r) + &
         D_5*(C_Rs_4 + C_Rs_3/r + 6*C_Rs_2/r**2 + 6*C_Rs_1/r**3)

end function

! B16
real(dp) function v_s_t_tau(r) result(vsttau)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsttau = -C_6 * (C_Rs_2 - C_Rs_1/r) + &
         D_6*(C_Rs_4 + C_Rs_3/r + 6*C_Rs_2/r**2 + 6*C_Rs_1/r**3)


end function

! B17
real(dp) function v_s_b(r) result(vsb)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsb = -C_7 * C_Rs_1/r + D_7*(C_Rs_3/r + 2*C_Rs_2/r**2 - 2*C_Rs_1/r**3)

end function

! B18
real(dp) function v_s_b_tau(r) result(vsbtau)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsbtau = D_8*(C_Rs_3/r + 2*C_Rs_2/r**2 - 2*C_Rs_1/r**3)

end function

! B19
real(dp) function v_s_b_b(r) result(vsbb)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsbb = -D_9*(C_Rs_2 - C_Rs_1/r)/r**2

end function

! B20
real(dp) function v_s_q(r) result(vsq)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)

    vsq = -D_10*(C_Rs_2 - C_Rs_1/r)/r**2

end function

! B21
real(dp) function v_s_q_sigma(r) result(vsqsigma)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)

    vsqsigma = -D_11*(C_Rs_2 - C_Rs_1/r)/r**2

end function

! B22
real(dp) function v_s_p(r) result(vsp)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsp = D_12*(C_Rs_2 - C_Rs_1/r)

end function

! B23
real(dp) function v_s_p_sigma(r) result(vspsigma)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vspsigma = D_13*(-C_Rs_2 - 2C_Rs_1/r)

end function

! B24
real(dp) function v_s_p_t(r) result(vspt)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vspt = -D_14*(C_Rs_2 - C_Rs_1/r)

end function

! B25
real(dp) function v_s_p_t_tau(r) result(vspttau)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vspttau = -D_15*(C_Rs_2 - C_Rs_1/r)

end function

! B26
real(dp) function v_s_upperT(r) result(vsupperT)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vsupperT = C_IV_0*C_Rs + C_IT_1*(-C_Rs_2 - 2*C_Rs_1/r)

end function

! B27
real(dp) function v_s_tau_z(r) result(vstauz)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vstauz = C_IV_0*C_Rs + C_IT_1*(-C_Rs_2 - 2*C_Rs_1/r)

end function

! B28
real(dp) function v_sigma_upperT_s(r) result(vssigmaupperT)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vssigmaupperT = C_IT_2*(-C_Rs_2 - 2*C_Rs_1/r)

end function

! B29
real(dp) function v_s_sigma_tau_z(r) result(vssigmatauz)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vssigmatauz = C_IV_2*(-C_Rs_2 - 2*C_Rs_1/r)

end function

! B30
real(dp) function v_s_t_upperT(r)
    implicit none
    real(dp), intent(in) :: r ! r = radius

    C_Rs = C_Rs(r, 0)
    c_Rs_1 = C_Rs(r, 1)
    c_Rs_2 = C_Rs(r, 2)
    c_Rs_3 = C_Rs(r, 3)
    c_Rs_4 = C_Rs(r, 4)

    vssigmatauz = -C_IT_3*(C_Rs_2 - C_Rs_1/r)

end function

! B31
real(dp) function v_s_t_tau_z(r) result(vsttauz)
    implicit none
    real(dp) intent(in) :: r
    integer :: N

    C_Rs_1 = C_Rs(r, 1)
    C_Rs_2 = C_Rs(r, 2)

    vsttauz = -C_IV_3*(C_Rs_2 - C_Rs_1/r)

end function

! B32
real(dp) function v_b_upperT_s(r) result(vsbupperT)
    implicit none
    real(dp) intent(in) :: r
    integer :: N

    C_Rs_1 = C_Rs(r, 1)

    !! NOTE: C_4 is one of the LECs... how do I implement this?
    vsbupperT = - C_IT_4 * C_Rs_1/r  

end function

! B33
real(dp) function v_b_tau_z_s(r) result(vsbtauz)
    implicit none
    real(dp) intent(in) :: r
    integer :: N

    C_Rs_1 = C_Rs(r, 1)

    !! NOTE: C_4 is one of the LECs... how do I implement this?
    vsbtauz = - C_IV_4 * C_Rs_1/r   

end function

end module short_range_chiral