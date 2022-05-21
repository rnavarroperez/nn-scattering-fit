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

    public :: short_range_potentials

contains

subroutine short_range_potentials(r, short_lecs, v_short, d_v_short)
    implicit none
    real(dp) :: CRs, d1cRs, d2cRs, d3cRs, d4cRs ! CRs and its derivatives
    real(dp) :: C_s, C_T, C_1, C_2, C_3, C_4, C_5, C_6, C_7, D_1, D_2, D_3, D_4, D_5, D_6, D_7, &
                D_8, D_9, D_10, D_11, C_0_IV, C_0_IT, C_1_IT, C_2_IT, C_3_IT, C_4_IT, R_S ! short-range LECS
    real(dp), intent(in) :: r
    integer, parameter :: n_potentials = 19
    real(dp), intent(out), allocatable, dimension(:) :: v_short
    real(dp), intent(in), dimension(:) :: short_lecs
    real(dp), intent(out), allocatable, dimension(:,:) :: d_v_short

    allocate(v_short(1:n_potentials))
    allocate(d_v_short(1:n_potentials, 1:size(short_lecs)))

    ! initialize v_short and d_v_short arrays to contain 0s
    d_v_short = 0
    v_short = 0
  
    !unpack short-range low energy constants (short_lecs)
    C_s = short_lecs(1)
    C_T = short_lecs(2)

    C_1 = short_lecs(3)
    C_2 = short_lecs(4)
    C_3 = short_lecs(5)
    C_4 = short_lecs(6)
    C_5 = short_lecs(7)
    C_6 = short_lecs(8)
    C_7 = short_lecs(9)

    D_1 = short_lecs(10)
    D_2 = short_lecs(11)
    D_3 = short_lecs(12)
    D_4 = short_lecs(13)
    D_5 = short_lecs(14)
    D_6 = short_lecs(15)
    D_7 = short_lecs(16)
    D_8 = short_lecs(17)
    D_9 = short_lecs(18)
    D_10 = short_lecs(19)
    D_11 = short_lecs(20)

    C_0_IV = short_lecs(21)
    C_0_IT = short_lecs(22)

    C_1_IT = short_lecs(23)
    C_2_IT = short_lecs(24)
    C_3_IT = short_lecs(25)
    C_4_IT = short_lecs(26)

    R_S = short_lecs(27)

    !C_Rs and its derivatives
    CRs = C_Rs(r, R_S)
    d1cRs = d1_c_Rs(r, R_S)
    d2cRs = d2_c_Rs(r, R_S)
    d3cRs = d3_c_Rs(r, R_S)
    d4cRs = d4_c_Rs(r, R_S)

    ! short-range potentials
    v_short( 1) = C_s*CRs + C_1*(-d2cRs - 2*d1cRs/r) + D_1*(d4cRs + 4*d3cRs/r)           ! v_c
    v_short( 2) = C_2 * (-d2cRs - 2*d1cRs/r) + D_2*(d4cRs + 4*d3cRs/r)                   ! v_tau
    v_short( 3) = C_T*CRs + C_3*(-d2cRs - 2*d1cRs/r) + D_3*(d4cRs + 4*d3cRs/r)           ! v_sigma
    v_short( 4) = C_4*(-d2cRs - 2*d1cRs/r) + D_4*(d4cRs + 4*d3cRs/r)                     ! v_sigma_tau
    v_short( 5) = -C_5*(d2cRs - d1cRs/r) + D_5*(d4cRs + d3cRs/r - 6*d2cRs/r**2._dp + 6*d1cRs/r**3._dp)   ! v_t
    v_short( 6) = -C_6*(d2cRs - d1cRs/r) + D_6*(d4cRs + d3cRs/r - 6*d2cRs/r**2._dp + 6*d1cRs/r**3._dp)   ! v_t_tau
    v_short( 7) = -C_7*d1cRs/r + D_7*(d3cRs/r + 2*d2cRs/r**2._dp - 2*d1cRs/r**3._dp)     ! v_b = v_ls
    v_short( 8) = D_8*(d3cRs/r + 2*d2cRs/r**2._dp - 2*d1cRs/r**3._dp)                    ! v_b_tau = v_ls_tau
    v_short( 9) = -D_10*(d2cRs - d1cRs/r)/r**2._dp                                       ! v_q = v_l2

    v_short(11) = -D_11*(d2cRs - d1cRs/r)/r**2._dp      ! v_q_sigma or v_l2_sigma
    v_short(13) = -D_9*(d2cRs - d1cRs/r)/r**2._dp       ! v_b_b = v_ls_ls = v_ls2

    v_short(15) = C_0_IT*CRs + C_1_IT*(-d2cRs - 2*d1cRs/r)      ! v_T
    v_short(16) = C_2_IT*(-d2cRs -2*d1cRs/r)                    ! v_sigma_T
    v_short(17) = -C_3_IT*(d2cRs - d1cRs/r)                     ! v_t_T
    v_short(18) = -C_4_IT*d1cRs/r                               ! v_b_T = v_ls_T
    v_short(19) = C_0_IV*CRs                                    ! v_tau_z

    ! dv_c_s (derivatives of B11)
    d_v_short(1, 1) =  CRs
    d_v_short(1, 3) =  -d2cRs - 2*d1cRs/r
    d_v_short(1, 10) =  d4cRs + 4*d3cRs/r

    ! dv_tau_s (derivatives of B12)
    d_v_short(2, 4) = -d2cRs - 2*d1cRs/r
    d_v_short(2, 11) = d4cRs + 4*d3cRs/r

    ! dv_sigma_s (derivatives of B13)
    d_v_short(3, 2) = CRs
    d_v_short(3, 5) = -d2cRs - 2*d1cRs/r
    d_v_short(3, 12) = d4cRs + 4*d3cRs/r

    ! dv_sigma_tau_s (derivatives of B14)
    d_v_short(4, 6) = -d2cRs - 2*d1cRs/r
    d_v_short(4, 13) = d4cRs + 4*d3cRs/r

    ! dv_t_s (derivatives of B15)
    d_v_short(5, 7) = - d2cRs + d1cRs/r
    d_v_short(5, 14) = d4cRs + d3cRs/r - 6*d2cRs/r**2._dp + 6*d1cRs/r**3._dp

    ! dv_t_tau_s (derivatives of B16)
    d_v_short(6, 8) = - d2cRs + d1cRs/r
    d_v_short(6, 15) = d4cRs + d3cRs/r - 6*d2cRs/r**2._dp + 6*d1cRs/r**3._dp

    ! dv_b_s (derivatives of B17)
    d_v_short(7, 9) = - d1cRs/r
    d_v_short(7, 16) = d3cRs/r + 2*d2cRs/r**2._dp - 2*d1cRs/r**3._dp

    ! dv_b_tau_s (derivatives of B18)
    d_v_short(8, 17) = d3cRs/r + 2*d2cRs/r**2._dp - 2*d1cRs/r**3._dp

    ! dv_q_s (derivatives of B20)
    d_v_short(9, 19) = -(d2cRs - d1cRs/r)/r**2._dp

    ! dv_q_sigma_s (derivatives of B21)
    d_v_short(11, 20) = -(d2cRs - d1cRs/r)/r**2._dp

    ! dv_b_b_s (derivatives of B19)
    d_v_short(13, 18) = -(d2cRs - d1cRs/r)/r**2._dp

    ! dv_T_s (derivatives of B26)
    d_v_short(15, 22) = CRs
    d_v_short(15, 23) = -d2cRs - 2*d1cRs/r

    ! dv_sigma_T_s (derivatives of B28)
    d_v_short(16, 24) = -d2cRs -2*d1cRs/r

    ! dv_t_T_s (derivatives of B30)
    d_v_short(17, 25) = - d2cRs + d1cRs/r

    ! dv_b_T_s (derivatives of B32)
    d_v_short(18, 26) = -d1cRs/r

    ! dv_tau_z_s (derivatives of B27)
    d_v_short(19, 21) = CRs

end subroutine short_range_potentials

real(dp) function C_Rs(r, R_S) result(CRs)
    implicit none
    real(dp), intent(in) :: r, R_S ! radius and short-range LEC

    CRs = hbar_c*exp(-(r/R_S)**2._dp)/(sqrt(pi)*R_S)**3._dp

end function C_Rs
    
real(dp) function d1_c_Rs(r, R_S) result(d1cRs)
    implicit none
    real(dp), intent(in) :: r, R_S ! radius and short-range LEC

    d1cRs = -2*r*hbar_c*exp(-(r/R_S)**2._dp)/(pi**(3._dp/2._dp)*R_S**5._dp)
    
end function

real(dp) function d2_c_Rs(r, R_S) result(d2cRs)
    implicit none
    real(dp), intent(in) :: r, R_S ! radius and short-range LEC

    d2cRs = hbar_c*exp(-(r/R_S)**2._dp)*(4*r**2._dp - 2*R_S**2._dp)/(pi**(3._dp/2._dp)*R_S**7._dp)

end function

real(dp) function d3_c_Rs(r, R_S) result(d3cRs)
    implicit none
    real(dp), intent(in) :: r, R_S ! radius and short-range LEC

    d3cRs = hbar_c*exp(-(r/R_S)**2._dp)*(12*r*R_S**2._dp - 8*r**3._dp)/(pi**(3._dp/2._dp)*R_S**9._dp)

end function

real(dp) function d4_c_Rs(r, R_S) result(d4cRs)
    implicit none
    real(dp), intent(in) :: r, R_S ! radius and short-range LEC

    d4cRs = hbar_c*exp(-(r/R_S)**2._dp)*(16*r**4._dp - 48*(r*R_S)**2._dp + 12*R_S**4._dp) &
            /(pi**(3._dp/2._dp)*R_S**11._dp)

end function

end module short_range_chiral_potentials