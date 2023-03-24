!!
!> @brief       short range chiral potentials
!!
!! Module to calculate short range potentials from:
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

!!
!> @brief shortrange potentials
!!
!! Linear sums of potential components and their derivatives (derivatives used for optimization).
!!
!! Short-range potential components:
!! \f[
!! v^{c}_{S} = C_s  C_{R_{S}} + C_1 \left( -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}}\right) + D_1\left( C^{(4)}_{R_{S}} + \frac{4}{r}C^{(3)}_{R_{S}}\right) \\
!! v^{\tau}_{S} = C_2 \left(-C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}}\right) + D_2\left(C^{(4)}_{R_{S}} + \frac{4}{r}C^{(3)}_{R_{S}}\right) \\
!! v^{\sigma}_{S} = C_T C_{R_{S}} + C_3 \left(-C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}}\right) + D_3\left(C^{(4)}_{R_{S}} + \frac{4}{r}C^{(3)}_{R_{S}}\right) \\
!! v^{\sigma \tau}_{S} = C_4 \left(-C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}}\right) + D_4\left(C^{(4)}_{R_{S}} + \frac{4}{r}C^{(3)}_{R_{S}}\right) \\
!! v^{t}_S = -C_5 \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right) + D_5\left(C^{(4)}_{R_{S}} + \frac{1}{r}C^{(3)}_{R_{S}} - \frac{6}{r^2}C^{(2)}_{R_{S}} + \frac{6}{r^3}C^{(1)}_{R_{S}}\right) \\
!! v^{t\tau}_S = -C_6 \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right) + D_6\left(C^{(4)}_{R_{S}} + \frac{1}{r}C^{(3)}_{R_{S}} - \frac{6}{r^2}C^{(2)}_{R_{S}} + \frac{6}{r^3}C^{(1)}_{R_{S}}\right) \\ 
!! v^{b}_S = - C_7 \frac{1}{r}C^{(1)}{R_{S}} + D_7\left(\frac{1}{r}C^{(3)}_{R_{S}} + \frac{2}{r^2}C^{(2)}_{R_{S}} - \frac{2}{r^3}C^{(1)}_{R_{S}}\right) \\
!! v^{b\tau}_S = D_8\left(\frac{1}{r}C^{(3)}_{R_{S}} + \frac{2}{r^2}C^{(2)}_{R_{S}} - \frac{2}{r^3}C^{(1)}_{R_{S}}\right) \\
!! v^{q}_S = -D_{10} \frac{1}{r^2} \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right) \\
!!
!! v^{q \sigma}_S = -D_{11} \frac{1}{r^2} \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right) \\
!! v^{bb}_S = -D_{9} \frac{1}{r^2} \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right) \\
!!
!! v^{T}_S = C_0^{IT} C_{R_{S}} + C_1^{IT} \left( -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}} \right) \right)\\
!! v^{\sigma T}_S = C^{IT}_2 \left(-C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}} \right) \\
!! v^{tT}_S = -C_{3}^{IT} \left( C^{(2)}_{R_{S}} - \frac{1}{r} C^{(1)}_{R_{S}} \right) \\
!! v^{bT}_S = -C_{4}^{IT} \frac{1}{r} C^{(1)}_{R_{S}} \\
!! v^{\tau z}_S = C_{0}^{IV} C_{R_{S}}
!! \f]
!!
!!
!! Derivatives of short-range potential components with respect to LECs:
!!
!! Derivatives of \f$ v^{c}_{S} f\$:
!! \f[
!!  \frac{\partial v^{c}_{S}}{C_{S}} = C_{R_{S}} \\
!   \frac{\partial v^{c}_{S}}{C_{1}} = -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}} \\
!   \frac{\partial v^{c}_{S}}{D_{1}} = C^{(4)}_{R_{S}} + \frac{4}{r}C^{(3)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{\tau}_{S} f\$:
!! \f[
!!  \frac{\partial v^{\tau}_{S}}{C_{2}} = -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R{S}} \\
!!  \frac{\partial v^{\tau}_{S}}{D_{2}} = C^{(4)}_{R{S}} + \frac{4}{r}C^{(3)}_{R{S}}
!! \f]
!!
!! Derivatives of \f$ v^{\sigma}_{S} f\$:
!! \f[
!!  \frac{\partial v^{\sigma}_{S}}{C_{T}} = C_{R_{S}} \\
!!  \frac{\partial v^{\sigma}_{S}}{C_{3}} = -C^{(2)}_{R{S}} - \frac{2}{r}C^{(1)}_{R{S}}
!!  \frac{\partial v^{\sigma}_{S}}{D_{3}} = C^{(4)}_{R{S}} + \frac{4}{r}C^{(3)}_{R{S}}
!! \f]
!!
!! Derivatives of \f$ v^{\sigma \tau}_{S} f\$:
!! \f[
!!  \frac{\partial v^{\sigma \tau}_{S}}{C_{4}} = -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}} \\
!!  \frac{\partial v^{\sigma \tau}_{S}}{D_{4}} = C^{(4)}_{R_{S}} + \frac{4}{r}C^{(3)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{t}_{S} f\$:
!! \f[
!!  \frac{\partial v^{t}_{S}}{C_{5}} = C^{(2)}_{R{S}} - \frac{1}{r}C^{(1)}_{R{S}} \\
!!  \frac{\partial v^{t}_{S}}{D_{5}} = C^{(4)}_{R{S}} + \frac{1}{r}C^{(3)}_{R{S}} - \frac{6}{r^2}C^{(2)}_{R{S}} + \frac{6}{r^3}C^{(1)}_{R{S}}
!! \f]
!!
!! Derivatives of \f$ v^{t \tau}_{S} f\$:
!! \f[
!!  \frac{\partial v^{t \tau}_{S}}{C_{6}} = -\left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right) \\
!!  \frac{\partial v^{t \tau}_{S}}{D_{6}} = C^{(4)}_{R_{S}} + \frac{1}{r}C^{(3)}_{R_{S}} - \frac{6}{r^2}C^{(2)}_{R_{S}} + \frac{6}{r^3}C^{(1)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{b}_{S} f\$:
!! \f[
!!  \frac{\partial v^{b}_{S}}{C_{7}} = \frac{1}{r}C^{(1)}{R_{S}} \\
!!  \frac{\partial v^{b}_{S}}{D_{7}} = \frac{1}{r}C^{(3)}_{R_{S}} + \frac{2}{r^2}C^{(2)}_{R_{S}} - \frac{2}{r^3}C^{(1)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{b \tau}_{S} f\$:
!! \f[
!!  \frac{\partial v^{b \tau}_{S}}{D_{8}} = \frac{1}{r}C^{(3)}_{R_{S}} + \frac{2}{r^2}C^{(2)}_{R_{S}} - \frac{2}{r^3}C^{(1)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{q}_{S} f\$:
!! \f[
!!  \frac{\partial v^{q}_{S}}{D_{10}} = -\frac{1}{r^2} \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right)
!! \f]
!!
!! Derivatives of \f$ v^{q \sigma}_{S} f\$:
!! \f[
!!  \frac{\partial v^{q \sigma}_{S}}{D_{11}} = -\frac{1}{r^2} \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right)
!! \f]
!!
!! Derivatives of \f$ v^{bb}_{S} f\$:
!! \f[
!!  \frac{\partial v^{bb}_{S}}{D_{9}} = -\frac{1}{r^2} \left(C^{(2)}_{R_{S}} - \frac{1}{r}C^{(1)}_{R_{S}}\right)
!! \f]
!!
!! Derivatives of \f$ v^{T}_{S} f\$:
!! \f[
!!  \frac{\partial v^{T}_{S}}{C_0^{IT}} = C_{R_{S}} \\
!!  \frac{\partial v^{T}_{S}}{C_1^{IT}} = -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{\sigma T}_{S} f\$:
!! \f[
!!  \frac{\partial v^{\sigma T}_{S}}{C^{IT}_2} = -C^{(2)}_{R_{S}} - \frac{2}{r}C^{(1)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{t T}_{S} f\$:
!! \f[
!!  \frac{\partial v^{t T}_{S}}{C_{3}^{IT}} = -\left( C^{(2)}_{R_{S}} - \frac{1}{r} C^{(1)}_{R_{S}} \right)
!! \f]
!!
!! Derivatives of \f$ v^{b T}_{S} f\$:
!! \f[
!!  frac{\partial v^{t T}_{S}}{C_{4}^{IT}} = -\frac{1}{r} C^{(1)}_{R_{S}}
!! \f]
!!
!! Derivatives of \f$ v^{\tau z}_{S} f\$:
!! \f[
!!  \frac{\partial v^{\tau z}_{S}}{C_{0}^{IV}} = C_{R_{S}}
!! \f]
!! 
!! @author      Ky Putnam
!!
subroutine short_range_potentials(r, short_lecs, v_short, d_v_short)
    implicit none
    real(dp) :: CRs         !< short-range regulator
    real(dp) :: d1cRsor     !< 1st derivative of short-range regulator
    real(dp) :: d2cRs       !< 2nd derivative of short-range regulator
    real(dp) :: d3cRsor     !< 3rd derivative of short-range regulator
    real(dp) :: d4cRs       !< 4th derivative of short-range regulator
    real(dp) :: coeff       !< coefficient term frequently used in linear sums of linear sums
    real(dp) :: C_s, C_T, C_1, C_2, C_3, C_4, C_5, C_6, C_7, D_1, D_2, D_3, D_4, D_5, D_6, D_7, &
                D_8, D_9, D_10, D_11, C_0_IV, C_0_IT, C_1_IT, C_2_IT, C_3_IT, C_4_IT, R_S ! short-range low-energy constants (LECs)
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    integer, parameter :: n_potentials = 19
    real(dp), intent(out), allocatable, dimension(:) :: v_short         ! array of short-range potential components
    real(dp), intent(in), dimension(:) :: short_lecs                    !< array of low-energy constants
    real(dp), intent(out), allocatable, dimension(:,:) :: d_v_short     !< array of derivatives of short-range potential components with respect to the LECs

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

    ! C_Rs and its derivatives, + the coefficient term
    CRs = C_Rs(r, R_S)
    d1cRsor = d1_c_Rs_or(r, R_S) ! or for "over r"
    d2cRs = d2_c_Rs(r, R_S)
    d3cRsor = d3_c_Rs_or(r, R_S) ! or for "over r"
    d4cRs = d4_c_Rs(r, R_S)
    coeff = coefficient(r, R_S)

    ! short-range potentials
    v_short( 1) = C_s*CRs + C_1*(-d2cRs - 2*d1cRsor) + D_1*(d4cRs + 4*d3cRsor)     ! v_c
    v_short( 2) = C_2 * (-d2cRs - 2*d1cRsor) + D_2*(d4cRs + 4*d3cRsor)             ! v_tau
    v_short( 3) = C_T*CRs + C_3*(-d2cRs - 2*d1cRsor) + D_3*(d4cRs + 4*d3cRsor)     ! v_sigma
    v_short( 4) = C_4*(-d2cRs - 2*d1cRsor) + D_4*(d4cRs + 4*d3cRsor)               ! v_sigma_tau
    v_short( 5) = -C_5*(d2cRs - d1cRsor) + D_5*(d4cRs + d3cRsor - 6*coeff)   ! v_t
    v_short( 6) = -C_6*(d2cRs - d1cRsor) + D_6*(d4cRs + d3cRsor - 6*coeff)   ! v_t_tau
    v_short( 7) = -C_7*d1cRsor + D_7*(d3cRsor + 2*coeff)     ! v_b = v_ls
    v_short( 8) = D_8*(d3cRsor + 2*coeff)      ! v_b_tau = v_ls_tau
    v_short( 9) = -D_10*coeff       ! v_q = v_l2

    v_short(11) = -D_11*coeff      ! v_q_sigma or v_l2_sigma
    v_short(13) = -D_9*coeff     ! v_b_b = v_ls_ls = v_ls2

    v_short(15) = C_0_IT*CRs + C_1_IT*(-d2cRs - 2*d1cRsor)      ! v_T
    v_short(16) = C_2_IT*(-d2cRs -2*d1cRsor)                    ! v_sigma_T
    v_short(17) = -C_3_IT*(d2cRs - d1cRsor)                     ! v_t_T
    v_short(18) = -C_4_IT*d1cRsor                               ! v_b_T = v_ls_T
    v_short(19) = C_0_IV*CRs                                    ! v_tau_z

    ! dv_c_s (derivatives of B11)
    d_v_short(1, 1) =  CRs
    d_v_short(1, 3) =  -d2cRs - 2*d1cRsor
    d_v_short(1, 10) =  d4cRs + 4*d3cRsor

    ! dv_tau_s (derivatives of B12)
    d_v_short(2, 4) = -d2cRs - 2*d1cRsor
    d_v_short(2, 11) = d4cRs + 4*d3cRsor

    ! dv_sigma_s (derivatives of B13)
    d_v_short(3, 2) = CRs
    d_v_short(3, 5) = -d2cRs - 2*d1cRsor
    d_v_short(3, 12) = d4cRs + 4*d3cRsor

    ! dv_sigma_tau_s (derivatives of B14)
    d_v_short(4, 6) = -d2cRs - 2*d1cRsor
    d_v_short(4, 13) = d4cRs + 4*d3cRsor

    ! dv_t_s (derivatives of B15)
    d_v_short(5, 7) = - d2cRs + d1cRsor
    d_v_short(5, 14) = d4cRs + d3cRsor - 6*coeff

    ! dv_t_tau_s (derivatives of B16)
    d_v_short(6, 8) = - d2cRs + d1cRsor
    d_v_short(6, 15) = d4cRs + d3cRsor - 6*coeff

    ! dv_b_s (derivatives of B17)
    d_v_short(7, 9) = - d1cRsor
    d_v_short(7, 16) = d3cRsor + 2*coeff

    ! dv_b_tau_s (derivatives of B18)
    d_v_short(8, 17) = d3cRsor + 2*coeff

    ! dv_q_s (derivatives of B20)
    d_v_short(9, 19) = -coeff

    ! dv_q_sigma_s (derivatives of B21)
    d_v_short(11, 20) = -coeff

    ! dv_b_b_s (derivatives of B19)
    d_v_short(13, 18) = -coeff

    ! dv_T_s (derivatives of B26)
    d_v_short(15, 22) = CRs
    d_v_short(15, 23) = -d2cRs - 2*d1cRsor

    ! dv_sigma_T_s (derivatives of B28)
    d_v_short(16, 24) = -d2cRs -2*d1cRsor

    ! dv_t_T_s (derivatives of B30)
    d_v_short(17, 25) = - d2cRs + d1cRsor

    ! dv_b_T_s (derivatives of B32)
    d_v_short(18, 26) = -d1cRsor

    ! dv_tau_z_s (derivatives of B27)
    d_v_short(19, 21) = CRs

end subroutine short_range_potentials

!!
!> @brief short-range regulator/cutoff: \f$ C_{R_{s}}(r) \f$
!!
!! \f[ C_{R_{s}}(r) = \frac{e^{-(r/R_{S})^2}}{(\sqrt{\pi}R_{S})^3}  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function C_Rs(r, R_S) result(CRs)
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_S !< short-range regulator (in fm)

    CRs = hbar_c*exp(-(r/R_S)**2._dp)/(sqrt(pi)*R_S)**3._dp

end function C_Rs
    
!!
!> @brief 1st derivative of the short-range regulator/cutoff: \f$ \frac{dC_{R_{s}}(r)}{dr} \f$
!!
!! \f[ \frac{dC_{R_{s}}(r)}{dr} = -\frac{2e^{-(r/R_{S})^2}}{\pi^{3/2}R_{S}^{5}}  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function d1_c_Rs_or(r, R_S) result(d1cRs) ! or is "over r" -- derivative has been divided by  1/r
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_S !< short-range regulator (in fm)

    d1cRs = -2*hbar_c*exp(-(r/R_S)**2._dp)/(pi**(3._dp/2._dp)*R_S**5._dp)
    
end function d1_c_Rs_or

!!
!> @brief 2nd derivative of the short-range regulator/cutoff: \frac{d^{2}C_{R_{s}}(r)}{dr^{2}} \f$
!!
!! \f[ \frac{d^{2}C_{R_{s}}(r)}{dr^{2}} = \frac{1}{\pi^{3/2}R_S^{7}}e^{-(r/R_{S})^2}\left(4r^2-2R_S^2\right)  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function d2_c_Rs(r, R_S) result(d2cRs)
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_S !< short-range regulator (in fm)

    d2cRs = hbar_c*exp(-(r/R_S)**2._dp)*(4*r**2._dp - 2*R_S**2._dp)/(pi**(3._dp/2._dp)*R_S**7._dp)

end function d2_c_Rs

!!
!> @brief 3rd derivative of the short-range regulator/cutoff: \frac{d^{3}C_{R_{s}}(r)}{dr^{3}} \f$
!!
!! \f[ \frac{d^{3}C_{R_{s}}(r)}{dr^{3}} = \frac{ e^{-\left(r/R_{S}\right)^2} \left(12 R_{S}^{2}-8 r^{2}\right)}{\pi^{3/2} R_{S}^9}  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function d3_c_Rs_or(r, R_S) result(d3cRs) ! or is "over r" -- derivative has been divided by  1/r
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_S !< short-range regulator (in fm)

    d3cRs = hbar_c*exp(-(r/R_S)**2._dp)*(12*R_S**2._dp - 8*r**2._dp)/(pi**(3._dp/2._dp)*R_S**9._dp)

end function d3_c_Rs_or


!!
!> @brief 4th derivative of the short-range regulator/cutoff: \frac{d^4 C_{R_s}(r)}{dr^4} \f$
!!
!! \f[ \frac{d^4 C_{R_s}(r)}{dr^4} = \frac{1}{\pi^{3/2} R_S^{11}} e^{-\left(r/R_{S}\right)^2} \left(16r^4 - 48(rR_S)^2 + 12R_S^4\right)  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function d4_c_Rs(r, R_S) result(d4cRs)
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_S !< short-range regulator (in fm)

    d4cRs = hbar_c*exp(-(r/R_S)**2._dp)*(16*r**4._dp - 48*(r*R_S)**2._dp + 12*R_S**4._dp) &
            /(pi**(3._dp/2._dp)*R_S**11._dp)

end function d4_c_Rs

!!
!> @brief coefficient of some short-range potential function terms
!!
!! factored to make code implementation easier
!! \f[ \mathrm{coeff} = \frac{4e^{-(r/R_{S})^2}}{\pi^{3/2}R_{S}^7}  \f]
!!
!! @author      Ky Putnam
!!
real(dp) function coefficient(r, R_S) result(coeff)
    implicit none
    real(dp), intent(in) :: r   !< point at which potential will be evaluated (in fm)
    real(dp), intent(in) :: R_S !< short-range regulator (in fm)

    coeff = 4*hbar_c*exp(-(r/R_S)**2._dp)/(pi**(3._dp/2._dp)*R_S**7._dp)

end function coefficient


end module short_range_chiral_potentials