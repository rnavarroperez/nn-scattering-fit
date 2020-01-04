module amplitudes
use precisions, only : dp
use nn_phaseshifts, only : eta_prime
use constants, only : i_, m_p => proton_mass, hbar_c, alpha, pi, m_e => electron_mass, &
    mu_p => mu_proton
use num_recipes, only : cmplx_log_gamma, spherical_harmonic, kronecker_delta
implicit none

private

contains

! subroutine saclay_amplitudes(k_cm, theta, reaction, phases, d_phases, a, b, c, d, e, d_a, d_b, &
!     d_c, d_d, d_e)
!     implicit none
!     real(dp), intent(in) :: k_cm
!     real(dp), intent(in) :: theta
!     character(len=2), intent(in) :: reaction
!     real(dp), intent(in) :: phases(:, :)
!     real(dp), intent(in) :: d_phases(:, :, :)
!     complex(dp), intent(out) :: a
!     complex(dp), intent(out) :: b
!     complex(dp), intent(out) :: c
!     complex(dp), intent(out) :: d
!     complex(dp), intent(out) :: e
!     complex(dp), intent(out), allocatable :: d_a(:)
!     complex(dp), intent(out), allocatable :: d_b(:)
!     complex(dp), intent(out), allocatable :: d_c(:)
!     complex(dp), intent(out), allocatable :: d_d(:)
!     complex(dp), intent(out), allocatable :: d_e(:)

!     complex(dp) :: m_000, m_100, m_110, m_101, m_111, m_11m1
!     complex(dp), allocatable, dimension(:) :: d_m_000, d_m_100, d_m_110, d_m_101, d_m_111, d_m_11m1
!     integer :: n_params

!     n_params = size(d_phases, 1)
!     allocate(d_a(1: n_params))
!     d_a = (0, 0)
!     allocate(d_b, d_c, d_d, d_e, source = d_a)
    
! end subroutine saclay_amplitudes

! !PWAMp(s,ms,mj,k,theta,eta,ps,dps,matrixS,reac,iderive,m,dm)

subroutine partial_wave_amplitude_sum(s, ms, mj, k_cm, theta, reaction, phases, d_phases, m, d_m)
    implicit none
    integer, intent(in) :: s
    integer, intent(in) :: ms
    integer, intent(in) :: mj
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: theta
    character(len=2), intent(in) :: reaction
    real(dp), intent(in) :: phases(:, :)
    real(dp), intent(in) :: d_phases(:, :, :)
    complex(dp), intent(out) :: m
    complex(dp), intent(out), allocatable  :: d_m(:)

    integer :: n_params, j_max, j, l, lp
    real(dp) :: etap, sigma_0, sigma_l, sigma_lp, Ylp
    complex(dp) :: sm
    complex(dp), allocatable :: d_sm(:)

    m = (0, 0)
    n_params = size(d_phases, 1)
    allocate(d_m(1: n_params))
    d_m = (0, 0)

    if (reaction == 'pp') then
        etap = eta_prime(k_cm)
        sigma_0 = coulomb_sigma_l(0._dp, etap)
    endif

    j_max = size(phases, 2) - 1

    do j = 0, j_max
        do l = j-1, j+1
            if (reaction == 'pp') sigma_l = coulomb_sigma_l(real(l, kind=dp), etap)
            do lp = j-1, j+1
                if (reaction == 'pp') sigma_lp = coulomb_sigma_l(real(lp, kind=dp), etap)
                if (lp >= abs(mj - ms) .and. lp >= 0 .and. l >= 0) then
                    Ylp = real(spherical_harmonic(lp, mj-ms, theta, 0._dp))
                    call s_matrix(lp, l, j, s, phases, d_phases, k_cm, etap, reaction, sm, d_sm)
                endif
            enddo
        enddo
    enddo

    ! do J=0,nj-1
    !    do L=J-1,J+1
    !       if(reac.eq.'pp') call coulSigmal(real(l,kind=dp),eta,Sigl)
    !       do Lp=J-1,J+1
    !          if(reac.eq.'pp') call coulSigmal(real(lp,kind=dp),eta,Siglp)
    !          if(lp.ge.abs(Mj-MS)) then
    !             call SphericalHarmonic(lp,MJ-MS,Theta,Ylp)
    !          else
    !             Ylp = 0._dp
    !          endIf
    !          if(lp.ge.0.and.l.ge.0) then
    !             call matrixs(lp,l,j,s,ps,dps,k,eta,reac,rm,drm)
    !             m = m + 2.d0*Ylp*CG(Lp,S,J,mS,mJ)*(ii**(l-lp))&
    !                  *exp(II*(Siglp-Sig0))*rm*exp(II*(Sigl-Sig0))&
    !                  *CG(L,S,J,mJ,mJ)*sqrt(2*L+1.d0)*sqrt(4*Pi)/(2*ii*k)
    !             if (iderive.eq.1) then 
    !                dm = dm+2.d0*Ylp*CG(Lp,S,J,mS,mJ)*(ii**(l-lp))&
    !                     *exp(II*(Siglp-Sig0))*drm*exp(II*(Sigl-Sig0))&
    !                     *CG(L,S,J,mJ,mJ)*sqrt(2*L+1.d0)*sqrt(4*Pi)/&
    !                     (2*ii*k)
    !             else
    !                dm = 0._dp
    !             endif
    !          endIf
    !       enddo
    !    enddo
    ! enddo
    ! if(reac.eq.'np') then
    !    m = 0.5d0*m
    !    dm = 0.5d0*dm
    ! endif
    
end subroutine partial_wave_amplitude_sum

!  s_matrix(lp, l, j, s, phases, d_phases, k_cm, eta, reaction, sm, d_sm)

subroutine s_matrix(lp, l, j, s, phases, d_phases, k_cm, eta, reaction, sm, d_sm)
    implicit none
    integer, intent(in) :: lp
    integer, intent(in) :: l
    integer, intent(in) :: j
    integer, intent(in) :: s
    real(dp), intent(in) :: phases(:,:)
    real(dp), intent(in) :: d_phases(:, :, :)
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: eta
    character(len=2), intent(in) :: reaction
    complex(dp), intent(out) :: sm
    complex(dp), intent(out), allocatable :: d_sm(:)

    integer :: n_params
    real(dp) :: alphap, lambda, sigma_l, sigma_lambda, rho_l, t_lab, nu, tau00, tau0
    real(dp), allocatable :: mm_phases(:, :)

    n_params = size(d_phases, 1)
    allocate(d_sm(1: n_params))
    d_sm = (0, 0)

    allocate(mm_phases, mold = phases)
    mm_phases = 0

    select case(reaction)
    case('pp')
        alphap = 2*k_cm*eta/m_p*hbar_c
        lambda = (-1 + sqrt(1 + 4*l*(l+1) - 4*alpha*alphap))/2._dp
        sigma_l = coulomb_sigma_l(real(l, kind=dp), eta)
        sigma_lambda = coulomb_sigma_l(lambda, eta)
        rho_l = sigma_lambda - sigma_l + (l - lambda)*pi/2._dp

        t_lab = 2/m_p*(k_cm*hbar_c)**2
        nu = 4*m_e**2/(m_p*t_lab)
        tau00 = -alpha*eta/(6*pi)*(0.5_dp*(log(2._dp/nu))**2 + 1.7615_dp - 0.2804_dp*log(2._dp/nu))
        tau0 = alpha*eta/(3*Pi)*(5.2_dp*eta - 3.83_dp*eta*sqrt(nu/2._dp) - 0.29_dp*eta**2 &
             + (1.202_dp*eta**2 - 2.36_dp*eta*sqrt(nu/2._dp))*log(2._dp/nu)) + tau00
        if (s == 1) then
            call mm_phaseshifts(k_cm, eta, mm_phases)
        endif
    case('np')
        rho_l = 0._dp
        tau0 = 0._dp
    case default
        stop 's_matrix only takes pp and np as reaction channel'
    end select



    if (s == 0) then
        if (lp == l .and. l == j) then
            sm = (exp(2*i_*phases(1, j+1)) - 1)*exp(2*i_*(rho_l + tau0))
            d_sm = 2*i_*exp(2*i_*phases(1, j+1))*d_phases(:, 1, j+1)*exp(2*i_*(rho_l + tau0))
        else
            sm = 1 - kronecker_delta(lp, l)
            d_sm = (0, 0)
        endif
    else if (s == 1) then
        if (lp == l) then

        else if (abs(l - lp) == 2) then
            
        else
            
        endif
    else
    endif

    ! If(S.eq.1) then
    !    if(reac.eq.'pp') then
    !       fLS = -alpha*(8._dp*mup-2._dp)/(4._dp*mp**2)
    !       Call MMPhaseShifts(jmax,k,eta,fLS,MMPS)
    !    elseif(reac.eq.'np') then
    !       MMPS = 0._dp
    !    else
    !       write(*,*) 'WARNING: wrong reaction type in rmatrix' 
    !    endif
    !    If(abs(l-lp).eq.0) Then
    !       If(l.eq.j)Then
    !          Rm = (exp(2*ii*APS(J+1,2))-1.d0)*EXP(2*ii*rhol)&
    !               *exp(2*ii*MMPS(j+1,2))
    !          do iw = 1,nw
    !             do il = 1,nl
    !                drm(iw,il) = 2*ii*exp(2*ii*APS(J+1,2))&
    !                     *dAPS(j+1,2,iw,il)*EXP(2*ii*rhol)&
    !                     *exp(2*ii*MMPS(j+1,2))
    !             enddo
    !          enddo
    !          return
    !       EndIf
    !       If(l.eq.j-1) Then
    !          Rm = (cos(2*APS(j+1,4))*exp(2*ii*APS(j+1,3)) - 1.d0)&
    !               *EXP(2*ii*rhol)*exp(2*ii*MMPS(j+1,3))
    !          do iw = 1,nw
    !             do il = 1,nl
    !                drm(iw,il) = (-2*sin(2*APS(j+1,4))*dAPS(j+1,4,iw,il)&
    !                     *exp(2*ii*APS(j+1,3))+cos(2*APS(j+1,4))*2*ii&
    !                     *exp(2*ii*APS(j+1,3))*dAPS(j+1,3,iw,il))&
    !                     *EXP(2*ii*rhol)*exp(2*ii*MMPS(j+1,3))
    !             enddo
    !          enddo
    !          return
    !       EndIf
    !       If(l.eq.j+1) Then
    !          Rm = (cos(2*APS(j+1,4))*exp(2*ii*APS(j+1,5)) - 1.d0)&
    !               *EXP(2*ii*rhol)*exp(2*ii*MMPS(j+1,5))  
    !          do iw = 1,nw
    !             do il = 1,nl
    !                drm(iw,il) = (-2*dsin(2*APS(j+1,4))*dAPS(j+1,4,iw,il)&
    !                     *exp(2*ii*APS(j+1,5))&
    !                     +cos(2*APS(j+1,4))*2*ii*exp(2*ii*APS(j+1,5))&
    !                     *dAPS(j+1,5,iw,il))&
    !                     *EXP(2*ii*rhol)*exp(2*ii*MMPS(j+1,5))  
    !             enddo
    !          enddo
    !          return
    !       EndIf
    !    EndIf
    !    If(abs(lp-l).eq.2) Then
    !       if(reac.eq.'pp') then
    !          lambdap = (-1+SQRT(1+4*lp*(lp+1)-4*alpha*alphap))/2._dp
    !          Call CoulSigmal(real(lp,kind=dp),eta,Siglp)
    !          Call CoulSigmal(lambdap,eta,Siglamp)
    !          rholp = Siglamp - Siglp + (lp-lambdap)*PI/2._dp
    !       elseif(reac.eq.'np') then
    !          rholp = 0._dp
    !       else
    !          write(*,*) 'WARNING: wrong reaction type in rmatrix' 
    !       endif
    !       Rm = ii*sin(2*APS(j+1,4))*exp(ii*(APS(j+1,3)+APS(j+1,5)))&
    !            *EXP(ii*(rhol+rholp))*exp(ii*(MMPS(j+1,3)+MMPS(j+1,5)))
    !       do iw = 1,nw
    !          do il = 1,nl
    !             drm(iw,il) = ii*(2*cos(2*APS(j+1,4))*dAPS(j+1,4,iw,il)&
    !                  *exp(ii*(APS(j+1,3)+APS(j+1,5)))&
    !                  +sin(2*APS(j+1,4))*ii&
    !                  *exp(ii*(APS(j+1,3)+APS(j+1,5)))&
    !                  *(dAPS(j+1,3,iw,il)+dAPS(j+1,5,iw,il)))&
    !                  *EXP(ii*(rhol+rholp))*exp(ii*(MMPS(j+1,3)&
    !                  +MMPS(j+1,5)))
    !          enddo
    !       enddo
    !       return
    !    EndIf
    ! Endif
    ! Rm = 0._dp - deltak(lp,l)
    ! drm = 0._dp
    
end subroutine s_matrix

subroutine mm_phaseshifts(k_cm, eta, mm_phases)
    implicit none
    real(dp), intent(in) :: k_cm
    real(dp), intent(in) :: eta
    real(dp), intent(out) :: mm_phases(:, :)

    real(dp) :: f_T, f_ls, I_ll, I_lp2lp2, I_llp2
    integer :: l, j_max, i

    j_max = size(mm_phases, 2)

    f_T = -alpha*mu_p**2/(4*m_p**2)
    f_ls = -alpha*(8*mu_p - 2)/(4*m_p**2)
    mm_phases = 0._dp
    l = -1
    I_lp2lp2 = mm_coulomb_Ill(l+2, eta)

    mm_phases(5, 1) = ((2*l + 6)/(2*l + 3._dp)*f_T + (l + 3)*f_ls)*I_lp2lp2

    do i = 2, j_max - 1, 2
        l = i - 1
        I_ll = I_lp2lp2
        mm_phases(2, i) = -(2*f_T - f_ls)*I_ll
        mm_phases(3, i+1) = -(-2*l/(2*l+3._dp)*f_T + l*f_ls)*I_ll
        I_llp2 = mm_coulomb_Illp2(l, eta)
        mm_phases(4, i+1) = -(6*sqrt((l + 1._dp)*(l + 2))/(2*l + 3._dp)*f_T)*I_llp2
        I_lp2lp2 = mm_coulomb_Ill(l+2, eta)
        mm_phases(5, i+1) = ((2*l + 6)/(2*l + 3._dp)*f_T + (l + 3)*f_ls)*I_lp2lp2
    enddo
    mm_phases = m_p*k_cm*mm_phases*hbar_c

    ! integer :: j,l,i
    ! real(dp) :: Ill, Ilp2lp2, Illp2, fT
    ! fT = -alpha*mup**2/(4._dp*mp**2)
    ! MMPS = 0._dp
    ! l = -1
    ! call MMCDWBAIll(l+2,eta,Ilp2lp2)
    ! MMPS(1,5)=-(-(2*l+6)/(2*l+3._dp)*fT-(l+3)*fLS)*Ilp2lp2
    ! Ill = Ilp2lp2
    ! do i = 2, Nj-1,2
    !    l = i-1
    !    MMPS(i,2)=-(2*fT-fLS)*Ill
    !    MMPS(i+1,3)=-(-2*l/(2*l+3._dp)*fT+l*fLS)*Ill
    !    call MMCDWBAIllp2(l,eta,Illp2)
    !    MMPS(i+1,4)=-(6*sqrt((l+1._dp)*(l+2))/(2*l+3._dp)*fT)*Illp2
    !    call MMCDWBAIll(l+2,eta,Ilp2lp2)
    !    MMPS(i+1,5)=-(-(2*l+6)/(2*l+3._dp)*fT-(l+3)*fLS)*Ilp2lp2
    !    Ill = Ilp2lp2
    ! endDo
    ! MMPS = mp*k*MMPS*hc
    
end subroutine mm_phaseshifts

real(dp) function mm_coulomb_Ill(l, eta) result(I_ll)
    implicit none
    integer, intent(in) :: l
    real(dp), intent(in) :: eta
    integer :: i
    real(dp) ::  sum
    sum = 0._dp
    do i = 0, l
       sum = sum + 1._dp/(i**2 + eta**2)
    endDo
    I_ll = 1._dp/(2._dp*l*(l+1)) + (1 - pi*eta + pi*eta/tanh(pi*eta) - 2*eta**2*sum)/&
        (2._dp*l*(l+1)*(2*l+1))
end function mm_coulomb_Ill

real(dp) function mm_coulomb_Illp2(l, eta) result(I_llp2)
    implicit none
    integer, intent(in) :: l
    real(dp), intent(in) :: eta
    I_llp2 = 1._dp/(6*abs(l + 1 +i_*eta)*abs(l + 2 + i_*eta))
end function mm_coulomb_Illp2

real(dp) function coulomb_sigma_l(l, eta) result(sig_l)
    implicit none
    real(dp), intent(in) :: l
    real(dp), intent(in) :: eta
    complex(dp) :: z
    z = 1 + l + i_*eta
    sig_l = aimag(cmplx_log_gamma(z))
end function coulomb_sigma_l

    
end module amplitudes