module observables
use av18, only: av18_all_partial_waves, default_params
use nn_phaseshifts, only: all_phaseshifts, momentum_cm
use amplitudes, only: saclay_amplitudes
use precisions, only: dp
use constants
implicit none
public observable
private

integer, parameter :: n_obs = 26
integer, parameter :: j_max = 20
character(len=4), dimension(1:n_obs), parameter :: &
       obs_types = ['DSG ','DT  ','AYY ','D   ','P   ','AZZ ','R   '&
         ,'RT  ','RPT ','AT  ','D0SK','NSKN','NSSN','NNKK','A   '&
         ,'AXX ','CKP ','RP  ','MSSN','MSKN','AZX ','AP  ','DTRT'&
         ,'SGT ','SGTT','SGTL']

contains

!!
!> @brief Calculates a NN scattering observable for a given laboratory energy
!! in MeV and scattering angle in d_egrees.
!!
!! The observable to be calculated is d_etermined by the type integer
!! and corresponds to the list of observable labels in the strObs
!! array.
!!
!! To avoid recalculating phase-shifts (the more time consuming part
!! of the calculation) the phases are only calculated if t and tprev
!! are different. At the end of the subroutine tprev is upd_ated to
!! the value of t
!!
!! The d_erivative of the observable with respect of the parameters is stored
!! on the ap array
!!
!! @author Raul L Bernal-Gonzalez
subroutine observable(t_lab, pre_t_lab, angle, type, reac, obs, d_obs)
    implicit none
    real(dp), parameter :: r_max = 12.5_dp , dr = 0.01_dp !< integration radius and step in fm
    real(dp), intent(in) :: t_lab !< laboratory energy
    real(dp), intent(inout) :: pre_t_lab !< previous value of t_lab
    real(dp), intent(in) :: angle !< scattering angle in d_egrees
    character(len=*), intent(in) :: type !< ind_ex to indicate the type of observable
    character(len=*), intent(in) :: reac !< reaction channel
    real(dp), intent(out) :: obs !< NN scattering observable
    real(dp), allocatable ,intent(out), dimension(:) :: d_obs !< d_erivitive of observbles
    real(dp), allocatable :: phases(:,:)
    real(dp), allocatable :: d_phases(:,:,:)
    complex(dp) :: a, b, c, d, e
    complex(dp), allocatable :: d_a(:), d_b(:), d_c(:), d_d(:), d_e(:) !< saclay parameters
    integer :: n_parameters !< number of parameters from the mod_el
    real(dp) :: k, theta, sg, num, denom !< place hold_er values
    real(dp), allocatable :: d_sg(:), d_num(:), d_denom(:) !< place hold_er values
    real(dp), allocatable :: dsg1(:), dsg2(:), dsg3(:), dsg4(:), dsg5(:), dsg6(:) !< place hold_er values
    integer :: i !< loop ind_ex
    save k

    ! Set number of parameters
    n_parameters = size(default_params)

    ! allocate all arrays
    allocate(d_obs(1:n_parameters))
    allocate(phases(1:5, 1:j_max))
    allocate(d_phases(1:5, 1:j_max, 1:n_parameters))
    allocate(d_sg(1:n_parameters))
    allocate(d_num(1:n_parameters))
    allocate(d_denom(1:n_parameters))
    allocate(dsg1, dsg2, dsg3, dsg4, dsg5, dsg6, mold=d_sg)


    if(t_lab /= pre_t_lab) then
        k = momentum_cm(t_lab, reac)
        call all_phaseshifts(av18_all_partial_waves, default_params, t_lab, reac, r_max, dr, k, phases, d_phases)
    end if
    theta = angle*pi/180.0_dp ! angle in d_egrees to radians
    call saclay_amplitudes(k, theta, reac, phases, d_phases, a, b, c, d, e, &
     d_a, d_b, d_c, d_d, d_e)

     ! Initialize values for calculation observable
     obs = 0.0_dp
     d_obs = 0.0_dp
     sg = (abs(a)**2 + abs(b)**2 + abs(c)**2 + abs(d)**2 + abs(e)**2)*0.5_dp
     ! seperate dsg into terms
     dsg1 = real(a)*real(d_a)
     dsg2 = aimag(a)*aimag(d_a)
     dsg3 = real(b)*real(d_b) + aimag(b)*aimag(d_b)
     dsg4 = real(c)*real(d_c) + aimag(c)*aimag(d_c)
     dsg5 = real(d)*real(d_d) + aimag(d)*aimag(d_d)
     dsg6 = real(e)*real(d_e) + aimag(e)*aimag(d_e)
     ! set d_sg
     d_sg = dsg1 + dsg2 + dsg3 + dsg4 + dsg5 + dsg6
     ! switch case for all observable
     select case (trim(type))
     case ('dsg')
         obs = sg*10.0_dp
         d_obs = d_sg*10.0_dp
     case ('dt')
         obs = 0.5_dp*(abs(a)**2-abs(b)**2+abs(c)**2-abs(d)**2+abs(e)**2)/sg
         do i = 1, n_parameters
             d_num(i) = real(a)*real(d_a(i))+aimag(a)*aimag(d_a(i)) &
             -real(b)*real(d_b(i))-aimag(b)*aimag(d_b(i)) &
             +real(c)*real(d_c(i))+aimag(c)*aimag(d_c(i)) &
             -real(d)*real(d_d(i))-aimag(d)*aimag(d_d(i)) &
             +real(e)*real(d_e(i))+aimag(e)*aimag(d_e(i))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('ayy')
         obs = 0.5_dp*(abs(a)**2-abs(b)**2-abs(c)**2+abs(d)**2+abs(e)**2)/sg
         do i = 1, n_parameters
             d_num(i) = real(a)*real(d_a(i))+aimag(a)*aimag(d_a(i)) &
             -real(b)*real(d_b(i))-aimag(b)*aimag(d_b(i)) &
             -real(c)*real(d_c(i))-aimag(c)*aimag(d_c(i)) &
             +real(d)*real(d_d(i))+aimag(d)*aimag(d_d(i)) &
             +real(e)*real(d_e(i))+aimag(e)*aimag(d_e(i))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('d')
         obs = 0.5_dp*(abs(a)**2+abs(b)**2-abs(c)**2-abs(d)**2+abs(e)**2)/sg
         do i = 1, n_parameters
             d_num(i) = real(a)*real(d_a(i))+aimag(a)*aimag(d_a(i)) &
             +real(b)*real(d_b(i))+aimag(b)*aimag(d_b(i)) &
             -real(c)*real(d_c(i))-aimag(c)*aimag(d_c(i)) &
             -real(d)*real(d_d(i))-aimag(d)*aimag(d_d(i)) &
             +real(e)*real(d_e(i))+aimag(e)*aimag(d_e(i))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('p')
         obs = real(conjg(a)*e)/sg
         do i = 1, n_parameters
             d_num(i) = real(conjg(d_a(i))*e+conjg(a)*d_e(i))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('azz')
         obs = (-real(conjg(a)*d)*cos(theta)+real(conjg(b)*c)+aimag(conjg(d)*e)*sin(theta))/sg
         do i = 1, n_parameters
             d_num(i) = -real(conjg(d_a(i))*d+conjg(a)*d_d(i)) &
             *cos(theta)+real(conjg(d_b(i))*c+conjg(b)*d_c(i)) &
             +aimag(conjg(d_d(i))*e+conjg(d)*d_e(i))*sin(theta)
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('r')
         obs = (cos(0.5_dp*theta)*(real(conjg(a)*b+conjg(c)*d)) &
        -sin(0.5_dp*theta)*(aimag(conjg(b)*e)))/sg
         do i = 1, n_parameters
             d_num(i) = cos(0.5_dp*theta)*(real(conjg(d_a(i))*b &
             +conjg(a)*d_b(i)+conjg(d_c(i))*d+conjg(c)*d_d(i))) &
             -sin(0.5_dp*theta)*(aimag(conjg(d_b(i))*e &
             +conjg(b)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('rt')
         obs = (-cos(0.5_dp*(pi+theta))*(real(conjg(a)*c)) &
         -cos(0.5_dp*(pi-theta))*(real(conjg(b)*d)) &
         +sin(0.5_dp*(pi+theta))*(aimag(conjg(c)*e)))/sg
         do i = 1, n_parameters
             d_num(i) = -cos(0.5_dp*(pi+theta))*(real(conjg(d_a(i))*c+conjg(a)*d_c(i))) &
             -cos(0.5_dp*(pi-theta))*(real(conjg(d_b(i))*d+conjg(b)*d_d(i))) &
             +sin(0.5_dp*(pi+theta))*(aimag(conjg(d_c(i))*e+conjg(c)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('rpt')
         obs = (sin(0.5_dp*(pi+theta))*(real(conjg(a)*c)) &
         +sin(0.5_dp*(pi-theta))*(real(conjg(b)*d)) &
         +cos(0.5_dp*(pi+theta))*(aimag(conjg(c)*e)))/sg
         do i = 1, n_parameters
             d_num(i) = sin(0.5_dp*(pi+theta))*(real(conjg(d_a(i))*c+conjg(a)*d_c(i))) &
             +sin(0.5_dp*(pi-theta))*(real(conjg(d_b(i))*d+conjg(b)*d_d(i))) &
             +cos(0.5_dp*(pi+theta))*(aimag(conjg(d_c(i))*e+conjg(c)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('at')
         obs = -(sin(0.5_dp*(pi+theta))*(real(conjg(a)*c)) &
         -sin(0.5_dp*(pi-theta))*(real(conjg(b)*d))&
         +cos(0.5_dp*(pi+theta))*(aimag(conjg(c)*e)))/sg
         do i = 1, n_parameters
             d_num(i) = -sin(0.5_dp*(pi+theta))*(real(conjg(d_a(i))*c+conjg(a)*d_c(i)))&
             +sin(0.5_dp*(pi-theta))*(real(conjg(d_b(i))*d+conjg(b)*d_d(i)))&
             -cos(0.5_dp*(pi+theta))*(aimag(conjg(d_c(i))*e+conjg(c)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('d0sk')
         obs = (sin(0.5_dp*(pi+theta))*(real(conjg(a)*b)) &
         -sin(0.5_dp*(pi-theta))*(real(conjg(c)*d)) &
         +cos(0.5_dp*(pi+theta))*(aimag(conjg(b)*e)))/sg
         do i = 1, n_parameters
             d_num(i) = sin(0.5_dp*(pi+theta))*(real(conjg(d_a(i))*b)+real(conjg(a)*d_b(i))) &
             -sin(0.5_dp*(pi-theta))*(real(conjg(d_c(i))*d)+real(conjg(c)*d_d(i))) &
             +cos(0.5_dp*(pi+theta))*(aimag(conjg(d_b(i))*e+conjg(b)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('nskn')
         obs = (sin(0.5_dp*(pi+theta))*(real(conjg(c)*e)) &
         -cos(0.5_dp*(pi+theta))*(aimag(conjg(a)*c)) &
         +cos(0.5_dp*(pi-theta))*(aimag(conjg(b)*d)))/sg
         do i = 1, n_parameters
         d_num(i) = sin(0.5_dp*(pi+theta))*(real(conjg(d_c(i))*e)+real(conjg(c)*d_e(i))) &
             -cos(0.5_dp*(pi+theta))*(aimag(conjg(d_a(i))*c)+aimag(conjg(a)*d_c(i))) &
             +cos(0.5_dp*(pi-theta))*(aimag(conjg(d_b(i))*d+conjg(b)*d_d(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('nssn')
         obs = (-cos(0.5_dp*(pi+theta))*(real(conjg(c)*e)) &
         -sin(0.5_dp*(pi+theta))*(aimag(conjg(a)*c)) &
         -sin(0.5_dp*(pi-theta))*(aimag(conjg(b)*d)))/sg
         do i = 1, n_parameters
             d_num(i) = -cos(0.5_dp*(pi+theta))*(real(conjg(d_c(i))*e)+real(conjg(c)*d_e(i))) &
             -sin(0.5_dp*(pi+theta))*(aimag(conjg(d_a(i))*c)+aimag(conjg(a)*d_c(i))) &
             -sin(0.5_dp*(pi-theta))*(aimag(conjg(d_b(i))*d+conjg(b)*d_d(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('nnkk')
         obs = (-cos(theta)*(real(conjg(d)*e))-sin(theta)*(aimag(conjg(a)*d)))/sg
         do i = 1, n_parameters
             d_num(i) = -cos(theta)*(real(conjg(d_d(i))*e &
             +conjg(d)*d_e(i))) &
             -sin(theta)*(aimag(conjg(d_a(i))*d+conjg(a)*d_d(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('a')
         obs = (-sin(0.5_dp*theta)*( real(conjg(a)*b+conjg(c)*d)) &
         -cos(0.5_dp*theta)*(aimag(conjg(b)*e)))/sg
         do i = 1, n_parameters
             d_num(i)=-sin(0.5_dp*theta)*(real(conjg(d_a(i))*b &
             +conjg(a)*d_b(i)+conjg(d_c(i))*d+conjg(c)*d_d(i))) &
             -cos(0.5_dp*theta)*(aimag(conjg(d_b(i))*e+conjg(b)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('axx')
         obs = (cos(theta)*real(conjg(a)*d)+real(conjg(b)*c) &
         -sin(theta)*aimag(conjg(d)*e))/sg
         do i = 1, n_parameters
             d_num(i) = cos(theta)*real(conjg(d_a(i))*d &
             +conjg(a)*d_d(i))+real(conjg(d_b(i))*c+conjg(b)*d_c(i)) &
             -sin(theta)*aimag(conjg(d_d(i))*e+conjg(d)*d_e(i))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('ckp')
         obs = aimag(conjg(d)*e)/sg
         do i = 1, n_parameters
             d_num(i) = aimag(conjg(d_d(i))*e+conjg(d)*d_e(i))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('rp')
         obs = (sin(0.5d0*theta)*( real(conjg(a)*b-conjg(c)*d)) &
         +cos(0.5d0*theta)*(aimag(conjg(b)*e)))/sg
         do i = 1, n_parameters
             d_num(i) = sin(0.5d0*theta)*( real(conjg(d_a(i))*b &
             +conjg(a)*d_b(i)-conjg(d_c(i))*d-conjg(c)*d_d(i))) &
             +cos(0.5d0*theta)*(aimag(conjg(d_b(i))*e+conjg(b)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('mssn')
         obs = (cos(0.5_dp*theta)*(real(conjg(b)*e)) &
         +sin(0.5_dp*theta)*(aimag(conjg(a)*b-conjg(c)*d)))/sg
         do i = 1,n_parameters
             d_num(i) = cos(0.5_dp*theta)*(real(conjg(d_b(i))*e+conjg(b)*d_e(i)))&
             +sin(0.5_dp*theta)*(aimag(conjg(d_a(i))*b &
             +conjg(a)*d_b(i)-conjg(d_c(i))*d-conjg(c)*d_d(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('mskn')
         obs = (-sin(0.5_dp*theta)*(real(conjg(b)*e)) &
         +cos(0.5_dp*theta)*(aimag(conjg(a)*b-conjg(c)*d)))/sg
         do i = 1, n_parameters
             d_num(i)= -sin(0.5_dp*theta)*(real(conjg(d_b(i))*e+conjg(b)*d_e(i))) &
             +cos(0.5_dp*theta)*(aimag(conjg(d_a(i))*b &
             +conjg(a)*d_b(i)-conjg(d_c(i))*d -conjg(c)*d_d(i)))
      end do
      num = obs*sg
      d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('azx')
         obs = (sin(theta)*(real(conjg(a)*d)) &
         +cos(theta)*(aimag(conjg(d)*e)))/sg
         do i = 1, n_parameters
             d_num(i)=+sin(theta)*( real(conjg(d_a(i))*d+conjg(a)*d_d(i))) &
             +cos(theta)*(aimag(conjg(d_d(i))*e+conjg(d)*d_e(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('ap')
         obs = (-sin(0.5_dp*theta)*(aimag(conjg(b)*e)) &
         +cos(0.5_dp*theta)*( real(conjg(a)*b-conjg(c)*d)))/sg
        do i = 1,n_parameters
             d_num(i)=-sin(0.5_dp*theta)*(aimag(conjg(d_b(i))*e+conjg(b)*d_e(i))) &
             +cos(0.5_dp*theta)*(real(conjg(d_a(i))*b+conjg(a)*d_b(i)-conjg(d_c(i))*d&
             -conjg(c)*d_d(i)))
         end do
         num = obs*sg
         d_obs = (d_num*sg-num*d_sg)/sg**2
     case ('dtrt')
         num = 0.5_dp*(abs(a)**2-abs(b)**2+abs(c)**2-abs(d)**2+abs(e)**2)
         denom = (-cos(0.5_dp*(pi+theta))*(real(conjg(a)*c)) &
         -cos(0.5_dp*(pi-theta))*( real(conjg(b)*d)) &
         +sin(0.5_dp*(pi+theta))*(aimag(conjg(c)*e)))
         do i = 1, n_parameters
             d_num(i)=(real(a)*real(d_a(i))+aimag(a)*aimag(d_a(i)) &
             -real(b)*real(d_b(i))-aimag(b)*aimag(d_b(i)) &
             +real(c)*real(d_c(i))+aimag(c)*aimag(d_c(i)) &
             -real(d)*real(d_d(i))-aimag(d)*aimag(d_d(i)) &
             +real(e)*real(d_e(i))+aimag(e)*aimag(d_e(i)))
             d_denom(i) = -cos(0.5_dp*(pi+theta))*(real(conjg(d_a(i))*c &
             +conjg(a)*d_c(i)))-cos(0.5_dp*(pi-theta))*(real(conjg(d_b(i))*d &
             +conjg(b)*d_d(i)))+sin(0.5_dp*(pi+theta))*(aimag(conjg(d_c(i))*e &
             +conjg(c)*d_e(i)))
         end do
         obs = num/denom
         d_obs = (d_num*denom-num*d_denom)/denom**2
     case ('sgt')
         obs = 20*pi*(aimag(a+b))/(k)
         d_obs= 20*pi*(aimag(d_a+d_b))/(k)
     case ('sgtt')
         obs =-40*pi*(aimag(a-b))/(k)
         d_obs=-40*pi*(aimag(d_a-d_b))/(k)
     case ('sgtl')
         obs =-40*pi*(aimag(c-d))/(k)
         d_obs=-40*pi*(aimag(d_c-d_d))/(k)
     case default
         stop 'WRONG TYPE OF OBSERVABLE'
     end select
end subroutine observable
end module observables
