!!
!> @brief      fundamental constants
!!
!! List of fundamental physical and mathematical constants. Unlike most other modules, all 
!! parameters are public.
!!
!! When possible, constants contain more digits than double precision, so that
!! they are rounded correctly. Single letter constants contain underscore so
!! that they do not clash with user variables ("e" and "i" are frequently used as
!! loop variables)
!!
!! @author     Rodrigo Navarro Perez
!!
module constants
use precisions, only: dp
implicit none

public

real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp !<  \f$\pi\f$
real(dp), parameter :: e_ = 2.7182818284590452353602874713527_dp !< base of the natural log, \f$ e \f$
real(dp), parameter :: euler_mascheroni = 0.577216_dp!< Euler Mascheroni constant, typically denoted by \f$ \gamma \f$ ! 0.5772156649015328606065120900824_dp 
complex(dp), parameter :: i_ = (0, 1) !< \f$ i = \sqrt{-1} \f$

! Values from National Institute of Standards and Technology (NIST)
! https://www.nist.gov/pml/fundamental-physical-constants
! Retrieved on June 6th 2019
real(dp), parameter :: hbar_c    = 197.327053_dp!< \f$\hbar c\f$ in units of MeV fm ! 197.3269804_dp  
real(dp), parameter :: proton_mass  = 938.27231_dp!< proton mass in units of MeV !938.27208816_dp 
real(dp), parameter :: neutron_mass = 939.56563_dp!< neutron mass in units of MeV !939.56542052_dp 
real(dp), parameter :: electron_mass = 0.510999_dp!< electron mass in units of MeV! 0.51099895000_dp 
real(dp), parameter :: mu_proton = 2.7928474_dp!< proton magnetic moment in units of the nuclear magneton \f$ \mu_n \f$ !2.79284734463_dp 
real(dp), parameter :: mu_neutron = -1.9130427_dp!< neutron magnetic moment in units of the nuclear magneton \f$ \mu_n \f$! -1.91304273_dp 
real(dp), parameter :: alpha = 1/137.035989_dp!< fine structure constant, dimensionless. !1/137.035999084_dp 


! Values from Particle Data Group (PDG)
! http://pdg.lbl.gov/
! M. Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018)
real(dp), parameter :: pion_c_mass =  139.5675_dp!< charged pion_mass in units of MeV !139.57061_dp 
real(dp), parameter :: pion_0_mass =  134.9739_dp!< neutral pion_mass in units of MeV !134.9770_dp  
real(dp), parameter :: pion_mass = 138.0363_dp!< (2*pion_c_mass + pion_0_mass)/3  !< average pion mass in units of MeV

! Historic, charge independent, recommended value. Might be modified later
real(dp), parameter :: f_pi_n_2 = 0.075_dp !< pion nucleon coupling constant \f$ f^2 \f$. Dimensionless
end module
