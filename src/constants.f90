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

! Values from National Institute of Standards and Technology (NIST CODATA)
! https://www.nist.gov/pml/fundamental-physical-constants
! Retrieved on April 8th 2021
real(dp), parameter :: hbar_c    = 197.3269804_dp!< \f$\hbar c\f$ in units of MeV fm
real(dp), parameter :: proton_mass  = 938.27208816_dp!< proton mass in units of MeV
real(dp), parameter :: neutron_mass = 939.56542052_dp!< neutron mass in units of MeV
real(dp), parameter :: electron_mass = 0.51099895000_dp!< electron mass in units of MeV
real(dp), parameter :: mu_proton = 2.79284734463_dp!< proton magnetic moment in units of the nuclear magneton \f$ \mu_n \f$
real(dp), parameter :: mu_neutron = -1.91304273_dp!< neutron magnetic moment in units of the nuclear magneton \f$ \mu_n \f$
real(dp), parameter :: alpha = 1/137.035999084_dp!< fine structure constant, dimensionless. 


! Values from Particle Data Group (PDG)
! http://pdg.lbl.gov/
! P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020, 083C01 (2020)
real(dp), parameter :: pion_c_mass =  139.57039_dp!< charged pion_mass in units of MeV
real(dp), parameter :: pion_0_mass =  134.9768_dp!< neutral pion_mass in units of MeV
real(dp), parameter :: pion_mass = (2*pion_c_mass + pion_0_mass)/3  !< average pion mass in units of MeV

! Historic, charge independent, recommended value.
real(dp), parameter :: f_pi_n_2 = 0.075_dp !< pion nucleon coupling constant \f$ f^2 \f$. Dimensionless

! Values from "Minimally nonlocal nucleon-nucleon potentials with chiral two-pion exchange including delta resonances".
! https://journals.aps.org/prc/abstract/10.1103/PhysRevC.91.024003
! ! Low Energy Constants
real(dp), parameter :: gA = 1.29_dp !< nucleon axial coupling constant, adimensional
real(dp), parameter :: hA = 2.74_dp !< N-to-delta axial coupling constant, adimensional
real(dp), parameter :: pion_decay_amplitude = 184.80_dp !< pion decay amplitude ("Fpi") in units of MeV
real(dp), parameter :: c1 = -0.57_dp !< necessary for "subleading N2LO terms", in units of GeV^-1
real(dp), parameter :: c2 = -0.25_dp !< necessary for "subleading N2LO terms", in units of GeV^-1
real(dp), parameter :: c3 = -0.79_dp !< necessary for "subleading N2LO terms", in units of GeV^-1
real(dp), parameter :: c4 = 1.33_dp !< necessary for "subleading N2LO terms", in units of GeV^-1
real(dp), parameter :: b3_b8 = 1.40_dp !< necessary for "subleading N2LO terms", in units of GeV^-1
! !
real(dp), parameter :: delta_nucleon_mass_difference = 293.1_dp !< delta nucleon mass difference, in units of MeV

! ! Original AV18 values
! real(dp), parameter :: hbar_c    = 197.327053_dp!< \f$\hbar c\f$ in units of MeV fm
! real(dp), parameter :: proton_mass  = 938.27231_dp!< proton mass in units of MeV
! real(dp), parameter :: neutron_mass = 939.56563_dp!< neutron mass in units of MeV
! real(dp), parameter :: electron_mass = 0.510999_dp!< electron mass in units of MeV
! real(dp), parameter :: mu_proton = 2.7928474_dp!< proton magnetic moment in units of the nuclear magneton \f$ \mu_n \f$
! real(dp), parameter :: mu_neutron = -1.9130427_dp!< neutron magnetic moment in units of the nuclear magneton \f$ \mu_n \f$
! real(dp), parameter :: alpha = 1/137.035989_dp!< fine structure constant, dimensionless.
! real(dp), parameter :: pion_c_mass =  139.5675_dp!< charged pion_mass in units of MeV
! real(dp), parameter :: pion_0_mass =  134.9739_dp!< neutral pion_mass in units of MeV
! real(dp), parameter :: pion_mass = 138.0363_dp!< average pion mass in units of MeV
! real(dp), parameter :: f_pi_n_2 = 0.075_dp !< pion nucleon coupling constant \f$ f^2 \f$. Dimensionless

end module
