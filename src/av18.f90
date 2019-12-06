

!!
!> @brief      av18 potential
!!
!! Module to calculate the AV18 potential in operator basis. See derivations and details in:
!!
!! Phys.Rev. C51 (1995) 38-51
!!
!! @author Rodrigo Navarro Perez
!!
module av18


use precisions, only : dp

implicit none
! use physical_constants, only : 

private 

integer, parameter :: n_parameters = 44 !< Number of phenomenological parameters
integer, parameter :: n_operators = 18  !< Number of operators in the AV18 basis

contains

!!
!> @brief      av18 potential in operator basis
!!
!! Given a set of parameters and a radius (in units of fm) returns the nuclear part of the 
!! AV18 potential (in units of MeV) in operator basis.
!! 
!! The role of every parameter and the derivation of the potential can be found
!!
!! Phys.Rev. C51 (1995) 38-51
!!
!! @author     Rodrigo Navarro Perez
!!
function av18_operator(parameters, radius) result(potential)
    implicit none
    real(dp), intent(in) :: parameters(1:n_parameters) !< Phenomenological parameters
    real(dp), intent(in) :: radius !< radius in units of fm
    real(dp) :: potential(1:n_operators) !< AV18 potential in operator basis, units of MeV
    
end function av18_operator


end module av18