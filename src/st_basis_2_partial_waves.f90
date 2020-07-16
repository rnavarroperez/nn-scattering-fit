module st_basis_2_partial_waves
use precisions, only : dp
implicit none
private

public :: n_st_terms, uncoupled_pot, coupled_pot


integer, parameter :: n_st_terms = 5 !< Number of terms in the spin-isospin basis

contains

!!
!> @brief      potential in uncoupled partial wave
!!
!! Given a potential in spin-isospin basis and a set of
!! angular momentum quantum numbers (l, s, j), calculates
!! the uncoupled potential in the corresponding partial wave
!!
!! @return     uncoupled potential in the given partial wave
!!
!! @author     Rodrigo Navarro Perez
!!
function uncoupled_pot(l, s, j, v_st) result(v_pw)
    implicit none
    integer, intent(in) :: l !< orbital angular momentum quantum number
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: j !< total angular momentum quantum number
    real(dp), intent(in) :: v_st(1:n_st_terms) !< potential in spin-isospin basis
    real(dp) :: v_pw !< uncoupled potential in the given partial wave
    integer :: ls, l2, s_12


    if (s == 0 .and. l == j) then ! singlets
        s_12 = 0
    elseif (s == 1 .and. l == j) then !triplet
        s_12 = 2
    elseif (s == 1 .and. l == 1 .and. j == 0) then !3p0
        s_12 = -4
    else
        stop 'incorrect quantum numbers in uncoupled_pot'
    endif

    ls = (j*(j + 1) - l*(l + 1) - s*(s + 1))/2 !numerator is always even, no need to worry about integer division
    l2 = l*(l + 1)
    v_pw =v_st(1) + s_12*v_st(2) + ls*v_st(3) + l2*v_st(4) + ls**2*v_st(5)
end function uncoupled_pot

!!
!> @brief      potential in uncoupled partial wave
!!
!! Given a potential in spin-isospin basis and the total angular momentum j, calculates the
!! coupled potential in the corresponding channel
!!
!! @return     coupled potential in the given j channel
!!
!! @author     Rodrigo Navarro Perez
!!
function coupled_pot(j, v_st) result(v_pw)
    implicit none
    integer, intent(in) :: j !< total angular momentum quantum number
    real(dp), intent(in) :: v_st(1:n_st_terms) !< potential in spin-isospin basis
    real(dp) :: v_pw(1:3) !< coupled potential in the j channel
    real(dp) :: s_12m, s_12, s_12p
    integer :: lm, lp, lm2, lp2, lsm, lsp

    lm = j - 1
    lp = j + 1
    lm2 = lm*(lm + 1)
    lp2 = lp*(lp + 1)
    lsm = j - 1    ! (j*(j + 1) - l*(l + 1) - s*(s + 1))/2 with l = j - 1 and s = 1
    lsp = -(j + 2) ! (j*(j + 1) - l*(l + 1) - s*(s + 1))/2 with l = j + 1 and s = 1
    s_12m = -2*(j - 1)/(2*j + 1._dp)
    s_12  = sqrt(36._dp*j*(j+1))/(2*j + 1)
    s_12p = -2*(j + 2)/(2*j + 1._dp)
    v_pw(1) = v_st(1) + s_12m*v_st(2) + lsm*v_st(3) + lm2*v_st(4) + lsm**2*v_st(5)
    v_pw(2) = s_12*v_st(2)
    v_pw(3) = v_st(1) + s_12p*v_st(2) + lsp*v_st(3) + lp2*v_st(4) + lsp**2*v_st(5)

end function coupled_pot
    
end module st_basis_2_partial_waves