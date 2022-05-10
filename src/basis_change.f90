module basis_change
use precisions, only : dp
implicit none
private

public :: n_st_terms, uncoupled_pot, coupled_pot, st_2_pw_basis, operator_2_partial_waves


integer, parameter :: n_st_terms = 5 !< Number of terms in the spin-isospin basis

contains

subroutine operator_2_partial_waves(reaction, v_nn, dv_nn, v_pw, dv_pw)
    implicit none
    character(len=2), intent(in) :: reaction
    real(dp), intent(in), dimension(:) :: v_nn
    real(dp), intent(in), dimension(:, :) :: dv_nn
    real(dp), intent(out), dimension(:, :) :: v_pw
    real(dp), intent(out), allocatable, dimension(:, :, :) :: dv_pw

    integer :: tz1, tz2
    real(dp) :: v_00(1:n_st_terms), v_01(1:n_st_terms), v_10(1:n_st_terms), v_11(1:n_st_terms)
    real(dp), allocatable, dimension(:, :) :: dv_00, dv_01, dv_10, dv_11

    allocate(dv_pw(1:size(dv_nn, 2), 1:size(v_pw, 1), 1:size(v_pw, 2)))
    allocate(dv_00(1:n_st_terms, 1:size(dv_pw,1)))
    dv_00 = 0
    allocate(dv_01, dv_10, dv_11, source=dv_00)

    v_pw = 0
    dv_pw = 0

    select case (trim(reaction))
    case ('pp')
        tz1 = 1
        tz2 = 1
    case ('np')
        tz1 = -1
        tz2 = 1
    case ('nn')
        tz1 = -1
        tz2 = -1
    case default
        stop 'incorrect reaction channel in av18_all_partial_waves'
    end select

    v_01 = operator_2_st_basis(tz1, tz2, 0, 1, v_nn)
    v_11 = operator_2_st_basis(tz1, tz2, 1, 1, v_nn)
    dv_01 = d_operator_2_st_basis(tz1, tz2, 0, 1, dv_nn)
    dv_11 = d_operator_2_st_basis(tz1, tz2, 1, 1, dv_nn)

    v_00 = 0._dp
    v_10 = 0._dp
    dv_00 = 0._dp
    dv_10 = 0._dp
    
    if (tz1*tz2 == -1) then
        v_00 = operator_2_st_basis(tz1, tz2, 0, 0, v_nn)
        v_10 = operator_2_st_basis(tz1, tz2, 1, 0, v_nn)
        dv_00 = d_operator_2_st_basis(tz1, tz2, 0, 0, dv_nn)
        dv_10 = d_operator_2_st_basis(tz1, tz2, 1, 0, dv_nn)
    endif

    call st_2_pw_basis(reaction, v_00, v_01, v_10, v_11, dv_00, dv_01, dv_10, dv_11, v_pw, dv_pw)

    
end subroutine operator_2_partial_waves

!!
!> @brief      transform potential from operator to st basis
!!
!! Given a potential in the AV18 operator basis and the corresponding spin and isospin quantum
!! numbers, returns the potential in the spin-isospin basis.
!!
!! The five terms in the basis are central, tensor, spin-orbit, l squared, and spin-orbit squared
!!
!! @return     potential in spin-isospin basis
!!
!! @author     Rodrigo Navarro Perez
!!
function operator_2_st_basis(tz1, tz2, s, t, v_op) result(v_st)
    implicit none
    integer, intent(in) :: tz1 !< isospin projected in the z direction for the first particle
    integer, intent(in) :: tz2 !< isospin projected in the z direction for the second particle
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: t !< isospin quantum number 
    real(dp), intent(in) :: v_op(:) !< potential in the AV18 operator basis
    real(dp) :: v_st(1:n_st_terms) !< potential in the spin-isospin basis
    integer :: s1ds2, t1dt2,t12
    s1ds2 = 4*s - 3
    t1dt2 = 4*t - 3
    t12 = 3*tz1*tz2 - t1dt2
    ! central term, c
    v_st(1) = v_op(1) + t1dt2*v_op(2) + s1ds2*v_op(3) + s1ds2*t1dt2*v_op(4) + t12*v_op(15) & 
        + s1ds2*t12*v_op(16) + (tz1+tz2)*v_op(19)
    ! tensor term, t
    v_st(2) = v_op(5) + t1dt2*v_op(6) + t12*v_op(17)
    ! spin-orbit term, ls
    v_st(3) = v_op(7) + t1dt2*v_op(8) + t12*v_op(18)
    ! l squared term, l2
    v_st(4) = v_op(9) + t1dt2*v_op(10) + s1ds2*v_op(11) + s1ds2*t1dt2*v_op(12)
    ! spin-orbit squared term, ls2
    v_st(5) = v_op(13) + t1dt2*v_op(14)
end function operator_2_st_basis

!!
!> @brief      transform derivatives potential from operator to st basis
!!
!! Given the derivatives of potential in the AV18 operator basis with respect to phenomenological
!! parameters and the corresponding spin and isospin quantum numbers, returns the derivatives in the
!! spin-isospin basis.
!!
!! The five terms in the basis are central, tensor, spin-orbit, l squared, and spin-orbit squared
!!
!! @return     derivatives of the potential in spin-isospin basis
!!
!! @author     Rodrigo Navarro Perez
!!
function d_operator_2_st_basis(tz1, tz2, s, t, dv_op) result(dv_st)
    implicit none
    integer, intent(in) :: tz1 !< isospin projected in the z direction for the first particle
    integer, intent(in) :: tz2 !< isospin projected in the z direction for the second particle
    integer, intent(in) :: s !< spin quantum number
    integer, intent(in) :: t !< isospin quantum number 
    real(dp), intent(in) :: dv_op(:, :) !< derivatives of the potential in the AV18 operator basis
    real(dp), allocatable :: dv_st(:, :) !< derivatives of the potential in the spin-isospin basis
    integer :: n_p, n_o, i
    n_p = size(dv_op,2)
    n_o = size(dv_op,1)
    ! if (n_o /= n_operators) stop 'incorrect number of operators in d_operator_2_st_basis'
    allocate(dv_st(1:n_st_terms, 1:n_p))
    do i = 1, n_p
        dv_st(:, i) = operator_2_st_basis(tz1, tz2, s, t, dv_op(:, i))
    enddo

end function d_operator_2_st_basis



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

subroutine st_2_pw_basis(reaction, v_00, v_01, v_10, v_11, dv_00, dv_01, dv_10, dv_11, v_pw, dv_pw)
    implicit none
    character(len=2), intent(in) :: reaction
    real(dp), intent(in), dimension(:) :: v_00
    real(dp), intent(in), dimension(:) :: v_01
    real(dp), intent(in), dimension(:) :: v_10
    real(dp), intent(in), dimension(:) :: v_11
    real(dp), intent(in), dimension(:, :) :: dv_00
    real(dp), intent(in), dimension(:, :) :: dv_01
    real(dp), intent(in), dimension(:, :) :: dv_10
    real(dp), intent(in), dimension(:, :) :: dv_11
    real(dp), intent(out), dimension(:, :) :: v_pw
    real(dp), intent(out), allocatable, dimension(:, :, :) :: dv_pw

    integer :: l, s, j, t, ip, ij
    integer :: n_waves, j_max, n_parameters

    n_waves = size(v_pw, 1)
    j_max = size(v_pw, 2)
    n_parameters = size(dv_00, 2)

    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max))
    dv_pw = 0

    ! 1s0
    l = 0
    s = 0
    j = 0
    v_pw(1, 1) = uncoupled_pot(l, s, j, v_01)
    do ip = 1, n_parameters
        dv_pw(ip, 1, 1) = uncoupled_pot(l, s, j, dv_01(:, ip))
    enddo

    ! 3p0
    l = 1
    s = 1
    j = 0
    v_pw(5, 1) = uncoupled_pot(l, s, j, v_11)
    do ip = 1, n_parameters
        dv_pw(ip, 5, 1) = uncoupled_pot(l, s, j, dv_11(:, ip))
    enddo

    ! everything with j >= 1
    do ij = 2, j_max
        j = ij - 1
        l = j
        ! singlets
        s = 0
        t = 1 - mod(l+s, 2)
        if (t == 1) then
            v_pw(1, ij) = uncoupled_pot(l, s, j, v_01)
            do ip = 1, n_parameters
                dv_pw(ip, 1, ij) = uncoupled_pot(l, s, j, dv_01(:, ip))
            enddo
        elseif (trim(reaction) == 'np') then ! only present in np
            v_pw(1, ij) = uncoupled_pot(l, s, j, v_00)
            do ip = 1, n_parameters
                dv_pw(ip, 1, ij) = uncoupled_pot(l, s, j, dv_00(:, ip))
            enddo
        endif
        ! triplets
        s = 1
        t = 1 - mod(l+s, 2)
        if (t == 1) then
            v_pw(2, ij) = uncoupled_pot(l, s, j, v_11)
            do ip = 1, n_parameters
                dv_pw(ip, 2, ij) = uncoupled_pot(l, s, j, dv_11(:, ip))
            enddo
        elseif (trim(reaction) == 'np') then !only present in np
            v_pw(2, ij) = uncoupled_pot(l, s, j, v_10)
            do ip = 1, n_parameters
                dv_pw(ip, 2, ij) = uncoupled_pot(l, s, j, dv_10(:, ip))
            enddo
        endif
        ! coupled channels
        t = 1 - mod(j,2)
        if (t == 1) then
            v_pw(3:5, ij) = coupled_pot(j, v_11)
            do ip = 1, n_parameters
                dv_pw(ip, 3:5, ij) = coupled_pot(j, dv_11(:, ip))
            enddo
        elseif (trim(reaction) == 'np') then !only present in np
            v_pw(3:5, ij) = coupled_pot(j, v_10)
            do ip = 1, n_parameters
                dv_pw(ip, 3:5, ij) = coupled_pot(j, dv_10(:, ip))
            enddo
        endif
    enddo
   
end subroutine st_2_pw_basis
    
end module basis_change