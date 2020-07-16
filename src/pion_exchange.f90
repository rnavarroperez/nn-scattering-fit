module pion_exchange
use precisions, only : dp
use constants, only : hbar_c, mpi0=>pion_0_mass, mpic=>pion_c_mass, mpi=>pion_mass, &
    pi
use em_nn_potential, only : n_em_terms, em_potential, add_em_potential
use st_basis_2_partial_waves, only : n_st_terms, uncoupled_pot, coupled_pot

implicit none

private

public :: ope_all_partial_waves

integer, parameter :: n_fpis = 3

contains

subroutine ope_all_partial_waves(parameters, r, reaction, v_pw, dv_pw)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: reaction
    real(dp), intent(out), dimension(:, :) :: v_pw
    real(dp), intent(out), allocatable, dimension(:, :, :) :: dv_pw

    real(dp), dimension(1:n_st_terms) :: v_00, v_10, v_01, v_11, v_01_p_em
    real(dp), allocatable, dimension(:, :) :: dv_00, dv_10, dv_01, dv_11

    integer :: n_waves, j_max, n_parameters, s, t, l, j, ij, ip
    real(dp) :: v_em(1:n_em_terms)


    n_waves = size(v_pw, 1)
    j_max = size(v_pw, 2)
    n_parameters = size(parameters)

    allocate(dv_pw(1:n_parameters, 1:n_waves, 1:j_max))

    s = 0
    t = 0
    call ope_st_basis(parameters, r, reaction, s, t, v_00, dv_00)
    s = 0
    t = 1
    call ope_st_basis(parameters, r, reaction, s, t, v_01, dv_01)
    s = 1
    t = 0
    call ope_st_basis(parameters, r, reaction, s, t, v_10, dv_10)
    s = 1
    t = 1
    call ope_st_basis(parameters, r, reaction, s, t, v_11, dv_11)

    v_em = em_potential(r)
    v_01_p_em = v_01
    call add_em_potential(reaction, 0, v_em, v_01_p_em)
    ! 1s0
    l = 0
    s = 0
    j = 0
    v_pw(1, 1) = uncoupled_pot(l, s, j, v_01_p_em)
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
end subroutine ope_all_partial_waves

subroutine ope_st_basis(parameters, r, reaction, s,t, v_st, dv_st)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: reaction
    integer, intent(in) :: s
    integer, intent(in) :: t
    real(dp), intent(out), dimension(1:n_st_terms) :: v_st
    real(dp), intent(out), allocatable, dimension(:, :) ::  dv_st
    real(dp) :: y_0, t_0, y_c, t_c, y_total, t_total
    real(dp) :: f2c, f2pp, f2nn, f2np
    real(dp), dimension(1:n_fpis) :: fpi, df2c, df2pp, df2nn, df2np, dy_total, dt_total
    
    integer :: s1ds2, n_parameters, i
    s1ds2 = 4*s - 3
    n_parameters = size(parameters)
    allocate(dv_st(1:n_st_terms, 1:n_parameters))
    v_st = 0._dp
    dv_st = 0._dp
    fpi = parameters(n_parameters-2:n_parameters)
    call set_pion_couplings(fpi, f2c, f2pp, f2nn, f2np, df2c, df2pp, df2nn, df2np)
    y_0 = v_ope_s(mpi0, r)
    t_0 = v_ope_t(mpi0, r)
    select case(trim(reaction))
    case('pp')
        y_total = f2pp*y_0
        t_total = f2pp*t_0
        dy_total = df2pp*y_0
        dt_total = df2pp*t_0
    case('nn')
        y_total = f2nn*y_0 
        t_total = f2nn*t_0
        dy_total = df2nn*y_0 
        dt_total = df2nn*t_0
    case('np')
        y_c = v_ope_s(mpic, r)
        t_c = v_ope_t(mpic, r)
        y_total = f2np*y_0 + (-1)**(t+1)*2*f2c*y_c
        t_total = f2np*t_0 + (-1)**(t+1)*2*f2c*t_c
        dy_total = df2np*y_0 + (-1)**(t+1)*2*df2c*y_c
        dt_total = df2np*t_0 + (-1)**(t+1)*2*df2c*t_c
    case default
        stop 'incorrect reaction channel in ope_st_basis'
    end select
    v_st(1)  = s1ds2*y_total
    v_st(2)  = t_total
    do i = 1, n_fpis
        dv_st(1, n_parameters - 3 + i) = s1ds2*dy_total(i)
        dv_st(2, n_parameters - 3 + i) = dt_total(i)
    enddo
end subroutine ope_st_basis

subroutine set_pion_couplings(fpi, f2c, f2pp, f2nn, f2np, df2c, df2pp, df2nn, df2np)
    implicit none
    real(dp), intent(in), dimension(1:n_fpis) :: fpi
    real(dp), intent(out) :: f2c
    real(dp), intent(out) :: f2pp
    real(dp), intent(out) :: f2nn
    real(dp), intent(out) :: f2np
    real(dp), intent(out), dimension(1:n_fpis) :: df2c
    real(dp), intent(out), dimension(1:n_fpis) :: df2pp
    real(dp), intent(out), dimension(1:n_fpis) :: df2nn
    real(dp), intent(out), dimension(1:n_fpis) :: df2np

    real(dp), parameter :: small = 1.e-12_dp
    integer :: i, n_active

    df2c = 0._dp
    df2pp = 0._dp
    df2nn = 0._dp
    df2np = 0._dp

    n_active = 0
    do i = 1, n_fpis
        if (abs(fpi(i)) > epsilon(fpi(i))) n_active = n_active + 1
    enddo

    select case(n_active)
    case(0)
        f2c  =  0.075_dp
        f2pp =  0.075_dp
        f2nn =  0.075_dp
        f2np = -0.75_dp
    case(1)
        f2c  =  fpi(1)**2
        f2pp =  f2c
        f2nn =  f2c
        f2np = -f2c
        df2c(1) = 2*fpi(1)
        df2pp =  df2c
        df2nn =  df2c
        df2np = -df2c
    case(2)
        f2c  =  fpi(1)**2
        f2pp =  fpi(2)**2
        f2nn =  f2pp
        f2np = -f2c
        df2c(1)  = 2*fpi(1)
        df2pp(2) = 2*fpi(2)
        df2nn =  df2pp
        df2np = -df2c
    case(3)
        f2c  = fpi(1)**2
        f2pp = fpi(2)**2
        f2nn = fpi(3)**2
        f2np = fpi(2)*fpi(3)
        df2c(1)  = 2*fpi(1)
        df2pp(2) = 2*fpi(2)
        df2nn(3) = 2*fpi(3)
        df2np(2) = fpi(3)
        df2np(3) = fpi(2)
    end select
end subroutine set_pion_couplings

real(dp) function v_ope_s(m, r) result(v)
    implicit none
    real(dp) :: m, r
    real(dp) :: xi
    xi = m/mpic
    v = xi**2*m*yukawa_y(m/hbar_c, r)/3
end function v_ope_s

real(dp) function v_ope_t(m, r) result(v)
    implicit none
    real(dp) :: m, r
    real(dp) :: xi
    xi = m/mpic
    v = xi**2*m*yukawa_t(m/hbar_c, r)/3
end function v_ope_t


real(dp) function yukawa_y(m, r) result(y)
    implicit none
    real(dp), intent(in) :: m
    real(dp), intent(in) :: r
    y = exp(-m*r)/(m*r)
end function yukawa_y

real(dp) function yukawa_t(m, r) result(t)
    implicit none
    real(dp), intent(in) :: m
    real(dp), intent(in) :: r
    t = yukawa_y(m, r)*(1 + 3/(m*r) + 3/(m*r)**2)
end function yukawa_t

end module pion_exchange
