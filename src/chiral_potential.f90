module chiral_potential

use precisions, only : dp
use delta_shell, only : nn_model
use long_range_chiral_potentials, only : long_range_potentials

implicit none

integer, parameter :: n_parameters = 28 !< Number of phenomenological parameters
integer, parameter :: n_operators = 19   !< Number of operators in the chiral potential
real(dp), parameter :: r_max = 12.5_dp !< Maximum integration radius for phase-shifts. In units of fm
real(dp), parameter :: delta_r = 1/128._dp !< Integration step for phases and the dueteron. In units of fm

real(dp), parameter, dimension(1:n_parameters) :: default_parameters = &
    [   2.936041_dp,     & !C_s
       -0.4933897_dp,    & !C_T
       -0.1013462_dp,    & !C_1
       -0.1444844_dp,    & !C_2
       -0.03647634_dp,   & !C_3
       -0.01630825_dp,   & !C_4
       -0.006658100_dp,  & !C_5
       -0.06176835_dp,   & !C_6
       -0.9578191_dp,    & !C_7
       -0.03102824_dp,   & !D_1
       -0.004438695_dp,  & !D_2
       -0.01351171_dp,   & !D_3
       -0.0007084459_dp, & !D_4
        0.01110108_dp,   & !D_5
       -0.008598857_dp,  & !D_6
       -0.05367908_dp,   & !D_7
        0.03119241_dp,   & !D_8
        0.03281636_dp,   & !D_9
       -0.08647128_dp,   & !D_10
       -0.01167788_dp,   & !D_11
        0.009575695_dp,  & !C_0_IV
        0.02194758_dp,   & !C_0_IT
       -0.001550501_dp,  & !C_1_IT
       -0.008354679_dp,  & !C_2_IT
       -0.006682746_dp,  & !C_3_IT
        0.01276971_dp,   & !C_4_IT
        0.6_dp,          & !R_s
        0.8_dp           & !R_L
    ] !< default parameters in the chiral potential

contains


subroutine set_chiral_potential(potential, parameters)
    implicit none
    type(nn_model), intent(out) :: potential
    real(dp), intent(out), allocatable, dimension(:) :: parameters

    allocate(parameters, source = default_parameters)
    potential%potential => chiral_potential_all_partial_waves
    potential%potential_components => chiral_potential_operator
    potential%display_subroutine => display_parameters
    potential%r_max = r_max
    potential%dr = delta_r
    potential%potential_type = 'local'
    potential%name = 'N2LO'
    potential%relativistic_deuteron = .False.
    potential%full_em_wave = .true.
    potential%n_components = n_operators

    ! Properties of a delta-shell potential. We set them to zero
    potential%n_lambdas = 0
    potential%dr_core = 0._dp
    potential%dr_tail = 0._dp
end subroutine set_chiral_potential

subroutine chiral_potential_all_partial_waves(parameters, r, reaction, v_pw, dv_pw)
    implicit none
    real(dp), intent(in) :: parameters(:) !< Phenomenological parameters
    real(dp), intent(in) :: r !< radius in units of fm
    character(len=2), intent(in) :: reaction !< reaction channel: np, pp, or nn
    real(dp), intent(out) :: v_pw(:, :) !< AV18 potential in all partial waves, units of MeV
    real(dp), allocatable, intent(out) :: dv_pw(:, :, :) !< derivatives of v_nn with respect to the parameters in parameters

    real(dp) :: v_nn(1:n_operators)
    real(dp), allocatable :: dv_nn(:, :)

    call chiral_potential_operator(parameters, r, v_nn, dv_nn)

    call operator_2_partial_waves(reaction, v_nn, dv_nn, v_pw, dv_pw)

    call add_em_potential_to_s_waves(r, reaction, v_pw)

end subroutine chiral_potential_all_partial_waves

subroutine chiral_potential_operator(parameters, r, v_nn, dv_nn)
    implicit none
    real(dp), intent(in)  :: parameters(:) !< Phenomenological parameters
    real(dp), intent(in)  :: r !< radius in units of fm
    real(dp), intent(out) :: v_nn(:) !< AV18 potential in operator basis, units of MeV
    real(dp), allocatable, intent(out) :: dv_nn(:, :) !< derivatives of v_nn with respect to the parameters in parameters


    ! call long_range_potentials(r, R_L, a_L, v_c, v_tau, v_sigma, v_sigma_tau, v_t, v_t_tau)
    ! call short_range_potentials(r, short_lecs, v_short, d_v_short)

end subroutine chiral_potential_operator

subroutine display_parameters(parameters, mask, unit, cv)
    implicit none
    real(dp), intent(in), dimension(:) :: parameters !< parameters for the AV18 potential
    logical, intent(in), dimension(:) :: mask !< Indicates which parameters are optimized
    integer, intent(in) :: unit !< Unit where the output is sent to. Either and already opened file or output_unit from iso_fortran_env
    real(dp), intent(in), optional, dimension(:, :) :: cv !< Covariance matrix of the parameters


end subroutine display_parameters

    
end module chiral_potential