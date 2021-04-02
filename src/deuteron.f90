!!
!> @brief      Calculating the deuteron
!!
!! Set of functions and subroutines to calculate the 
!! binding energy of the deuteron with a given local nn interaction
!!
!! The derivatives of the binding energy with respect of the potential
!! parameters are also calculated
!!
!! @author     Rodrigo Navarro Perez
!!
module deuteron
use precisions, only : dp
use delta_shell, only : nn_model, all_delta_shells
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi, alpha
use ode_solver, only : solve_runge_kutta_4
use quadrature, only : booles_quadrature
implicit none

private

public :: binding_energy

character(len=2), parameter :: channel = 'np' !< reaction channel. For calling all_delta_shells
real(dp), parameter :: k_cm = 0._dp !< center of mass momentum. For calling all_delta_shells
integer, parameter :: j_max = 2 !< maximum angular momentum index number. For calling all_delta_shells
integer, parameter :: n_waves = 5 !< number of waves per angular momentum quantum number 
real(dp), parameter :: reduced_mass = m_p*m_n/(m_p + m_n) !< reduced mass of a proton neutron pair

!!
!> @brief      Piece wise deuteron wave function
!!
!! All the information necessary to characterize
!! a piecewise bound solution of the Schrödinger equation in
!! in the Deuteron channel with a delta shell potential.
!!
!! @author     Rodrigo Navarro Perez
!!
type :: ds_wave_function
    integer :: n_radii !< number of radii where a delta shell is located
    real(dp) :: gamma !< wave number. In fm\f$^{-1}\f$
    real(dp), allocatable, dimension(:) :: a_s !< regular components in the S channel
    real(dp), allocatable, dimension(:) :: b_s !< irregular components in the S channel
    real(dp), allocatable, dimension(:) :: a_d !< regular components in the D channel
    real(dp), allocatable, dimension(:) :: b_d !< irregular components in the D channel
    real(dp), allocatable, dimension(:) :: radii !< radii where the delta shells a are located. In fm
    real(dp) :: dr !< distance between the delta shells. In fm
end type ds_wave_function

!!
!> @brief      Evaluated deuteron wave function
!!
!! A wave function obtained by solving Schrödinger equation in
!! in the Deuteron channel with a local potential
!!
!! @author     Rodrigo Navarro Perez
!!
type :: local_wave_function
    real(dp) :: gamma !< wave number. In fm\f$^{-1}\f$
    real(dp), allocatable, dimension(:) :: radii !< radii where the wave function is evaluated. In fm
    real(dp), allocatable, dimension(:) :: s_wave !< Wave function in the S channel
    real(dp), allocatable, dimension(:) :: d_wave !< Wave function in the D channel
end type local_wave_function

contains

!!
!> @brief      Binding energy of the deuteron
!!
!! Given a nn local model (local potential, maximum integration radius,
!! integration step, and potential parameters), calculates the deuteron
!! binding energy and its derivatives with respect to the potential parameters
!!
!! Once the wave number \f$k\f$ is determined, the binding energy can be calculated
!! using relativistic kinematics (see Phys. Rev. C 48 (1993), 792-815)
!! \f[ B = M_p + M_n - \sqrt{M_p^2 - k^2} - \sqrt{M_n^2 - k^2}  \f],
!! or non-relativistic kinematics (see Phys. Rev. C 51 (1995), 38-51)
!! \f[ B = \frac{k^2}{2 M_{np}}  \f],
!! where \f$ k \f$ has to be in units of MeV
!!
!! The partial derivatives of the binding energy with respect of the potential
!! parameters are calculated using the Feynman-Hellman theorem (see 
!! Phys. Rev. 56, 340 ) where if \f$ \hat{H} | \psi \rangle = -k^2 | \psi \rangle \f$,
!! then the partial derivative of \f$-k^2\f$ with respect of one of the potential
!! parameters \f$p_j\f$ is given by
!! \f[ \frac{\partial ( -k^2)}{\partial p_j} =  \left\langle \psi \left| \frac{\partial \hat{H}}{\partial p_j }
!!   \right| \psi \right\rangle. \f] 
!! In the case of the \f$^3S_1\f$-\f$^3D_1\f$ coupled channel, the Hamiltonian \f$\hat{H}\f$
!! is a \f$2 \times 2\f$ matrix with the \f$V_{\epsilon_1}\f$ potential in the off 
!! diagonal elements, and the \f$ | \psi \rangle \f$ ket is a column vector with the 
!! S and D wave functions. 
!!
!! When a sum of delta shells representation is used for the NN potential, the 
!! integral for the expectation value \f$  \left\langle \psi \left| \frac{\partial \hat{H}}{\partial p_j }
!!   \right| \psi \right\rangle \f$ becomes straightforward
!!  
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine binding_energy(model, parameters, be, dbe)
    implicit none
    type(nn_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) ::  parameters !< parameters of the nn interaction
    real(dp), intent(out) :: be !< deuteron binding energy. In MeV
    real(dp), intent(out), allocatable, dimension(:) :: dbe !< derivatives of the deuteron binding energy with respect of the potential parameters

    real(dp) :: k
    
    k = wave_number(model, parameters)
    if (model%relativistic_deuteron) then
        be = m_p + m_n - sqrt(m_p**2 - (k*hbar_c)**2) - sqrt(m_n**2 - (k*hbar_c)**2)
    else
        be = (k*hbar_c)**2/(2*reduced_mass)
    endif

    select case(trim(model%potential_type))
    case('local')
        dbe = local_deuteron_derivatives(k, model, parameters)
    case('delta_shell')
        dbe = ds_deuteron_derivatives(k, model, parameters)
    case default
        stop 'unrecognized potential type in binding_energy'
    end select
    
end subroutine binding_energy

!!
!> @brief      Determines the wave number of the deuteron
!!
!! Given a nn local model (local potential, maximum integration radius,
!! integration step, and potential parameters), uses a simple secant
!! method to determine the value of wave number \f$k\f$ that results
!! in a wave function where the irregular component of the S wave vanishes
!! near \f$r=0\f$. 
!!
!! The secant method stops when the change in the trial wave number
!! is smaller than the optional argument tolerance. The default value
!! for the tolerance is \f$10^{-10}\f$
!!
!! @return     wave number of the deuteron
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function wave_number(model, parameters, tolerance) result(r)
    implicit none
    type(nn_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) :: parameters !< parameters of the nn interaction
    real(dp), intent(in), optional :: tolerance !< tolerance in the secant method

    real(dp), parameter :: default_tol = 1.e-10_dp
    real(dp), parameter :: k_init_l = 0.18_dp
    real(dp), parameter :: k_init_h = 0.30_dp
    real(dp) :: k_low
    real(dp) :: k_high
    real(dp) :: bs_low, bs_high, k_0, bs_0, slope, dk
    real(dp) :: local_tolerance

    if (present(tolerance)) then
        local_tolerance = tolerance
    else
        local_tolerance = default_tol
    end if
    k_low = k_init_l
    k_high = k_init_h
    bs_low = irregularity(k_low, model, parameters)
    bs_high = irregularity(k_high, model, parameters)

    if (sign(1._dp, bs_low) == sign(1._dp, bs_high)) then
        print*, k_low, bs_low,  k_high, bs_high     
        stop 'deuteron wave function irregularity is not bracketed in wave_number function'
    endif
    do
        slope = (bs_low - bs_high)/(k_low - k_high)
        dk = -bs_low/slope
        k_0 = -(bs_low - slope*k_low)/slope
        bs_0 = irregularity(k_0, model, parameters)
        if ((abs(k_low-k_0) < local_tolerance) .or. (abs(k_high-k_0) < local_tolerance)) exit
        if (sign(1._dp, bs_low) == sign(1._dp, bs_0)) then
            k_low = k_0
            bs_low = bs_0
        else
            k_high = k_0
            bs_high = bs_0
        endif
    enddo
    r = k_0
end function wave_number

!!
!> @brief      Irregular component of the S wave solution
!!
!! After integrating the deuteron wave function, returns the irregular component
!! of the S wave solution when \f$r \rightarrow 0\f$
!!
!! If the NN model consists of a local potential, the S and D waves are integrated
!! directly with a 4th order Runge Kutta and the irregularity is simply the 
!! S wave at \f$r = 0\f$
!!
!! If the NN model is a delta shell potential, the wave function is integrated piecewise
!! where each piece consists of the linear combinations
!!
!! \f[ v_i(r) = A_{S,i} \tilde{j}_0(k r) +  B_{S,i} \tilde{y}_0(k r) \f]
!! \f[ w_i(r) = A_{D,i} \tilde{j}_2(k r) +  B_{D,i} \tilde{y}_2(k r) \f]
!!
!! In this case the irregularity in the S wave when \f$r = 0\f$ is given by 
!! \f$ B_{S,0} \f$
!!
!!
!! @return     Irregular component of the S wave when \f$r \rightarrow 0\f$
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function irregularity(k, model, parameters) result(r)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    type(nn_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) :: parameters !< parameters of the nn interaction

    type(ds_wave_function) :: ds_deuteron_wf
    type(local_wave_function) :: local_deuteron_wf

    r = 0
    select case(trim(model%potential_type))
    case('local')
        local_deuteron_wf = wave_function_local(k, model, parameters)
        r = local_deuteron_wf%s_wave(1)
    case('delta_shell')
        ds_deuteron_wf = wave_function_ds(k, model, parameters)
        r = ds_deuteron_wf%b_s(0)
    case default
        stop 'unrecognized potential type in irregularity'
    end select
end function irregularity


!!
!> @brief      Direct integration of the deuteron with a local potential
!!
!! Given an object of type nn_model (a local potential,
!! a maximum radius of integration and an integration step), determines
!! the deuteron wave function for a specific wave number \f$k\f$
!!
!! The bound Schrödinger equation is integrated directly using a 4th order
!! Runge Kutta algorithm. In order to guarantee a bound wave function we
!! integrate backwards starting at a radius \f$r_{\rm max} \f$ where the 
!! potential is numerically zero and  plug in the asymptotic wave functions 
!! \f$ v(r) \propto e^{-k r}  \f$ and \f$ w(r) \propto e^{-k r}(1 + 3/(kr) + 3/(kr)^2)\f$
!! 
!! The final wave function will be a linear combination of what we call
!! an \f$ \alpha \f$ solution, resulting from integrating backwards starting with
!! \f$v(r_{\rm max}) = e^{-k r_{\rm max}}, \quad w(r_{\rm max}) = 0 \f$, and a \f$ \beta \f$ solution, 
!! resulting from integrating backwards with
!! \f$v(r_{\rm max}) = 0, \quad w(r_{\rm max}) = e^{-k r_{\rm max}}(1 + 3/(k r_{\rm max}) + 3/(k r_{\rm max})^2) \f$.
!!
!! The \f$\beta\f$ parameter in the final linear combination
!! \f[v(r) = v_\alpha(r) + \beta v_\beta(r) \\
!!     w(r) = w_\alpha(r) + \beta w_\beta(r) \f]
!! is determined by making sure that the D wave \f$ w(r)\f$ 
!! vanishes at \f$ r=0 \f$
!! (.i.e \f$w_\alpha(0) + \beta w_\beta(0) = 0\f$).
!!
!! While technically there should be an \f$ \alpha \f$ parameter in the final
!! linear combination, it can be arbitrarily set to unity since the wave function
!! can be normalized later.
!!
!! @returns    An object of type local_wave_function 
!!
!! @author Rodrigo Navarro Perez
!!
type(local_wave_function) function wave_function_local(k, model, parameters) result(wf)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    type(nn_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) :: parameters !< parameters of the nn interaction

    integer, parameter :: n_variables = 4
    real(dp), dimension(1:n_variables) :: boundary_condition_alpha, boundary_condition_beta
    real(dp), allocatable, dimension(:) :: work
    real(dp), allocatable, dimension(:) :: radii
    real(dp), allocatable, dimension(:, :) :: alpha_solution, beta_solution
    real(dp) :: r_max, r_min, beta
    integer :: n_points
    logical :: first_warning = .True.
    save first_warning

    wf%gamma = k

    allocate(work(1: size(parameters) + 1))
    work(1) = k
    work(2:) = parameters

    r_max = model%r_max
    r_min = 0._dp
    n_points = int((r_max-r_min)/model%dr)

    ! We want to ensure that the number of point will be a multiple of 4 so that 
    ! we can later use Boole's quadrature
    if (mod(n_points,4) /= 0) then
        if (first_warning) then
            print*, "r_max and dr from the model in deuteron_wf_local don't result in a"
            print*, "number of points that is a multiple of 4."
            print*, "adjusting from:", n_points, "to:", n_points + 4 - mod(n_points,4), "points"
            print*, "to guarantee proper Boole's rule integration for normalization and derivatives"
            first_warning = .False.
        endif
        n_points = n_points + 4 - mod(n_points,4)
    else
        first_warning = .True.
    endif

    boundary_condition_alpha = 0._dp
    boundary_condition_beta  = 0._dp
    boundary_condition_alpha(1) = v_asymptotic(k, r_max)
    boundary_condition_beta( 2) = w_asymptotic(k, r_max)
    boundary_condition_alpha(3) = vprime_asymptotic(k, r_max)
    boundary_condition_beta( 4) = wprime_asymptotic(k, r_max)

    call solve_runge_kutta_4(deuteron_ode, work, model, r_max, r_min, n_points, boundary_condition_alpha, radii, alpha_solution)
    call solve_runge_kutta_4(deuteron_ode, work, model, r_max, r_min, n_points, boundary_condition_beta,  radii, beta_solution)

    !since we integrated Runge Kutta backwards (from r_max to r_min) the last element in the solution is 
    !the corresponding wave function evaluated at r_min = 0

    beta = -alpha_solution(2, n_points)/beta_solution(2, n_points)

    allocate(wf%radii(1:n_points+1))
    allocate(wf%s_wave, wf%d_wave, mold = wf%radii)
    wf%radii(1:n_points) = radii(n_points:1:-1)
    wf%s_wave(1:n_points) = alpha_solution(1, n_points:1:-1) + beta*beta_solution(1, n_points:1:-1)
    wf%d_wave(1:n_points) = alpha_solution(2, n_points:1:-1) + beta*beta_solution(2, n_points:1:-1)

    wf%radii(n_points+1) = r_max
    wf%s_wave(n_points+1) = boundary_condition_alpha(1)
    wf%d_wave(n_points+1) = beta*boundary_condition_beta(2)
    
end function wave_function_local

!!
!> @brief      Asymptotic S wave
!! 
!! Asymptotic S wave in the deuteron wave function
!!
!! @return     \f$ e^{- k r} \f$
!!
!! @author Rodrigo Navarro Perez
!!
real(dp) function v_asymptotic(k, r) result(v)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    real(dp), intent(in) :: r !< radius at which the wave will be evaluated

    v = exp(-k*r)
end function v_asymptotic

!!
!> @brief      Derivative of the asymptotic S wave
!! 
!! Derivative of the asymptotic S wave in the deuteron wave function
!!
!! @return     \f$ -k e^{- k r} \f$
!!
!! @author Rodrigo Navarro Perez
!!
real(dp) function vprime_asymptotic(k, r) result(vp)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    real(dp), intent(in) :: r !< radius at which the wave will be evaluated

    vp = -k*exp(-k*r)
end function vprime_asymptotic

!!
!> @brief      Asymptotic D wave
!! 
!! Asymptotic D wave in the deuteron wave function
!!
!! @return     \f$ e^{- k r}\left( 1 + \frac{3}{k r} + \frac{3}{(k r)^2} \right) \f$
!!
!! @author Rodrigo Navarro Perez
!!
real(dp) function w_asymptotic(k, r) result(w)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    real(dp), intent(in) :: r !< radius at which the wave will be evaluated

    w = exp(-k*r)*(1 + 3/(k*r) + 3/(k*r)**2)
end function w_asymptotic

!!
!> @brief      Derivative of the asymptotic D wave
!! 
!! Derivative of the asymptotic D wave in the deuteron wave function
!!
!! @return     \f$ -e^{- k r}\left( \frac{6 + 6 k r + 3 (k r)^2 + (k r)^3}{k^2 r^3} \right) \f$
!!
!! @author Rodrigo Navarro Perez
!!
real(dp) function wprime_asymptotic(k, r) result(wp)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    real(dp), intent(in) :: r !< radius at which the wave will be evaluated

    wp = -exp(-k*r)*(6 + 6*k*r + 3*(k*r)**2 + (k*r)**3)/(k**2*r**3)
end function wprime_asymptotic

!!
!> @brief      Deuteron Schrödinger in ODE form for Runge Kutta
!!
!! The Schrödinger equation for the deuteron channel 
!!
!! \f[ -v''(k, r) + U_{^3S_1}(r) v(k, r) + U_{\epsilon_1}(r) w(k,r) = - k^2 v(k, r) \f]
!! \f[ -w''(k, r) + [U_{^3D_1}(r) + 6/r^2 ]w(k, r) + U_{\epsilon_1}(r) v(k,r) = - k^2 w(k, r) \f]
!!
!! is recast as a set of four coupled ordinary differential equations 
!!
!! \f[v'(k, r) = \nu(k,r) \f]
!! \f[w'(k, r) = \omega(k,r) \f]
!! \f[\nu'(k, r) = [k^2 + U_{^3S_1}(r)] v(k, r) + U_{\epsilon_1}(r) w(k,r) \f]
!! \f[\omega'(k, r) = [k^2 + U_{^3D_1}(r) + 6/r^2 ]w(k, r) + U_{\epsilon_1}(r) v(k,r) \f]
!!
!! so that a Runge Kutta algorithm can integrate them.
!! 
!!
!! @return     Derivatives of the deuteron wave functions and their derivatives
!!
function deuteron_ode(wave_functions, r, work, model) result(f)
    implicit none
    real(dp), intent(in), dimension(:) :: wave_functions !< Array with the S and D waves and their derivatives
    real(dp), intent(in) :: r !< Radius at which the functions will be evaluated. In fm
    real(dp), intent(in) :: work(:) !< Array with all reals necessary to calculate derivatives, mainly the wave number and the potential parameters
    type(nn_model), intent(in) :: model !< NN potential 
    real(dp), allocatable :: f(:) !< First and second derivatives of the S and D waves

    real(dp) :: k, U_3S1, U_3D1, U_EP1, mu
    real(dp), allocatable, dimension(:) :: parameters
    real(dp) :: v_pw(n_waves, j_max)
    real(dp), allocatable, dimension(:, :, :) :: dv_pw

    allocate(f, mold = wave_functions)
    allocate(parameters(1 : size(work) - 1))
    k = work(1)
    parameters = work(2:)

    call model%potential(parameters, r, channel, v_pw, dv_pw)
    mu = 2*reduced_mass
    U_3S1 = v_pw(3,2)*mu/hbar_c**2
    U_EP1 = v_pw(4,2)*mu/hbar_c**2
    U_3D1 = v_pw(5,2)*mu/hbar_c**2

    f(1) = wave_functions(3)
    f(2) = wave_functions(4)
    f(3) = (k**2 + u_3s1)*wave_functions(1) + u_ep1*wave_functions(2)
    f(4) = (k**2 + u_3d1 + 6/r**2)*wave_functions(2) + u_ep1*wave_functions(1)

end function deuteron_ode

!!
!> @brief      Piecewise integration of the deuteron with a local potential
!!
!! Given an object of type nn_model (a local potential,
!! a maximum radius of integration and an integration step), determines
!! the deuteron wave function for an specific wave number \f$k\f$
!!
!! Constructing a delta shell representation of a local potential,
!! the deuteron wave function can be calculated using piecewise 
!! linear combinations of "free" wave solutions that solve the 
!! deuteron (bound) Schrödinger equation when the potential is zero.
!!
!! The piecewise solution in between the delta shells (located at \f$r_i\f$)
!! can be expressed as 
!!
!! \f[ v(r) = A_i \tilde{j}_0(k r) +  B_i \tilde{y}_0(k r) \f]
!! \f[ w(r) = C_i \tilde{j}_2(k r) +  D_i \tilde{y}_2(k r) \f]
!!
!! Plugging in these "free" solutions into the Schrödinger equation with
!! the delta shell potential results in a set of boundary conditions that 
!! must hold at every integration radius \f$r_i\f$ and allow to obtain
!! the \f$ABCD_{i-1} \f$ parameters in terms of the \f$ABCD_i \f$ ones
!! (see documentation in reverse_variable_wave)
!! 
!! The wave function must fulfill the boundary condition of \f$v(r)\f$ and \f$w(r)\f$
!! to vanish as \f$r \rightarrow \infty \f$. This is accomplished with \f$A_N = B_N\f$
!! and \f$ C_N = D_N \f$ so that the linear combinations result in decaying
!! exponentials.
!!
!! At \f$r=0\f$ the boundary condition that the wave function be \f$L^2\f$ integrable
!! is met with \f$ D_1 = 0 \f$
!!
!! The final wave function will be a linear combination of what we call
!! an \f$ \alpha \f$ solution, resulting from integrating backwards starting with
!! \f$A_N = B_N = 1, \quad C_N = D_N = 0 \f$, and a \f$ \beta \f$ solution, 
!! resulting from integrating backwards with \f$A_N = B_N = 0, \quad C_N = D_N = 1 \f$.
!!
!! The \f$\beta\f$ parameter in the final linear combination
!! \f[v(r) = v_\alpha(r) + \beta v_\beta(r) \\
!!     w(r) = w_\alpha(r) + \beta w_\beta(r) \f]
!! is determined by making sure that the \f$r=0\f$ boundary condition is met
!! (.i.e \f$D_{1, \alpha} + \beta D_{1, \beta} = 0\f$).
!!
!! While technically there should be an \f$ \alpha \f$ parameter in the final
!! linear combination, it can be arbitrarily set to unity since the wave function
!! can be normalized later.
!!
!! @returns    An object of type ds_wave_function 
!!
!! @author Rodrigo Navarro Perez
!!
type(ds_wave_function) function wave_function_ds(k, model, parameters) result(r)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    type(nn_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) :: parameters !< parameters of the nn interaction
    
    real(dp) ::  r_i
    real(dp), allocatable, dimension(:) :: radii
    real(dp), allocatable, dimension(:, :, :) :: v_pw
    real(dp), allocatable, dimension(:, :, :, :) :: dv_pw
    real(dp), allocatable, dimension(:) :: a_alpha, b_alpha, c_alpha, d_alpha
    real(dp), allocatable, dimension(:) :: a_beta, b_beta, c_beta, d_beta
    integer :: i
    real(dp) :: beta

    call all_delta_shells(model, parameters, channel, k_cm, j_max, radii, v_pw, dv_pw)

    r%n_radii = size(radii)
    r%gamma = k
    allocate(r%radii(1:r%n_radii))
    allocate(r%a_s(0:r%n_radii))
    allocate(r%b_s, r%a_d, r%b_d, a_alpha, b_alpha, c_alpha, d_alpha, a_beta, b_beta, c_beta, d_beta, mold=r%a_s)
    r%dr = model%dr
    a_alpha(r%n_radii) = -1._dp
    b_alpha(r%n_radii) = -1._dp
    c_alpha(r%n_radii) = 0._dp
    d_alpha(r%n_radii) = 0._dp

    a_beta(r%n_radii) = 0._dp
    b_beta(r%n_radii) = 0._dp
    c_beta(r%n_radii) = -1._dp
    d_beta(r%n_radii) = -1._dp

    do i = r%n_radii , 1, -1
        a_alpha(i-1) = a_alpha(i)
        b_alpha(i-1) = b_alpha(i)
        c_alpha(i-1) = c_alpha(i)
        d_alpha(i-1) = d_alpha(i)

        a_beta(i-1) = a_beta(i)
        b_beta(i-1) = b_beta(i)
        c_beta(i-1) = c_beta(i)
        d_beta(i-1) = d_beta(i)

        r_i = radii(i)
        call reverse_variable_wave(k, r_i, v_pw(3:5, 2, i), a_alpha(i-1), b_alpha(i-1), c_alpha(i-1), d_alpha(i-1))
        call reverse_variable_wave(k, r_i, v_pw(3:5, 2, i),  a_beta(i-1),  b_beta(i-1),  c_beta(i-1),  d_beta(i-1))
        r%radii(i) = r_i
    enddo

    beta = -d_alpha(0)/d_beta(0)
    r%a_s = (a_alpha + beta*a_beta)
    r%b_s = (b_alpha + beta*b_beta)
    r%a_d = (c_alpha + beta*c_beta)
    r%b_d = (d_alpha + beta*d_beta)
    
end function wave_function_ds

!!
!> @brief      Delta shell backward integration for the deuteron
!!
!! Given a set of delta shell strength parameters in the \f${^3}S_1\f$-\f${^3}D_1\f$
!! channel, makes a single step of a backwards integration in the linear combination 
!! parameters that characterize the wave function
!!
!! The deuteron wave function for a delta shell potential
!! can be solved piecewise. In the neighborhood of delta shell radius
!! \f$ r_i \f$ the S and D waves are given by 
!! 
!! \f[ v(r) = A_i \tilde{j}_0(k r) +  B_i \tilde{y}_0(k r) \f]
!! \f[ w(r) = C_i \tilde{j}_2(k r) +  D_i \tilde{y}_2(k r) \f]
!!
!! The backwards integration of the linear combination parameters is 
!! given by
!!
!! \f[ B_{i-1} =  B_i - \frac{\tilde{j}_0}{k} \left[ 
!!     \lambda_i^{^3S_1} (A_i \tilde{j}_0 + B_i \tilde{y}_0) + 
!!     \lambda_i^{\epsilon_1} (C_i \tilde{j}_2 + D_i \tilde{y}_2) \right] \\ 
!!     A_{i-1} = A_i - \frac{\tilde{y}_0}{\tilde{j}_0} (B_{i-1} - B_i ) \\
!!     D_{i-1} =  D_i - \frac{\tilde{j}_2}{k} \left[ 
!!     \lambda_i^{^3D_1} (C_i \tilde{j}_2 + D_i \tilde{y}_2) + 
!!     \lambda_i^{\epsilon_1} (A_i \tilde{j}_0 + B_i \tilde{y}_0) \right] \\
!!     C_{i-1} = C_i - \frac{\tilde{y}_0}{\tilde{j}_0} (D_{i-1} - D_i ) \f]
!!
!! Where the evaluation of the \f$ \tilde{j} \f$ and \f$ \tilde{y} \f$ functions has been 
!! obviated for simplicity. 
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine reverse_variable_wave(k, r, lambdas, a, b, c, d)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< radius of integration (where the delta shell is). In fm
    real(dp), intent(in), dimension(1:3) :: lambdas !< delta shell strength parameters. In fm\f$^{-1}\f$
    real(dp), intent(inout) :: a !< linear combination parameter
    real(dp), intent(inout) :: b !< linear combination parameter
    real(dp), intent(inout) :: c !< linear combination parameter
    real(dp), intent(inout) :: d !< linear combination parameter

    real(dp) :: j0, y0, j2, y2
    real(dp) :: l_3s1, l_ep1, l_3d1
    real(dp) :: diff_b, diff_d
    j0 = j0_bound(k*r)
    y0 = y0_bound(k*r)
    j2 = j2_bound(k*r)
    y2 = y2_bound(k*r)
    l_3s1 = lambdas(1)
    l_ep1 = lambdas(2)
    l_3d1 = lambdas(3)
    
    diff_b = -j0*(l_3s1*(a*j0 + b*y0) + l_ep1*(c*j2 + d*y2))/k
    diff_d = -j2*(l_3d1*(c*j2 + d*y2) + l_ep1*(a*j0 + b*y0))/k

    a = a - diff_b*y0/j0
    b = b + diff_b
    c = c - diff_d*y2/j2
    d = d + diff_d
end subroutine reverse_variable_wave

!!
!> @brief      regular "free" wave function for the deuteron
!!
!! Bound version of the regular "free" wave function in the S 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function j0_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluated
    r = sinh(x)
end function j0_bound

!!
!> @brief      irregular "free" wave function for the deuteron
!!
!! Bound version of the irregular "free" wave function in the S 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function y0_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluated
    r = -cosh(x)
end function y0_bound

!!
!> @brief      regular "free" wave function for the deuteron
!!
!! Bound version of the regular "free" wave function in the D 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function j2_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluated
    r = sinh(x) - 3*cosh(x)/x + 3*sinh(x)/(x**2)
end function j2_bound

!!
!> @brief      irregular "free" wave function for the deuteron
!!
!! Bound version of the irregular "free" wave function in the D 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function y2_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluated
    r = -cosh(x) + 3*sinh(x)/x -3*cosh(x)/(x**2)
end function y2_bound

!!
!> @brief      Derivatives of the deuteron for a local potential
!!
!! Given a wave number, a NN potential, and the potential parameters,
!! uses the Feynman-Hellman theorem to calculate the derivatives of 
!! the deuteron binding energy with respect of the potential parameters
!!
!! The partial derivatives of the binding energy with respect of the potential
!! parameters are calculated using the Feynman-Hellman theorem (see 
!! Phys. Rev. 56, 340 ) where if \f$ \hat{H} | \psi \rangle = -k^2 | \psi \rangle \f$,
!! then the partial derivative of \f$-k^2\f$ with respect of one of the potential
!! parameters \f$p_j\f$ is given by
!! \f[ \frac{\partial ( -k^2)}{\partial p_j} =  \left\langle \psi \left| \frac{\partial \hat{H}}{\partial p_j }
!!   \right| \psi \right\rangle. \f]
!!
!! In the case of the \f$^3S_1\f$-\f$^3D_1\f$ coupled channel, the Hamiltonian \f$\hat{H}\f$
!! is a \f$2 \times 2\f$ matrix with the \f$V_{\epsilon_1}\f$ potential in the off 
!! diagonal elements, and the \f$ | \psi \rangle \f$ ket is a column vector with the 
!! normalized S and D wave functions. 
!!
!!
!! @return     Array with derivatives with respect of potential parameters
!!
!! @author     Rodrigo Navarro Perez
!!
function local_deuteron_derivatives(k, model, parameters) result(dbe)
    implicit none
    real(dp), intent(in) :: k !< Deuteron wave number. In fm\f$^{-1}\f$
    type(nn_model), intent(in) :: model !< NN potential
    real(dp), intent(in), dimension(:) :: parameters !< potential parameters
    real(dp), allocatable, dimension(:) :: dbe !< Derivatives of the deuteron binding energy with respect of the potential parameters

    type(local_wave_function) :: wavefunc
    integer :: i
    real(dp) :: v_pw(5,2), mu, r, delta_r
    real(dp), allocatable, dimension(:, :, :) :: dv_pw
    real(dp), allocatable, dimension(:, :) :: du_3s1, du_ep1, du_3d1
    real(dp), allocatable, dimension(:) :: integrand

    wavefunc = wave_function_local(k, model, parameters)
    call normalize_local(wavefunc)

    allocate(dbe, mold=parameters)
    allocate(du_3s1(1:size(parameters), 1:size(wavefunc%radii)))
    allocate(du_ep1, du_3d1, mold=du_3s1)
    dbe = 0._dp

    ! Evaluating derivatives of the 3S1-3D1 potential as a function of r
    mu = 2*reduced_mass
    do i=1, size(wavefunc%radii)
        r = wavefunc%radii(i)
        call model%potential(parameters, r, channel, v_pw, dv_pw)
        du_3s1(:, i) = dv_pw(:, 3,2)*mu/hbar_c**2
        du_ep1(:, i) = dv_pw(:, 4,2)*mu/hbar_c**2
        du_3d1(:, i) = dv_pw(:, 5,2)*mu/hbar_c**2
    enddo
    
    ! Feynman-Hellman theorem to get the derivative of -k**2 with
    ! respect of the potential parameters
    allocate(integrand, mold=wavefunc%radii)
    delta_r = abs(wavefunc%radii(2) - wavefunc%radii(1))
    do i=1, size(parameters)
        integrand = du_3s1(i, :)*(wavefunc%s_wave)**2 + 2*du_ep1(i, :)*wavefunc%s_wave*wavefunc%d_wave &
            + du_3d1(i, :)*(wavefunc%d_wave)**2
        dbe(i) = booles_quadrature(integrand, delta_r)
    enddo

    ! Going from derivative of -k**2 to derivative of the binding energy
    call scale_derivatives(model%relativistic_deuteron, k, dbe)

end function local_deuteron_derivatives

!!
!> @brief      Normalize a deuteron wave function
!!
!! Given an already evaluated (from 0 to \f$r_{\rm max}\f$) solution of
!! the deuteron wave function, uses Boole's quadrature to calculate 
!! the integral of the probability function in the evaluated range.
!!
!! The integral from 0 to \f$r_{\rm max}\f$ to infinity
!! is calculated analytically based on the exponential decay
!! of the wave function that has been imposed by construction
!! by integrating backwards (see documentation of wave_function_local
!! for details)
!!
!! The full integral is then used to calculate the normalization
!! constant and to normalize the wave function
!!
!! The analytic expressions for the integrals to infinity are
!! given by
!!
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine normalize_local(wavefunc)
    implicit none
    type(local_wave_function), intent(inout) :: wavefunc !< evaluated deuteron wave function

    real(dp), allocatable, dimension(:) :: probability_density
    real(dp) :: delta_r, s, a_s, a_d, r_max, k, a_norm
    integer :: n_points

    n_points = size(wavefunc%radii)

    allocate(probability_density, mold=wavefunc%radii)
    probability_density = (wavefunc%s_wave)**2 + (wavefunc%d_wave)**2

    delta_r = abs(wavefunc%radii(2) - wavefunc%radii(1))
    s = booles_quadrature(probability_density, delta_r)
    r_max = wavefunc%radii(n_points)
    k = wavefunc%gamma
    a_s = wavefunc%s_wave(n_points)/v_asymptotic(k, r_max)
    a_d = wavefunc%d_wave(n_points)/w_asymptotic(k, r_max)
    s = s + asymptotic_integral(a_s, a_d, k, r_max)
    a_norm = sqrt(1._dp/s)
    wavefunc%s_wave = a_norm*wavefunc%s_wave
    wavefunc%d_wave = a_norm*wavefunc%d_wave
end subroutine normalize_local

!!
!> @brief     Integral of the asymptotic deuteron wave function
!!
!! Given the S and D wave amplitudes in the asymptotic region, the wave number,
!! and the radius where the asymptotic function starts, uses an analytic expression
!! to calculate the integral of probability density from the given radius to infinity
!!
!! \f[ \int_{r_N}^\infty  A_S^2 e^{-2 k r} dr = A_S^2 \frac{e^{-2 k r_N}}{2 k}  \\
!!     \int_{r_N}^\infty  A_D^2 e^{-2 k r} \left( 1 + \frac{3}{k r} + \frac{3}{(k r)^2} \right)^2 dr =
!!     A_D^2 e^{-2 k r_N} \frac{ 6 + 12 k r_N + 6 (k r_N)^2 + (k r_N)^3}{2 k^4 r_N^3}
!! \f]
!!
!! @return     Integral of probability density from the given radius to infinity
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function asymptotic_integral(a_s, a_d, k, r) result(s)
    implicit none
    real(dp), intent(in) :: a_s !< S wave asymptotic amplitude
    real(dp), intent(in) :: a_d !< D wave asymptotic amplitude
    real(dp), intent(in) :: k !< Wave number. In fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< Lower limit of the asymptotic integral. In fm
    s = exp(-2*k*r)*(a_s**2/(2*k) + a_d**2*(6 + k*r*(12 + k*r*(6 + k*r)))/(2*k**4*r**3))
end function asymptotic_integral


!!
!> @brief      Convert derivatives of \f$-k^2\f$ to derivatives of \f$B(k)\f$
!!
!! The formula used depends on whether or not relativistic kinematics where used to 
!! calculate the binding energy
!!
!! @return     Derivatives of the deuteron binding energy
!!
subroutine scale_derivatives(is_relativistic, k, dbe)
    implicit none
    logical, intent(in) :: is_relativistic !< Are relativistic kinematics being used?
    real(dp), intent(in) :: k !< Wave number. In fm\f$^{-1}\f$
    real(dp), intent(inout), dimension(:) :: dbe !< In input, derivatives of \f$-k^2\f$. On output, derivatives of \f$B(k)\f$

    if (is_relativistic) then
        dbe = -0.5_dp*(1/sqrt(m_p**2 - (k*hbar_c)**2) + 1/sqrt(m_p**2 - (k*hbar_c)**2))*hbar_c**2*dbe
    else
        dbe = -(hbar_c)**2/(2*reduced_mass)*dbe
    endif
    
end subroutine scale_derivatives

!!
!> @brief      Derivatives of the deuteron for a Delta-Shell potential
!!
!! Given a wave number, a NN potential, and the potential parameters,
!! uses the Feynman-Hellman theorem to calculate the derivatives of 
!! the deuteron binding energy with respect of the Delta-Shell potential parameters
!!
!! The partial derivatives of the binding energy with respect of the potential
!! parameters are calculated using the Feynman-Hellman theorem (see 
!! Phys. Rev. 56, 340 ) where if \f$ \hat{H} | \psi \rangle = -k^2 | \psi \rangle \f$,
!! then the partial derivative of \f$-k^2\f$ with respect of one of the potential
!! parameters \f$p_j\f$ is given by
!! \f[ \frac{\partial ( -k^2)}{\partial p_j} =  \left\langle \psi \left| \frac{\partial \hat{H}}{\partial p_j }
!!   \right| \psi \right\rangle. \f]
!!
!! In the case of the \f$^3S_1\f$-\f$^3D_1\f$ coupled channel, the Hamiltonian \f$\hat{H}\f$
!! is a \f$2 \times 2\f$ matrix with the \f$V_{\epsilon_1}\f$ potential in the off 
!! diagonal elements, and the \f$ | \psi \rangle \f$ ket is a column vector with the 
!! normalized S and D wave functions. 
!!
!!
!! @return     Array with derivatives with respect of potential parameters
!!
!! @author     Rodrigo Navarro Perez
!!
function ds_deuteron_derivatives(k, model, parameters) result(dbe)
    implicit none
    real(dp), intent(in) :: k !< Wave number. In fm\f$^{-1}\f$
    type(nn_model), intent(in) :: model !< Delta Shell NN potential
    real(dp), intent(in), dimension(:) :: parameters !< Parameters of the DS potential
    real(dp), allocatable, dimension(:) :: dbe !< Derivatives of the binding energy with respect of the potential parameters

    real(dp) :: r, v, w
    type(ds_wave_function) :: wavefunc
    integer :: i
    real(dp), allocatable, dimension(:) :: radii
    real(dp), allocatable, dimension(:, :, :) :: v_pw
    real(dp), allocatable, dimension(:, :, :, :) :: dv_pw

    wavefunc = wave_function_ds(k, model, parameters)
    call normalize_ds(wavefunc)

    allocate(dbe, mold=parameters)
    dbe = 0._dp
    call all_delta_shells(model, parameters, channel, k_cm, j_max, radii, v_pw, dv_pw)
    ! Feynman-Hellman theorem to get the derivative of -k**2 with
    ! respect of the potential parameters
    do i = 1, wavefunc%n_radii
        r = wavefunc%radii(i)
        v = s_wave_function(wavefunc, i)
        w = d_wave_function(wavefunc, i)
        dbe = dbe + dv_pw(:, 3, 2, i)*v**2 + 2*dv_pw(:, 4, 2, i)*v*w + dv_pw(:, 5, 2, i)*w**2
    enddo

    ! Going from derivative of -k**2 to derivative of the binding energy
    call scale_derivatives(model%relativistic_deuteron, k, dbe)

    
end function ds_deuteron_derivatives

!!
!> @brief      Normalize a piecewise wave function from a delta shell potential
!!
!! Given a piecewise solution of the deuteron wave function, 
!! uses piecewise definite integrals to calculate the integral
!! of the probability density between each pair of consecutive
!! integration radii.
!!
!! The integral from the last integration radius to infinity
!! is calculated analytically based on the exponential decay
!! of the wave function that has been imposed by construction
!! by integrating backwards (see documentation of wave_function
!! for details)
!!
!! The analytic expressions for the integrals to infinity are
!! given by
!!
!! \f[ \int_{r_N}^\infty  A_S^2 \left[ \tilde{j}_0(kr) + \tilde{y}_0(kr) \right]^2 dr = 
!!      A_S^2 \frac{e^{-2 k r_N}}{2 k}  \\
!!     \int_{r_N}^\infty  A_D^2 \left[ \tilde{j}_2(kr) + \tilde{y}_2(kr) \right]^2 dr =
!!     A_D^2 e^{-2 k r_N} \frac{ 6 + 12 k r_N + 6 (k r_N)^2 + (k r_N)^3}{2 k^4 r_N^3}
!! \f]
!!
!! @author     Rodrigo Navarro Perez
subroutine normalize_ds(wavefunc)
    implicit none
    type(ds_wave_function), intent(inout) :: wavefunc !< piecewise deuteron wave function

    real(dp) :: s, a_s, a_d, r, k, a_norm
    integer :: i, n_points

    n_points = wavefunc%n_radii
    s = 0._dp
    do i = 1, n_points
        s = s + definite_density_integral(wavefunc, i)
    enddo
    a_s = wavefunc%a_s(n_points)
    a_d = wavefunc%a_d(n_points)
    r = wavefunc%radii(n_points)
    k = wavefunc%gamma
    ! Integral from last grid point to infinity where the wave function has exponential decay
    s = s + asymptotic_integral(a_s, a_d, k, r)
    a_norm = sqrt(1._dp/s)
    wavefunc%a_s = a_norm*wavefunc%a_s
    wavefunc%b_s = a_norm*wavefunc%b_s
    wavefunc%a_d = a_norm*wavefunc%a_d
    wavefunc%b_d = a_norm*wavefunc%b_d
end subroutine normalize_ds

!!
!> @brief      Analytic definite integral of the deuteron probability density
!!
!! Given a piecewise solution to the deuteron wave function and the position
!! index of a concentration radius, calculates the definite integral of the
!! deuteron probability density between the previous radius and the current one
!!
!! @return     Analytic integral of deuteron probability density between 2
!!             consecutive interaction radii
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function definite_density_integral(wavefunc, i) result(s)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc !< Piecewise deuteron solution
    integer, intent(in) :: i !< concentration radius index

    real(dp) :: k, a_s, b_s, a_d, b_d, r_i, r_f

    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in definite_density_integral is incompatible with wavefunc grid'
    endif
    k = wavefunc%gamma
    r_f = wavefunc%radii(i)
    a_s = wavefunc%a_s(i - 1)
    b_s = wavefunc%b_s(i - 1)
    a_d = wavefunc%a_d(i - 1)
    b_d = wavefunc%b_d(i - 1)
    s = s_density_integral(a_s, b_s, k, r_f) + d_density_integral(a_d, b_d, k, r_f)
    if (i == 1) then
        s = s + a_s*b_s/(2*k)
    else
        r_i = wavefunc%radii(i - 1)
        s = s - s_density_integral(a_s, b_s, k, r_i) - d_density_integral(a_d, b_d, k, r_i)
    endif
    
end function definite_density_integral

!!
!> @brief      Analytic integral of S wave probability density
!!
!! Given the linear combination parameters that determine the
!! deuteron S wave function in a certain interval, returns the 
!! integral of the corresponding probability density evaluated
!! at the given radius
!!
!! The analytic integral is given by
!!
!! \f[  \int \left[A \tilde{j}_0(kr) + B \tilde{y}_0(kr) \right]^2 dr = 
!!        \frac{2(B^2 - A^2)kr - 2AB\cosh(2kr) +(B^2+A^2)\sinh(2kr)}{4k} \f]
!!
!!
!! @return     Analytic integral of S wave probability density
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function s_density_integral(a, b, k, r) result(s)
    implicit none
    real(dp), intent(in) :: a !< regular component in the S channel
    real(dp), intent(in) :: b !< irregular component in the S channel
    real(dp), intent(in) :: k !< wavenumber in fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< radius in fm
    s = (2*(b**2 - a**2)*k*r - 2*a*b*cosh(2*k*r) + (b**2 + a**2)*sinh(2*k*r))/(4*k)
end function s_density_integral

!!
!> @brief      Analytic integral of D wave probability density
!!
!! Given the linear combination parameters that determine the
!! deuteron D wave function in a certain interval, returns the 
!! integral of the corresponding probability density evaluated
!! at the given radius
!!
!! The analytic integral is given by
!!
!! \f[  \int \left[A \tilde{j}_2(kr) + B \tilde{y}_2(kr) \right]^2 dr = 
!!        \frac{2(A^2 - B^2)(3 - 3 k^2 r^2 - k^4 r^4) 
!!              -2(3(A^2+B^2)(1 + k^2r^2) + ABkr(12 + k^2r^2))\cosh(2kr) 
!!              +(kr(A^2+B^2)(12 + k^2r^2) + 12 AB(1 + k^2r^2))\sinh(2kr)}{4k^4r^3} \f]
!!
!!
!! @return     Analytic integral of D wave probability density
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function d_density_integral(a, b, k, r) result(s)
    implicit none
    real(dp), intent(in) :: a !< regular component in the D channel
    real(dp), intent(in) :: b !< irregular component in the D channel
    real(dp), intent(in) :: k !< wavenumber in fm\f$^{-1}\f$
    real(dp), intent(in) :: r !< radius in fm

    s =  2*(a**2 - b**2)*(3 - 3*(k*r)**2 - (k*r)**4) &
        -2*(3*(a**2 + b**2)*(1 + (k*r)**2) + a*b*k*r*(12 + (k*r)**2))*cosh(2*k*r) &
        +(k*r*(a**2 + b**2)*(12 + (k*r)**2) + 12*a*b*(1 + (k*r)**2) )*sinh(2*k*r)
    s = s/(4*k**4*r**3)    
end function d_density_integral

!!
!> @brief      S wave at a particular integration radius
!!
!! Given a piecewise solution of the deuteron, evaluates the 
!! S wave at a particular integration radius characterized 
!! by the index i
!!
!! @return     S wave at a particular integration radius
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function s_wave_function(wavefunc, i) result(v)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc !< Piecewise deuteron solution
    integer, intent(in) :: i !< integration radius index

    real(dp) :: k ,r
    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in s_wave_function is incompatible with wavefunc grid'
    endif
    k = wavefunc%gamma
    r = wavefunc%radii(i)
    v = wavefunc%a_s(i)*j0_bound(k*r) + wavefunc%b_s(i)*y0_bound(k*r)
end function s_wave_function

!!
!> @brief      D wave at a particular integration radius
!!
!! Given a piecewise solution of the deuteron, evaluates the 
!! D wave at a particular integration radius characterized 
!! by the index i
!!
!! @return     D wave at a particular integration radius
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function d_wave_function(wavefunc, i) result(w)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc !< Piecewise deuteron solution
    integer, intent(in) :: i !< integration radius index

    real(dp) :: k ,r
    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in d_wave_function is incompatible with wavefunc grid'
    endif
    k = wavefunc%gamma
    r = wavefunc%radii(i)
    w = wavefunc%a_d(i)*j2_bound(k*r) + wavefunc%b_d(i)*y2_bound(k*r)
end function d_wave_function


end module deuteron