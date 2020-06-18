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
use nn_phaseshifts, only : nn_local_model
use constants, only : hbar_c, m_p=>proton_mass, m_n=>neutron_mass, pi, alpha
implicit none

private

public :: binding_energy

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

contains

!!
!> @brief      Binding energy of the deuteron
!! Given a nn local model (local potential, maximum integration radius,
!! integration step, and potential parameters), calculates the deuteron
!! binding energy and its derivatives with respect to the potential parameters
!!
!! Once the wave number \f$k\f$ is determined, the binding energy can be calculated
!! using relativistic kinematics (see Phys. Rev. C 48 (1993), 792-815)
!! \f[ B = M_p + M_n - \sqrt{M_p^2 - k^2} - \sqrt{M_n^2 - k^2}  \f]
!! where \f$ k \f$ has to be in units of MeV.
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
    type(nn_local_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) ::  parameters !< parameters of the nn interaction
    real(dp), intent(out) :: be !< deuteron binding energy. In MeV
    real(dp), intent(out), allocatable, dimension(:) :: dbe !< derivatives of the deuteron binding energy with respect of the potential parameters

    real(dp) :: k, r, v, w
    type(ds_wave_function) :: wavefunc
    integer :: i
    real(dp), parameter :: mu = 2*m_p*m_n/(m_p + m_n)
    real(dp), dimension(1:5, 1:2) :: v_pw
    real(dp), allocatable, dimension(:, :, :) :: dv_pw

    k = wave_number(model, parameters)
    be = m_p + m_n - sqrt(m_p**2 - (k*hbar_c)**2) - sqrt(m_n**2 - (k*hbar_c)**2)
    wavefunc = wave_function(k, model, parameters)
    call normalize(wavefunc)
    allocate(dbe, mold=parameters)
    dbe = 0._dp
    ! Feynman-Hellman theorem to get the derivative of -k**2 with
    ! respect of the potential parameters
    do i = 1, wavefunc%n_radii
        r = wavefunc%radii(i)
        call model%potential(parameters, r, 'np', v_pw, dv_pw)
        dv_pw = dv_pw*mu*model%dr/(hbar_c**2)
        v = s_wave_function(wavefunc, i)
        w = d_wave_function(wavefunc, i)
        dbe = dbe + dv_pw(:, 3, 2)*v**2 + 2*dv_pw(:, 4, 2)*v*w + dv_pw(:, 5, 2)*w**2
    enddo
    ! Going from derivative of -k**2 to derivative of the binding energy
    dbe = -0.5_dp*(1/sqrt(m_p**2 - (k*hbar_c)**2) + 1/sqrt(m_p**2 - (k*hbar_c)**2))*hbar_c**2*dbe
end subroutine binding_energy

!!
!> @brief      Normalize a piecewise wave function
!!
!! Given a piecewise solution of the deuteron wave function, 
!! uses a simple trapezoid rue to integrate the probability
!! density, evaluating it at the integration radii, to get 
!! a normalization constant and normalize the wave function
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
!!
subroutine normalize(wavefunc)
    implicit none
    type(ds_wave_function), intent(inout) :: wavefunc !< piecewise deuteron wave function

    real(dp) :: r, k, s
    real(dp) :: a_s, a_d, a_norm
    integer :: i, n_points
    n_points = wavefunc%n_radii
    ! one half from the trapezoid rule plus one quarter from the triangle from r = 0 to r = r_1 = 0.5*dr
    s = 0.75_dp*probability_density(wavefunc, 1)
    ! trapezoid rule for all intermediate points
    do i = 2, n_points - 1
        s = s + probability_density(wavefunc, i)
    enddo
    ! one half of the last point from trapezoid rule
    s = s + 0.5_dp*probability_density(wavefunc, n_points)
    ! overall dr factor from the trapezoid rule
    s = s*wavefunc%dr
    a_s = wavefunc%a_s(n_points)
    a_d = wavefunc%a_d(n_points)
    r = wavefunc%radii(n_points)
    k = wavefunc%gamma
    ! Integral from last grid point to infinity where the wave function has exponential decay
    s = s + exp(-2*k*r)*(a_s**2/(2*k) + a_d**2*(6 + k*r*(12 + k*r*(6 + k*r)))/(2*k**4*r**3))
    a_norm = sqrt(1._dp/s)
    wavefunc%a_s = a_norm*wavefunc%a_s
    wavefunc%b_s = a_norm*wavefunc%b_s
    wavefunc%a_d = a_norm*wavefunc%a_d
    wavefunc%b_d = a_norm*wavefunc%b_d
end subroutine normalize

!!
!> @brief      Probability density at a particular integration radius
!!
!! Given a piecewise solution of the deuteron, evaluates the 
!! probability density \f$ v(r)^2 + w(r)^2 \f$ at a particular integration
!! radius characterized by the index i
!!
!! @return     Probability density at a particular integration radius
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function probability_density(wavefunc, i) result(p)
    implicit none
    type(ds_wave_function), intent(in) :: wavefunc !< Piecewise deuteron solution
    integer, intent(in) :: i !< integration radius index

    if (i <= 0 .or. i > wavefunc%n_radii) then
        print*, i, wavefunc%n_radii
        stop 'index i in probability_density is incompatible with wavefunc grid'
    endif
    p = s_wave_function(wavefunc, i)**2 + d_wave_function(wavefunc, i)**2
end function probability_density

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
        stop 'index i in probability_density is incompatible with wavefunc grid'
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
        stop 'index i in probability_density is incompatible with wavefunc grid'
    endif
    k = wavefunc%gamma
    r = wavefunc%radii(i)
    w = wavefunc%a_d(i)*j2_bound(k*r) + wavefunc%b_d(i)*y2_bound(k*r)
end function d_wave_function

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
    type(nn_local_model), intent(in) :: model !< local model for nn the interaction
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
!! After integrating the piecewise deuteron wave function
!! of the type
!! \f[ v_i(r) = A_{S,i} \tilde{j}_0(k r) +  B_{S,i} \tilde{y}_0(k r) \f]
!! \f[ w_i(r) = A_{D,i} \tilde{j}_2(k r) +  B_{D,i} \tilde{y}_2(k r) \f]
!!
!! Returns the irregular component of the S wave solution \f$ B_{S,i} \f$
!! when \f$r \rightarrow 0\f$
!!
!! @return     Irregular component of the S wave when \f$r \rightarrow 0\f$
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function irregularity(k, model, parameters) result(r)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    type(nn_local_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) :: parameters !< parameters of the nn interaction


    type(ds_wave_function) :: deuteron_wf

    deuteron_wf = wave_function(k, model, parameters)
    r = deuteron_wf%b_s(0)
end function irregularity

!!
!> @brief      Piecewise integration of the deuteron with a local potential
!!
!! Given an object of type nn_local_model (a local potential,
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
!!     v(r) = v_\alpha(r) + \beta v_\beta(r) \f]
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
type(ds_wave_function) function wave_function(k, model, parameters) result(r)
    implicit none
    real(dp), intent(in) :: k !< deuteron wave number. In fm\f$^{-1}\f$ 
    type(nn_local_model), intent(in) :: model !< local model for nn the interaction
    real(dp), intent(in), dimension(:) :: parameters !< parameters of the nn interaction
    
    real(dp) ::  r_i
    real(dp), parameter :: mu = 2*m_p*m_n/(m_p + m_n)
    real(dp), dimension(1:5, 1:2) :: v_pw
    real(dp), allocatable, dimension(:, :, :) :: dv_pw
    real(dp), allocatable, dimension(:) :: a_alpha, b_alpha, c_alpha, d_alpha
    real(dp), allocatable, dimension(:) :: a_beta, b_beta, c_beta, d_beta
    integer :: i
    real(dp) :: beta

    r%n_radii = int(model%r_max/model%dr)
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

        r_i = model%dr*(i - 0.5_dp)
        call model%potential(parameters, r_i, 'np', v_pw, dv_pw)
        v_pw = v_pw*mu*model%dr/(hbar_c**2)
        call reverse_variable_wave(k, r_i, v_pw(3:5, 2), a_alpha(i-1), b_alpha(i-1), c_alpha(i-1), d_alpha(i-1))
        call reverse_variable_wave(k, r_i, v_pw(3:5, 2),  a_beta(i-1),  b_beta(i-1),  c_beta(i-1),  d_beta(i-1))
        r%radii(i) = r_i
    enddo

    beta = -d_alpha(0)/d_beta(0)
    r%a_s = (a_alpha + beta*a_beta)
    r%b_s = (b_alpha + beta*b_beta)
    r%a_d = (c_alpha + beta*c_beta)
    r%b_d = (d_alpha + beta*d_beta)
    
end function wave_function

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
!> @brief      regular "free" wave funcntion for the deuteron
!!
!! Bound version of the regular "free" wave function in the S 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function j0_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluted
    r = sinh(x)
end function j0_bound

!!
!> @brief      irregular "free" wave funcntion for the deuteron
!!
!! Bound version of the irregular "free" wave function in the S 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function y0_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluted
    r = -cosh(x)
end function y0_bound

!!
!> @brief      regular "free" wave funcntion for the deuteron
!!
!! Bound version of the regular "free" wave function in the D 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function j2_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluted
    r = sinh(x) - 3*cosh(x)/x + 3*sinh(x)/(x**2)
end function j2_bound

!!
!> @brief      irregular "free" wave funcntion for the deuteron
!!
!! Bound version of the irregular "free" wave function in the D 
!! channel
!!
!! @author     Rodrigo Navarro Perez
!!
real(dp) function y2_bound(x) result(r)
    implicit none
    real(dp), intent(in) :: x !< point at which the function is evaluted
    r = -cosh(x) + 3*sinh(x)/x -3*cosh(x)/(x**2)
end function y2_bound

end module deuteron