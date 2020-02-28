!!
!> @brief fits a NN potential to scattering data
!!
!! Uses the Levenberg Marquardt algorithm to adjust the parameters of a NN interaction and
!! reproduce experimental data collected since 1954.
!!
!! @author Rodrigo Navarro Perez
!!
program nn_fit

use read_write, only: print_em_np_amplitudes

implicit none

call print_em_np_amplitudes()


end program nn_fit
