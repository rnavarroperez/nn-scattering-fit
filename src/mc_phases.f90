!!
!> @brief      Writes phases from a MC sample of potential parameters.
!!
!! Reads a set of MC potential paramters (determined from a legacy code)
!! and writes the corresponding np phase shifts from 1 to 350 MeV
!!
!! @author     Rodrigo Navarro Perez
!!
program write_mc_phases
use precisions, only : dp
use delta_shell, only : set_ds_potential, nn_model
use read_write, only : read_montecarlo_parameters, write_montecarlo_phases

implicit none

type(nn_model) :: ds_model
real(dp), allocatable, dimension(:, :) :: mc_lambas
real(dp), allocatable, dimension(:) :: parameters
integer :: itlab

call set_ds_potential('ds_ope30_fff', ds_model, parameters)
call read_montecarlo_parameters(abs(parameters) > 0, 'mc_parameters/mc_lambdas_fixedZ_0001_1000ope30fff.dat', mc_lambas)
do itlab = 1, 350
    call write_montecarlo_phases(ds_model, mc_lambas, 'np', itlab)    
enddo
    
end program write_mc_phases