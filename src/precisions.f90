!!
!> @brief precisions to declare variables
!! 
!! Defines the integer parameters to be used as kinds to define real variables as single precision 
!! (32-bit), double precision (64-bit) and quadruple precision (128-bit) (depending on the machine 
!! and compiler, this last one may not always be available).
!!
!! The intrinsic module 'iso_fortran_env' (Fortran 2008 and later) is
!! used.
!!
!! More precisions may be defined later (i.e. larger integers)
!!
!! @author Rodrigo Navarro Perez
!!
module precisions

use iso_fortran_env

implicit none

integer, parameter :: sp = REAL32 !< single precision kind
integer, parameter :: dp = REAL64 !< double precision kind
integer, parameter :: qp = REAL128!< quadruple precision kind

end module precisions
