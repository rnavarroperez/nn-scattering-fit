module av18_compatibility

use precisions, only : dp
implicit none

private

public :: read_marias_format, write_marias_format

contains

subroutine write_marias_format(filename, parameters)
    implicit none
    character(len=*), intent(in) :: filename
    real(dp), intent(in), dimension(:)  :: parameters

    integer :: unit

    open(newunit=unit, file=trim(filename))
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.075_dp, 'cfsq0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.075_dp, 'cfsqc'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(55), 'cc'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(56), 'cr'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(57), 'ca'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 1), 'ccI01pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 2), 'ccP01pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 3), 'ccR01pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 1) + parameters(43), 'ccI01np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 2) + parameters(44), 'ccP01np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 3) + parameters(45), 'ccR01np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 1) + parameters(46), 'ccI01nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 2) + parameters(47), 'ccP01nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 3) + parameters(48), 'ccR01nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 4), 'cl2I01'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 5), 'cl2P01'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 6), 'cl2R01'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 7), 'ccI00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 8), 'ccP00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 9), 'ccR00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(10), 'cl2I00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(11), 'cl2P00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(12), 'cl2R00'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(13), 'ccI11pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(14), 'ccP11pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(15), 'ccR11pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(13) + parameters(46)/2._dp, 'ccI11np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(14) + parameters(47)/2._dp, 'ccP11np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(15) + parameters(48)/2._dp, 'ccR11np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(13) + parameters(46), 'ccI11nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(14) + parameters(47), 'ccP11nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(15) + parameters(48), 'ccR11nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(16), 'cl2I11'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(17), 'cl2P11'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(18), 'cl2R11'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(19), 'ctI1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(20), 'ctQ1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(21), 'ctR1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(19) + parameters(49), 'ctI1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(20) + parameters(50), 'ctQ1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(21) + parameters(51), 'ctR1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(19), 'ctI1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(20), 'ctQ1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(21), 'ctR1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(22), 'clsI1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(23), 'clsP1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(24), 'clsR1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(22) + parameters(52), 'clsI1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(23) + parameters(53), 'clsP1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(24) + parameters(54), 'clsR1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(22), 'clsI1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(23), 'clsP1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(24), 'clsR1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(25), 'cls2I1'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(26), 'cls2P1'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(27), 'cls2R1'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(28), 'ccI10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(29), 'ccP10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(30), 'ccR10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(31), 'cl2I10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(32), 'cl2P10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(33), 'cl2R10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(34), 'ctI0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(35), 'ctQ0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(36), 'ctR0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(37), 'clsI0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(38), 'clsP0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(39), 'clsR0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(40), 'cls2I0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(41), 'cls2P0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(42), 'cls2R0'
    close(unit)
    
end subroutine write_marias_format

subroutine read_marias_format(filename, parameters)
    implicit none
    character(len=*), intent(in) :: filename
    real(dp), intent(out), dimension(:)  :: parameters

    logical :: file_exists, drop
    integer :: unit, i
    real(dp) :: ccI01np, ccP01np, ccR01np, ccI01nn, ccP01nn, ccR01nn, ccI11np, ccP11np, ccR11np, ccI11nn, &
        ccP11nn, ccR11nn, ctI1np, ctQ1np, ctR1np, ctI1nn, ctQ1nn, ctR1nn, clsI1np, clsP1np, clsR1np, clsI1nn, &
        clsP1nn, clsR1nn

    inquire(file=trim(filename), exist = file_exists)
    if (file_exists) then
        open(newunit= unit, file=trim(filename))
        do i=1,2
            read(unit, *)
        enddo
        read(unit, *) drop, parameters(55)
        read(unit, *) drop, parameters(56)
        read(unit, *) drop, parameters(57)
        read(unit, *) drop, parameters( 1)
        read(unit, *) drop, parameters( 2)
        read(unit, *) drop, parameters( 3)
        read(unit, *) drop, ccI01np
        read(unit, *) drop, ccP01np
        read(unit, *) drop, ccR01np
        read(unit, *) drop, ccI01nn
        read(unit, *) drop, ccP01nn
        read(unit, *) drop, ccR01nn
        read(unit, *) drop, parameters( 4)
        read(unit, *) drop, parameters( 5)
        read(unit, *) drop, parameters( 6)
        read(unit, *) drop, parameters( 7)
        read(unit, *) drop, parameters( 8)
        read(unit, *) drop, parameters( 9)
        read(unit, *) drop, parameters(10)
        read(unit, *) drop, parameters(11)
        read(unit, *) drop, parameters(12)
        read(unit, *) drop, parameters(13)
        read(unit, *) drop, parameters(14)
        read(unit, *) drop, parameters(15)
        read(unit, *) drop, ccI11np
        read(unit, *) drop, ccP11np
        read(unit, *) drop, ccR11np
        read(unit, *) drop, ccI11nn
        read(unit, *) drop, ccP11nn
        read(unit, *) drop, ccR11nn
        read(unit, *) drop, parameters(16)
        read(unit, *) drop, parameters(17)
        read(unit, *) drop, parameters(18)
        read(unit, *) drop, parameters(19)
        read(unit, *) drop, parameters(20)
        read(unit, *) drop, parameters(21)
        read(unit, *) drop, ctI1np
        read(unit, *) drop, ctQ1np
        read(unit, *) drop, ctR1np
        read(unit, *) drop, ctI1nn
        read(unit, *) drop, ctQ1nn
        read(unit, *) drop, ctR1nn
        read(unit, *) drop, parameters(22)
        read(unit, *) drop, parameters(23)
        read(unit, *) drop, parameters(24)
        read(unit, *) drop, clsI1np
        read(unit, *) drop, clsP1np
        read(unit, *) drop, clsR1np
        read(unit, *) drop, clsI1nn
        read(unit, *) drop, clsP1nn
        read(unit, *) drop, clsR1nn
        read(unit, *) drop, parameters(25)
        read(unit, *) drop, parameters(26)
        read(unit, *) drop, parameters(27)
        read(unit, *) drop, parameters(28)
        read(unit, *) drop, parameters(29)
        read(unit, *) drop, parameters(30)
        read(unit, *) drop, parameters(31)
        read(unit, *) drop, parameters(32)
        read(unit, *) drop, parameters(33)
        read(unit, *) drop, parameters(34)
        read(unit, *) drop, parameters(35)
        read(unit, *) drop, parameters(36)
        read(unit, *) drop, parameters(37)
        read(unit, *) drop, parameters(38)
        read(unit, *) drop, parameters(39)
        read(unit, *) drop, parameters(40)
        read(unit, *) drop, parameters(41)
        read(unit, *) drop, parameters(42)
        close(unit)
    else
        print*, trim(filename), ' does not exist. stopping the program'
        stop
    endif

    parameters(43) = parameters( 1) - ccI01np
    parameters(44) = parameters( 2) - ccP01np
    parameters(45) = parameters( 3) - ccR01np
    parameters(46) = parameters( 1) - ccI01nn
    parameters(47) = parameters( 2) - ccP01nn
    parameters(48) = parameters( 3) - ccR01nn
    parameters(49) = parameters(19) - ctI1np
    parameters(50) = parameters(20) - ctQ1np
    parameters(51) = parameters(21) - ctR1np
    parameters(52) = parameters(22) - clsI1np
    parameters(53) = parameters(23) - clsP1np
    parameters(54) = parameters(24) - clsR1np
    
    ! call verify_symmetries(ccI01pp, ccP01pp, ccR01pp, ccI01nn, ccP01nn, ccR01nn, &
    !                        ccI11pp, ccP11pp, ccR11pp, ccI11np, ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, &
    !                         ctI1pp,  ctQ1pp,  ctR1pp,  ctI1np,  ctQ1np,  ctR1np,  ctI1nn,  ctQ1nn,  ctR1nn, &
    !                        clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, clsP1nn, clsR1nn, &
    !                        filename)

end subroutine read_marias_format

! subroutine verify_symmetries(ccI01pp, ccP01pp, ccR01pp, ccI01nn, ccP01nn, ccR01nn, &
!                              ccI11pp, ccP11pp, ccR11pp, ccI11np, ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, &
!                               ctI1pp,  ctQ1pp,  ctR1pp,  ctI1np,  ctQ1np,  ctR1np,  ctI1nn,  ctQ1nn,  ctR1nn, &
!                              clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, clsP1nn, clsR1nn, &
!                              filename)
!     implicit none
!     real(dp), intent(in) :: ccI01pp, ccP01pp, ccR01pp, ccI01nn, ccP01nn, ccR01nn, &
!                             ccI11pp, ccP11pp, ccR11pp, ccI11np, ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, &
!                              ctI1pp,  ctQ1pp,  ctR1pp,  ctI1np,  ctQ1np,  ctR1np,  ctI1nn,  ctQ1nn,  ctR1nn, &
!                             clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, clsP1nn, clsR1nn
!     character(len=*), intent(in) :: filename

!     logical :: violation
!     real(dp), parameter :: delta = 1.0e-8_dp
!     violation = .false.

!     if(abs(ccI01pp - ccI01nn) > delta) then
!         print*, 'CSB in ccI01!'
!         print*, 'I_01 c pp:', ccI01pp
!         print*, 'I_01 c nn:', ccI01nn
!         violation = .true.
!     endif

!     if(abs(ccR01pp - ccR01nn) > delta) then
!         print*, 'CSB in ccR01!'
!         print*, 'R_01 c pp:', ccR01pp
!         print*, 'R_01 c nn:', ccR01nn
!         violation = .true.
!     endif

!     if(abs(ctI1pp - ctI1np) > delta) then
!         print*, 'CSB in ctI1!'
!         print*, 'I_11 t pp:', ctI1pp
!         print*, 'I_11 t np:', ctI1np
!         violation = .true.
!     endif

!     if(abs(ctI1pp - ctI1nn) > delta) then
!         print*, 'CSB in ctI1!'
!         print*, 'I_11 t pp:', ctI1pp
!         print*, 'I_11 t nn:', ctI1nn
!         violation = .true.
!     endif

!     if(abs(ctQ1pp - ctQ1np) > delta) then
!         print*, 'CSB in ctQ1!'
!         print*, 'Q_11 t pp:', ctQ1pp
!         print*, 'Q_11 t np:', ctQ1np
!         violation = .true.
!     endif

!     if(abs(ctQ1pp - ctQ1nn) > delta) then
!         print*, 'CSB in ctQ1!'
!         print*, 'Q_11 t pp:', ctQ1pp
!         print*, 'Q_11 t nn:', ctQ1nn
!         violation = .true.
!     endif

!     if(abs(ctR1pp - ctR1np) > delta) then
!         print*, 'CSB in ctR1!'
!         print*, 'R_11 t pp:', ctR1pp
!         print*, 'R_11 t np:', ctR1np
!         violation = .true.
!     endif

!     if(abs(ctR1pp - ctR1nn) > delta) then
!         print*, 'CSB in ctR1!'
!         print*, 'R_11 t pp:', ctR1pp
!         print*, 'R_11 t nn:', ctR1nn
!         violation = .true.
!     endif


!     if(abs(clsI1pp - clsI1nn) > delta) then
!         print*, 'CSB in clsI1!'
!         print*, 'I_11 ls pp:', clsI1pp
!         print*, 'I_11 ls nn:', clsI1nn
!         violation = .true.
!     endif


!     if(abs(clsP1pp - clsP1nn) > delta) then
!         print*, 'CSB in clsP1!'
!         print*, 'P_11 ls pp:', clsP1pp
!         print*, 'P_11 ls nn:', clsP1nn
!         violation = .true.
!     endif


!     if(abs(clsR1pp - clsR1nn) > delta) then
!         print*, 'CSB in clsR1!'
!         print*, 'R_11 ls pp:', clsR1pp
!         print*, 'R_11 ls nn:', clsR1nn
!         violation = .true.
!     endif

!     if(abs(ccI01nn - ccI01pp - (ccI11nn - ccI11pp)) > delta) then
!         print*, 'CSB between I00 and I11'
!         print*, 'I_01_nn - I_01_pp:', ccI01nn - ccI01pp
!         print*, 'I_11_nn - I_11_pp:', ccI11nn - ccI11pp
!         violation = .true.        
!     endif

!     if(abs(ccP01nn - ccP01pp - (ccP11nn - ccP11pp)) > delta) then
!         print*, 'CSB between P00 and P11'
!         print*, 'P_01_nn - P_01_pp:', ccP01nn - ccP01pp
!         print*, 'P_11_nn - P_11_pp:', ccP11nn - ccP11pp
!         violation = .true.        
!     endif

!     if(abs(ccR01nn - ccR01pp - (ccR11nn - ccR11pp)) > delta) then
!         print*, 'CSB between R00 and R11'
!         print*, 'R_01_nn - R_01_pp:', ccR01nn - ccR01pp
!         print*, 'R_11_nn - R_11_pp:', ccR11nn - ccR11pp
!         violation = .true.        
!     endif

!     if(abs(ccI11nn - ccI11pp - 2*(ccI11np - ccI11pp)) > delta) then
!         print*, 'CSB between I11_pp, I11_np and I11_nn'
!         print*, '  I_11_nn - I_11_pp :', ccI11nn - ccI11pp
!         print*, '2(I_11_np - I_11_pp):', 2*(ccI11np - ccI11pp)
!         violation = .true.        
!     endif

!     if(abs(ccP11nn - ccP11pp - 2*(ccP11np - ccP11pp)) > delta) then
!         print*, 'CSB between P11_pp, P11_np and P11_nn'
!         print*, '  P_11_nn - P_11_pp :', ccP11nn - ccP11pp
!         print*, '2(P_11_np - P_11_pp):', 2*(ccP11np - ccP11pp)
!         violation = .true.        
!     endif

!     if(abs(ccR11nn - ccR11pp - 2*(ccR11np - ccR11pp)) > delta) then
!         print*, 'CSB between R11_pp, R11_np and R11_nn'
!         print*, '  R_11_nn - R_11_pp :', ccR11nn - ccR11pp
!         print*, '2(R_11_np - R_11_pp):', 2*(ccR11np - ccR11pp)
!         violation = .true.        
!     endif

!     if (violation) then
!         print*, 'Unexpected charge symmetry breaking in ', trim(filename)
!         print*, 'Check that those symmetry breaks are small enough for you'
!     endif

! end subroutine verify_symmetries
    
end module av18_compatibility