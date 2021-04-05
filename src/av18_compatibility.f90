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
    write(unit, '(1l3,1E24.9,1x,a)') .false., 2.1_dp,   'cc'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.5_dp,   'cr'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.2_dp,   'ca'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(31), 'ccI01pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(32), 'ccP01pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  0.0_dp,         'ccR01pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(33), 'ccI01np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(34), 'ccP01np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.0_dp,         'ccR01np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(31), 'ccI01nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(41) + parameters(32), 'ccP01nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.0_dp, 'ccR01nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(35), 'cl2I01'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(36), 'cl2P01'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  0.0_dp,         'cl2R01'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(37), 'ccI00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(38), 'ccP00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.0_dp,         'ccR00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(39), 'cl2I00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(40), 'cl2P00'
    write(unit, '(1l3,1E24.9,1x,a)') .false., 0.0_dp,         'cl2R00'
    write(unit, '(1l3,1E24.9,1x,a)') .true., parameters( 1), 'ccI11pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true., parameters( 2), 'ccP11pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true., parameters( 3), 'ccR11pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 1), 'ccI11np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(41)*0.5_dp + parameters(2), 'ccP11np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 3), 'ccR11np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 1), 'ccI11nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(41) + parameters(2), 'ccP11nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 3), 'ccR11nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(10), 'cl2I11'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(11), 'cl2P11'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(12), 'cl2R11'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 4), 'ctI1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 5), 'ctQ1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 6), 'ctR1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 4), 'ctI1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 5), 'ctQ1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 6), 'ctR1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 4), 'ctI1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 5), 'ctQ1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 6), 'ctR1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 7), 'clsI1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 8), 'clsP1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters( 9), 'clsR1pp'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 7), 'clsI1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 8), 'clsP1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 9), 'clsR1np'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 7), 'clsI1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 8), 'clsP1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters( 9), 'clsR1nn'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(13), 'cls2I1'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(14), 'cls2P1'
    write(unit, '(1l3,1E24.9,1x,a)') .true.,  parameters(15), 'cls2R1'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(16), 'ccI10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(17), 'ccP10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(18), 'ccR10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(25), 'cl2I10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(26), 'cl2P10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(27), 'cl2R10'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(19), 'ctI0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(20), 'ctQ0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(21), 'ctR0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(22), 'clsI0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(23), 'clsP0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(24), 'clsR0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(28), 'cls2I0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(29), 'cls2P0'
    write(unit, '(1l3,1E24.9,1x,a)') .false., parameters(30), 'cls2R0'
    close(unit)
    
end subroutine write_marias_format

subroutine read_marias_format(filename, parameters)
    implicit none
    character(len=*), intent(in) :: filename
    real(dp), intent(out), dimension(:)  :: parameters

    logical :: file_exists, drop
    integer :: unit, i
    real(dp) :: ccI01pp, ccP01pp, ccR01pp, ccI01np, ccP01np, ccR01np, ccI01nn, ccP01nn, ccR01nn, cl2I01, &
        cl2P01, cl2R01, ccI00, ccP00, ccR00, cl2I00, cl2P00, cl2R00, ccI11pp, ccP11pp, ccR11pp, ccI11np, &
        ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, cl2I11, cl2P11, cl2R11, ctI1pp, ctQ1pp, ctR1pp, ctI1np, &
        ctQ1np, ctR1np, ctI1nn, ctQ1nn, ctR1nn, clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, &
        clsP1nn, clsR1nn, cls2I1, cls2P1, cls2R1, ccI10, ccP10, ccR10, cl2I10, cl2P10, cl2R10, ctI0, ctQ0, &
        ctR0, clsI0, clsP0, clsR0, cls2I0, cls2P0, cls2R0

    inquire(file=trim(filename), exist = file_exists)
    if (file_exists) then
        open(newunit= unit, file=trim(filename))
        do i=1,5
            read(unit, *)
        enddo
        read(unit, *) drop, ccI01pp
        read(unit, *) drop, ccP01pp
        read(unit, *) drop, ccR01pp
        read(unit, *) drop, ccI01np
        read(unit, *) drop, ccP01np
        read(unit, *) drop, ccR01np
        read(unit, *) drop, ccI01nn
        read(unit, *) drop, ccP01nn
        read(unit, *) drop, ccR01nn
        read(unit, *) drop, cl2I01
        read(unit, *) drop, cl2P01
        read(unit, *) drop, cl2R01
        read(unit, *) drop, ccI00
        read(unit, *) drop, ccP00
        read(unit, *) drop, ccR00
        read(unit, *) drop, cl2I00
        read(unit, *) drop, cl2P00
        read(unit, *) drop, cl2R00
        read(unit, *) drop, ccI11pp
        read(unit, *) drop, ccP11pp
        read(unit, *) drop, ccR11pp
        read(unit, *) drop, ccI11np
        read(unit, *) drop, ccP11np
        read(unit, *) drop, ccR11np
        read(unit, *) drop, ccI11nn
        read(unit, *) drop, ccP11nn
        read(unit, *) drop, ccR11nn
        read(unit, *) drop, cl2I11
        read(unit, *) drop, cl2P11
        read(unit, *) drop, cl2R11
        read(unit, *) drop, ctI1pp
        read(unit, *) drop, ctQ1pp
        read(unit, *) drop, ctR1pp
        read(unit, *) drop, ctI1np
        read(unit, *) drop, ctQ1np
        read(unit, *) drop, ctR1np
        read(unit, *) drop, ctI1nn
        read(unit, *) drop, ctQ1nn
        read(unit, *) drop, ctR1nn
        read(unit, *) drop, clsI1pp
        read(unit, *) drop, clsP1pp
        read(unit, *) drop, clsR1pp
        read(unit, *) drop, clsI1np
        read(unit, *) drop, clsP1np
        read(unit, *) drop, clsR1np
        read(unit, *) drop, clsI1nn
        read(unit, *) drop, clsP1nn
        read(unit, *) drop, clsR1nn
        read(unit, *) drop, cls2I1
        read(unit, *) drop, cls2P1
        read(unit, *) drop, cls2R1
        read(unit, *) drop, ccI10
        read(unit, *) drop, ccP10
        read(unit, *) drop, ccR10
        read(unit, *) drop, cl2I10
        read(unit, *) drop, cl2P10
        read(unit, *) drop, cl2R10
        read(unit, *) drop, ctI0
        read(unit, *) drop, ctQ0
        read(unit, *) drop, ctR0
        read(unit, *) drop, clsI0
        read(unit, *) drop, clsP0
        read(unit, *) drop, clsR0
        read(unit, *) drop, cls2I0
        read(unit, *) drop, cls2P0
        read(unit, *) drop, cls2R0
        close(unit)
    else
        print*, trim(filename), ' does not exist. stopping the program'
        stop
    endif

    parameters( 1) = ccI11pp
    parameters( 2) = ccP11pp
    parameters( 3) = ccR11pp
    parameters( 4) = ctI1pp
    parameters( 5) = ctQ1pp
    parameters( 6) = ctR1pp
    parameters( 7) = clsI1pp
    parameters( 8) = clsP1pp
    parameters( 9) = clsR1pp
    parameters(10) = cl2I11
    parameters(11) = cl2P11
    parameters(12) = cl2R11
    parameters(13) = cls2I1
    parameters(14) = cls2P1
    parameters(15) = cls2R1
    parameters(16) = ccI10
    parameters(17) = ccP10
    parameters(18) = ccR10
    parameters(19) = ctI0
    parameters(20) = ctQ0
    parameters(21) = ctR0
    parameters(22) = clsI0
    parameters(23) = clsP0
    parameters(24) = clsR0
    parameters(25) = cl2I10
    parameters(26) = cl2P10
    parameters(27) = cl2R10
    parameters(28) = cls2I0
    parameters(29) = cls2P0
    parameters(30) = cls2R0
    parameters(31) = ccI01pp
    parameters(32) = ccP01pp
    parameters(33) = ccI01np
    parameters(34) = ccP01np
    parameters(35) = cl2I01
    parameters(36) = cl2P01
    parameters(37) = ccI00
    parameters(38) = ccP00
    parameters(39) = cl2I00
    parameters(40) = cl2P00
    parameters(41) = ccP01nn - ccP01pp
    
    call verify_symmetries(ccI01pp, ccP01pp, ccR01pp, ccI01nn, ccP01nn, ccR01nn, &
                           ccI11pp, ccP11pp, ccR11pp, ccI11np, ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, &
                            ctI1pp,  ctQ1pp,  ctR1pp,  ctI1np,  ctQ1np,  ctR1np,  ctI1nn,  ctQ1nn,  ctR1nn, &
                           clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, clsP1nn, clsR1nn, &
                           filename)

end subroutine read_marias_format

subroutine verify_symmetries(ccI01pp, ccP01pp, ccR01pp, ccI01nn, ccP01nn, ccR01nn, &
                             ccI11pp, ccP11pp, ccR11pp, ccI11np, ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, &
                              ctI1pp,  ctQ1pp,  ctR1pp,  ctI1np,  ctQ1np,  ctR1np,  ctI1nn,  ctQ1nn,  ctR1nn, &
                             clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, clsP1nn, clsR1nn, &
                             filename)
    implicit none
    real(dp), intent(in) :: ccI01pp, ccP01pp, ccR01pp, ccI01nn, ccP01nn, ccR01nn, &
                            ccI11pp, ccP11pp, ccR11pp, ccI11np, ccP11np, ccR11np, ccI11nn, ccP11nn, ccR11nn, &
                             ctI1pp,  ctQ1pp,  ctR1pp,  ctI1np,  ctQ1np,  ctR1np,  ctI1nn,  ctQ1nn,  ctR1nn, &
                            clsI1pp, clsP1pp, clsR1pp, clsI1np, clsP1np, clsR1np, clsI1nn, clsP1nn, clsR1nn
    character(len=*), intent(in) :: filename

    logical :: violation
    real(dp), parameter :: delta = 1.0e-8_dp
    violation = .false.

    if(abs(ccI01pp - ccI01nn) > delta) then
        print*, 'CSB in ccI01!'
        print*, 'I_01 c pp:', ccI01pp
        print*, 'I_01 c nn:', ccI01nn
        violation = .true.
    endif

    if(abs(ccR01pp - ccR01nn) > delta) then
        print*, 'CSB in ccR01!'
        print*, 'R_01 c pp:', ccR01pp
        print*, 'R_01 c nn:', ccR01nn
        violation = .true.
    endif

    if(abs(ccI11pp - ccI11np) > delta) then
        print*, 'CSB in ccI11!'
        print*, 'I_11 c pp:', ccI11pp
        print*, 'I_11 c np:', ccI11np
        violation = .true.
    endif

    if(abs(ccI11pp - ccI11nn) > delta) then
        print*, 'CSB in ccI11!'
        print*, 'I_11 c pp:', ccI11pp
        print*, 'I_11 c nn:', ccI11nn
        violation = .true.
    endif

    if(abs(ccR11pp - ccR11np) > delta) then
        print*, 'CSB in ccR11!'
        print*, 'R_11 c pp:', ccR11pp
        print*, 'R_11 c np:', ccR11np
        violation = .true.
    endif

    if(abs(ccR11pp - ccR11nn) > delta) then
        print*, 'CSB in ccR11!'
        print*, 'R_11 c pp:', ccR11pp
        print*, 'R_11 c nn:', ccR11nn
        violation = .true.
    endif

    if(abs(ctI1pp - ctI1np) > delta) then
        print*, 'CSB in ctI1!'
        print*, 'I_11 t pp:', ctI1pp
        print*, 'I_11 t np:', ctI1np
        violation = .true.
    endif

    if(abs(ctI1pp - ctI1nn) > delta) then
        print*, 'CSB in ctI1!'
        print*, 'I_11 t pp:', ctI1pp
        print*, 'I_11 t nn:', ctI1nn
        violation = .true.
    endif

    if(abs(ctQ1pp - ctQ1np) > delta) then
        print*, 'CSB in ctQ1!'
        print*, 'Q_11 t pp:', ctQ1pp
        print*, 'Q_11 t np:', ctQ1np
        violation = .true.
    endif

    if(abs(ctQ1pp - ctQ1nn) > delta) then
        print*, 'CSB in ctQ1!'
        print*, 'Q_11 t pp:', ctQ1pp
        print*, 'Q_11 t nn:', ctQ1nn
        violation = .true.
    endif

    if(abs(ctR1pp - ctR1np) > delta) then
        print*, 'CSB in ctR1!'
        print*, 'R_11 t pp:', ctR1pp
        print*, 'R_11 t np:', ctR1np
        violation = .true.
    endif

    if(abs(ctR1pp - ctR1nn) > delta) then
        print*, 'CSB in ctR1!'
        print*, 'R_11 t pp:', ctR1pp
        print*, 'R_11 t nn:', ctR1nn
        violation = .true.
    endif

    if(abs(clsI1pp - clsI1np) > delta) then
        print*, 'CSB in clsI1!'
        print*, 'I_11 ls pp:', clsI1pp
        print*, 'I_11 ls np:', clsI1np
        violation = .true.
    endif

    if(abs(clsI1pp - clsI1nn) > delta) then
        print*, 'CSB in clsI1!'
        print*, 'I_11 ls pp:', clsI1pp
        print*, 'I_11 ls nn:', clsI1nn
        violation = .true.
    endif

    if(abs(clsP1pp - clsP1np) > delta) then
        print*, 'CSB in clsP1!'
        print*, 'P_11 ls pp:', clsP1pp
        print*, 'P_11 ls np:', clsP1np
        violation = .true.
    endif

    if(abs(clsP1pp - clsP1nn) > delta) then
        print*, 'CSB in clsP1!'
        print*, 'P_11 ls pp:', clsP1pp
        print*, 'P_11 ls nn:', clsP1nn
        violation = .true.
    endif

    if(abs(clsR1pp - clsR1np) > delta) then
        print*, 'CSB in clsR1!'
        print*, 'R_11 ls pp:', clsR1pp
        print*, 'R_11 ls np:', clsR1np
        violation = .true.
    endif

    if(abs(clsR1pp - clsR1nn) > delta) then
        print*, 'CSB in clsR1!'
        print*, 'R_11 ls pp:', clsR1pp
        print*, 'R_11 ls nn:', clsR1nn
        violation = .true.
    endif

    if(abs(ccP01nn - ccP01pp - (ccP11nn - ccP11pp)) > delta) then
        print*, 'CSB between P00 and P11'
        print*, 'P_01_nn - P_01_pp:', ccP01nn - ccP01pp
        print*, 'P_11_nn - P_11_pp:', ccP11nn - ccP11pp
        violation = .true.        
    endif

    if(abs(ccP11nn - ccP11pp - 2*(ccP11np - ccP11pp)) > delta) then
        print*, 'CSB between P11_pp, P11_np and P11_nn'
        print*, '  P_11_nn - P_11_pp :', ccP11nn - ccP11pp
        print*, '2(P_11_np - P_11_pp):', 2*(ccP11np - ccP11pp)
        violation = .true.        
    endif

    if (violation) then
        print*, 'stopping due to unexpected charge symmetry breaking in ', trim(filename)
        stop
    endif

end subroutine verify_symmetries
    
end module av18_compatibility