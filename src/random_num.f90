        module random_num
                IMPLICIT NONE
                PRIVATE 
                PUBLIC box_muller_num
                CONTAINS
                subroutine box_muller_num(z_1)
                use precisions, only: dp
                use constants, only: pi
                        IMPLICIT NONE
                        REAL(dp), intent(out) :: z_1
                        REAL(dp) :: x, y
                        REAL(dp), save :: z_0
                        LOGICAL, save :: k = .TRUE.
                        IF (k) THEN
                                k = .FALSE.
                        call random_number(x)
                        call random_number(y)
                        z_0 = SQRT(-2.0*LOG(x)) * COS(2.0*pi*y)
                        z_1 = SQRT(-2.0*LOG(x)) * SIN(2.0*pi*y)
                        ELSE
                                z_1 = z_0
                                k = .TRUE.
                        END IF 
                end subroutine box_muller_num
                subroutine verify_box_muller_num(n_samples)
                use precisions, only: dp
                        IMPLICIT NONE
                        REAL (dp) :: mean, stan_dev, num
                        Integer, intent(in) :: n_samples
                        INTEGER :: i
                        mean = 0.0_dp
                        DO i = 1, n_samples
                        call box_muller_num(num)
                                mean = mean + num
                        END DO
                        mean = mean/n_samples
                end subroutine verify_box_muller_num       
        end module random_num  
