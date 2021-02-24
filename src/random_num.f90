        module random_num
                use precisions, only: dp
                use constants, only: pi
                IMPLICIT NONE
                PRIVATE 
                PUBLIC box_muller_num, generator_100_num, &
                        verify_box_muller_num
                CONTAINS
                subroutine box_muller_num(z_1)
                        IMPLICIT NONE
                        REAL(dp), intent(out) :: z_1
                        REAL(dp) :: x, y
                        REAL(dp), save :: z_0
                        LOGICAL, save :: not_available = .TRUE.
                        IF (not_available) THEN
                                not_available = .FALSE.
                        call random_number(x)
                        call random_number(y)
                        z_0 = SQRT(-2.0*LOG(x)) * COS(2.0*pi*y)
                        z_1 = SQRT(-2.0*LOG(x)) * SIN(2.0*pi*y)
                        ELSE
                                z_1 = z_0
                                not_available = .TRUE.
                        END IF 
                end subroutine box_muller_num
                subroutine generator_100_num
                        IMPLICIT NONE
                        real(dp) :: x
                        integer :: i
                        call box_muller_num(x)
                        PRINT*, x
                        open(7, file = 'generator_info.dat', status =&
                                'unknown')
                        DO i = 1,100
                                call box_muller_num(x)
                                WRITE(7,*) x
                        END DO
                        close(7)
                end subroutine generator_100_num
                subroutine verify_box_muller_num(n_samples)
                        IMPLICIT NONE
                        real(dp) :: mean, num, variance, stan_dev
                        real, dimension(n_samples) :: random_numbers 
                        integer, intent(in) :: n_samples
                        integer :: i
                        DO i = 1, n_samples
                        call box_muller_num(num)
                                random_numbers(i) = num
                        END DO
                        mean = 0.0_dp
                        DO i = 1, n_samples
                                mean = mean + random_numbers(i)
                        END DO
                        mean = mean/n_samples
                        variance = 0.0_dp
                        DO i = 1, n_samples
                                variance = variance + (&
                                        random_numbers(i) - mean)&
                                        **2.0_dp
                        END DO
                        variance = variance / (n_samples - 1)
                        stan_dev = SQRT(variance)
                        print*, mean
                        print*, stan_dev
                end subroutine verify_box_muller_num       
        end module random_num  
