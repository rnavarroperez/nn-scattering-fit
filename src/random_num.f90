!!
!> @brief       Random Number Generator
!! 
!! Module to generate random numbers that follow a Gaussian distribution
!! and the Box-Muller form. In addtion, there is a 100 number generator that
!! writes the numbers in an external file and a subroutine to verify the
!! Box-Muller number generator.
!!
!! @author      Marielle Elizabeth Duran
!!
      module random_num
                use precisions, only: dp
                use constants, only: pi
                IMPLICIT NONE
                PRIVATE 
                PUBLIC box_muller_num
                CONTAINS

!!
!> @brief       Box-Muller random number generator
!!
!! Generates random numbers that fit a Gaussian distribution and the
!! Box-Muller transform. Each value is saved so that a repeat value is
!! rejected.
!!
!! @author      Marielle Elizabeth Duran
!!
                subroutine box_muller_num(z_1)
                        IMPLICIT NONE
                        REAL(dp), intent(out) :: z_1
                        REAL(dp) :: x, y
                        REAL(dp), save :: z_0 !< saved Box-Muller number
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

!!
!> @brief       Generates 100 numbers with Box-Muller form
!!
!! Calls a random number using the box_muller_num subroutine. This
!! subroutine generates 100 of these numbers and writes them to the
!! external file named 'generator_info.dat'.
!!
!! @author      Marielle Elizabeth Duran
!!
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

!!
!> @brief       Verifies the box_muller_num subroutine results
!!
!! Creates an array to put in random numbers from box_muller_num. This
!! array is used to calculate the mean and standard deviation of the
!! random numbers. These steps are completed to verify that the
!! box_muller_num subroutine follows a Gaussian distribution.
!!
!! @author      Marielle Elizabeth Duran
!!
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
