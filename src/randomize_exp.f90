!!
!> @brief       Randomizes data from a single experiment
!!
!! Module to randomize experimental point values to create new
!! experimental points. These points are built to follow the same
!! pattern as old experimental values from the Granada database.
!!
!! @author      Marielle Elizabeth Duran
!! 
      module randomize_exp
                use exp_data
                use precisions, only: dp
                use random_num, only: box_muller_num
                implicit none
                private
                public randomize_experiment, calc_mean, calc_stan_dev
                contains

!!
!> @brief       Creates new experimental values based on old points
!!
!! Function to form new eperimental point values based on old
!experimental point values and the random number generator subroutine.
! ADD MORE HERE 
                function randomize_experiment(old_exp) result(new_exp)
                        implicit none
                        type(nn_experiment), intent(in) :: old_exp
                        type(nn_experiment) :: new_exp
                        integer :: i, length
                        real(dp) :: num, mean, stan_dev
                        real(dp), dimension(SIZE(old_exp%data_points)) &
                                :: point_values
                        do i = 1, length
                           point_values(i) = old_exp%data_points(i)%value
                        end do
                        mean = calc_mean(point_values)
                        stan_dev = calc_stan_dev(point_values)
                        new_exp = old_exp
                        do i = 1, old_exp%n_data 
                        call box_muller_num(num)
                           new_exp%data_points(i)%value = &
                                        old_exp%data_points(i)%value&
                                        + num*stan_dev
                        end do
                        new_exp = new_exp 
                end function randomize_experiment

!!
!> @brief       Calculates mean value
!!
!! Function to calculate the mean of a generic array. This function is
!! used above with a more specific array in order to help find new
!! experimental points that are similar to points in the Granada
!! database.
!!
!! @author      Marielle Elizabeth Duran
!!  
                function calc_mean(numbers) result(mean)
                        implicit none
                        real(dp) :: mean
                        real(dp), dimension( : ) :: numbers !< generic array
                        integer :: i, length
                        length = SIZE(numbers) !< size of array
                        mean = 0.0_dp
                        do i = 1, length
                                mean = mean + numbers(i)
                        end do
                        mean = mean/length
                end function calc_mean

!!
!> @brief       Calculates standard deviation
!!
!! Function to calculate the standard deviation of a generic array. This
!! function uses a function that calculates the mean of a generic array
!! so that it is possible to calculate the standard deviation.
!!
!! @author      Marielle Elizabeth Duran
!!
                function calc_stan_dev(numbers) result(stan_dev)
                        implicit none
                        real(dp) :: variance, stan_dev, mean
                        real(dp), dimension( : ) :: numbers !< generic array
                        integer :: i, length
                        length = SIZE(numbers) !< size of array
                        variance = 0.0_dp
                        mean = calc_mean(numbers)
                        do i = 1, length
                                variance = variance + (&
                                        numbers(i) - mean)&
                                        **2.0_dp
                        end do
                        variance = variance / (length - 1)
                        stan_dev = SQRT(variance)
                end function calc_stan_dev
        
        end module randomize_exp
