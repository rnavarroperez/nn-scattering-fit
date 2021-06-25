!!
!> @brief       Randomizes data from a single experiment
!!
!! Module to randomize experimental point values to create new
!! experimental points. These points are built to follow the same
!! pattern as old experimental values from the Granada database.
!! This module also randomizes databases using the function written
!! to randomize single experiments.      
!!      
!! @author      Marielle Elizabeth Duran
!! 
      module randomize_exp
                use exp_data
                use precisions, only: dp
                use random_num, only: box_muller_num
                implicit none
                private
                public randomize_experiment, randomize_database, calc_stan_dev
                contains

!!
!> @brief       Creates new experimental values based on old points
!!
!! Function to form new experimental point values based on old
!! experimental point values and the random number generator subroutine.
!! The values for the old experimental point values are stored in the
!! array, point_values. The subroutine for the random number generator
!! is called within the do loop. New experimental point values are
!! calculated using the mean, the randomly generated number, and
!! standard deviation.
!!
!! @author      Marielle Elizabeth Duran
!!
                function randomize_experiment(old_exp) result(new_exp)
                        implicit none
                        type(nn_experiment), intent(in) :: old_exp
                        type(nn_experiment) :: new_exp
                        integer :: i
                        real(dp) :: num, mean, stan_dev
                        real(dp), dimension(SIZE(old_exp%data_points)) &
                                :: point_values !< array for old values
                        new_exp = old_exp
                        do i = 1, old_exp%n_data
                                mean = old_exp%data_points(i)%value
                             stan_dev =old_exp%data_points(i)%stat_error
                           point_values(i) = old_exp%data_points(i)%value
                        call box_muller_num(num)
                           new_exp%data_points(i)%value = &
                                        mean + (num*stan_dev)
                        end do
                end function randomize_experiment

!!
!> @brief       Creates a database with randomly generated numbers
!!
!! Function to randomize a data base using the function,
!! randomize_experiment. That function randomizes a single experiment
!! and randomize_database puts each experiment into an array that 
!! represents a database.
!!
!! @author      Marielle Elizabeth Duran
!! 
                function randomize_database(old_data) result(new_data)
                        implicit none
               type(nn_experiment), dimension(:), intent(in) :: old_data
              type(nn_experiment), allocatable, dimension(:) :: new_data
                        integer :: i
                        allocate(new_data(1:SIZE(old_data)))
                        do i = 1, SIZE(old_data)
                         new_data(i) = randomize_experiment(old_data(i))
                        end do
                end function randomize_database

!!
!> @brief       Calculates mean value
!!
!! Function to calculate the mean of a generic array. This function was
!! not used but was saved.
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
