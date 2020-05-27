module exp_data
use precisions, only: dp
implicit none
private
public read_old_data_base, nn_experiment

type :: exp_point
    real(dp) :: t_lab
    real(dp) :: theta
    real(dp) :: value
    real(dp) :: stat_error
    complex(dp) :: saclay_a
    complex(dp) :: saclay_b
    complex(dp) :: saclay_c
    complex(dp) :: saclay_d
    complex(dp) :: saclay_e
end type exp_point

type :: nn_experiment
    integer :: n_data
    real(dp) :: sys_error
    real(dp) :: normalization
    character(len=4) :: obs_type
    logical :: rejected
    integer :: year
    character(len=120) :: reference
    type(exp_point), allocatable :: data_points(:)
end type nn_experiment
contains

subroutine read_old_data_base(data_file, reaction, experiments)
    implicit none
    character(len=*), intent(in) :: data_file
    character(len=2), intent(in) :: reaction
    type(nn_experiment), intent(out), allocatable :: experiments(:)

    integer :: unit, counter, n, i_data
    ! type(nn_experiment), allocatable, dimension(:) :: temp

    real(dp) :: t_lab, sys_error, normalization
    integer :: n_data, n_keep, year
    character(len=1) :: keepN, keepL, KeepS, rejecN, rejecL, rejecS, space, a_char
    character(len=2) :: a_string, rejecA, rejecB, ref1
    character(len=4) :: reac, string_obs, published
    character(len=9) :: some_string
    character(len=115) :: reference

    logical :: in_data_base

    n = 1
    counter = 1
    allocate(experiments(1:325))
    open(newunit=unit, file=trim(data_file), status='OLD')
    do
        ! if (i>n) then
        !     allocate(temp(1:2*n))
        !     temp(1:n) = r
        !     call move_alloc(temp, r)
        !     n = 2*n
        ! endif
        read(unit,102) t_lab, n_data, sys_error, normalization, n_keep, keepN, keepL, keepS, &
            a_string, rejecA, rejecN, rejecL, rejecS, rejecB, reac, &
            space, string_obs, a_char, ref1, space, year, some_string, published
        if(t_lab.gt.350.0D0) exit
        read(unit,'(A115)') reference
        select case (reaction)
        case ('pp')
            in_data_base = published.eq.'PUBL'
        case ('np')
            in_data_base = (RejecB.ne.'c '.and.RejecB.ne.'U ')
        case default
            stop 'incorrect reaction channel in all_phaseshifts'
        end select
        allocate(experiments(counter)%data_points(1:n_data))
        do i_data=1, n_data
            read(unit, *) experiments(counter)%data_points(i_data)%t_lab, &
                experiments(counter)%data_points(i_data)%theta, &
                experiments(counter)%data_points(i_data)%value, &
                experiments(counter)%data_points(i_data)%stat_error
        enddo
        counter = counter + 1 
    enddo
    close(unit)

102   format(f10.6,I4,2f7.3,I3,3A1,A2,A2,3A1,A2,A4,A1,A4,A1,A2,A1,I2,A9,A4)


end subroutine read_old_data_base

end module exp_data