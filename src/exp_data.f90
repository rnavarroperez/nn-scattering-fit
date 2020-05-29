!!
!> @brief      Experimental nn scattering data
!!
!! Module to read NN experimental data and load it into a database
!!
!! @author     Rodrigo Navarro Perez
!!
module exp_data
use precisions, only: dp
use string_functions, only: lower
use num_recipes, only: int_to_logical
implicit none
private
public read_old_data_base, nn_experiment, write_database, read_database

!!
!> @brief      a single experimental point
!!
!! Type to define a single experimental NN scattering point.
!! It contains the variables kinematic (laboratory energy and scattering angle),
!! the experimental value, statistical uncertainty, and the corresponding
!! electromagnetic amplitude
!!
!! @author     Rodrigo Navarro Perez
!!
type :: exp_point
    real(dp) :: t_lab !< laboratory energy in MeV
    real(dp) :: theta !< scattering angle in degrees
    real(dp) :: value !< experimental value
    real(dp) :: stat_error !< statistical error
    complex(dp), dimension(1:5) :: em_amplitude !< electromagnetic amplitude
end type exp_point

!!
!> @brief      A single experiment
!!
!! type to define a NN scattering experiment, it consists of an array of 
!! experimental points and data common to all experimental points; namely 
!! number of data points, systematic error, type of observable, reaction 
!! channel, whether or not the experiment is to be rejected from the 
!! analysis, the year, and a bibliographical reference.
!!
!! @author     Rodrigo Navarro Perez
!!
type :: nn_experiment
    integer :: n_data !< number of data points in the experiment
    real(dp) :: sys_error !< systematic error
    character(len=4) :: obs_type !< type of observable
    character(len=2) :: channel !< reaction channel
    logical :: rejected !< should the experiment be rejected from the analysis
    integer :: year !< year the experiment was published
    character(len=120) :: reference !< bibliographical reference
    type(exp_point), allocatable :: data_points(:) !< data points
end type nn_experiment
contains

!!
!> @brief      Reads a database.
!!
!! Reads a database from a text file and stores it into a 
!! database. See the 'granada_database.dat' file for 
!! the specific format
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine read_database(data_file, experiments)
    use iso_fortran_env, only : iostat_end
    implicit none
    character(len=*), intent(in) :: data_file !< Text file with the database
    type(nn_experiment), intent(out), allocatable, dimension(:) :: experiments !< database

    integer :: n_data, year
    real(dp) :: sys_error, t_lab, theta, value, stat_error
    character(len=4) :: obs_type
    character(len=2) :: channel
    logical :: rejected
    character(len=120) :: reference

    integer :: counter, limit, unit, io_status, i
    character(len=300) :: line
    allocate(experiments(1:1))
    counter = 0
    limit = 1
    open(newunit=unit, file=trim(data_file), status='OLD')
    do
        read(unit, '(a)', iostat=io_status) line
        if (io_status == iostat_end) exit
        counter = counter + 1
        if (counter > limit) then
            call double_allocation(experiments)
            limit = 2*limit
        endif
        read(line, *) n_data, sys_error, rejected, channel, obs_type, year, reference
        experiments(counter)%n_data = n_data
        experiments(counter)%sys_error = sys_error
        experiments(counter)%rejected = rejected
        experiments(counter)%channel = trim(channel)
        experiments(counter)%obs_type = trim(obs_type)
        experiments(counter)%year = year
        experiments(counter)%reference = trim(reference)
        allocate(experiments(counter)%data_points(1:n_data))
        do i=1, n_data
            read(unit, *) t_lab, theta, value, stat_error
            experiments(counter)%data_points(i)%t_lab = t_lab
            experiments(counter)%data_points(i)%theta = theta
            experiments(counter)%data_points(i)%value = value
            experiments(counter)%data_points(i)%stat_error = stat_error
        enddo
    enddo
    close(unit)
    call trim_database(counter, experiments)
    
end subroutine read_database

!!
!> @brief      Writes a database.
!!
!! Given two databases (one for pp data and one for np data),
!! the subroutine combines them in a single array and writes 
!! the full database into a text file
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine write_database(pp_data, np_data, file_name)
    implicit none
    type(nn_experiment), intent(in), dimension(:) :: pp_data !< pp database
    type(nn_experiment), intent(in), dimension(:) :: np_data !< np database
    character(len=*), intent(in) :: file_name !< file name to write down the full database

    type(nn_experiment), allocatable, dimension(:) :: all_data
    integer :: unit, pp_size, np_size, i, j

    pp_size = size(pp_data)
    np_size = size(np_data)
    allocate(all_data(1:pp_size + np_size))
    all_data(1:pp_size) = pp_data
    all_data(pp_size + 1:pp_size + np_size) = np_data
    open(newunit=unit, file=trim(file_name))
    do i=1, size(all_data)
        write(unit,'(i4,f15.8,l5,2(3x,a4),i5,a)') all_data(i)%n_data, all_data(i)%sys_error, all_data(i)%rejected, &
            all_data(i)%channel, all_data(i)%obs_type, all_data(i)%year, trim(all_data(i)%reference)
            do j=1, all_data(i)%n_data
                write(unit, '(4f18.8)') all_data(i)%data_points(j)%t_lab, all_data(i)%data_points(j)%theta, &
                    all_data(i)%data_points(j)%value, all_data(i)%data_points(j)%stat_error
            enddo
    enddo
    close(unit)
    
end subroutine write_database

!!
!> @brief      Reads a database in the old format.
!!
!! Given a text file with database in the old format, and the
!! corresponding reaction channel, the subroutine reads the database
!! and stores it in the experiments array
!!
!! The old format contains a lot of additional information that is 
!! not relevant. This makes the reading of the file a much more
!! convoluted process. The new format is preferred when reading 
!! and writing any database
!!
!! @author    Rodrigo Navarro Perez
!!
subroutine read_old_data_base(data_file, reaction, experiments)
    implicit none
    character(len=*), intent(in) :: data_file !< Text file with the database in the old format
    character(len=2), intent(in) :: reaction !< reaction channel
    type(nn_experiment), intent(out), allocatable :: experiments(:) !< experimental database

    integer :: unit, counter, n, i_data, rejec_unit
    real(dp) :: t_lab, sys_error, normalization
    integer :: n_data, n_keep, year, n_final
    character(len=1) :: keepN, keepL, KeepS, rejecN, rejecL, rejecS, space, a_char
    character(len=2) :: a_string, rejecA, rejecB, ref1
    character(len=4) :: reac, string_obs, published
    character(len=9) :: some_string
    character(len=115) :: reference
    real(dp), allocatable, dimension(:) :: tlabs, angles, values, deltas

    logical :: in_data_base
    integer, dimension(1:2) :: nn_reject
    integer :: ireject

    integer, parameter :: n_np_obs = 26, n_pp_obs = 23
    character(len=4), dimension(1:n_np_obs), parameter :: &
       all_types = ['DSG ','DT  ','AYY ','D   ','P   ','AZZ ','R   '&
         ,'RT  ','RPT ','AT  ','D0SK','NSKN','NSSN','NNKK','A   '&
         ,'AXX ','CKP ','RP  ','MSSN','MSKN','AZX ','AP  ','DTRT'&
         ,'SGT ','SGTT','SGTL']

    n = 1
    counter = 0
    allocate(experiments(1:1))
    open(newunit=unit, file=trim(data_file), status='OLD')
    open(newunit=rejec_unit, file='database/old_format/rejected.dat', status='OLD')
    do
        read(unit,102) t_lab, n_data, sys_error, normalization, n_keep, keepN, keepL, keepS, &
            a_string, rejecA, rejecN, rejecL, rejecS, rejecB, reac, &
            space, string_obs, a_char, ref1, space, year, some_string, published
        if(t_lab.gt.350.0D0) exit
        read(unit,'(A115)') reference
        select case (reaction)
        case ('pp')
            in_data_base = published.eq.'PUBL' .and. any(all_types(1:n_pp_obs) == string_obs)
            ireject = 1
        case ('np')
            in_data_base = (RejecB.ne.'c '.and.RejecB.ne.'U ') .and. any(all_types == string_obs)
            ireject = 2
        case default
            stop 'incorrect reaction channel in all_phaseshifts'
        end select
        if (allocated(tlabs)) deallocate(tlabs, angles, values, deltas)
        allocate(tlabs(1:n_data))
        allocate(angles, values, deltas, mold=tlabs)
        read(unit, 164) (tlabs(i_data), angles(i_data), values(i_data), deltas(i_data), i_data=1, n_data)
        if (in_data_base) then
            counter = counter + 1
            if (counter>n) then
                call double_allocation(experiments)
                n = 2*n
            endif
            n_final = 0
            do i_data=1, n_data
                if (deltas(i_data)>0) n_final = n_final + 1
            enddo
            experiments(counter)%n_data = n_final
            experiments(counter)%sys_error = abs(sys_error)
            experiments(counter)%obs_type = lower(trim(string_obs))
            experiments(counter)%reference = trim(reference)
            experiments(counter)%channel = lower(trim(reac))
            read(rejec_unit,*) nn_reject(1:2)
            experiments(counter)%rejected = int_to_logical(nn_reject(ireject))
            if(year.le.14) then
                year = year + 2000
            else
                year = year + 1900
            endif
            experiments(counter)%year = year
            allocate(experiments(counter)%data_points(1:n_final))
            n_final = 0
            do i_data=1, n_data
                if (deltas(i_data)>0) then
                    n_final = n_final + 1
                    experiments(counter)%data_points(n_final)%t_lab = tlabs(i_data)
                    experiments(counter)%data_points(n_final)%theta = angles(i_data)
                    experiments(counter)%data_points(n_final)%value = values(i_data)
                    experiments(counter)%data_points(n_final)%stat_error = deltas(i_data)
                endif
            enddo
            
        endif
    enddo
    close(unit)
    close(rejec_unit)
    call trim_database(counter, experiments)

102   format(f10.6,I4,2f7.3,I3,3A1,A2,A2,3A1,A2,A4,A1,A4,A1,A2,A1,I2,A9,A4)
164   format(2f13.6, f13.5, f13.5)

end subroutine read_old_data_base

!!
!> @brief      trims a database
!!
!! Given an array of type nn_experiment, returns a reallocated
!! array of a size determined by the cut_off
!!
!! This subroutine is necessary due to the fact that the size
!! of the database is not know when the database is about to be 
!! read from a file and the size of the array needs to be 
!! constantly increased. Once the reading is done, it is possible
!! that the allocated array is larger than the actual number of experiments.
!! This subroutine removes the extra elements of the array that are note used
!!
!! @author Rodrigo Navarro Perez
!!
subroutine trim_database(cut_off, array)
    implicit none
    integer, intent(in) :: cut_off !< Size of the database to keep
    type(nn_experiment), intent(inout), allocatable, dimension(:) :: array !< database to trim

    type(nn_experiment), allocatable, dimension(:) :: temp
    if (cut_off > size(array)) then
        stop 'cut_off has to be smaller than array size in trim_database'
    endif
    allocate(temp(1:cut_off))
    temp = array(1:cut_off)
    call move_alloc(temp, array)
    
end subroutine trim_database

!!
!> @brief      Doubles the allocation of a database
!!
!! Given an array of type nn_experiment, returns a reallocated
!! array with double the size. The original elements are kept
!! in the first half of the new array
!!
!! This subroutine is necessary due to the fact that the size
!! of the database is not know when the database is about to be 
!! read from a file and the size of the array needs to be 
!! constantly doubled.
!!
!! @author     Rodrigo Navarro Perez
!!
subroutine double_allocation(array)
    implicit none
    type(nn_experiment), intent(inout), allocatable, dimension(:) :: array !< database to expand

    type(nn_experiment), allocatable, dimension(:) :: temp
    integer :: array_size
    array_size = size(array)
    allocate(temp(1:2*array_size))
    temp(1:array_size) = array
    call move_alloc(temp, array)
    
end subroutine double_allocation

end module exp_data