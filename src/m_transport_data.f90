!> Module that provides routines for reading in transport data
module m_transport_data

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  !> The maximum number of rows of transport data for an entry
  integer, parameter :: max_num_rows = 2000
  integer, parameter :: lineLen = 200

  public :: TD_get_td_from_file

contains

  !> Routine to read in tabulated data from a file
  subroutine TD_get_td_from_file(file_name, data_name, x_data, y_data, data_found)
    character(len=*), intent(in)       :: file_name, data_name
    real(dp), allocatable, intent(out) :: x_data(:), y_data(:)
    logical, optional, intent(out)     :: data_found

    ! Temporary variables
    integer                :: ioState, nL
    integer                :: n_rows
    integer                :: my_unit
    character(LEN=40)      :: line_fmt
    character(LEN=lineLen) :: line
    real(dp)               :: temp_table(2, max_num_rows)
    real(dp)               :: factor

    nL = 0 ! Set the number of lines to 0

    ! Set the line format to read, only depends on string_len currently
    write(line_fmt, FMT = "(I6)") lineLen
    line_fmt = "(A" // trim(adjustl(line_fmt)) // ")"

    ! Open 'file_name' (with error checking)
    open(newunit=my_unit, file = trim(file_name), action = "read", &
         err = 999, iostat = ioState, status="old")

    ! Table format

    !     table_name
    !     FACTOR: 1.0                   [optional: multiply with this factor]
    !     [other lines]
    !     ------------------            [at least 5 dashes]
    !     xxx       xxx                 [data in two column format]
    !     ...       ...
    !     xxx       xxx
    !     ------------------

    ! The outer DO loop, running until the end of the file is reached
    do
       ! Search for 'data_name' in the file
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 888) line; nL = nL+1
          if (line == data_name) exit
       end do

       factor = 1.0_dp

       ! Now we can check whether there is a comment, while scanning lines until
       ! dashes are found, which indicate the start of the data
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit
          else if (line(1:7) == "FACTOR:") then
             read(line(8:), *) factor
          else if (line(1:8) == "COMMENT:") then
             continue
          else if (len_trim(line) < 5) then
             ! Assume this is the name of a gas mixture
             print *, "In file ", trim(file_name), " at line", nL
             print *, "Warning: gas names are now deprecated"
             continue
          else
             print *, "In file ", trim(file_name), " at line", nL
             print *, trim(line)
             error stop "Unknown statement in input file"
          end if
       end do

       ! Read the data into a temporary array
       n_rows = 0
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit  ! Dashes mark the end of the data
          else if (trim(line) == "" .or. line(1:1) == "#") then
             cycle ! Ignore whitespace or comments
          else if (n_rows < max_num_rows) then
             n_rows = n_rows + 1
             read(line, FMT = *, ERR = 999, end = 777) temp_table(:, n_rows)
          else
             print *, "Too many rows in ", file_name, " at line ", nL
          end if
       end do

       ! Store the data in the actual table
       if (allocated(x_data)) deallocate(x_data)
       if (allocated(y_data)) deallocate(y_data)
       allocate(x_data(n_rows))
       allocate(y_data(n_rows))

       x_data = temp_table(1, 1:n_rows)
       y_data = factor * temp_table(2, 1:n_rows)

       exit                   ! Done
    end do

    close(my_unit)
    if (present(data_found)) data_found = .true.

    ! Normal exit here
    return

777 continue ! If the end of the file is reached after finding data
    print *, "TD_get_td_from_file reached end of file while reading data"
    print *, "searching '" // trim(data_name) // "'"
    stop

888 continue ! If the end of the file is reached without finding data
    if (present(data_found)) then
       data_found = .false.
       return
    else
       print *, "TD_get_td_from_file reached end of file without finding data"
       print *, "searching '" // trim(data_name) // "'"
       stop
    end if

999 continue ! If there was an input error, the routine will end here
    print *, "TD_get_td_from_file error at line", nL
    print *, "ioState = ", ioState, " in ", file_name
    print *, "searching '" // trim(data_name) // "'"
    stop

  end subroutine TD_get_td_from_file

end module m_transport_data
