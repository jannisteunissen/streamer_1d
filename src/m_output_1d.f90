module m_output_1d
  use m_model_choice

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer :: OUT_num_bins, OUT_num_fields
  real(dp) :: OUT_head_density
  real(dp) :: OUT_eedf_eV_range(2)
  real(dp), allocatable :: OUT_eedf_min_max_fields(:)

  public :: OUT_write_vars
  public :: OUT_write_coeffs
  public :: OUT_init

contains

  subroutine OUT_init(cfg)
    use m_config
    type(CFG_t), intent(in) :: cfg

    call CFG_get(cfg, "output_head_density", OUT_head_density)
    call CFG_get(cfg, "output_num_bins", OUT_num_bins)
    call CFG_get_size(cfg, "output_eedf_min_max_fields", OUT_num_fields)
    OUT_num_fields = OUT_num_fields / 2
    allocate(OUT_eedf_min_max_fields(2 * OUT_num_fields))
    call CFG_get(cfg, "output_eedf_min_max_fields", OUT_eedf_min_max_fields)
    call CFG_get(cfg, "output_eedf_eV_range", OUT_eedf_eV_range)
  end subroutine OUT_init

  subroutine OUT_write_vars(sim_name, sim_type, pc, cntr, time)
    use m_fluid_dd_1d, only: FL_get_output
    use m_particle_core

    type(PC_t), intent(in) :: pc
    character(len=*), intent(in)   :: sim_name
    integer, intent(in)            :: sim_type, cntr
    real(dp), intent(in)           :: time

    integer                        :: ix, n_pos, n_sca
    character(len=100)             :: filename, eedf_names(2)

    real(dp), allocatable          :: pos_data(:,:), eedf(:,:), sca_data(:)
    character(len=20), allocatable :: data_names(:)

    select case (sim_type)
    case (MODEL_part)
       write(filename, fmt="(A,I0,A)"), "output/" // trim(sim_name) // "_part_", cntr, ".txt"

       call PM_get_output(pc, pos_data, sca_data, data_names, n_pos, n_sca, time, OUT_head_density)

       if (n_sca > 0) then
          call write_line("output/" // trim(sim_name) // "_part_scalar.txt", sca_data, cntr == 1)
       end if

       if (n_pos > 0) then
          call write_data_2d(filename, pos_data, data_names(n_sca+1:n_sca+n_pos), 30)
       end if

       eedf_names(1) = "eV"
       eedf_names(2) = "eedf"

       do ix = 1, OUT_num_fields
          call PM_get_eedf(eedf, OUT_eedf_eV_range, OUT_eedf_min_max_fields(2*ix-1:2*ix), OUT_num_bins)
          write(filename, fmt="(A,I0,A,I0,A)"), "output/" // trim(sim_name) // "_part_", &
               cntr, "_eedf_", ix, ".txt"
          call write_data_2d(filename, eedf, eedf_names, 30, do_transpose = .true.)
          deallocate(eedf)
       end do

    case (MODEL_fluid_lfa)
       write(filename, fmt="(A,I0,A)"), "output/" // trim(sim_name) // "_fluid_lfa_", cntr, ".txt"
       call FL_get_output(pos_data, sca_data, data_names, n_pos, n_sca, time, OUT_head_density)

       if (n_sca > 0) then
          call write_line("output/" // trim(sim_name) // "_fluid_lfa_scalar.txt", sca_data, cntr == 1)
       end if

       if (n_pos > 0) then
          call write_data_2d(filename, pos_data, data_names(n_sca+1:n_sca+n_pos), 30)
       end if

    case (MODEL_fluid_ee)
       write(filename, fmt="(A,I0,A)"), "output/" // trim(sim_name) // "_fluid_ee_", cntr, ".txt"
       call FL_get_output(pos_data, sca_data, data_names, n_pos, n_sca, time, OUT_head_density)

       if (n_sca > 0) then
          call write_line("output/" // trim(sim_name) // "_fluid_ee_scalar.txt", sca_data, cntr == 1)
       end if

       if (n_pos > 0) then
          call write_data_2d(filename, pos_data, data_names(n_sca+1:n_sca+n_pos), 30)
       end if

    end select

    print *, "Written output " // trim(filename) // " at t = ", time
  end subroutine OUT_write_vars

  subroutine OUT_write_coeffs(pc, sim_name, sim_type)
    use m_particle_core
    use m_fluid_dd_1d, only: FL_get_coeffs
    type(PC_t), intent(in) :: pc

    character(len=*), intent(in)   :: sim_name
    integer, intent(in)            :: sim_type

    integer                        :: n_coeffs
    character(len=100)             :: filename

    real(dp), allocatable          :: coeff_data(:,:)
    character(len=20), allocatable :: coeff_names(:)

    select case (sim_type)
    case (MODEL_part)
       filename = "output/" // trim(sim_name) // "_part_coeffs.txt"
       call pc%get_coeffs(coeff_data, coeff_names, n_coeffs)

       if (n_coeffs > 0) then
          call write_data_2d(filename, coeff_data, coeff_names, 30, do_transpose = .true.)
       end if
    case (MODEL_fluid_lfa)
       filename = "output/" // trim(sim_name) // "_fluid_lfa_coeffs.txt"
       call FL_get_coeffs(coeff_data, coeff_names, n_coeffs)

       if (n_coeffs > 0) then
          call write_data_2d(filename, coeff_data, coeff_names, 30, do_transpose = .true.)
       end if
    case (MODEL_fluid_ee)
       filename = "output/" // trim(sim_name) // "_fluid_ee_coeffs.txt"
       call FL_get_coeffs(coeff_data, coeff_names, n_coeffs)

       if (n_coeffs > 0) then
          call write_data_2d(filename, coeff_data, coeff_names, 30, do_transpose = .true.)
       end if
    end select

    print *, "Written output " // trim(filename)
  end subroutine OUT_write_coeffs

  subroutine write_line(filename, line_data, new_file)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: line_data(:)
    logical, intent(in) :: new_file

    integer :: my_unit = 333

    if (new_file) then
       open(my_unit, file = trim(filename))
    else
       open(my_unit, file = trim(filename), position="APPEND")
    end if

    write(my_unit, *) line_data
    close(my_unit)
  end subroutine write_line

  subroutine write_data_2d(filename, data_2d, col_names, col_width, do_transpose)
    character(len=*), intent(in)        :: filename
    real(dp), intent(in)                :: data_2d(:,:)
    integer, intent(in)                 :: col_width
    character(len=*), intent(in) :: col_names(:)
    logical, intent(in), optional       :: do_transpose

    integer                             :: my_unit, n, n_rows, n_cols, io_state
    character(len=20)             :: fmt_string
    real(dp), allocatable               :: copy_of_data(:,:)
    logical                             :: transpose_data

    my_unit = 333
    if (present(do_transpose)) then
       transpose_data = do_transpose
    else
       transpose_data = .false.
    end if

    if (transpose_data) then
       n_rows = size(data_2d, 2)
       n_cols = size(data_2d, 1)
       allocate(copy_of_data(n_rows, n_cols))
       copy_of_data = transpose(data_2d)
    else
       n_rows = size(data_2d, 1)
       n_cols = size(data_2d, 2)
       allocate(copy_of_data(n_rows, n_cols))
       copy_of_data = data_2d
    end if

    if (size(col_names) /= n_cols) then
       print *, "write_data_2d: incompatible argument sizes"
       stop
    end if

    open(unit=my_unit, file=filename, iostat=io_state)
    if (io_state /= 0) then
       print *, "write_data_2d: error writing " // filename
       stop
    end if

    ! Create format string for the header
    write(fmt_string, fmt="(A,I0,A)") "(A", col_width, ")"

    ! Write header
    do n = 1, n_cols
       write(my_unit, advance="NO", FMT=fmt_string) "  # " // col_names(n)
    end do

    write(my_unit, *) ""

    ! Create format string for data
    write(fmt_string, fmt="(A,I0,A,I0,A,I0,A)") "(", n_cols, "E", col_width, ".", col_width - 9, "e3)"

    ! Write data
    do n = 1, n_rows
       write(my_unit, fmt_string) copy_of_data(n, :)
    end do

    close(my_unit)

  end subroutine write_data_2d

end module m_output_1d
