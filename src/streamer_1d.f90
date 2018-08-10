!> Main program to simulate discharges in 1D, using either a particle or a fluid
!> model
program streamer_1d
  use m_config
  use m_generic
  use m_fluid_1d
  use m_particle_1d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  character(len=200) :: output_name

  integer :: it, end_iteration
  integer :: output_ix
  integer :: info_cntr

  real(dp) :: time, end_time
  real(dp) :: dt, dt_next, dt_max, dt_min, dt_limit
  real(dp) :: dt_output
  real(dp) :: time_elapsed
  integer  :: time_now, time_start, count_rate

  logical     :: write_output
  type(CFG_t) :: cfg

  call CFG_update_from_arguments(cfg)

  end_time = 3.0e-9_dp
  call CFG_add_get(cfg, "end%time", end_time, "End time (s)")

  end_iteration = huge(1)
  call CFG_add_get(cfg, "end%iteration", end_iteration, "End iteration")

  dt_next = 1e-12_dp
  call CFG_add_get(cfg, "dt%initial", dt_next, "Initial time step (s)")

  dt_max = 1e-11_dp
  call CFG_add_get(cfg, "dt%max", dt_max, "Maximal time step (s)")

  dt_min = 1e-13_dp
  call CFG_add_get(cfg, "dt%min", dt_max, "Minimal time step (s)")

  output_name = "output/sim"
  call CFG_add_get(cfg, "output%filename", output_name, &
       "Base file name for output")

  dt_output = 1.0e-10_dp
  call CFG_add_get(cfg, "output%dt", dt_output, &
         "The time step for writing output")

  call generic_initialize(cfg)
  call fluid_initialize(cfg)
  call particle_initialize(cfg)

  call check_output_folder(trim(output_name))

  if (model_type == model_particle) then
     call CFG_write(cfg, trim(output_name) // "_particle.cfg", .true.)
  else
     call CFG_write(cfg, trim(output_name) // "_fluid.cfg", .true.)
  end if

  ! Initialize variables
  time      = 0.0_dp
  output_ix = 0
  info_cntr = 1
  it        = 0

  call system_clock(time_start, count_rate)

  ! Here the simulation starts
  do while (time < end_time .and. it < end_iteration)
     call system_clock(time_now)
     time_elapsed = (time_now - time_start) / real(count_rate, dp)

     if (time_elapsed > 10 * info_cntr) then
        info_cntr = info_cntr + 1
        write(*, "(F6.2,A,E9.2)") 100.0_dp * time/end_time, &
             "% done, dt: ", dt_next
     end if

     dt = dt_next
     write_output = (time + dt >= output_ix * dt_output)

     if (write_output) then
        dt        = output_ix * dt_output - time
        output_ix = output_ix + 1
     end if

     if (model_type == model_fluid) then
        call fluid_advance(dt, time, dt_limit)
        dt_next = get_new_dt(dt_next, dt_limit)

        if (write_output) then
           call fluid_write_output(output_name, time + dt, &
                dt_next, output_ix)
        end if
     else
        call particle_advance(dt, time, dt_limit)
        dt_next = get_new_dt(dt_next, dt_limit)

        if (write_output) then
           call particle_write_output(output_name, time + dt, &
                dt_next, output_ix)
        end if
     end if

     it   = it + 1
     time = time + dt
  end do

  write(*, "(A,E10.4,A)") "Simulation ended after ", 1.0d9 * time, " ns"

contains

  !> Check whether the output folder exists and is writable
  subroutine check_output_folder(filename)
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, iostate

    open(newunit=my_unit, file=trim(filename)//"DUMMY", iostat=iostate)

    if (iostate /= 0) then
       print *, "Output base file name: ", trim(filename)
       error stop "Cannot write to output folder"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_output_folder

  !> Get the next time step
  real(dp) function get_new_dt(dt, dt_limit)
    real(dp), intent(in) :: dt, dt_limit

    if (dt > dt_limit) then
       get_new_dt = dt_limit
    else
       ! Increase time step at most by 10%
       get_new_dt = min(1.1_dp * dt, dt_limit)
    end if
  end function get_new_dt

end program
