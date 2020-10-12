program test_outside_check
  use m_particle_core
  use m_cross_sec
  use m_units_constants

  implicit none

  integer, parameter          :: dp       = kind(0.0d0)
  character(len=*), parameter :: cs_file  = "cs_example.txt"
  character(len=*), parameter :: gas_name = "N2"

  integer, parameter      :: max_num_part  = 1000*1000
  integer, parameter      :: init_num_part = 10*1000
  integer, parameter      :: max_num_steps = 100
  integer, parameter      :: lkp_tbl_size  = 1000
  real(dp), parameter     :: delta_t       = 1.0e-11_dp
  real(dp), parameter     :: max_en_eV     = 500.0_dp
  real(dp), parameter     :: pressure      = 1.0_dp
  real(dp), parameter     :: temperature   = 300.0_dp
  real(dp), parameter     :: field_z       = -1e7_dp
  real(dp), parameter     :: part_mass     = UC_elec_mass
  real(dp)                :: init_accel(3)
  real(dp)                :: pos(3), vel(3), accel(3), weight
  real(dp)                :: neutral_dens
  integer                 :: ll, step
  type(CS_t), allocatable :: cross_secs(:)
  type(PC_t)              :: pc

  print *, "Testing m_particle_core.f90 implementation"

  ! Use constant momentum transfer cross section, so that we get a Druyvesteyn
  ! distribution
  neutral_dens = pressure * 1e5_dp / (UC_boltzmann_const * temperature)

  print *, "Reading in cross sections from ", trim(cs_file)
  call CS_add_from_file(cs_file, gas_name, neutral_dens, max_en_eV, &
       cross_secs)

  print *, "Initializing particle module"
  call pc%initialize(part_mass, max_num_part)
  call pc%use_cross_secs(max_en_eV, lkp_tbl_size, cross_secs)

  pc%outside_check => outside_check

  where (pc%colls(:)%type == CS_ionize_t)
     pc%coll_is_event(:) = .true.
  end where

  print *, "Creating initial particles"
  init_accel = [0.0_dp, 0.0_dp, field_z * UC_elec_q_over_m]

  do ll = 1, init_num_part
     pos    = 0.0_dp
     vel    = 0.0_dp
     accel  = init_accel
     weight = 1.0_dp
     call pc%create_part(pos, vel, accel, weight, 0.0_dp)
  end do

  do step = 1, max_num_steps
     call pc%advance_openmp(delta_t)
     print *, step, pc%n_part, &
          count(pc%event_list(1:pc%n_events)%ctype == PC_particle_went_out), &
          count(pc%event_list(1:pc%n_events)%ctype == CS_ionize_t), &
          minval(pc%event_list(1:pc%n_events)%part%w), &
          maxval(pc%event_list(1:pc%n_events)%part%w)
  end do

contains

  integer function outside_check(my_part)
    type(PC_part_t), intent(inout) :: my_part

    outside_check = 0

    if (norm2(my_part%x) > 5.0e-5_dp) then
       outside_check = 1
    end if
  end function outside_check

end program
