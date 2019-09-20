program test_m_particle_core
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
  real(dp), parameter     :: field_z       = -5e6_dp
  real(dp), parameter     :: part_mass     = UC_elec_mass
  real(dp)                :: init_accel(3)
  real(dp)                :: pos(3), vel(3), accel(3), weight
  real(dp)                :: neutral_dens
  integer                 :: ll, step, rng_seed(4)
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
  rng_seed = get_random_seed()
  call pc%initialize(part_mass, max_num_part, rng_seed=rng_seed)
  call pc%use_cross_secs(max_en_eV, lkp_tbl_size, cross_secs)

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
     write(*, '(F10.2,A,I10,A)') (step * 100.0_dp)/max_num_steps, '%, ', &
          pc%get_num_sim_part(), ' particles'
     call pc%advance_openmp(delta_t)
     print *, step, "total number of ionizations: ", pc%n_events
  end do

  call print_stats()

contains

  subroutine part_stats(part, vec)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(out) :: vec(:)
    vec(1:3) = part%w * part%x
    vec(4:6) = part%w * part%v
    vec(7:9) = part%w * part%a
    vec(10) = part%w * PC_v_to_en(part%v, part_mass)
    vec(11) = part%w
  end subroutine part_stats

  subroutine print_stats()
    integer :: n_part
    real(dp) :: sum_x(3), sum_v(3), sum_a(3), sum_en
    real(dp) :: sum_weight
    real(dp) :: sum_vec(11)

    real(dp), parameter :: x3_stored        = 1.9793569827856038e-4_dp
    real(dp), parameter :: n_stored         = 64106.0_dp
    real(dp), parameter :: max_deviation(2) = [2.0e-3_dp, 2.0e-2_dp]
    real(dp)            :: rel_deviation(2)

    n_part = 0

    n_part = n_part + pc%get_num_sim_part()
    call pc%compute_vector_sum(part_stats, sum_vec)

    sum_weight = sum_vec(11)
    sum_x = sum_vec(1:3)
    sum_v = sum_vec(4:6)
    sum_a = sum_vec(7:9)
    sum_en =  sum_vec(10)

    print *, ""
    write(*, '(A,3E12.4)') " mean(x):", sum_x / sum_weight
    write(*, '(A,E12.4)') " n_particles", sum_weight
    rel_deviation = [abs(x3_stored - sum_x(3) / sum_weight) / x3_stored, &
         abs(n_stored - sum_weight) / n_stored]
    write(*, '(A,2E12.4)') " Rel. diff. with stored run: ", rel_deviation
    if (all(rel_deviation < max_deviation)) then
       print *, "PASS"
    else
       print *, "Deviation too large"
       error stop "FAIL"
    end if
  end subroutine print_stats

  !> Get a random seed based on the current time
  function get_random_seed() result(seed)
    integer :: seed(4)
    integer :: time, i

    call system_clock(time)
    do i = 1, 4
       seed(i) = ishftc(time, i*8)
    end do
  end function get_random_seed

end program test_m_particle_core
