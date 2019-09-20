program ion_drift
  use m_particle_core
  use m_units_constants
  use m_cross_sec

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  real(dp), parameter :: part_mass     = 32 * UC_atomic_mass
  integer, parameter  :: max_num_part  = 100*1000
  integer, parameter  :: init_num_part = 10000
  integer, parameter  :: max_num_steps = 100
  real(dp), parameter :: delta_t       = 1e-10
  real(dp), parameter :: gas_temp      = 300.0_dp
  real(dp), parameter :: gas_N0        = 2.5e25_dp
  real(dp)            :: max_en_eV, mobility
  integer             :: i
  type(rate_func_t)   :: rate_funcs(1)
  type(PC_t)          :: pc
  type(PC_part_t)     :: my_part

  mobility  = 7.1e21_dp / gas_N0
  max_en_eV = PC_speed_to_en(1e7 * mobility, part_mass) / UC_elec_volt

  rate_funcs(1)%coll%type = CS_attach_t
  rate_funcs(1)%ptr => my_coll_rate

  pc%particle_mover => PC_tracer_advance_midpoint
  pc%tracer_velocity => velocity

  call pc%initialize(part_mass, max_num_part)
  call pc%use_rate_funcs(max_en_eV, 100, rate_funcs)

  pc%coll_is_event(:) = .true.

  do i = 1, init_num_part
     my_part%x = 0.0_dp
     my_part%v = velocity(my_part)
     my_part%a = 0.0_dp
     my_part%w = 1.0_dp
     call pc%add_part(my_part)
  end do

  do i = 1, 100
     call pc%advance_openmp(delta_t)
     print *, i, "events/n_part: ", pc%n_events, pc%n_part
  end do

contains

  real(dp) function my_coll_rate(v)
    real(dp), intent(in) :: v
    real(dp)             :: theta, k0, delta_e

    theta        = 0.5_dp * acos(-1.0_dp) * part_mass * v**2 + &
         UC_boltzmann_const * gas_temp
    k0           = 1.22e-11_dp * 1e-6_dp
    delta_e      = 0.78_dp * UC_elec_volt
    my_coll_rate = k0 * gas_N0 * exp(-delta_e / theta)
  end function my_coll_rate

  function velocity(part) result(v)
    type(PC_part_t), intent(inout) :: part
    real(dp)                       :: v(3)

    v = [1500.0_dp, 0.0_dp, 0.0_dp]
  end function velocity

end program ion_drift
