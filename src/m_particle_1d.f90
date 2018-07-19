module m_particle_1d
  use m_generic
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp          = kind(0.0d0)
  integer, parameter :: PM_num_vars = 3
  integer, parameter :: iv_elec  = 1, iv_pion = 2, iv_en = 3

  integer, parameter :: num_scalars = 2
  integer, parameter :: i_lbound_elec = 1
  integer, parameter :: i_rbound_elec = 2

  real(dp), allocatable :: PM_vars(:, :)
  real(dp), allocatable :: PM_scalars(:)

  real(dp) :: PM_transverse_area
  real(dp) :: PM_grid_volume, PM_inv_grid_volume
  real(dp) :: PM_part_per_cell, PM_v_merge_weight
  real(dp) :: merge_factor
  real(dp) :: merge_prev_npart

  character(len=10), allocatable :: gas_components(:)
  real(dp), allocatable          :: gas_fractions(:)

  type(PC_t)        :: pc
  type(PC_events_t) :: events

  ! Public routines
  public :: particle_initialize
  public :: particle_advance
  public :: particle_write_output

contains

  subroutine particle_initialize(cfg)
    use m_config
    use m_random
    use m_units_constants
    use m_cross_sec

    type(CFG_t), intent(inout) :: cfg

    character(len=200)      :: cs_file
    type(CS_t), allocatable :: cross_secs(:)
    integer                 :: n_gas_comp
    integer                 :: n_part_init, n_part_max
    integer                 :: n, i, num_part, tbl_size
    real(dp)                :: x, sum_elec_dens
    real(dp)                :: elec_dens, ion_dens
    real(dp)                :: pos(3), vel(3), accel(3), max_eV

    call CFG_add(cfg, "gas%components", ["N2"], &
         "Gas components", .true.)
    call CFG_add(cfg, "gas%fractions", [1.0_dp], &
         "Gas components fractions", .true.)

    call CFG_get_size(cfg, "gas%components", n_gas_comp)
    allocate(gas_components(n_gas_comp))
    allocate(gas_fractions(n_gas_comp))

    call CFG_get(cfg, "gas%components", gas_components)
    call CFG_get(cfg, "gas%fractions", gas_fractions)

    cs_file = "input/n2_cross_sections_siglo.txt"
    call CFG_add_get(cfg, "gas%cross_sections", cs_file, &
         "File with cross section data")

    max_eV = 1e3_dp
    call CFG_add_get(cfg, "particle%max_eV", max_eV, &
         "Maximum energy (eV) for particles")

    do n = 1, n_gas_comp
       call CS_add_from_file(trim(cs_file), trim(gas_components(n)), &
            gas_fractions(n) * gas_number_dens, max_eV, cross_secs)
    end do

    if (size(cross_secs) < 1) &
         error stop "No cross sections given to particle module!"

    ! call CS_write_summary(cross_secs, trim(output_name) // "_cs.txt")

    n_part_init = 10*1000
    call CFG_add_get(cfg, "particle%initial_number", n_part_init, &
         "The number of initial simulation particles")

    n_part_max = 10*1000*1000
    call CFG_add_get(cfg, "particle%_max_number", n_part_max, &
         "The maximum number of particles")

    sum_elec_dens = 0.0_dp
    do n = 1, domain_nx
       x             = (n - 0.5_dp) * domain_dx
       elec_dens     = init_elec_dens(x)
       sum_elec_dens = sum_elec_dens + elec_dens
    end do

    if (sum_elec_dens < epsilon(1.0_dp)) &
         error stop "Initial electron density seems to be zero!"

    PM_grid_volume     = n_part_init / sum_elec_dens
    PM_transverse_area = PM_grid_volume / domain_dx
    PM_inv_grid_volume = 1 / PM_grid_volume

    PM_part_per_cell = 2.0e2_dp
    call CFG_add_get(cfg, "particle%part_per_cell", PM_part_per_cell, &
         "The desired number of particles per cell")

    PM_v_merge_weight = 1.0e-12_dp
    call CFG_add_get(cfg, "particle%v_merge_weight", PM_v_merge_weight, &
         "Velocity weight factor (s) for particle merging")

    merge_factor = 1.4_dp
    call CFG_add_get(cfg, "particle%merge_increase_factor", &
         merge_factor, &
         "Adapt weights if particle count increases by this factor")


    allocate(PM_vars(0:domain_nx+1, PM_num_vars))
    PM_vars(:, :) = 0.0_dp
    allocate(PM_scalars(num_scalars))
    PM_scalars(:) = 0.0_dp

    tbl_size = 1000
    call CFG_add_get(cfg, "particle%table_size", tbl_size, &
         "Size of the lookup table for collision rates")

    pc%accel_function => accel_func
    print *, "init", tbl_size, max_eV, n_part_max
    call pc%initialize(UC_elec_mass, cross_secs, tbl_size, max_eV, n_part_max)

    pc%outside_check => outside_check

    where (pc%colls(:)%type == CS_ionize_t .or. &
         pc%colls(:)%type == CS_attach_t)
       pc%coll_is_event(:) = .true.
    end where

    ! Initialization of electron and ion density
    do n = 1, domain_nx
       x         = (n - 0.5_dp) * domain_dx
       elec_dens = init_elec_dens(x)
       ion_dens  = init_pos_ion_dens(x)

       num_part = nint(elec_dens * PM_grid_volume)

       ! Uniform position in the cell
       pos(2:3) = [0.0_dp, 0.0_dp]
       vel(:)   = 0
       ! Accel is set after computing electric field
       accel(:) = 0.0_dp

       do i = 1, num_part
          pos(1) = x + (i-0.5_dp) * domain_dx / num_part
          call pc%create_part(pos, vel, accel, 1.0_dp, 0.0_dp)
       end do

       pm_vars(n, iv_pion) = num_part / PM_grid_volume
    end do

    call PM_update_efield(0.0_dp)
    call pc%set_accel()

    merge_prev_npart = pc%get_num_sim_part()
  end subroutine particle_initialize

  subroutine set_elec_density()
    PM_vars(:, iv_elec) = 0.0_dp
    call pc%loop_ipart(add_elec_to_dens)
  end subroutine set_elec_density

  subroutine set_elec_en_density()
    PM_vars(:, iv_en) = 0.0_dp
    call pc%loop_ipart(add_elec_en_to_dens)
  end subroutine set_elec_en_density

  subroutine add_to_dens_cic(amount, x, iv)
    real(dp), intent(in) :: amount, x
    integer, intent(in)  :: iv

    integer                 :: low_ix
    real(dp)                :: tmp, dens_dif

    ! Do linear interpolation, note that dens(1) is defined at 0.5 * dx
    tmp    = x * domain_inv_dx - 0.5_dp
    low_ix = floor(tmp) + 1
    ! Coefficient for upper point between 0 and 1
    tmp    = low_ix - tmp

    dens_dif = amount * PM_inv_grid_volume

    if (low_ix < 0 .or. low_ix > domain_nx) then
       error stop "a particle is outside the computational domain"
    end if

    if (low_ix == 0) then
       PM_vars(1, iv) = PM_vars(1, iv) + dens_dif
    else if (low_ix == domain_nx) then
       PM_vars(domain_nx, iv) = PM_vars(domain_nx, iv) + dens_dif
    else
       PM_vars(low_ix, iv)   = PM_vars(low_ix, iv) + (1.0_dp - tmp) * dens_dif
       PM_vars(low_ix+1, iv) = PM_vars(low_ix+1, iv) + tmp * dens_dif
    end if
  end subroutine add_to_dens_cic

  subroutine add_elec_to_dens(my_part)
    type(PC_part_t), intent(in) :: my_part
    call add_to_dens_cic(my_part%w, my_part%x(1), iv_elec)
  end subroutine add_elec_to_dens

  subroutine add_elec_en_to_dens(my_part)
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: energy
    energy = 0.5_dp * my_part%w * UC_elec_mass * sum(my_part%v**2)
    call add_to_dens_cic(energy, my_part%x(1), iv_en)
  end subroutine add_elec_en_to_dens

  subroutine PM_update_efield(time)
    use m_units_constants
    real(dp), intent(in) :: time
    real(dp) :: source(1:domain_nx)
    real(dp) :: surface_charge(2)

    call set_elec_density()

    source = UC_elem_charge * (PM_vars(1:domain_nx, iv_pion) - &
         PM_vars(1:domain_nx, iv_elec))
    surface_charge = -UC_elem_charge * PM_scalars([i_lbound_elec, i_rbound_elec])

    call compute_field(source, surface_charge, time)
  end subroutine PM_update_efield

  function accel_func(my_part) result(accel)
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: accel(3)
    real(dp), parameter         :: accel_fac = UC_elec_charge / UC_elec_mass

    ! Only acceleration in the z-direction
    accel(1) = accel_fac * get_field_at(my_part%x(1))
    accel(2:3) = 0.0_dp
  end function accel_func

  subroutine particle_advance(dt, time)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: time

    call pc%advance_openmp(dt, events)
    call handle_events(events)
    call PM_update_efield(time)
    call pc%after_mover(dt)

    if (pc%get_num_sim_part() > merge_factor * merge_prev_npart) then
       call PM_adapt_weights()
       merge_prev_npart = pc%get_num_sim_part()
       call PM_update_efield(time)
       call pc%set_accel()
    end if
  end subroutine particle_advance

  subroutine handle_events(events)
    use m_cross_sec
    type(PC_events_t), intent(inout) :: events
    integer                          :: n

    do n = 1, events%n_stored
       select case (events%list(n)%ctype)
       case (CS_ionize_t)
          call add_to_dens_cic(events%list(n)%part%w, &
               events%list(n)%part%x(1), iv_pion)
       case (CS_attach_t)
          call add_to_dens_cic(-events%list(n)%part%w, &
               events%list(n)%part%x(1), iv_pion)
       case (PC_particle_went_out)
          if (events%list(n)%part%x(1) < 0.0_dp) then
             PM_scalars(i_lbound_elec) = PM_scalars(i_lbound_elec) + &
                  events%list(n)%part%w / PM_transverse_area
          else
             PM_scalars(i_rbound_elec) = PM_scalars(i_rbound_elec) + &
                  events%list(n)%part%w / PM_transverse_area
          end if
       case default
          print *, events%list(n)%ctype
          error stop "Unknown event"
       end select
    end do

    events%n_stored = 0
  end subroutine handle_events

  subroutine PM_adapt_weights()
    real(dp) :: v_fac
    v_fac = PM_v_merge_weight * domain_dx / PM_part_per_cell

    call pc%merge_and_split([.true., .false., .false.], v_fac, .true., &
         get_desired_weight, domain_dx, PC_merge_part_rxv, PC_split_part)
    ! write(*,'(A, E12.4, A, I0)') " After merging real/sim part: ", &
    !      pc%get_num_real_part(), " / ", pc%get_num_sim_part()
  end subroutine PM_adapt_weights

  real(dp) function get_desired_weight(my_part)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: elec_dens
    integer                     :: ix

    ix = floor(my_part%x(1) * domain_inv_dx)
    elec_dens = PM_vars(ix, iv_elec)
    get_desired_weight = elec_dens * PM_grid_volume / PM_part_per_cell
    get_desired_weight = max(1.0_dp, get_desired_weight)
  end function get_desired_weight

  subroutine particle_write_output(base_fname, time, ix)
    use m_units_constants
    character(len=*), intent(in) :: base_fname
    real(dp), intent(in)         :: time
    integer, intent(in)          :: ix
    integer                      :: n, nx, my_unit
    character(len=200)           :: fname
    real(dp)                     :: total_charge

    nx = domain_nx
    call set_elec_en_density()

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_particle_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    do n = 1, domain_nx
       write(my_unit, *) domain_dx * (n-0.5_dp), &
            field_cc(n), &
            PM_vars(n, iv_elec), &
            PM_vars(n, iv_pion), &
            PM_vars(n, iv_en)
    end do
    close(my_unit)

    write(*, "(A,E9.2,A,A)") " t = ", time, " wrote ", trim(fname)

    write(fname, "(A,A,I0.6,A)") trim(base_fname), "_particle_scalars.txt"
    if (ix == 1) then
       open(newunit=my_unit, file=trim(fname))
       write(my_unit, *) "# Header TODO"
    else
       open(newunit=my_unit, file=trim(fname), access='append')
    end if

    total_charge = domain_dx * sum(PM_vars(1:nx, iv_pion) &
         - PM_vars(1:nx, iv_elec)) &
         - PM_scalars(i_lbound_elec) &
         - PM_scalars(i_rbound_elec)

    write(my_unit, *) time, PM_scalars(i_lbound_elec), &
         PM_scalars(i_rbound_elec), total_charge
    close(my_unit)

  end subroutine particle_write_output

  ! subroutine PM_get_eedf(eedf, energy_range, field_range, n_bins)
  !   use m_units_constants
  !   real(dp), intent(inout), allocatable :: eedf(:,:)
  !   real(dp), intent(in) :: energy_range(2), field_range(2)
  !   integer, intent(in) :: n_bins

  !   integer :: ix
  !   real(dp) :: accel_range(2)
  !   allocate(eedf(2, n_bins))
  !   do ix = 1, n_bins
  !      eedf(1, ix) = energy_range(1) + (ix-0.5_dp) * &
  !           (energy_range(2) - energy_range(1)) / n_bins
  !   end do

  !   ! Convert to range of velocity**2
  !   eedf(1,:) = eedf(1,:) * UC_elec_volt / (0.5_dp * pc%get_mass())
  !   accel_range = field_range * UC_elem_charge / UC_elec_mass

  !   call pc%histogram(get_vel_squared, select_part_by_accel, &
  !        accel_range, eedf(1,:), eedf(2,:))

  !   ! Convert to energy range
  !   eedf(1,:) = eedf(1,:) * 0.5_dp * pc%get_mass() / UC_elec_volt
  ! end subroutine PM_get_eedf

  ! real(dp) function get_vel_squared(my_part)
  !   type(PC_part_t), intent(in) :: my_part
  !   get_vel_squared = sum(my_part%v**2)
  ! end function get_vel_squared

  ! logical function select_part_by_accel(my_part, accel_range)
  !   type(PC_part_t), intent(in) :: my_part
  !   real(dp), intent(in)        :: accel_range(:)
  !   real(dp)                    :: accel
  !   accel                = norm2(my_part%a)
  !   select_part_by_accel = accel >= accel_range(1) .and. accel <= accel_range(2)
  ! end function select_part_by_accel

  pure integer function outside_check(my_part)
    type(PC_part_t), intent(in) :: my_part

    if (my_part%x(1) < 0.0_dp .or. my_part%x(1) > domain_length) then
       outside_check = 1
    else
       outside_check = 0
    end if
  end function outside_check

end module m_particle_1d
